############################ SETUP #############################################

rm(list = ls())

library(tidyverse)
library(readr)
library(sf)
library(rWCVP)
#install.packages("rWCVPdata", repos=c("https://matildabrown.github.io/drat", "https://cloud.r-project.org"))
library(rWCVPdata)
library(CoordinateCleaner)
library(conflicted)
library(patchwork)
library(here)
library(RANN)

source(here("rscripts/Cotoneaster review functions.R"))
source(here("rscripts/nn_dist_to_arc.R"))

theme_set(theme_classic())

conflict_scout()
conflicts_prefer(dplyr::filter,
                 dplyr::lag)

######################## LOAD WCVP DATA ########################################

wcvp_distributions<-rWCVPdata::wcvp_distributions

wcvp_names<-rWCVPdata::wcvp_names

country_shapefiles<-wgsrpd3%>%
  st_set_crs(4326)

###################### LOAD DATA FROM GBIF AND GLONAF ##########################

gbif_raw <- read_delim(here("data/GBIF downloads/0075125-260226173443078.csv"), 
                       delim = "\t", 
                       escape_double = FALSE,
                       trim_ws = TRUE)


#GloNAF taxon by region data

glonaf_flora2<-
  read.csv2(here("data/GloNAF 2024/glonaf_flora2.csv"))

glonaf_taxon<-
  read.csv2(here("data/GloNAF 2024/glonaf_taxon_wcvp.csv"))

glonaf_region<-
  read.csv2(here("data/GloNAF 2024/glonaf_region.csv"))

glonaf_taxon_x_region<-
  left_join(glonaf_flora2,
            glonaf_taxon%>%
              select(id,accepted_plant_name_id,taxa_accepted)%>%
              rename(taxon_wcvp_id = id))%>%
  filter(startsWith(taxa_accepted,"Cotoneaster"))%>%
  left_join(glonaf_region%>%
              select(id,code,name,tdwg4_id)%>%
              rename(region_id = id))%>%
  mutate(tdwg3 = str_sub(tdwg4_id,1,3)
         ,IDregion = code)%>%
  mutate(tdwg3 = case_when(tdwg3 != "000" ~ tdwg3,
                           name == "Italy" ~ "ITA",
                           name == "South Africa" ~ "CPP,NAT,OFS,TVL",
                           name == "New Zealand" ~ "ATP,CTM,KER,NZN,NZS",
                           startsWith(name,"United States") ~ "ALA,ALU,ARI,ARK,CAL,CNT,COL,DEL,FLA,GEO,IDA,ILL,INI,IOW,KAN,KTY,LOU,MAI,MAS,MIC,MIN,MNT,MRY,MSI,MSO,NCA,NDA,NEB,NEV,NWH,NWJ,NWM,NWY,OHI,OKL,ORE,PEN,RHO,SCA,SDA,TEN,TEX,UTA,VER,VRG,WAS,WDC,WIS,WVA,WYO",
                           startsWith(name,"France") == T ~ "FRA",
                           startsWith(name,"United Kingdom") == T ~ "GRB",
                           name == "Turkey" ~ "TUE,TUR",
                           name == "Papua New Guinea" ~ "NWG",
                           name == "Australia" ~ "NFK,NSW,NTA,QLD,SOA,TAS,VIC,WAU",
                           startsWith(name,"Russia") == T ~ "ALT,AMU,RUC,RUN,RUS,RUW,CTA,KAM,KHA,KRA,KUR,MAG,NCS,PRM,RUE,SAK,TVA,WSB,YAK,BRY,IRK",
                           TRUE ~ NA))%>%
  separate_rows(tdwg3,sep = ",")


# WCVP_native distributions 
wcvp_native<-
  wcvp_names%>%
  filter(genus == "Cotoneaster")%>%
  left_join(wcvp_distributions)#%>%view()





###################### CLEAN GBIF DATA #########################################


gbif_clean<-
  gbif_raw%>%
  select(gbifID,
         genus,species,
         taxonRank,
         decimalLatitude,decimalLongitude,
         day,month,year,
         establishmentMeans
  )%>%
  clean_coordinates()%>%
  filter(.summary)%>%
  select(-(.val:.summary))





######################### CREATE TRANSLATION TABLE #############################


translation_table<-
  gbif_clean%>%
  mutate(taxon = ifelse(taxonRank == "GENUS",genus,species))%>%
  group_by(taxon)%>%
  summarise(num_recs = length(species)
  )%>%
  wcvp_match_names(.,
                   name_col="taxon",
                   author_col=NULL,
                   fuzzy=TRUE,
                   progress_bar=TRUE
  )%>%
  multi_correct()%>%
  mutate(plant_name_id = wcvp_accepted_id
  )%>%
  left_join(wcvp_names%>%
              select(genus,species,taxon_name, taxon_authors,plant_name_id)
  )%>%
  select(-plant_name_id)%>%
  rename(wcvp_accepted_name = taxon_name
  )%>%
  filter(multi_correct == "keeper" | multi_correct == "keeper_checked")%>%
  select(taxon,num_recs,wcvp_accepted_id,wcvp_accepted_name,genus,species,taxon_authors)%>%
  rename(taxon_stand = taxon,
         num_recs_stand = num_recs,
         plant_name_id_stand = wcvp_accepted_id,
         genus_stand = genus,
         species_stand = species,
         authors_stand = taxon_authors)%>%
  filter(!is.na(wcvp_accepted_name))%>%
  distinct()

# any names not in the translation table will be removed from further analysis


# Checks that every entry in the translation table is unique
length(unique(translation_table$taxon_stand)) == nrow(translation_table)





######################## STANDARDIZE NOMENCLATURE ##############################


gbif_geo<-
  gbif_clean%>%
  mutate(taxon_stand = ifelse(taxonRank == "GENUS",genus,species))%>%
  left_join(translation_table)%>%
  
  
  
  
  
######################### MAKE GBIF INTO SHAPEFILE ##############################


st_as_sf(coords = c("decimalLongitude","decimalLatitude"),crs = 4326)%>%
  st_join(country_shapefiles) %>%
  
  filter(!is.na(LEVEL3_COD))





############################ FINAL PRODUCTS ####################################


gbif_final<-
  gbif_geo%>%
  
  mutate(tdwg3 = LEVEL3_COD)%>%
  
  left_join(glonaf_taxon_x_region%>%
              group_by(status,tdwg3,taxa_accepted)%>%
              summarise(obs = length(status))%>%
              group_by(tdwg3,taxa_accepted)%>%
              filter(!(status == "Naturalized" & "Invasive" %in% status))%>%
              rename(wcvp_accepted_name = taxa_accepted)%>%
              ungroup())%>%
  left_join(wcvp_native%>%
              mutate(plant_name_id_stand = accepted_plant_name_id,
                     tdwg3 = area_code_l3)%>%
              select(plant_name_id_stand,tdwg3,introduced)%>%
              filter(!is.na(tdwg3))%>%
              unique())%>%
  mutate(glonaf_status = status,
         wcvp_status = case_when(introduced == 1 ~ "introduced",
                                 introduced == 0 ~ "native",
                                 TRUE ~ NA))%>%
  select(-status,-introduced,-taxon_stand,-obs)%>%
  cbind(st_coordinates(.))%>%
  st_drop_geometry()%>%
  rename(decimalLongitude = X,
         decimalLatitude = Y)%>%
  filter(!is.na(wcvp_accepted_name),
         startsWith(wcvp_accepted_name,"Cotoneaster"))%>%
  
  
  
  
  
########################  ALLOCATING NATIVENESS  ###############################


  mutate(nativeness = case_when((glonaf_status == "Naturalized"|
                                 glonaf_status == "Invasive") &
                                (wcvp_status == "introduced") ~ "non-native",
                              is.na(glonaf_status) & wcvp_status == "native" ~ "native",
                              is.na(glonaf_status) & wcvp_status == "introduced"~ "non-native"))%>%
  mutate(nn_nativeness = NA,
         nn_id = NA,
         nn_dist = NA)%>%
  mutate(x = sphere(lat = decimalLatitude,lon = decimalLongitude)[[1]],
         y = sphere(lat = decimalLatitude,lon = decimalLongitude)[[2]],
         z = sphere(lat = decimalLatitude,lon = decimalLongitude)[[3]])  




# This algorithm from the RANN package uses euclidean distance, therefore coordinates are transformed into 3d coordinates and distances are then back-transformed using arc_distance

for (i in seq_len(nrow(gbif_final))) {
  
  if (is.na(gbif_final$nativeness[i])) {
    
    subset_df <- gbif_final %>%
      filter(!is.na(nativeness),
             species_stand == gbif_final$species_stand[i])
    
    if (nrow(subset_df) == 0) next
    
    coords <- subset_df %>%
      select(x, y, z) %>%
      as.matrix()
    
    query <- matrix(c(gbif_final$x[i],
                      gbif_final$y[i],
                      gbif_final$z[i]),
                    nrow = 1)
    
    nn <- nn2(data = coords, query = query, k = 1)
    
    nearest_point <- subset_df[nn$nn.idx, ]
    
    gbif_final$nn_nativeness[i] <- nearest_point$nativeness
    gbif_final$nn_id[i] <- nn$nn.idx
    gbif_final$nn_dist[i] <- nn_dist_to_arc(nn$nn.dists)
  }
}

############################# EXPORTING ########################################

write_csv(gbif_final,here("data","Cotoneaster review output","Final gbif_records_090426.csv"))



########################## CODE GRAVEYARD ######################################


## Plotting records of different species ####

# species_list<-result%>%
#   group_by(species_stand)%>%
#   summarise(n = length(species_stand))

# results<-list()

# for(i in 1:length(species_list$species_stand)){
#   
#   plot<-
#     ggplot()+
#     annotation_borders()+
#     #geom_sf(data = country_shapefiles)+
#     geom_sf(data = result%>%
#               mutate(nativeness = ifelse(!is.na(nativeness),nativeness,nn_nativeness))%>%
#               mutate(nativeness = ifelse(nativeness == 0,"unknown",nativeness))%>%
#               filter(#nn_dist<10,
#                      species_stand == species_list$species_stand[i])%>%
#               st_as_sf(coords = c("decimalLongitude","decimalLatitude"),crs = 4326),
#             aes(colour = nativeness))+
#     labs(title = paste("Cotoneaster",species_list$species_stand[i]))
#   
#   results[[i]]<-plot
#   
# }

# ggplot()+
#   #geom_sf(data = country_shapefiles)+
#   annotation_borders()+
#   geom_sf(data = result%>%
#             mutate(nativeness = ifelse(!is.na(nativeness),nativeness,nn_nativeness))%>%
#             mutate(nativeness = ifelse(nativeness == 0 ,"unknown",nativeness))%>%
#             st_as_sf(coords = c("decimalLongitude","decimalLatitude"),crs = 4326),
#           aes(colour = nativeness))



# results[[runif(1,1,151)]]+
#   theme(legend.position = "none")+
#   results[[runif(1,1,151)]]+
#   theme(legend.position = "none")+
#   results[[runif(1,1,151)]]+
#   theme(legend.position = "none")+
#   results[[runif(1,1,151)]]+
#   theme(legend.position = "none")+
#   results[[runif(1,1,151)]]+
#   theme(legend.position = "none")+
#   results[[runif(1,1,151)]]+
#   theme(legend.position = "none")



## Final spreadsheet ###########################################################


# final_spreadsheet<-
#   gbif_final%>%
#   group_by(wcvp_accepted_name,tdwg3,glonaf_status,wcvp_status)%>%
#   summarise(earliestYear = custom_min(year),
#             latestYear = custom_max(year),
#             obs = length(wcvp_accepted_name))%>%
#   left_join(country_shapefiles%>%select(LEVEL3_NAM,LEVEL3_COD)%>%
#               rename(tdwg3 = LEVEL3_COD,
#                      tdwg3_name = LEVEL3_NAM)%>%
#               st_drop_geometry())%>%
#   st_drop_geometry()


## Exporting ###################################################################

#write.csv2(st_drop_geometry(final_spreadsheet),"Cotoneaster review output/Final summary spreadsheet_050225.csv")
#write.csv2(st_drop_geometry(gbif_final),"Cotoneaster review output/Final gbif records_050225.csv")
# write.csv2(translation_table,"translation_table.csv")

# ggplot(data = gbif_final,aes(wcvp_status,glonaf_status))+ 
#   geom_count(show.legend = T)+ 
#   facet_wrap(~wcvp_accepted_name)+
#   geom_text(aes(label = after_stat(n)),
#             stat = "sum",
#             vjust = -0.5,
#             size = 1)
# 


# gbif_final%>%
#   rgl::plot3d(x = .$x,y = .$y,z = .$z,colour = .$species_stand)