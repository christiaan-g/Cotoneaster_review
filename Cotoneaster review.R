
## Setup #######################################################################

rm(list = ls())
#setwd("~/Universiteit werk/Masters/Cotoneaster review")
#setwd("~/Cotoneaster_review")
library(tidyverse)
library(readr)
library(sf)
library(rWCVP)
#install.packages("rWCVPdata", repos=c("https://matildabrown.github.io/drat", "https://cloud.r-project.org"))
library(rWCVPdata)
library(CoordinateCleaner)
library(conflicted)
#library(tictoc)
library(patchwork)
#library(readxl)
library(here)
source("sphere.R")
theme_set(theme_classic())

conflict_scout()
conflicts_prefer(dplyr::filter,
                 dplyr::lag)


## Load WCVP data ##############################################################

wcvp_distributions<-rWCVPdata::wcvp_distributions

wcvp_names<-rWCVPdata::wcvp_names

country_shapefiles<-wgsrpd3%>%
  st_set_crs(4326)

# ggplot(country_shapefiles)+
#   geom_sf()


## Load data from gbif and GloNAF #########################################################

#gbif data
gbif <- read_delim("0006576-241126133413365.csv", 
                   delim = "\t", 
                   escape_double = FALSE,
                   trim_ws = TRUE)

#GloNAF taxon by region data

#glonaf_taxon_x_region_orig<-
#  read.csv2("257_9_257_2_GLONAF/GLONAF/Taxon_x_List_GloNAF_vanKleunenetal2018Ecology.csv")

glonaf_flora2<-
  read.csv2("GloNAF 2024/glonaf_flora2.csv")

glonaf_taxon<-
  read.csv2("GloNAF 2024/glonaf_taxon_wcvp.csv")

glonaf_region<-
  read.csv2("GloNAF 2024/glonaf_region.csv")

# glonaf_shapefiles<-
#   read_sf("GloNAF 2024/glonaf_2024_nol/glonaf_2024_nol/glonaf_2024_nol.shp")

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


# wcvp_native%>%
#   filter(taxon_status == "Accepted")%>%
#   mutate(wcvp_status = ifelse(introduced == 0,"native","introduced"))%>%
#   select(accepted_plant_name_id,area_code_l3,wcvp_status)%>%
#   rename(tdwg3 = area_code_l3)%>%view()


## Define multi_correct function ###############################################

multi_correct<-function (df){
  
  #synonym_codes <- c("Synonym", "Orthographic", "Artificial Hybrid", "Unplaced")
  
  synonym_codes <- c("Synonym")
  
  accepted_matches <-
    df %>%
    group_by(taxon)%>%
    mutate(n_matches = n())%>%
    mutate(multi_correct = case_when(multiple_matches == FALSE ~ "keeper",
                                     multiple_matches == TRUE & wcvp_status == "Accepted" ~ "keeper",
                                     TRUE ~ "remove"))%>%
    mutate(multi_correct = case_when(multiple_matches == TRUE & 
                                       wcvp_status%in%synonym_codes == TRUE & 
                                       "keeper" %in% multi_correct == FALSE ~ "keeper",
                                     TRUE ~ multi_correct))%>%
    ungroup()
  
  #Shows species for which no legitimate name was found
  no_legit_match<-
    df%>%
    distinct(taxon)%>%
    mutate(keeper_found = (df%>%distinct(taxon))$taxon %in% (accepted_matches%>%filter(multi_correct == "keeper"))$taxon)%>%
    filter(keeper_found == F)
  
  #Shows species for which multiple legitimate names were found, `dups` is the number of multiples found
  multi_legit_match<-
    accepted_matches%>%
    filter(multi_correct == "keeper")%>%
    group_by(taxon)%>%
    summarise(dups = length(taxon))%>%
    filter(dups>1)
  
  accepted_matches2<-
    accepted_matches%>%
    mutate(multi_correct = ifelse(taxon%in%no_legit_match$taxon | taxon%in%multi_legit_match$taxon,"manual_check_needed",multi_correct))
  
  return(accepted_matches2)
  
}

## Define custom min and max functions #########################################

custom_min <- function(x) {
  if (all(is.na(x))) {
    return(NA)
  } else {
    return(min(x, na.rm = TRUE))
  }
}

custom_max <- function(x) {
  if (all(is.na(x))) {
    return(NA)
  } else {
    return(max(x, na.rm = TRUE))
  }
}

## Clean gbif data #############################################################

str(gbif)


# WARNING: this function can take a while to run
# Here's a question... should the gbif cleaning happen before or after the matching to the WCVP?


gbif_clean<-
  gbif%>%
  cc_aohi()%>%
  cc_cen()%>%
  cc_dupl()%>%
  cc_equ()%>%
  cc_inst()%>%
  cc_val()%>%
  cc_zero()

# Visualising all the removed records #####

# gbif%>%
#   filter(! gbifID %in% gbif_clean$gbifID)%>%
#   ggplot(aes(decimalLongitude,decimalLatitude))+
#   geom_point()+
# gbif_clean%>%
#   ggplot(aes(decimalLongitude,decimalLatitude))+
#   geom_point()

# Cleaned the coordinates, refer to citation for details
# citation("CoordinateCleaner")



## Create translation table ######################################################

gbif_list<-
  gbif_clean%>%
  mutate(taxon = ifelse(taxonRank == "GENUS",genus,species))%>%
  group_by(taxon)%>%
  summarise(num_recs = length(species))#%>%view()


# Note: this block takes about a minute to run
matched_list<-
  #full_join(gbif_list,GloNAF_list)%>%
  gbif_list%>%
  wcvp_match_names(.,
                   name_col="taxon",
                   author_col=NULL,
                   fuzzy=TRUE,
                   progress_bar=TRUE)%>%
  multi_correct()%>%
  mutate(plant_name_id = wcvp_accepted_id)%>%
  left_join(wcvp_names%>%
             select(genus,species,taxon_name, taxon_authors,plant_name_id))%>%
  select(-plant_name_id)%>%
  rename(wcvp_accepted_name = taxon_name)


matched_list_checked<-
  matched_list%>%
  mutate(multi_correct = 
           case_when(taxon == "Cotoneaster dammeri x integrifolius" & wcvp_accepted_name == "Cotoneaster integrifolius" ~ "remove_checked",
                     multi_correct == "manual_check_needed" & wcvp_accepted_name == "Cotoneaster leveillei" ~ "remove_checked",
                     multi_correct == "manual_check_needed" & wcvp_accepted_name == "Cotoneaster nan-shan" ~ "keeper_checked",
                     multi_correct == "manual_check_needed" & wcvp_accepted_name == "Cotoneaster obscurus" ~ "keeper_checked",
                     multi_correct == "manual_check_needed" & wcvp_accepted_name == "Cotoneaster hebephyllus var. hebephyllus" ~ "remove_checked",
                     multi_correct == "manual_check_needed" & wcvp_accepted_name == "Cotoneaster soonogoricus" ~ "keeper_checked",
                     multi_correct == "manual_check_needed" & wcvp_accepted_name == "Cotoneaster affinis" ~ "remove_checked",
                     multi_correct == "manual_check_needed" & wcvp_accepted_name == "Cotoneaster horizontalis" ~ "keeper_checked",
                     taxon == "Cotoneaster parnassicus" & wcvp_accepted_name == "Eriobotrya elliptica" ~ "remove_checked",
                     taxon == "Cotoneaster parnassicus" & wcvp_accepted_name == "Cotoneaster lindleyi" ~ "keeper_checked",
                     TRUE ~ multi_correct))#%>%view()

translation_table<-
  matched_list_checked%>%
  filter(multi_correct == "keeper" | multi_correct == "keeper_checked")%>%
  select(taxon,num_recs,wcvp_accepted_id,wcvp_accepted_name,genus,species,taxon_authors)%>%
  rename(taxon_stand = taxon,
         num_recs_stand = num_recs,
         plant_name_id_stand = wcvp_accepted_id,
         genus_stand = genus,
         species_stand = species,
         authors_stand = taxon_authors)%>%
  filter(!is.na(wcvp_accepted_name))%>%
  distinct()#%>%view()


# any names not in the translation table will be removed from further analysis


# Checks that every entry in the translation table is unique
length(unique(translation_table$taxon_stand)) == nrow(translation_table)

## Standardize nomenclature ####################################################

gbif_stand<-
  gbif_clean%>%
  mutate(taxon_stand = ifelse(taxonRank == "GENUS",genus,species))%>%
  left_join(translation_table)#%>%view()


## Make gbif into shapefile ###################################################

gbif_stand_geo<-
  gbif_stand%>%
  st_as_sf(coords = c("decimalLongitude","decimalLatitude"),crs = 4326)%>%
  st_join(country_shapefiles)


# # This shows an example in Norway of some records not within the country boundaries, 
# # much of this is likely due to rounded coordinates 
# ggplot(gbif_stand_geo%>%filter(is.na(LEVEL3_NAM)))+
#   geom_sf(data = country_shapefiles)+
#   geom_sf()+
#   coord_sf(xlim = c(5,10), ylim = c(60, 58), expand = FALSE)


gbif_stand_geo%>%filter(is.na(LEVEL3_NAM))%>%
  nrow()

gbif_stand_geo_clean<-
  gbif_stand_geo%>%
  filter(!is.na(LEVEL3_COD))#%>%view()

## Final_products ##############################################################

gbif_final<-
  gbif_stand_geo_clean%>%
  mutate(tdwg3 = LEVEL3_COD)%>%
  #left_join(GloNAF_contoneaster%>%select(status,tdwg3,taxon_stand))%>%
  left_join(glonaf_taxon_x_region%>%
              group_by(status,tdwg3,taxa_accepted)%>%
              summarise(obs = length(status))%>%
              #select(status,tdwg3,taxa_accepted)%>%
              group_by(tdwg3,taxa_accepted)%>%
              filter(!(status == "Naturalized" & "Invasive" %in% status))%>%
              rename(wcvp_accepted_name = taxa_accepted)%>%
              ungroup())%>%
              #unique())#%>%
  left_join(wcvp_native%>%
              mutate(plant_name_id_stand = accepted_plant_name_id,
                     tdwg3 = area_code_l3)%>%
              select(plant_name_id_stand,tdwg3,introduced)%>%
              filter(!is.na(tdwg3))%>%
              unique())%>%
  mutate(glonaf_status = status,
         wcvp_status = case_when(introduced == 1 ~ "introduced",
                                 is.na(introduced) ~ "likely_introduced",
                                 introduced == 0 & taxonRank != "GENUS" ~ "native",
                                 
                                 #taxonRank == "GENUS" ~ NA,
                                 TRUE ~ NA))%>%
  select(-status,-introduced,-taxon_stand,-obs)%>%
  cbind(st_coordinates(.))%>%
  st_drop_geometry()%>%
  rename(decimalLongitude = X,
         decimalLatitude = Y)%>%
  filter(!is.na(wcvp_accepted_name),
         startsWith(wcvp_accepted_name,"Cotoneaster"))%>%
  select(gbifID,
         taxonRank,
         coordinateUncertaintyInMeters,
         eventDate,
         day,
         month,
         year,
         speciesKey,
         basisOfRecord,
         license,
         recordedBy,
         plant_name_id_stand,
         wcvp_accepted_name,
         genus_stand,
         species_stand,
         authors_stand,
         tdwg3,
         glonaf_status,
         wcvp_status,
         decimalLongitude,
         decimalLatitude)%>%


## Allocating nativeness #######################################################


# gbif_final<-
#   readr::read_csv2("Cotoneaster review output/Final gbif records_050225.csv")

  gbif_final%>%
  mutate(nativeness = case_when((glonaf_status == "Naturalized"|
                                  glonaf_status == "Invasive") &
                                  (wcvp_status == "introduced"|
                                  wcvp_status == "likely_introduced")~ "non-native",
                                is.na(glonaf_status) & wcvp_status == "native" ~ "native",
                                is.na(glonaf_status) & wcvp_status == "introduced"~ "non-native"))%>%
  mutate(nn_nativeness = "0",
         nn_id = "0",
         nn_dist = "0")%>%
  mutate(x = sphere(lat = decimalLatitude,lon = decimalLongitude)[[1]],
         y = sphere(lat = decimalLatitude,lon = decimalLongitude)[[2]],
         z = sphere(lat = decimalLatitude,lon = decimalLongitude)[[3]])  



gbif_final%>%
  rgl::plot3d(x = .$x,y = .$y,z = .$z,colour = .$species_stand)

# htmlwidgets::saveWidget(rglwidget(width = 520, height = 520), 
#                         file = here("3dscatter.html"),
#                         libdir = "libs",
#                         selfcontained = FALSE
# )


library(RANN)
# This algorithm uses euclidean distance which may be inaccurate but I'll leave it as is for now...


for (i in seq_len(nrow(gbif_final))) {
  
  if (is.na(gbif_final$nativeness[i])) {
    
    subset_df <- gbif_final %>%
      filter(!is.na(nativeness),
             species_stand == gbif_final$species_stand[i])
    
    if (nrow(subset_df) == 0) next
    
    coords <- subset_df %>%
      select(decimalLongitude, decimalLatitude) %>%
      as.matrix()
    
    
    query <- matrix(c(gbif_final$decimalLongitude[i],
                      gbif_final$decimalLatitude[i]),
                    nrow = 1)
    
    nn <- nn2(data = coords, query = query, k = 1)
    
    nearest_point <- subset_df[nn$nn.idx, ]
    
    gbif_final$nn_nativeness[i] <- nearest_point$nativeness
    gbif_final$nn_id[i] <- nn$nn.idx
    gbif_final$nn_dist[i] <- nn$nn.dists
  }
}

# gbif_final%>%
#   st_as_sf(coords = c("decimalLongitude","decimalLatitude"),crs = 4326)


## Plotting records of different species ####

species_list<-gbif_final%>%
  group_by(species_stand)%>%
  summarise(n = length(species_stand))

results<-list()

for(i in 1:length(species_list$species_stand)){
  
  plot<-
    ggplot()+
    geom_sf(data = country_shapefiles)+
    geom_sf(data = gbif_final%>%
              mutate(nativeness = ifelse(!is.na(nativeness),nativeness,nn_nativeness))%>%
              mutate(nativeness = ifelse(nativeness == 0,"unknown",nativeness))%>%
              filter(nn_dist<10,
                     species_stand == species_list$species_stand[i])%>%
              st_as_sf(coords = c("decimalLongitude","decimalLatitude"),crs = 4326),
            aes(colour = nativeness))+
    labs(title = paste("Cotoneaster",species_list$species_stand[i]))
  
  results[[i]]<-plot
  
}

ggplot()+
  geom_sf(data = country_shapefiles)+
  geom_sf(data = gbif_final%>%
            mutate(nativeness = ifelse(!is.na(nativeness),nativeness,nn_nativeness))%>%
            mutate(nativeness = ifelse(nativeness == 0 | nn_dist>20,"unknown",nativeness))%>%
            st_as_sf(coords = c("decimalLongitude","decimalLatitude"),crs = 4326),
          aes(colour = nativeness))

m <- runif(6,0,139)%>%round()


plots_random<-list()

for (j in 1:6) {

  plots_random[[j]]<-results[[m[[j]]]]

}

plots_random[[2]]

plots_random[1:6]



## Final spreadsheet ###########################################################


final_spreadsheet<-
  gbif_final%>%
  group_by(wcvp_accepted_name,tdwg3,glonaf_status,wcvp_status)%>%
  summarise(earliestYear = custom_min(year),
            latestYear = custom_max(year),
            obs = length(wcvp_accepted_name))%>%
  left_join(country_shapefiles%>%select(LEVEL3_NAM,LEVEL3_COD)%>%
              rename(tdwg3 = LEVEL3_COD,
                     tdwg3_name = LEVEL3_NAM)%>%
              st_drop_geometry())%>%
  st_drop_geometry()



## Exporting ###################################################################

#write.csv2(st_drop_geometry(final_spreadsheet),"Cotoneaster review output/Final summary spreadsheet_050225.csv")
#write.csv2(st_drop_geometry(gbif_final),"Cotoneaster review output/Final gbif records_050225.csv")
# write.csv2(translation_table,"translation_table.csv")



## Check the distributions for a few species ########

# r<-wcvp_distribution(taxon = "Cotoneaster washanensis",
#                      introduced = F)
# wcvp_distribution_map(r)
# view(r)
# 
# 
# gbif_final%>%
#   filter(wcvp_accepted_name == "Cotoneaster horizontalis")%>%
#   st_as_sf(coords = c("decimalLongitude","decimalLatitude"),crs = 4326)%>%
#   ggplot()+
#   geom_sf(data = country_shapefiles,fill = "white")+
#   geom_sf(aes(colour = glonaf_status,shape = wcvp_status))
# 
# 
# 
# get_wgsrpd3_codes("Russia")

