#Overview

# This R script was made for assessing the Cotoneaster data

# Note: does need an internet connection, some of the packages like CoordinateCleaner and 
# rWCVP require connection to the internet though it does not consume a lot of data


############################ SETUP #############################################

rm(list = ls())

library(tidyverse)
library(readr)
library(sf)
library(rWCVP)
#install.packages("rWCVPdata", repos=c("https://matildabrown.github.io/drat", "https://cloud.r-project.org"))
library(rWCVPdata)
library(CoordinateCleaner)
library(here)
library(RANN)
library(conflicted)

source(here("rscripts/Cotoneaster review functions.R"))

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


glonaf_taxon_x_region<-
  
  left_join(read_glonaf("glonaf_flora2"),
            read_glonaf("glonaf_taxon")%>%
              select(id,accepted_plant_name_id,taxa_accepted)%>%
              rename(taxon_wcvp_id = id))%>%
  filter(startsWith(taxa_accepted,"Cotoneaster"))%>%
  
  left_join(read_glonaf("glonaf_region")%>%
              select(id,code,name,tdwg4_id)%>%
              rename(region_id = id))%>%
  
  mutate(tdwg3 = str_sub(tdwg4_id,1,3),
         IDregion = code)%>%
  
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
  left_join(wcvp_distributions)



###################### CLEAN GBIF DATA #########################################


gbif_clean<-
  gbif_raw%>%
  mutate(taxon = ifelse(taxonRank == "GENUS",genus,species))%>%
  select(gbifID,
         taxon,
         taxonRank,
         decimalLatitude,decimalLongitude,
         day,month,year,
         establishmentMeans,
         countryCode
  )%>%
  clean_coordinates(
    countries = "countryCode",
    species = "taxon",
    tests = c("capitals", "centroids", "equal", "gbif", "institutions", "seas", "zeros")
  )%>%
  filter(.summary)%>%
  select(-(.val:.summary))


######################### CREATE TRANSLATION TABLE #############################


translation_table<-
  
  gbif_clean%>%
  group_by(taxon)%>%
  summarise(num_recs = length(taxon))%>%
  
  wcvp_match_names(.,
                   name_col="taxon",
                   author_col=NULL,
                   fuzzy=TRUE,
                   progress_bar=FALSE)%>%
  
  multi_correct()%>%
  
  filter(match_similarity > 0.9)%>%
  
  mutate(wcvp_accepted_id = ifelse(!is.na(wcvp_accepted_id),wcvp_accepted_id,wcvp_id))%>%
  
  left_join(wcvp_names%>%
              select(genus,species,taxon_name, taxon_authors,plant_name_id)%>%
              rename(wcvp_accepted_id = plant_name_id))%>%
  
  rename(wcvp_accepted_name = taxon_name)%>%
  # 
  group_by(taxon)%>%
  
  filter(multi_correct == "keeper" | multi_correct == "keeper_checked")%>%
  #select(taxon,num_recs,wcvp_accepted_id,wcvp_accepted_name,genus,species,taxon_authors,wcvp_status)%>%
  
  rename(wcvp_species = species,
         wcvp_genus = genus)%>%
  filter(!is.na(wcvp_accepted_name))%>%
  ungroup()%>%
  distinct()

# any names not in the translation table will be removed from further analysis


# Checks that every entry in the translation table is unique
length(unique(translation_table$taxon)) == nrow(translation_table)


######################## CREATE SHAPEFILES ####################################



wcvp_glonaf_sf<-
  right_join(country_shapefiles%>%mutate(tdwg3 = LEVEL3_COD),
             
             full_join(glonaf_taxon_x_region%>%
                         select(status,accepted_plant_name_id,tdwg3)%>%
                         unique(),
                       
                       wcvp_native%>%
                         filter(taxon_name %in% translation_table$wcvp_accepted_name)%>%
                         select(accepted_plant_name_id,
                                area_code_l3,
                                introduced)%>%
                         
                         unique(),
                       by = join_by(tdwg3 == area_code_l3,accepted_plant_name_id))%>%
               
               
               mutate(glonaf_status = status,
                      wcvp_status = case_when(introduced == 1 ~ "introduced",
                                              introduced == 0 ~ "native",
                                              TRUE ~ NA),
                      nativeness = case_when(wcvp_status == "native" & !is.na(glonaf_status) ~ "disputed",
                                             wcvp_status == "native" & is.na(glonaf_status) ~ "native",
                                             wcvp_status == "introduced" | !is.na(glonaf_status) ~ "nonnative"
                      ))%>%
               
               select(accepted_plant_name_id,tdwg3,glonaf_status,wcvp_status,nativeness)%>%
               
               left_join(translation_table%>%
                           select(wcvp_accepted_id,wcvp_accepted_name)%>%
                           unique(),
                         by = join_by(accepted_plant_name_id==wcvp_accepted_id))
  )


#plotting some example data

wcvp_glonaf_sf%>%
  filter(wcvp_accepted_name == "Cotoneaster pannosus")%>%
  ggplot()+
  geom_sf(data = rWCVPdata::coast)+
  geom_sf(aes(fill = nativeness))+
  theme_minimal()+
  theme(legend.position = "inside",
        legend.position.inside = c(0.1,0.4))



######################## STANDARDIZE NOMENCLATURE ##############################


gbif_geo<-
  gbif_clean%>%
  #mutate(taxon_stand = ifelse(taxonRank == "GENUS",genus,species))%>%
  left_join(translation_table)%>%
  
  
  ######################### MAKE GBIF INTO SHAPEFILE ##############################


st_as_sf(coords = c("decimalLongitude","decimalLatitude"),crs = 4326)%>%
  st_join(country_shapefiles) %>%
  filter(!is.na(LEVEL3_COD))


############################ FINAL PRODUCTS ####################################


gbif_final<-
  gbif_geo%>%
  
  mutate(tdwg3 = LEVEL3_COD,
         accepted_plant_name_id = wcvp_accepted_id)%>%
  
  left_join(st_drop_geometry(wcvp_glonaf_sf)%>%
              select(accepted_plant_name_id,tdwg3,#glonaf_status, wcvp_status,
                     nativeness
              )%>%
              unique(),
            by = c("accepted_plant_name_id","tdwg3"))%>%
  
  
  group_by(wcvp_accepted_name) %>%
  group_modify(~{
    
    group_name <- .y$wcvp_accepted_name
    
    known <- .x %>% filter(!is.na(nativeness))
    unknown <- .x %>% filter(is.na(nativeness))
    
    if (nrow(unknown) == 0) return(.x)
    
    candidates <- wcvp_glonaf_sf %>%
      filter(wcvp_accepted_name == group_name)
    
    if (nrow(candidates) == 0) return(.x)
    
    idx <- st_nearest_feature(unknown, candidates)
    dists <- st_distance(unknown, candidates[idx, ], by_element = TRUE)
    
    unknown$nn_nativeness <- candidates$nativeness[idx]
    unknown$nn_id <- idx
    unknown$nn_dist <- dists
    
    bind_rows(known, unknown)
    
  }) %>%
  ungroup()%>%
  
  st_as_sf()%>%
  
  mutate(
    longitude = st_coordinates(.)[,1],
    latitude  = st_coordinates(.)[,2]
  ) %>%
  
  st_drop_geometry()%>%
  units::drop_units()


############################# EXPORTING ########################################


write_csv(gbif_final,here("data","Cotoneaster review output","Final gbif_records_090426.csv"))


