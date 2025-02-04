rm(list = ls())
#setwd("~/Universiteit werk/Masters/Cotoneaster review")
setwd("~/Cotoneaster_review")
library(tidyverse)
library(readr)
library(sf)
library(rWCVP)
#install.packages("rWCVPdata", repos=c("https://matildabrown.github.io/drat", "https://cloud.r-project.org"))
library(rWCVPdata)
library(CoordinateCleaner)
library(conflicted)
library(tictoc)
library(patchwork)
library(readxl)

conflict_scout()
conflicts_prefer(dplyr::filter,
                 dplyr::lag)


## Load WCVP data ##############################################################

wcvp_distributions<-rWCVPdata::wcvp_distributions

wcvp_names<-rWCVPdata::wcvp_names

country_shapefiles<-wgsrpd3%>%
  st_set_crs(4326)

ggplot(country_shapefiles)+
  geom_sf()


## Load data from gbif and GloNAF #########################################################



#gbif data
gbif <- read_delim("0006576-241126133413365.csv", delim = "\t", escape_double = FALSE,trim_ws = TRUE)

#GloNAF taxon by region data

GloNAF_taxon_x_region_orig<-
  read.csv2("257_9_257_2_GLONAF/GLONAF/Taxon_x_List_GloNAF_vanKleunenetal2018Ecology.csv")



GloNAF_taxon_x_region<-
  read.csv2("GloNAF 2024/glonaf_flora2.csv")

Glonaf_region<-
  read.csv2("GloNAF 2024/glonaf_region.csv")

glonaf_TxR<-read.csv2("GloNAF 2024/glonaf_TxR.csv")

glonaf_shapefiles<-
  read_sf("GloNAF 2024/glonaf_2024_nol/glonaf_2024_nol/glonaf_2024_nol.shp")










# WCVP_native distributions 
WCVP_native<-
  wcvp_names%>%
  filter(genus == "Cotoneaster")%>%
  left_join(wcvp_distributions)#%>%view()

#GloNAF 2024 download ######




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

# Visualiszing all the removed records #####

gbif%>%
  filter(! gbifID %in% gbif_clean$gbifID)%>%
  ggplot(aes(decimalLongitude,decimalLatitude))+
  geom_point()+
gbif_clean%>%
  ggplot(aes(decimalLongitude,decimalLatitude))+
  geom_point()

# Cleaned the coordinates, refer to citation for details
citation("CoordinateCleaner")



## Create translation table ######################################################

gbif_list<-
  gbif_clean%>%
  #filter(taxonRank == "SPECIES")%>%
  mutate(taxon = ifelse(taxonRank == "GENUS",genus,species))%>%
  group_by(taxon)%>%
  summarise(num_recs = length(species))#%>%view()

# GloNAF_list<-
# GloNAF_taxon_x_region%>%
#   filter(startsWith(standardized_name, "Cotoneaster"))%>%
#   left_join(GloNAF_regions)%>%
#   group_by(standardized_name)%>%
#   summarise()%>%
#   rename(taxon = standardized_name)#%>%view()


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
  #rename(plant_name_id = wcvp_accepted_id)%>%
  #left_join(wcvp_names%>%
  #            select(genus,species, taxon_authors,plant_name_id))%>%
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

## Standardized glonaf #########################################################
# GloNAF_contoneaster<-
# GloNAF_taxon_x_region_orig%>%
#   filter(startsWith(standardized_name, "Cotoneaster"))%>%
#   mutate(taxon_stand = standardized_name)%>%
#   left_join(translation_table)%>%
#   #left_join(GloNAF_regions)%>%
#   left_join(Glonaf_region)%>%
#   group_by(genus_stand,species_stand,tdwg3,tdwg2,tdwg1,country,status)%>%
#   summarise(obs = length(species_stand))%>%
#   mutate(taxon_stand = paste(genus_stand,species_stand))%>%
#   ungroup()%>%
#   select(-country)%>%
#   unique()

## Make gbif into shapefile ###################################################

gbif_stand_geo<-
gbif_stand%>%
  st_as_sf(coords = c("decimalLongitude","decimalLatitude"),crs = 4326)%>%
  st_join(country_shapefiles)#%>%view()

# This shows an example in Norway of some records not within the country boundaries, 
# much of this is likely due to rounded coordinates 
ggplot(gbif_stand_geo%>%filter(is.na(LEVEL3_NAM)))+
  geom_sf(data = country_shapefiles)+
  geom_sf()+
  coord_sf(xlim = c(5,10), ylim = c(60, 58), expand = FALSE)




gbif_stand_geo%>%filter(is.na(LEVEL3_NAM))%>%
  nrow()

gbif_stand_geo_clean<-
  gbif_stand_geo%>%
  filter(!is.na(LEVEL3_COD))#%>%view()

## Final_products ##############################################################

gbif_final<-
  gbif_stand_geo_clean%>%
  mutate(tdwg3 = LEVEL3_COD)%>%
  left_join(GloNAF_contoneaster%>%select(status,tdwg3,taxon_stand))%>%
  left_join(WCVP_native%>%mutate(plant_name_id_stand = accepted_plant_name_id,
                                 tdwg3 = area_code_l3)%>%
              select(plant_name_id_stand,tdwg3,introduced)%>%
              unique())%>%
  mutate(glonaf_status = status,
         wcvp_status = case_when(introduced == 1 ~ "introduced",
                                 is.na(introduced) ~ "likely_introduced",
                                 introduced == 0 & taxonRank != "GENUS" ~ "native",
                                 
                                 #taxonRank == "GENUS" ~ NA,
                                 TRUE ~ NA))%>%
  select(-status,-introduced,-taxon_stand)%>%
  cbind(st_coordinates(.))%>%
  st_drop_geometry()%>%
  rename(decimalLongitude = X,
         decimalLatitude = Y)%>%
  filter(!is.na(wcvp_accepted_name),
         startsWith(wcvp_accepted_name,"Cotoneaster"))

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

# write.csv2(st_drop_geometry(final_spreadsheet),"Final summary spreadsheet.csv")
# write.csv2(st_drop_geometry(gbif_final),"Final gbif records.csv")
# write.csv2(translation_table,"translation_table.csv")



## Check the distributions for a few species ########

r<-wcvp_distribution(taxon = "Cotoneaster integerrimus",
                     introduced = F)
wcvp_distribution_map(r)
view(r)


## Playing around with GeoRange ################################################
# 
# a<-
# gbif_clean%>%
#   filter(species == "Cotoneaster franchetii")%>%
#   mutate(decimalLongitude = as.numeric(decimalLongitude),
#          decimalLatitude = as.numeric(decimalLatitude))%>%
#   st_as_sf(coords = c("decimalLongitude","decimalLatitude"),crs = 4326,remove = F)%>%
#   sf_make_grid(cellsize = 10,what = "polygons")
#   
#     #view()
#   
#   
#   grid <- st_make_grid(sf_data, cellsize = 0.01, what = "centers")  # Adjust cellsize
# 
# 
#       
# a%>%
#   group_by(countryCode)%>%
#   st_convex_hull()%>%
#   plot()
# 
# 
# a%>%
#   mutate(cluster = kmeans(select(decimalLongitude,decimalLatitude),centers = 6)$cluster)%>%
#   as.factor%>%view()
# 
# km_result<-kmeans(x = select(a,decimalLongitude,decimalLatitude),centers = 6)$cluster
# 
# 
# a%>%
#   mutate(cluster = kmeans(x = select(a,decimalLongitude,decimalLatitude),centers = 10)$cluster)%>%
#   st_as_sf(coords = c("decimalLongitude","decimalLatitude"),crs = 4326,remove = F)%>%
#   group_by(cluster)%>%
#   summarise(geometry = st_convex_hull(st_union(geometry)))%>%
#   ggplot()+
#   geom_sf()
#     
# PlotConvexHull(xcoord = as.numeric(a$decimalLongitude), ycoord = as.numeric(a$decimalLatitude))
#     plot(    )
#     
#     CellCount(longs = as.numeric(a$decimalLongitude), lats = as.numeric(a$decimalLatitude), CellSize = 10,
#               latBounds = c(-90, 90))%>%
#       view()->b
#     
#       ggplot(b)+
#       geom_point(aes(Coordinates.Longitude,Coordinates.Latitude))

## Code graveyard ########

# gbif_matched_2<-
#   gbif_matched%>%
#   mutate(species = taxon)%>%
#   multi_correct()%>%
#   filter(multi_correct == "keeper")%>%
#   select(taxon,species,wcvp_name,wcvp_authors,wcvp_accepted_id,obs)#%>%view()

# joined_gbif<-
# left_join(gbif,gbif_matched_2)%>%view()

# C_list_wcvp<-
#   joined_gbif%>%
#   group_by(species,countryCode)%>%
#   summarise(num_recs = length(scientificName))%>%view()

# species_list_matched_check%>%
#   select(taxon,num_recs,wcvp_accepted_id)%>%
#   view()
# 
# GloNAF_WCVP_match_check%>%
#   select(taxon,wcvp_accepted_id)%>%
#   view()
# new_gbif<-
#   gbif_clean%>%
#   mutate(taxon = species)%>%
#   mutate(taxon = ifelse(taxon %in% gbif_list_matched$taxon,taxon,"Cotoneaster"))%>%
#   left_join(species_list_matched_check)%>%
#   mutate(taxon = ifelse(taxon == "","Cotoneaster",taxon))%>%view()

# First match names from GBIF
# gbif_list<-
#   gbif_clean%>%
#   filter(taxonRank == "SPECIES")%>%
#   group_by(species)%>%
#   summarise(num_recs = length(scientificName))#%>%view()
# 
# gbif_list_matched<-
# gbif_list%>%
#   rename(taxon = species)%>%
#   wcvp_match_names(.,
#                    name_col="taxon",
#                    author_col=NULL,
#                    fuzzy=TRUE,
#                    progress_bar=TRUE)
# 
# gbif_list_matched_check<-
#   gbif_list_matched%>%
#   multi_correct()%>%
#   filter(multi_correct == "keeper")%>%
#   filter(taxon != 'Cotoneaster dammeri x integrifolius',
#          taxon != 'Cotoneaster amphigenus')

## Notes
# Removed 'dammeri x integrifolius' (hybrid not known to WCVP) and 'amphigenus' (unplaced name)
# Found in gbif only 97 of the 268 accepted species (accoring to plants of the world online)
# Including herbarium records is likely to increase this number though

# Second, match names from GloNAF

# GloNAF_list<-
# GloNAF_taxon_x_region%>%
#   filter(startsWith(standardized_name, "Cotoneaster"))%>%
#   left_join(GloNAF_regions)%>%
#   group_by(standardized_name)%>%
#   summarise()%>%
#   rename(taxon = standardized_name)#%>%view()
# 
# 
# GloNAF_WCVP_match<-
#   GloNAF_list%>%
#   group_by(taxon)%>%
#   summarise(obs = length(taxon))%>%
#   wcvp_match_names(.,
#                    name_col="taxon",
#                    author_col=NULL,
#                    fuzzy=TRUE,
#                    progress_bar=TRUE)
#   
# GloNAF_WCVP_match_check<-
#   GloNAF_WCVP_match%>%
#   multi_correct()%>%
#   filter(multi_correct != "remove")#%>%view()

# There is one species that needed manual checking but they are practically identical, 
# so it doen't make a difference for the rest of the workflow as it is the same accepted species
# maybe it's worth removing one of them 

# Species list translation table 

# translation_table<-
# full_join(gbif_list_matched_check%>%
#             select(taxon,num_recs,wcvp_accepted_id),
#           GloNAF_WCVP_match_check%>%
#             select(taxon,wcvp_accepted_id))%>%
#   rename(plant_name_id = wcvp_accepted_id)%>%view()
#   left_join(wcvp_names)
# 
# unique_species<-
#   translation_table%>%
#   group_by(taxon_name,plant_name_id)%>%
#   summarise(num_synonyms = length(plant_name_id))%>%view()

# translation_table<-
#   full_join(gbif_list,GloNAF_list)%>%
#   wcvp_match_names(.,
#                    name_col="taxon",
#                    author_col=NULL,
#                    fuzzy=TRUE,
#                    progress_bar=TRUE)
# 
# translation_table2<-
# translation_table%>%
#   multi_correct()%>%view()
# 
# translation_table3<-
#   translation_table2%>%
#   mutate(multi_correct = ifelse(multi_correct == "manual_check_needed" &
#                                   wcvp_id == "2944452","keeper_checked",multi_correct))%>%
#   filter(multi_correct == "keeper" | multi_correct == "keeper_checked")%>%
#   filter(!is.na(wcvp_accepted_id))
# 
# length(unique(translation_table3$taxon)) == nrow(translation_table3)
# 
# translation_table_proper<-
#   translation_table3%>%
#   select(taxon,num_recs,wcvp_accepted_id)%>%
#   rename(plant_name_id = wcvp_accepted_id)%>%
#   left_join(wcvp_names%>%
#               select(genus,species, taxon_authors,plant_name_id))%>%view()
# gbif_clean2%>%
#   rename(taxon = species)%>%
#   left_join(translation_table%>%
#               select(taxon,plant_name_id,genus,species,taxon_authors),
#             by = join_by(taxon))%>%view()
## Create country and species list ##

# gbif_geo<-
# gbif_clean%>%
#   filter(!is.na(decimalLongitude) | decimalLongitude != "",
#          !is.na(decimalLongitude) | decimalLatitude != "")%>%
#   st_as_sf(coords = c("decimalLongitude","decimalLatitude"),crs = 4326)
# 
# ggplot()+
#   #geom_sf()+
#   geom_sf(data = country_shapefiles)
# 
# ggplot()+
#   geom_sf(data = country_shapefiles)+
#   geom_sf(data = gbif_geo)
# gbif_country<-
#   gbif_geo%>%
#   st_join(country_shapefiles)
# 
# country_list<-
#   gbif_country%>%
#   filter(!is.na(LEVEL3_NAM))%>%
#   group_by(taxon,wcvp_name,wcvp_authors,wcvp_accepted_id,LEVEL3_NAM,LEVEL3_COD,
#            LEVEL2_COD,LEVEL1_COD)%>%
#   summarise(obs = length(taxon))#%>%view()

## Get native ranges for each species ##
# wcvp_coton_dist<-
#   wcvp_distributions%>%
#   filter(plant_name_id%in%species_list_matched_check$wcvp_accepted_id)%>%
#   filter(introduced != 1)%>%
#   select(plant_name_id,area_code_l3,area)%>%
#   mutate(wcvp_native = 1)
# 
# country_list_native<-
# country_list%>%
#   rename(plant_name_id = wcvp_accepted_id,
#          area = LEVEL3_NAM,
#          area_code_l3 = LEVEL3_COD)%>%
#   left_join(wcvp_coton_dist)%>%
#   mutate(wcvp_native = ifelse(is.na(wcvp_native),0,wcvp_native))%>%view()
#  
# Final_spreadsheet<-
# country_list_native%>%
#   left_join(GloNAF)%>%
#   mutate(status = ifelse(is.na(status),"unknown",status))%>%
#   rename(glonaf_status = status)%>%view()

## Country colours ##

# Final_spreadsheet%>%
#   st_drop_geometry()%>%
#   group_by(area_code_l3)%>%
#   summarise(obs = length(taxon))%>%
#   rename(LEVEL3_NAM = area_code_l3)%>%
#   left_join(country_shapefiles)%>%view()
#   ggplot()+
#   geom_sf(aes(colour = obs))
# country_shapefiles%>%
#   left_join(Final_spreadsheet%>%
#               st_drop_geometry()%>%
#               group_by(area_code_l3)%>%
#               summarise(n.spp = length(taxon))%>%
#               rename(LEVEL3_COD = area_code_l3))%>%
#   mutate(n.spp = ifelse(is.na(n.spp),0,n.spp))%>%
#   ggplot()+
#   geom_sf(aes(fill = n.spp))
# 
# country_shapefiles%>%
#   left_join(Final_spreadsheet%>%
#               st_drop_geometry()%>%
#               group_by(area_code_l3)%>%
#               summarise(abund = sum(obs))%>%
#               rename(LEVEL3_COD = area_code_l3))%>%
#   mutate(abund = ifelse(is.na(abund),0,abund))%>%
#   ggplot()+
#   geom_sf(aes(fill = abund))
# 
# country_shapefiles%>%
#   left_join(Final_spreadsheet%>%
#               st_drop_geometry()%>%
#               group_by(area_code_l3)%>%
#               summarise(abund = sum(obs))%>%
#               rename(LEVEL3_COD = area_code_l3))%>%
#   mutate(abund = ifelse(is.na(abund),0,abund))%>%
#   ggplot()+
#   geom_sf(aes(fill = log(abund)))
# gbif_stand_geo%>%
#   mutate(tdwg3 = LEVEL3_COD)%>%
#   left_join(GloNAF_contoneaster%>%
#               select(tdwg3,status))%>%
#   left_join(wcvp_distributions%>%
#               select())
# gbif<-read.csv2("0013728-241107131044228.csv")#%>%
#   #sample_n(10000)
# gbif_new<-
#   read.csv2("0006576-241126133413365.csv",sep = "\t")
# institutions<-CoordinateCleaner::institutions%>%view()
# 
# ggplot(data = institutions,aes(decimalLongitude,decimalLatitude))+
#   geom_point()
# filter(decimalLongitude != "" & !is.na(decimalLongitude),
  #        decimalLatitude != "" & !is.na(decimalLatitude))%>%
  # mutate(decimalLongitude = as.numeric(decimalLongitude),
  #        decimalLatitude = as.numeric(decimalLatitude))%>%
# gbif%>%
#   select(species,decimalLatitude,decimalLongitude)%>%
#   unique()%>%
#   ggplot(aes(decimalLongitude,decimalLatitude,colour = species))+
#   geom_point()+
#   theme(legend.position = "none")
# gbif_clean2<-
#   gbif_clean%>%
#   filter(.summary == TRUE)%>%
#   select(-.val,-.equ,-.zer,-.cap,-.cen,-.gbf,-.inst,-.summary)#%>%view()
# vec1 <- c(3, 7, NA, 2, 5) # Contains non-NA values
# vec2 <- c(NA, NA, NA)      # Contains only NA values
# 
# # Using min() with na.rm = TRUE
# 
# min(vec1)
# min(vec2)
# min(vec1,na.rm = T)
# min(vec2,na.rm = T)
# 
# custom_min(vec1)
# custom_min(vec2)
# 
# custom_min <- function(x) {
#   if (all(is.na(x))) {
#     return(NA) # Return NA if all values are NA
#   } else {
#     return(min(x, na.rm = TRUE)) # Ignore NA when other values exist
#   }
#}
# ?min
