
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

sphere <- function(lat, lon, R = 6371000) {
  # Convert to radians
  phi <- lat * pi / 180
  lambda <- lon * pi / 180
  
  # Calculate coordinates
  x <- R * cos(phi) * cos(lambda)
  y <- R * cos(phi) * sin(lambda)
  z <- R * sin(phi)
  
  return(list(x, y, z))
}

read_glonaf <- function(which) {
  files <- c(
    glonaf_flora2 = "data/GloNAF 2024/glonaf_flora2.csv",
    glonaf_taxon  = "data/GloNAF 2024/glonaf_taxon_wcvp.csv",
    glonaf_region = "data/GloNAF 2024/glonaf_region.csv"
  )
  
  if (!which %in% names(files)) {
    stop("Invalid input")
  }
  
  read.csv2(here::here(files[[which]]))
}

nn_dist_to_arc<-
  function (x, r = 6371000){
    d <- acos( 1 - ( (x^2) / ( 2*(r^2)) ) ) *r
    round(d)
  }

