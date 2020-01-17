library(tidyverse)

convDecDate <- function(input){
  # This function converts the ForestPlots decimal date to a POSIXct
  year  = floor(input)
  doy   = round(365*(input-year))
  date  = lubridate::parse_date_time(paste(year,doy),"yj")
  return(date)
}

get_time_diffs <- function(date_vec){
  # This function calculates difference between dates, returning the number of days
  date_diff_vec <- numeric(length(date_vec)-1)
  for(i in 2:length(date_vec)){
    date_diff_vec[i-1] <- lubridate::int_length(lubridate::interval(date_vec[i-1], date_vec[i]))/86400
  }
  return(date_diff_vec);
}

#--- START GET WOOD DENSITY ----------------------------------------------------------------------------------#
find_wd <- function(census_sp_list){
  library(tidyverse);
  if(exists("zanne")==F) stop("load zanne wood density database")
  zanne_t <- zanne
  zanne_t$tropical <- F
  zanne_t$tropical[grepl("\\btropical\\b", zanne_t$region)] <- T
  zanne_t <- zanne_t %>% filter(tropical==T)
  zanne_t$region <- recode(zanne_t$region, 
                           `South America (tropical)`="South America", 
                           `Africa (tropical)`="Africa",
                           `South-East Asia (tropical)`="Asia")
  
  gwd_sps1 <- zanne_t %>% group_by(genus, species, region) %>% 
    summarize(wd_median=median(wd,na.rm=T),wd_mean=mean(wd,na.rm=T), nobs=n(),sd=sd(wd,na.rm=T)); 
  gwd_sps2 <- zanne_t %>% group_by(genus, species) %>% 
    summarize(wd_median=median(wd,na.rm=T),wd_mean=mean(wd,na.rm=T), nobs=n(),sd=sd(wd,na.rm=T)); 
  gwd_genus1 <- zanne_t %>% group_by(genus, region) %>% 
    summarize(wd_median=median(wd,na.rm=T),wd_mean=mean(wd,na.rm=T), nobs=n(),sd=sd(wd,na.rm=T)); 
  gwd_genus2 <- zanne_t %>% group_by(genus) %>% 
    summarize(wd_median=median(wd,na.rm=T),wd_mean=mean(wd,na.rm=T), nobs=n(),sd=sd(wd,na.rm=T)); 
  gwd_family1 <- zanne_t %>% group_by(family, region) %>% 
    summarize(wd_median=median(wd,na.rm=T),wd_mean=mean(wd,na.rm=T), nobs=n(),sd=sd(wd,na.rm=T)); 
  
  tmp <- census_sp_list
  tmp$region <- tmp$continent
  tmp$species <- sapply(strsplit(tmp$species," "), "[",2)
  # 1-) by Genus+Sps in SAME Region:
  tmp$wd_sps1 <- gwd_sps1[match(with(tmp, paste(genus, species, region)), with(gwd_sps1, paste(genus, species, region))),]$wd_median
  # 2-) by Genus+Sps in ANY Region:
  tmp$wd_sps2 <- gwd_sps2[match(with(tmp, paste(genus, species)), with(gwd_sps2, paste(genus, species))),]$wd_median
  # 3-) by Genus in SAME Region:
  tmp$wd_genus1 <- gwd_genus1[match(with(tmp, paste(genus, region)), with(gwd_genus1, paste(genus, region))),]$wd_median
  # 4-) by Genus in ANY Region:
  tmp$wd_genus2 <- gwd_genus2[match(with(tmp, paste(genus)), with(gwd_genus2, paste(genus))),]$wd_median
  # 5-) by Family in SAME Region:
  tmp$wd_family <- gwd_family1[match(with(tmp, paste(family, region)), with(gwd_family1, paste(family, region))),]$wd_median 
  
  # Assigning the wd in order of preference to the species with common names+scientific names:
  tmp$wd <- apply(tmp[,c("wd_sps1","wd_sps2","wd_genus1","wd_genus2","wd_family")], 1,function(x)na.omit(x)[1])
  print(paste("number of taxa without a match: ", sum(is.na(tmp$wd)))) #are all trees identified? 
  print("assigning species mean to unknowns"); 
  tmp$wd[is.na(tmp$wd)==T] <- mean(tmp$wd, na.rm=T)
  tmp$density <- tmp$wd
  tmp <- tmp %>% select(continent, family, genus, species,region,wd)
  return(tmp)
}
#--- END GET WOOD DENSITY ----------------------------------------------------------------------------------#


#--- FELDPAUCH TREE HEIGHT ESTIMATION SESSION NEEDS WORK, WILL BREAK EASILY --------------------------------#
## Feldpauch correction procedure for heights, diameters and densitys:

feldpauch_height <- function(dbh_cm){
  # Brazilian shield
  Bo    = 0.6373
  B1    = 0.4647  # E.C. Amazonia
  
  So1   = 0.012  # E.C. Amazonia
  Abar  = 20.4  # mean cenetered basal area m-2 ha-1
  
  n01   = 0.0034 # E.C. Amazonia
  Pvbar = 0.68 # mean centered precipitation coefficient of variation
  
  n02   = -0.0449 # E.C. Amazonia
  Sdbar = 5     # mean centered dry season length no months less than 100mm
  
  n03   = 0.0191  # E.C. Amazonia
  Tabar = 25.0  # mean centered annual averge temperature
  
  pred_height <- 10^(Bo + B1*log10(dbh_cm/10) + Abar*So1 + n01*Pvbar + n02*Sdbar + n03*Tabar); 
  return(pred_height)
}
#--- END FELDPAUCH -----------------------------------------------------------------------------------------------#


#--- START ALLOMETRIC EQUATIONS ----------------------------------------------------------------------------------#

## Chave et al, 2005 gives different allometric models to estimate aboveground biomass (AGB), three of which will be used in the R-Code:
# p. 91, Chave et al, 2005: model I.3 dry: 
Chave2005_dry <- function(diax, density, height, wd) {
  density=wd; 
  AGB_est <- 0.112*(density*((diax)^2)*height)^0.916
  return(AGB_est)
}

# p. 90, Chave et al, 2005: model I.6 moist:
Chave2005_moist <- function(diax, density, height, wd) {
  density=wd; 
  AGB_est <- 0.0509*density*((diax)^2)*height
  return(AGB_est)
}

# p. 95, Chave et al, 2005: model I.6 wet:
Chave2005_wet <- function(diax, density, height, wd) {
  density=wd; 
  AGB_est <- 0.0776*(density*((diax)^2)*height)^0.940
  return(AGB_est)
}

# Chave et al, 2014:
Chave2014 <- function(diax, density, height, wd) {
  density=wd; 
  AGB_est <- 0.0673*(density*((diax)^2)*height)^0.976 
  return(AGB_est)
}

# Chave et al, 2014:
Chave2014_eq4 <- function(diam_cm, density, height_m, wd) {
  density=wd; 
  AGB_est <- 0.0673*(density*((diam_cm)^2)*height_m)^0.976 
  return(AGB_est)
}


# Tree ferns Farfan et al., in prep:
# Farfan2015 <- function() {}

## Chambers et al, 2004: estimate surface area of a tree:
# see also RAINFOR manual, p. 52; diameter in cm 
# Note: A. Shenkin is developing our own surface area equation. Check if this is the most up to date equation that we use.
Chambers2004_surfaceArea <- function(diameter) {
  surface_area = 10^(-0.015-0.686*log10(diameter)+(2.208*log10(diameter)^2)-(0.627*log10(diameter)^3))
  return(surface_area)
}
#--- END ALLOMETRIC EQUATIONS ----------------------------------------------------------------------------------#

#--- MISC CENSUS CLEANING FUNCTIONS -----------------------------------------------------------------------#

get_dbh <- function(census, min_dend_date, cens_dend_dist=2){
  #--- This function:
  #--- 1) finds the max diameter of dead trees from previous census measurement dates 
  #--- 2) replaces the 0s in the original forestPlots census object with the prev max dbh
  #--- 3) filters the census file for being within 2 years of the first dendrometer measurment date
  
  tmp <- census 
  
  # ! There are a LOT of 0s, which I think are meant to be NAs! All these trees are dead.
  tmp <- tmp %>% mutate(d0=replace(d0, d0==0, NA), 
                        d1=replace(d1, d1==0, NA), 
                        d2=replace(d2, d2==0, NA), 
                        d3=replace(d3, d3==0, NA), 
                        d4=replace(d4, d4==0, NA)) 
  
  dead_tags <- tmp %>% filter(f1=="0") %>% pull(tree_tag) %>% unique()
  
  tmp_dead <- tmp %>% filter(tree_tag %in% dead_tags)
  
  tmp_dead <- tmp_dead %>%
    rowwise() %>% mutate(dbh_dead=max(c(d0,d1,d2,d3,d4),na.rm=T)) %>%
    ungroup() %>%  filter(dbh_dead<5000) %>% 
    group_by(date, tree_tag) %>% 
    summarize(dbh_dead=max(dbh_dead))
  
  # generate the dbh from the mean of the measurements
  tmp <- tmp %>%
    rowwise() %>% mutate(#dbh_mm_sd=sd(c(d0,d1,d2,d3,d4),na.rm=T), 
      dbh=max(c(d0,d1,d2,d3,d4),na.rm=T)) %>%
    ungroup() %>%  filter(dbh<5000)
  
  tmp <- left_join(tmp,tmp_dead,by=c("tree_tag","date")) %>% 
    rowwise() %>% 
    mutate(dbh=max(c(dbh, dbh_dead),na.rm=T)) %>% 
    select(-dbh_dead) %>% ungroup()
  
  
  # filter census for dates within $cens_dend_dist years of min dend date
  min_relevant_census_date <- min_dend_date-lubridate::years(cens_dend_dist)
  tmp <- tmp %>% dplyr::filter(date > min_relevant_census_date)
  
  return(tmp)
}

est_dbh_from_dendro <- function(dbh, dendrometer_reading_mm){
  if(is.na(dendrometer_reading_mm)==T){
    dbh_new = dbh;
  }else{
  circumf_old <- dbh*pi; 
  circumf_new <- circumf_old + dendrometer_reading_mm
  dbh_new = circumf_new/pi
  }
  return(dbh_new)
}
est_dbh_from_dendro <- Vectorize(est_dbh_from_dendro)

ppa <- function(census, plot_area_ha=1){
  tmp_df <- tibble(
    height = census$height_m,
    # crown_area = (0.036*(census$dbh**1.28)),
    crown_area = 0.06955046*(census$dbh**1.305),
    tree_tag = census$tree_tag, 
    canopy_layer = NA
  )
  n_canopyLayers <- ceiling(sum(tmp_df$crown_area)/(plot_area_ha*10000))
  tmp_df <- tmp_df %>% arrange(desc(height)) 
  
  for(i in 1:dim(tmp_df)[1]){
    if(i==1){
      idx_cl <- 1
      zstar <- 1}
    cap <- sum(tmp_df$crown_area[zstar:i])
    ifelse(cap<=(plot_area_ha*10000),
           {tmp_df$canopy_layer[i] <- idx_cl},
           {zstar <- i
           idx_cl <- idx_cl+1
           tmp_df$canopy_layer[i] <- idx_cl
           }
    )
  }
  census <- left_join(census, tmp_df %>% select(tree_tag,canopy_layer), by="tree_tag"); 
  return(census)
}

#--- END MISC CENSUS CLEANING FUNCTIONS -------------------------------------------------------------------#
