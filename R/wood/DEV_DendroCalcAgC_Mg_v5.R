## Scaling dendrometer data to monthly stem NPP
# original: Cecil Girardin [2014]
# methodological revisions: Sami Rifai [2017-2018]

# Assumptions:
# dendrometers measured every 3 months
# assumes all trees over 10 cm measured at 1.3 meters
# dbh is reported in mm

# Notes: 
# I've written in a number of data quality checks, but anyone (YOU!) running
# this script needs to look carefully at what it is doing. Don't assume that
# it is working correctly until you understand what it is doing and
# can verify that the growth increments between dendrometer measurements
# match up with the biomass incrementes between the census intervals.
# Note that this is only estimating a LIVE tree growth rate on a
# monthly basis, and does not account for tree mortality between
# census intervals.
# It's likely there are mistakes so use with caution!!! 
# STILL PROBLEMS WITH UPDATING BASELINE DBH AND CURRENT DBH !!!!!!!!!!!!!!!!!!!!!!!!!!!!

dendro_calc_agC_Mg_dev <- function(census, 
                                   mod_height, 
                                   dendrometer, 
                                   plot_code, 
                                   allometric_option="Default",
                                   ret="nppacw.perday.pertree", 
                                   census_year, 
                                   are_dend_units_correct=TRUE) { #census_year 1 
  
  library(tidyverse); library(mgcv); 
  dend <- dendrometer
  dend <- dend %>% distinct() # GET RID OF THE AWFUL DAMAGING DUPLICATES
  ## get data for all trees that are in the plot selected from census & dendrometer files
  # cen  <- subset(census, plot_code==plot_code) 
  dend$decDate <- decimal_date(dend$date)
  min_dend_date <- dend$decDate %>% min()
  start_census_date <- census[which.min(abs(census$census_date-min_dend_date)),]$census_date
  # cen <- cen %>% filter(near(census_date, start_census_date))
  
  # filter the census to plausible dates when dendro started to get the first DBH record from the census object
  census <- census %>% filter(census_date>= start_census_date) %>% arrange(census_date, tree_tag)
  
  # get list of tree_tags in dend object with no corresponding record in census
  vec_dtree_noCensus <- unique(dend$tree_tag)[!(unique(dend$tree_tag) %in% unique(census$tree_tag))]
  
  # filter out tree_tags with no census record
  dend <- dend %>% filter(!tree_tag %in% vec_dtree_noCensus) 
  
  # Get the date of the first dendrometer reading for each tree in the dend object
  baseline_dend <- dend %>% group_by(tree_tag) %>% filter(date==min(date)) %>%
    arrange(date) %>% select(tree_tag, date) %>%
    mutate(decDate=decimal_date(date))
  
  # Get the corresponding census date to start the dendro measures
  vec_census_dates <- census$census_date %>% unique() %>% sort()
  baseline_dend$census_date <- NA
  for(i in 1:dim(baseline_dend)[1]){
    baseline_dend$census_date[i] <- (vec_census_dates[which.min(abs(baseline_dend$decDate[i] - vec_census_dates))])
  }
  
  baseline_j <- left_join(baseline_dend, 
                          census %>% select(tree_tag, census_date, dbh, height_m,wd, canopy_layer),
                          by=c("tree_tag","census_date"))
  
  jnk_A <- baseline_j %>% filter(is.na(dbh)==T) # these jerks don't have dbh
  
  # Here we are getting the first dbh reported in the subsetted census for all those trees in jnk_A missing dbh
  jnk_B <- census %>% filter(tree_tag %in% jnk_A$tree_tag) %>% 
    arrange(census_date) %>% 
    select(tree_tag, census_date, dbh, height_m,wd, canopy_layer) %>% 
    group_by(tree_tag) %>% 
    filter(dbh==first(dbh)) %>% 
    group_by(tree_tag) %>% 
    filter(census_date==first(census_date))
  
  jnk_C <- baseline_j %>% filter(is.na(dbh)==F)
  baseline_j <- bind_rows(jnk_B, jnk_C)
  
  # --- ESTIMATE TREE CANOPY POSITION FOR EACH DEND MEASURMENT DATE ---
  vec_tree_tags <- dend %>% pull(tree_tag) %>% unique() %>%  sort()
  dend$canopy_layer_pred <- NA
  for(i in 1:length(vec_tree_tags)){
    tmp_fit <- census %>% filter(tree_tag == vec_tree_tags[i]) %>% 
      lm(canopy_layer~date, data=.)
    dend$canopy_layer_pred[dend$tree_tag==vec_tree_tags[i]] <- 
      predict(tmp_fit, newdata=data.frame(date=dend$date[dend$tree_tag==vec_tree_tags[i]]))
  }
  
  # Check tree_tags are the same in dendrometer and census datasets
  unique(dend$tree_tag) %in% unique(census$tree_tag) %>% table # --- REQUIRES ATTENTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!--- #
  print(paste("!!! dendrometer tree_tag not in the census:",dend$tree_tag[!(dend$tree_tag %in% census$tree_tag)] %>% unique))
  
  ## get the data you need from the census file into the dendrometer data frame: wd, height, first dbh measurement, date of first dbh measurement
  # str(dend) # check you have data in here. If not, make sure dend1 and cen are in the right formats, e.g. using sapply(cen, class).
  
  dend <- left_join(dend,
                    baseline_j, 
                    by="tree_tag")
  dend <- dend %>% filter(is.na(dbh)==F) # filter out the NAs
  
  ## Allometric equation option. Set of allometric equations after Chave et al. 2005
  ## and Chave et al. 2014 are defined in allometricEquations.R. Options defined here:
  if(allometric_option == 2 | allometric_option == "dry") {
    allometrix <- 2
    print("dry equation is used for estimating AGB, model I.3 (see Chave et al., 2005)")
  } else if (allometric_option == 3 | allometric_option == "moist" | allometric_option == "Default" | allometric_option == 1) {
    allometrix <- 3
    print("moist equation is used for estimating AGB, model I.6 (see Chave et al., 2005)")
  } else if (allometric_option == 4 | allometric_option == "wet") {
    allometrix <- 4
    print("wet equation is used for estimating AGB, model I.3 (see Chave et al., 2005)")
  } else if (allometric_option == 5 | allometric_option == "Chave2014") {
    allometrix <- 5
    print("pantropical equation is used for estimating AGB, model (4) (see Chave et al., 2014)")
  } else {
    print("Please specify a valid allometric_option!")
    return()
  }
  
  # data cleaning
  dend$dendrometer_reading_mm[which(dend$dendrometer_reading_mm > 2000)] <- NaN
  
  dend <- dend %>% mutate(dbh_first_date = census_date, 
                          date=lubridate::parse_date_time(paste(year,month,day), "ymd"))
  count_obs_by_date <- dend %>% group_by(date) %>% summarize(n_obs=n())
  dates_too_few_obs <- count_obs_by_date %>% filter(n_obs<5) %>% pull(date)
  dend <- dend %>% filter(!date %in% dates_too_few_obs) # This gets ride of the date(s) where too few observations were made. Arbitrarily set at 5 for now!!! 
  
  # !!! IMPORTANT, the dendrometer band readings are assumed to be reported in cm if the min reading is not <0 #
  if(min(dend$dendrometer_reading_mm, na.rm=T) <0 | 
     max(dend$dendrometer_reading_mm,na.rm=T)>100 | 
     # unique(dend$plot_code)=="NXV-01" | 
     unique(dend$plot_code) !="NXV-02" |
     unique(dend$plot_code)=="ANK-01") {
    are_dend_units_correct = T
  }else{
    are_dend_units_correct = F
    dend$dendrometer_reading_mm <- dend$dendrometer_reading_mm*10
  }
  
  if(unique(dend$plot_code) == "MNG-03" | unique(dend$plot_code)=="MNG-04" | unique(dend$plot_code)=="NXV-02"){
    dend$dendrometer_reading_mm <- dend$dendrometer_reading_mm*10
  }
  
  
  ################################################################################################
  # --- begin new stff --- THIS PART IS CRUCIAL FOR CATCHING THE ERRORS AND NEW DBANDS
  ################################################################################################
  dend$baseline_dbh <- dend$dbh
  dend <- dend %>% mutate(dendrometer_reading_mm = 
                            ifelse(is.na(dendrometer_reading_replaced_mm)==F &
                                     dendrometer_reading_replaced_mm<dendrometer_reading_mm, 
                                   dendrometer_reading_replaced_mm, dendrometer_reading_mm))
  dend$current_dbh <- est_dbh_from_dendro(dbh=dend$baseline_dbh, dendrometer_reading_mm = as.numeric(dend$dendrometer_reading_mm))
  
  pass1 <- dend %>% 
    group_by(tree_tag) %>% 
    arrange(date) %>% 
    mutate(delta1 = current_dbh-lag(current_dbh,order_by = date)) %>% 
    mutate(delta1_std = as.double(delta1/current_dbh)) %>%
    mutate(baseline_dbh = if_else(delta1_std < -0.0075, lag(current_dbh), as.double(baseline_dbh))) %>%  # if the tree shrinks too much 
    mutate(baseline_dbh = if_else(is.na(baseline_dbh)==T, lead(baseline_dbh), baseline_dbh)) %>% 
    mutate(baseline_dbh1 = cummax(baseline_dbh)) %>%
    mutate(current_dbh = est_dbh_from_dendro(dbh=baseline_dbh1, dendrometer_reading_mm = dendrometer_reading_mm)) 
  
  pass2 <- pass1 %>% 
    group_by(tree_tag) %>% 
    arrange(date) %>% 
    mutate(delta1 = current_dbh-lag(current_dbh,order_by = date)) %>% 
    mutate(delta1_std = as.double(delta1/current_dbh)) %>%
    mutate(baseline_dbh = if_else(delta1_std < -0.0075, lag(current_dbh), as.double(baseline_dbh))) %>% 
    mutate(baseline_dbh = if_else(is.na(baseline_dbh)==T, lead(baseline_dbh), baseline_dbh)) %>% 
    mutate(baseline_dbh1 = cummax(baseline_dbh)) %>%
    mutate(current_dbh = est_dbh_from_dendro(dbh=baseline_dbh1, dendrometer_reading_mm = dendrometer_reading_mm)) 

  dend <- pass2
  dend$thisdbh_cm <- dend$current_dbh/10
  dend$agC_Mg <- NA
  dend %>% names()
  dend$height_pred <- predict(mod_height, newdata=dend, type="response")
  rm(pass1, pass2)
  ################################################################################################
  # --- END ERROR CHECK!!!
  ################################################################################################
  
  #--- Unusually Large DBH-increment check ---------------------------------------------------#
  library(purrr)
  vec_gte3 <- dend %>% group_by(tree_tag) %>% summarize(nobs=n()) %>% 
    filter(nobs>=3) %>% pull(tree_tag)
  
  vec_max_resids <- dend %>% filter(tree_tag %in% vec_gte3) %>%
    split(.$tree_tag) %>% 
    map(~lm(current_dbh~date, data=.)) %>% 
    map(residuals) %>% 
    map_dbl(max)
  
  df_resids <- tibble(tree_tag=vec_max_resids %>% names() %>% as.character(), 
                      max_resid=vec_max_resids)
  
  dend <- left_join(dend, df_resids, by="tree_tag")
  dend <- dend %>% filter(max_resid<25); #judgement call, not a great fix
  #-------------------------------------------------------------------------------------------#  
  
  # estimate biomass of each tree for each new thisdbh_mm
  # loop through each tree to estimate biomass (bm) and convert to above ground carbon (agC)
  for(ii in 1:length(unique(dend$tree_tag))){  
    thistree <- which(dend$tree_tag==unique(dend$tree_tag)[ii])
    dbh_tree <- dend$thisdbh_cm[thistree]
    den_tree <- dend$wd[thistree]
    h_tree   <- dend$height_m[thistree]
    
    # this uses allometric equations from allometricEquations.R
    if (allometrix == 2) {
      bm <- Chave2005_dry(diax=dbh_tree, wd=den_tree, height=h_tree)
    } else if (allometrix == 3) {
      bm <- Chave2005_moist(diax=dbh_tree, wd=den_tree, height=h_tree)
    } else if (allometrix == 4) {
      bm <- Chave2005_wet(diax=dbh_tree, wd=den_tree, height=h_tree)
    } else if (allometrix == 5) {
      bm <- Chave2014(diax=dbh_tree, wd=den_tree, height=h_tree)
    }
    
    # Unit conversions 
    dend$agC_Mg[thistree] <- (bm)*(1/(2.1097*1000)) # convert kg to Mg=1/1000=10 and convert to carbon = 47.8% (ADD REF!! Txx et al?)
  }
  
  dend <- dend %>% mutate(date=parse_date_time(paste(year,month,15),"ymd"))
  
  # NPPacw per tree: substact bm(t) - bm(t+1) / (t)-(t+1)
  # 1. Tree above ground Carbon stock difference
  ################################################################################################
    dend <- dend %>% arrange(date) %>% filter(is.na(dendrometer_reading_mm)==F)
  
  uid             <- unique(dend$tree_tag) %>% sort
  aa              <- c()
  bb              <- c()
  cc              <- c()
  dd              <- c()
  ee              <- c()
  ff              <- c()
  vec_date_diffs              <- c()
  vec_current_dbh <- c() 
  vec_pmer        <- c()
  vec_tag <- c()
  vec_canopyLayer <- c()
  vec_height <- c()
  vec_wd <- c()
  vec_dates <- c()
  
  for (ii in 1:length(uid)){  
    thistree  <- which(dend$tree_tag == uid[ii]) #dend$tree_tag[ii])
    if(length(thistree)==1){next}
    agC_Mg       <- dend$agC_Mg[thistree]
    tag       <- dend$tree_tag[thistree]
    #agCdiff   <- dend$agCdiff[thistree]
    agC_Mg_diff   <- ave(dend$agC_Mg[thistree], FUN = function(x) c(NA, diff(x)))
    tmpyear      <- dend$year[thistree]
    tmpmonth     <- dend$month[thistree]
    tmpplot_code <- dend$plot_code[thistree]
    tmpDates <- dend$date[thistree]
    # tmpDayDiffs <- get_time_diffs(tmpDates) # This returns 'days, but it's much much slower than diff
    # But diff is "type unstable" meaning that it tries to be clever and return either days or seconds. 
    # That causes problems, so if we force dates to be seconds with as.numeric(), it removes the ambiguity and we get back seconds. 
    tmpDayDiffs <- diff(as.numeric(tmpDates))/(24*60*60) # THIS IS MUCH FASTER THAN get_time_diffs(). Force convert to date to seconds then transform back to days.
    tmpDayDiffs <- c(NA,tmpDayDiffs) # need to tack an NA to the front of the dayDiff vector
    print(tmpDayDiffs %>% length)
    dbhs <- dend$current_dbh[thistree]
    tree_tags <- as.character(dend$tree_tag[thistree])
    # this_pos_meas_errors <- dend$pos_meas_err[thistree]
    this_canopy_layer <- dend$canopy_layer_pred[thistree]
    this_height <- dend$height_pred[thistree]
    this_wd <- dend$wd[thistree]
    this_date <- dend$date[thistree]
    
    aa            <- c(aa, tmpplot_code)
    bb            <- c(bb, tag)
    cc            <- c(cc, tmpyear)
    dd            <- c(dd, tmpmonth)
    ee            <- c(ee, agC_Mg)
    ff            <- c(ff, agC_Mg_diff)
    vec_date_diffs            <- c(vec_date_diffs, tmpDayDiffs)
    vec_current_dbh            <- c(vec_current_dbh, dbhs)
    vec_tag       <- c(vec_tag, tree_tags)
    # vec_pmer       <- c(vec_pmer, this_pos_meas_errors)
    vec_canopyLayer <- c(vec_canopyLayer, this_canopy_layer); 
    vec_height <- c(vec_height, this_height); 
    vec_wd        <- c(vec_wd, this_wd)
    vec_dates     <- c(vec_dates, this_date);
    
    # if(ii%%100 == 0){print(ii)}
  }
  
  npp_tree <- tibble(plot_code=aa, tag=bb, year=cc, month=dd, 
                     agC_Mg=ee, agC_Mg_diff=ff, dateDiff=vec_date_diffs, current_dbh=vec_current_dbh,
                     tree_tag=as.character(vec_tag), 
                     # pos_meas_error=vec_pmer, 
                     # canopyLayer_pred = vec_canopyLayer,
                     height_pred = vec_height,
                     wd = vec_wd, 
                     date=vec_dates, 
                     date_diff=vec_date_diffs) %>% 
    mutate(date=as.POSIXct(date, origin="1970-01-01", tz = "UTC"))
  npp_tree <- npp_tree %>% filter(dateDiff>0) #%>% # otherwise we're dividing by 0!
  # filter(pos_meas_error != T) # filtering out possible measurement errors 
  
  # 3. NPP: MgC per tree per day #!!! SOMETHING IS CAUSING OCCASIONAL VERY LOW VALS - NEEDS CHECKING!!!
  npp_tree <- npp_tree %>% mutate(agC_Mg_diff_std=(agC_Mg_diff/agC_Mg))
  npp_tree <- npp_tree %>% filter(agC_Mg_diff_std>-0.05) # allowing for some shrinkage and measurement error. Still some unusually low vals!!!
  npp_tree$nppacw_Mg_tree_day  <- npp_tree$agC_Mg_diff/npp_tree$dateDiff
  
  # "This section is for estimating agC_Mg with multiple censuses"  
  #   tmp_npp_census_trees <- census %>% filter(dbh>=90) %>% 
  #     # filter(f2=="1") %>%  # Comment this out if you want to include dead trees
  #     rowwise() %>% 
  #     mutate(biomass_tree_kg=sum(Chave2014_eq4(diam_cm=dbh/10, wd=wd, height_m=height_m), na.rm=T), 
  #            basal_area_tree_cm=sum(pi*(0.1*dbh/2)^2, na.rm=T)) %>% 
  #     mutate(agC_Mg = (1/(2.1097*1000))*biomass_tree_kg) %>% ungroup()  
  #   
  #   tmp_npp_census <- tmp_npp_census_trees %>% group_by(year,month) %>% 
  #     summarize(n_trees=n(), 
  #               basal_area_cm2 = sum(basal_area_tree_cm,na.rm=T), 
  #               biomass_kg = sum(biomass_tree_kg, na.rm=T)) %>% 
  #     mutate(agC_Mg = (1/(2.1097*1000))*biomass_kg, 
  #            date=parse_date_time(paste(year,month,15), "ymd"))
  #   
  #   tmp_npp_census$day_interval <- c(NA,get_time_diffs(tmp_npp_census$date))
  #   tmp_npp_census$agC_Mg_diff <- c(NA, diff(tmp_npp_census$agC_Mg))
  #   tmp_npp_census$agC_Mg_day <- tmp_npp_census$agC_Mg_diff/tmp_npp_census$day_interval
  
  #-------------------------------------------------------------------------------------------------------------------
  # NPP of the plot from the dendrometer bands --------------------------------------------------------------------------------
  # THIS NEEDS TO BE REALITY CHECKED AGAINST THE CENSUS VALS !!!
  # Steps: 
  npp_tree <- npp_tree %>% mutate(decDate=decimal_date(parse_date_time(paste(year,month,15),"ymd")), 
                                  basal_area_cm2 = pi*((0.1*current_dbh)/2)^2)
  vec_census_date <- census %>% pull(census_date) %>% unique %>% sort
  census$basal_area_cm2 <- pi*(0.1*census$dbh/2)^2; 
  npp_tree$closestCensusDate <- NA; 
  npp_tree$weight <- NA
  
  # calculate which census date is the closest to the dendrometer measurement
  for(i in 1:dim(npp_tree)[1]){
    npp_tree$closestCensusDate[i] <- which(abs(npp_tree$decDate[i] - vec_census_date)==min(abs(npp_tree$decDate[i] - vec_census_date)))
  }
  
  # Calculate the basal area weighted plot mean growth rate in MgC ha-1 month-1
  tmp_npp_dend <- npp_tree %>% group_by(plot_code, year, month) %>% 
    summarize(npp_avgtrees_day_dend=mean(nppacw_Mg_tree_day,na.rm=T), 
              npp_avgtrees_day_dend_sd=sd(nppacw_Mg_tree_day, na.rm=T), 
              npp_weightedAvgTrees_day_dend = weighted.mean(nppacw_Mg_tree_day, w=weight, na.rm=T),
              nobs=n()) %>% 
    mutate(date=parse_date_time(paste(year,month,15), "ymd"), 
           npp_weightedAvgTrees_month_dend = npp_weightedAvgTrees_day_dend*30.4, 
           npp_weightedAvgTrees_yr_dend=npp_weightedAvgTrees_month_dend*12)
  
  
  #!!!!! EXPERIMENTAL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # Predict growth for each tree in census with a linear model
  # Then sum all predictions to avoid the issues with sampling bias. 
  # Should solve dendrometer sampling bias towards either large or small trees. 
  
  # tmp_preds <- npp_tree %>% group_by(plot_code, year, month) %>% 
  #   summarize(nobs=n(), closestCensusDate=unique(closestCensusDate)) %>% 
  #   mutate(date=parse_date_time(paste(year,month,15), "ymd")) %>% 
  #   ungroup() %>% group_by(date) %>% 
  #   mutate(n_trees_census=tmp_npp_census[which.min(abs(tmp_npp_census$date-date)),]$n_trees)
  
  # filter out trees that "appear" to have lost more than 5% of their carbon
  npp_tree <- npp_tree %>% 
    mutate(agC_Mg_diff_std = as.double(agC_Mg_diff/agC_Mg)) %>% 
    filter(agC_Mg_diff_std > -0.05)
  
  vec_dend_date <- npp_tree %>% group_by(date) %>% 
    summarize(nobs=n()) %>% 
    filter(nobs>=25) %>% pull(date) %>% unique() %>% sort()
  # tmp_preds$npp_month_dend <- NA
  
  census_max_basal_area_date <- census %>% group_by(census_date) %>% 
    summarize(ba=sum(basal_area_cm2)) %>% 
    filter(ba==max(ba)) %>% pull(census_date)
  census_max_basal_area <- census %>% 
    filter(dbh>=90) %>% 
    filter(near(census_date,max(census_date))) %>% 
    rename(canopyLayer_pred = canopy_layer) %>% 
    mutate(height_pred = height_m)
  
  ## --- just for NXV-02 --- #
  # census_max_basal_area <- census %>% filter(census_date==2011.255) %>% 
  # filter(dbh>=100) %>% 
  #   filter(near(census_date,max(census_date))) %>% 
  #   rename(canopyLayer_pred = canopy_layer) %>% 
  #   mutate(height_pred = height_m)
  
  npp_tree <- npp_tree %>% mutate(dbh=current_dbh)
  
  growth_lm_mgc_vec <- numeric()
  # growth_nl_mgc_vec <- numeric()
  growth_gam_mgc_vec <- numeric()
  growth_hybrid_gam_mgc_vec <- numeric()
  growth_hybrid_bestMod_mgc_vec <- numeric()
  growth_u_mgc_vec <- numeric()
  growth_wu_mgc_vec <- numeric()
  growth_sd_mgc_vec <- numeric()
  growth_bestEst_vec <- numeric()
  nobs_vec <- numeric()
  vec_best_dev <- numeric()
  vec_date_diffs <- numeric()
  vec_medtree_ratio <- numeric(); 
  vec_which_est <- c()
  
  for(idx in 1:length(vec_dend_date)){
    print(idx)
    tmp_dat <- npp_tree %>% filter(date==vec_dend_date[idx]) %>%
      filter(dbh>=90) %>%
      mutate(rat_growth = nppacw_Mg_tree_day/agC_Mg) %>% filter(rat_growth>=-0.0015)
    tmp_tree_tags <- tmp_dat$tree_tag %>% unique() %>% sort()
 
    medtree_ratio <- (tmp_dat %>% pull(current_dbh) %>% median())/(census_max_basal_area$dbh %>% median())
    vec_medtree_ratio <- c(vec_medtree_ratio, medtree_ratio)
    
    ntrees_census <- dim(census_max_basal_area)[1]
    nobs <- dim(tmp_dat)[1]; print(nobs); 
    nobs_vec <- c(nobs_vec, nobs)
    
    u_date_diff <- tmp_dat$date_diff %>% mean()
    vec_date_diffs <- c(vec_date_diffs, u_date_diff)
    
    growth_u_mgc <- tmp_dat %>% 
      summarize(u=mean(nppacw_Mg_tree_day, na.rm=T)) %>% pull(u)
    growth_u_mgc <- growth_u_mgc*30.4*dim(census_max_basal_area)[1]
    growth_u_mgc_vec <- c(growth_u_mgc_vec, growth_u_mgc)
    
    tmp_wfun <- approxfun(density(census_max_basal_area$dbh))
    wu_weights <- tmp_wfun(ifelse(tmp_dat$dbh>max(census_max_basal_area$dbh),max(census_max_basal_area$dbh),tmp_dat$dbh))
    growth_wu_mgc <- weighted.mean(tmp_dat$nppacw_Mg_tree_day,
                                   w = wu_weights,na.rm=T)*30.4*dim(census_max_basal_area)[1]
    growth_wu_mgc_vec <- c(growth_wu_mgc_vec, growth_wu_mgc)
    

    tmp_fit_lm <- gam(nppacw_Mg_tree_day~ height_pred+dbh*wd+I(dbh**2),
                      data=tmp_dat, weights = wu_weights, method='REML')
    tmp_fit_lm2 <- gam(nppacw_Mg_tree_day~ height_pred+dbh+wd,
                       data=tmp_dat, weights = wu_weights, method='REML')
    
    growth_lm_mgc <- predict(tmp_fit_lm, type = "response",
                             newdata = census_max_basal_area) %>% sum()
    growth_lm_mgc_vec <- c(growth_lm_mgc_vec, growth_lm_mgc*30.4)
    
    
    tmp_fit_gam1 <- gam(nppacw_Mg_tree_day~ti(wd,dbh, k=3)+s(height_pred,k=2),
                        data=tmp_dat,
                        select = T,
                        method='REML', weights = wu_weights)
    
    tmp_fit_gam4 <- gam(nppacw_Mg_tree_day~wd+s(dbh)+s(height_pred),
                        data=tmp_dat,
                        select = T, method='REML', 
                        weights = wu_weights)
    
    vec_mod_r2 <- c(summary(tmp_fit_lm)$r.sq, 
                    summary(tmp_fit_lm2)$r.sq, 
                    summary(tmp_fit_gam1)$r.sq, 
                    summary(tmp_fit_gam4)$r.sq); 
    max_r2 <- vec_mod_r2 %>% which.max()
    
    vec_mod_dev <- c(summary(tmp_fit_lm)$dev.expl, 
                     summary(tmp_fit_lm2)$dev.expl, 
                     summary(tmp_fit_gam1)$dev.expl, 
                     summary(tmp_fit_gam4)$dev.expl); 
    max_dev <- vec_mod_dev %>% which.max()
    
    mod_list <- list(tmp_fit_lm, tmp_fit_lm2, tmp_fit_gam1, tmp_fit_gam4)
    
    best_mod_dev <- summary(mod_list[[max_dev]])$dev.expl
    vec_best_dev <- c(vec_best_dev, best_mod_dev)
    
    growth_gam_mgc <- predict(tmp_fit_gam1, type = "response",
                              newdata = census_max_basal_area) %>% sum(na.rm=T)
    growth_gam_mgc_vec <- c(growth_gam_mgc_vec, growth_gam_mgc*30.4)
    
    
    growth_hybrid_gam_mgc <- sum(predict(tmp_fit_gam1, 
                                         type="response", 
                                         newdata = census_max_basal_area %>% 
                                           filter(!tree_tag %in% tmp_tree_tags)),na.rm=T)+sum(tmp_dat$nppacw_Mg_tree_day,na.rm=T)
    growth_hybrid_gam_mgc <- growth_hybrid_gam_mgc*30.4
    growth_hybrid_gam_mgc_vec <- c(growth_hybrid_gam_mgc_vec, growth_hybrid_gam_mgc)
    
    print(summary(mod_list[[max_dev]]))
    
    growth_hybrid_bestMod_mgc <- sum(predict(mod_list[[max_dev]], 
                                             type="response", 
                                             newdata = census_max_basal_area %>% 
                                               filter(!tree_tag %in% tmp_tree_tags)),na.rm=T)+sum(tmp_dat$nppacw_Mg_tree_day,na.rm=T)
    growth_hybrid_bestMod_mgc <- growth_hybrid_bestMod_mgc*30.4
    growth_hybrid_bestMod_mgc_vec <- c(growth_hybrid_bestMod_mgc_vec, growth_hybrid_bestMod_mgc)
    
    
    growth_sd_mgc <- npp_tree %>% filter(date==vec_dend_date[idx]) %>%
      filter(dbh>=90) %>%
      filter(nppacw_Mg_tree_day>=-0.0005) %>%
      summarize(sd=mean(nppacw_Mg_tree_day, na.rm=T)) %>% pull(sd)
    growth_sd_mgc_vec <- c(growth_sd_mgc_vec, growth_sd_mgc*30.4*dim(census_max_basal_area)[1])
    
    if(nobs>=(ntrees_census*0.5)&(near(medtree_ratio, 1.0, tol=0.15))){
      growth_bestEst_vec <- c(growth_bestEst_vec, growth_u_mgc)
      vec_which_est <- c(vec_which_est, "growth_u_mgc")
    }else if( (nobs>(ntrees_census*0.3)) & 
              (near(medtree_ratio, 1.0, tol=0.1))){
      growth_bestEst_vec <- c(growth_bestEst_vec, growth_u_mgc)
      vec_which_est <- c(vec_which_est, "growth_u_mgc")
    }else if( (nobs<(ntrees_census*0.5)) & 
              (near(medtree_ratio, 1.0, tol=0.1)==F &
              (best_mod_dev < 0.2))){
      growth_bestEst_vec <- c(growth_bestEst_vec, growth_wu_mgc)
      vec_which_est <- c(vec_which_est, "growth_wu_mgc")
    }else{
      tmp_x <- c(growth_wu_mgc, growth_hybrid_bestMod_mgc)
      tmp_x1 <- c("growth_wu_mgc", "growth_hybrid_bestMod_mgc")
      growth_bestEst_vec <- c(growth_bestEst_vec, tmp_x[which.max(tmp_x)])
      vec_which_est <- c(vec_which_est, tmp_x1[which.max(tmp_x)])
    }
    

      
    # if(nobs>=(ntrees_census*0.5)){
    #   growth_bestEst_vec <- c(growth_bestEst_vec, growth_u_mgc)
    #   vec_which_est <- c(vec_which_est, "growth_u_mgc")
    #   print("best estimate is mean")
    # }else if(nobs>=(ntrees_census*0.25) & best_mod_dev<0.2){
    #   growth_bestEst_vec <- c(growth_bestEst_vec, growth_u_mgc)
    #   vec_which_est <- c(vec_which_est, "growth_u_mgc")
    #   print("best estimate is mean")
    # }else if(nobs<(ntrees_census*0.25) & best_mod_dev<0.2){
    #   growth_bestEst_vec <- c(growth_bestEst_vec, growth_wu_mgc)
    #   vec_which_est <- c(vec_which_est, "growth_wu_mgc")
    #   print("best estimate is weighted mean")
    # }else if((medtree_ratio>1.1)&
    #          (nobs<=(ntrees_census*0.5))&
    #          (best_mod_dev>0.3)&
    #          growth_hybrid_bestMod_mgc>0){
    #   growth_bestEst_vec <- c(growth_bestEst_vec, growth_hybrid_bestMod_mgc)
    #   vec_which_est <- c(vec_which_est, "growth_hybrid_bestMod_mgc")
    #   print("best estimate is hybrid model")
    # }else if(nobs<=(ntrees_census*0.5) & 
    #          (best_mod_dev>0.3)&
    #          (medtree_ratio<1.1)&
    #          ((growth_hybrid_bestMod_mgc/growth_u_mgc)>1.5 |
    #          (growth_hybrid_bestMod_mgc/growth_u_mgc)<0.5)){
    #   growth_bestEst_vec <- c(growth_bestEst_vec, growth_u_mgc)
    #   vec_which_est <- c(vec_which_est, "growth_u_mgc")
    # }else{
    #   growth_bestEst_vec <- c(growth_bestEst_vec, growth_hybrid_bestMod_mgc)
    #   vec_which_est <- c(vec_which_est, "growth_hybrid_bestMod_mgc")
    #   print("best estimate is hybrid model")
    # }
    
  }  
  #!!! END REGRESSION SCALING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  out_df <- tibble(date=vec_dend_date, 
                   npp_lm_agC_month=growth_lm_mgc_vec, 
                   npp_u_agC_month=growth_u_mgc_vec,
                   npp_wu_agC_month=growth_wu_mgc_vec,
                   npp_gam_agC_month=growth_gam_mgc_vec,
                   npp_hybrid_gam_agC_month = growth_hybrid_gam_mgc_vec,
                   npp_hybrid_bestMod_agC_month = growth_hybrid_bestMod_mgc_vec,
                   npp_sd_agC_month=growth_sd_mgc_vec, 
                   npp_bestEst_agC_month=growth_bestEst_vec, 
                   which_bestEst = vec_which_est,
                   nobs=nobs_vec, 
                   dev=vec_best_dev,
                   medtree_ratio = vec_medtree_ratio,
                   date_diff = vec_date_diffs,
                   plot_code=this_plot_code, 
                   site=substr(this_plot_code,1,3)) 
  
  out_df <- left_join(out_df, 
                      npp_tree %>% group_by(date) %>% summarize(nobs=n()), 
                      by="date")
  
  out_df <- out_df %>% mutate(plot_code=dend$plot_code %>% unique())
  # out_df <- left_join(tmp_npp_dend, n_trees_df, by="date")
  
  switch(ret,
         nppacw.permonth.perha = {return(out_df)},
         nppacw.perday.pertree = {return(npp_tree)}
  )
  
}



###############################################################################
##---- SCRATCH 
# try(tmp_fit_nl <- nls(nppacw_Mg_tree_day~b2 +b0*exp((dbh**b1)/1500),
#     data=tmp_dat,
#     start = list(b0=0.00003, b1=1.15, b2=0),
#     control=list(maxiter = 100000, warnOnly=T),
#     lower=c(0,0,-1e-04), upper=c(0.001, 3,5e-04), algorithm = "port"), silent = T)
# 
# if(exists("tmp_fit_nl")==F){
#   tmp_fit_nl <- nls(nppacw_Mg_tree_day~ b2*height_pred+b0*exp(((dbh)^b1)/1500), 
#                     data=tmp_dat,
#                     start = list(b0=0.00003, b1=1.15, b2=0),
#                     control=list(maxiter = 100000, warnOnly=T),
#                     lower=c(0,0,-1e-04), upper=c(0.001, 3,5e-04), algorithm = "port")
# }

# growth_nl_mgc <- predict(tmp_fit_nl, type = "response",
#                          newdata = census_max_basal_area) %>% sum()
# growth_nl_mgc_vec <- c(growth_nl_mgc_vec, growth_nl_mgc*30.4)

# tmp_fit_gam <- gam(nppacw_Mg_tree_day~s(scale(dbh),bs="gp")+wd+s(height_pred,bs="gp")+s(scale(I(dbh**2)),bs="gp"),
#                    select=T,gamma=3,data=tmp_dat)
# tmp_fit_gam <- gam(nppacw_Mg_tree_day~s(scale(dbh),bs="gp")+wd+s(height_pred,bs="gp"),
#                    select=T,gamma=3,data=tmp_dat) 
