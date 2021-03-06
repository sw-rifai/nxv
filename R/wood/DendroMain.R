##############################################################################################
### NPPdendrometers REFACTORED FOR EASIER DEBUGGING WHEN THINGS BREAK. THEY WILL BREAK OFTEN.#
##############################################################################################
#--- NOTES ----------------------------------------------------------------------------------#
# ! EACH PLOT REQUIRES EXTENSIVE CHECKING TO ENSURE THE CALCULATIONS ARE WORKING CORRECTLY ! #
#DO:                                                                          CorrectUnits?
# "NXV-01" "NXV-02"                                                                   T
#--- END NOTES ------------------------------------------------------------------------------#

# rm(list=ls()) #--- CLEAR THE ENVIRONMENT (?)---

try(dev.off(),silent=T);
library(tidyverse); library(lubridate)
source("R/wood/DendroFunctions.R"); 
source("R/wood/DEV_DendroCalcAgC_Mg_v5.R") # DEVELOPMENT VERSION

# LOAD DATA ---------------------------------------------------------------------------------#
dend_all <- read_csv("data/wood/dend_nxv.csv", na=c("NA", "NaN", ""), guess_max = 100000)
census_all  <- read_csv("data/wood/census_nxv.csv", na=c("NA", "NaN", ""), guess_max = 100000) # Census_Santarem_2014_2016.csv
names(census_all) <- tolower(names(census_all))
names(census_all) <- gsub(pattern=" ",replacement="_",names(census_all))
census_all <- census_all %>% rename(tree_tag = tag_number) %>% mutate(tree_tag = as.character(tree_tag))
zanne <- read_csv("data/wood/GlobalWoodDensityDatabase.csv")
glimpse(census_all) 
# END LOAD DATA -----------------------------------------------------------------------------#

# USER OPTIONS ------------------------------------------------------------------------------#
# - Allometric options: [2-Chave2005,dry eqI.3][3-Chave2005,moist eqI.6][4-Chave2014,wet-eq4]
# - Height correction options: 
this_plot_code           = "NXV-02"
allometric_option        = "5" 
height_correction_option = "Default"
are_dend_units_correct   = F
# USER OPTIONS ------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------
#--- PREPROCESS DENDROMETER DATA 
# Steps: 1) guess whether units are correct [and then verify with census biomass increments diagnostic plot]
# 
# Site ---- dendrometer_reading_mm(mmOrCm?) [4.5 threshold?, check w/census]
# 
dendrometer <- dend_all %>% filter(plot_code==this_plot_code)
dendrometer <- dendrometer %>% 
  mutate(tree_tag = as.character(tree_tag), 
         date = parse_date_time(paste(year,month,day), "ymd")) #%>% 
dendrometer %>% pull(date) %>% is.na() %>% table() # number of NA dates
min_dend_date <- dendrometer %>% pull(date) %>% min(na.rm=T)
glimpse(dendrometer)


# site specific checks for idiosyncratic errors --- THIS IS A CRAP WAY TO DEAL WITH FUNDAMENTAL PROBLEMS, THE DATA CSV SHOULD BE CORRECTED.
if(substr(this_plot_code, 1,3)=="TAM"){
 dendrometer <- dendrometer %>% mutate(dendrometer_reading_mm = ifelse(year>=2015, dendrometer_reading_mm*10, dendrometer_reading_mm)) 
}
if(substr(this_plot_code, 1,3)=="CAX"){
  dendrometer <- dendrometer %>% mutate(dendrometer_reading_mm = dendrometer_reading_mm*100) 
}

#--- END PREPROCESS DENDROMETER DATA 
#--------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------
# PREPROCESS CENSUS DATA --------------------------------------------------------------------
# Steps: 
# 1) fix col names and dates
# 2) Get species list and wood densities from Zanne, 
# 3) Identify problem tree_tags
census <- census_all %>% filter(plot_code==this_plot_code)
census <- census %>% select(continent,country,plot_code,tree_tag,treeid,family,genus,species,census_date,
                            d0,d1,d2,d3,d4,f1,f2,f3,f4,f5,height)

# convert dates
census$date <- convDecDate(census$census_date)
census <- census %>% mutate(year=year(date), month=month(date), day=day(date))
census_year = min(census$year)

# fill in dbh values from forestPlots d0-d4
if((is.na(c(census$d0,census$d1,census$d2,census$d3,census$d4)) %>% sum())<100){
census <- get_dbh(census=census, min_dend_date=min_dend_date) # cens_dend_dist=5 for NXV-02 ... hack solution atm
}else{
  census$dbh <- census$d0;
}



# missing dbh value check
print(paste("are there missing dbh in the census?"))
census$dbh %>% is.na() %>% table()
census$dbh %>% is.infinite() %>% table()

# filter them out
census <- census %>% filter(is.na(dbh)==F & is.infinite(dbh)==F)

# get wood wd from Zanne global wood wd database
census$wd <- NA
census_sp_list <- census %>% distinct(family, genus, species, continent, tree_tag)

# identify the trees with multiple species ids 
# The GEM tags seems to be duplicated because these trees have different ForestPlots id numbers
bad_tags <- census_sp_list %>% group_by(tree_tag) %>% summarize(n_sp = length(unique(species)))
bad_tags <- bad_tags %>% filter(n_sp>1)
bad_tags_census <- census %>% filter(tree_tag %in% bad_tags$tree_tag)

# for now I'm just throwing those trees with duplicate tags out because it'll be 
# some work to figure out how to deal with them... 
census <- census %>% filter(!tree_tag %in% bad_tags$tree_tag)

# recast census species list
census_sp_list <- census %>% distinct(family, genus, species, continent)

census_sp_list$wd <- NA
census_sp_list <- find_wd(census_sp_list = census_sp_list)
census_sp_list <- census_sp_list %>% mutate(species=paste(genus,species))

# census_sp_list$species %in% census$species
# census$species %in% census_sp_list$species %>% table
census <- census %>% select(-wd)
census_sp_list <- census_sp_list %>% select(species, wd) %>% group_by(species) %>% 
  summarize(wd=mean(wd,na.rm=T)) %>% ungroup()

census <- left_join(census,census_sp_list, by="species")


#--- HEIGHT ESTIMATION -----------------------------------------------------------
# STEPS: filter census for trees with height measurements, filter heights from method of
# measurement, 
# f5 from the census tells how the tree was measured
# 1= Estimated by eye.
# 2= Manually by trigonometry (clinometer).
# 3= Manually by trigonometry (clinometer), carefully trained.
# 4= Laser or ultrasonic distance to tree, electronic tilt sensor for angle.
# 5= Laser hypsometer from directly below crown, “last return” filter function.
# 6= Directly (e.g. climbing, cutting, adjacent tower).

tmp_height_inv <- census_all %>% filter(substr(plot_code,1,3)==substr(this_plot_code,1,3)) %>% 
  filter(is.na(height)==F)
names(tmp_height_inv) <- tolower(names(tmp_height_inv))

if(substr(this_plot_code,1,3)=="CAX"){
  caxHeight <- read_csv("data/CAX_heights_cleaned.csv")
  log_fit <- lm(height~log(dbh), data=caxHeight %>% filter(is.na(wd)==F));
  lm_fit <- lm(height~dbh, data=caxHeight %>% filter(is.na(wd)==F))
  nl_fit <- nls(height~ b0*((1/wd)^b1)*dbh^b2, data=caxHeight,
                start = list(b0=2.643, b1=-0.3, b2=0.5))
  bbmle::AICctab(log_fit, lm_fit, nl_fit)
  mod_list <- list(log_fit, lm_fit, nl_fit)
  best_mod <- AIC(log_fit, lm_fit, nl_fit) %>% min
  
  tmp <- find_wd(census)
  census$wd <- tmp$wd
  census$height_m <- predict(mod_list[[best_mod]], 
                             newdata=data.frame(wd=census$wd, dbh=census$dbh))
  mod_height <- mod_list[[best_mod]]
  
}else if(substr(this_plot_code,1,3)=="KEN"){
  print("Special KEN routine for height")
  kenHeight <- read_csv("data/KEN_tree_heights.csv")
  names(kenHeight) <- tolower(names(kenHeight))
  kenHeight <- kenHeight %>% 
    rename(family=familia, genus=genero, species = especie, 
           dbh=`dap 20/12/2015`, height=`altura total (medida)`) %>% 
    mutate(continent="South America", species=paste(genus, species)) %>% 
    filter(dbh!="Muerto" & is.na(height)==F & height!="Liana" & height!="hemiepifito") %>% 
    mutate(dbh=as.double(dbh), height=as.double(height)) %>% 
    filter(is.na(height)==F) %>% 
    mutate(dbh=dbh*10)
  
  kenHeight_2 <- find_wd(census_sp_list = kenHeight)
  kenHeight_2 <- kenHeight_2 %>% rename(wd=wd) %>% 
    mutate(species=paste(genus, species))
  kenHeight <- left_join(kenHeight, kenHeight_2)

  log_fit <- lm(height~log(dbh), data=kenHeight %>% filter(is.na(wd)==F));
  lm_fit <- lm(height~dbh, data=kenHeight %>% filter(is.na(wd)==F))
  nl_fit <- nls(height~ b0*((1/wd)^b1)*dbh^b2, data=kenHeight,
                start = list(b0=2.643, b1=-0.3, b2=0.5))
  bbmle::AICctab(log_fit, lm_fit, nl_fit)
  mod_list <- list(log_fit, lm_fit, nl_fit)
  best_mod <- AIC(log_fit, lm_fit, nl_fit) %>% min
  
  census$wd <- census$wd
  census$height_m <- predict(mod_list[[best_mod]], 
                             newdata=data.frame(wd=census$wd, dbh=census$dbh))
  census <- census[census$f1!="0",] # REMOVE THE DEAD TREES
  mod_height <- mod_list[[best_mod]]
}else if(substr(this_plot_code,1,3)=="BOB"){
  print("Special Bobiri routine")
# --- UNCOMMENT FOR BOBIRI PLOTS ---
# this step is just for BOB because it's unclear what plot the heights come from 
bobHeight <- read_csv("data/bobHeights.csv")
log_fit <- lm(height~log(dbh), data=bobHeight %>% filter(is.na(wd)==F));
lm_fit <- lm(height~dbh, data=bobHeight %>% filter(is.na(wd)==F))
nl_fit <- nls(height~ b0*((1/wd)^b1)*dbh^b2, data=bobHeight,
              start = list(b0=2.643, b1=-0.3, b2=0.5))
AIC(log_fit, lm_fit, nl_fit)
census$height_m <- predict(nl_fit, newdata=data.frame(wd=census$wd, dbh=census$dbh))
census <- census[census$f1!="0",] # REMOVE THE DEAD TREES
mod_list <- list(log_fit, lm_fit, nl_fit)
best_mod <- AIC(log_fit, lm_fit, nl_fit) %>% min
# --- END BOBIRI SPECIAL SECTION ---
}else if(sum(is.na(tmp_height_inv$height)==F) < 50){
  print("Using Feldpauch diameter approximation")
  #--- IF MISSING ANY FIELD OBSERVED HEIGHT, USE FELDPAUCH EQ
  census$height_m <- feldpauch_height(dbh_cm =census$dbh/10)
  census <- census[census$f1!="0",] # REMOVE THE DEAD TREES
  census <- census %>% filter(is.na(dbh)==F)
} else{
  print("Using best of linear/log/non-linear allometric fit")
  # USE best fit model 
tmpHeight <- tmp_height_inv %>% filter(is.na(height)==F & f5>1)
tmpHeight <- tmpHeight %>% rowwise() %>% mutate(dbh=mean(c(d0,d1,d2,d3,d4),na.rm=T))
tmpHeight <- tmpHeight %>% filter(dbh>0) %>% ungroup()
tmpHeight$f5 %>% table

tmpHeight_sp_list <- tmpHeight %>% distinct(family, genus, species, continent)
tmpHeight_sp_list$wd <- NA
tmpHeight_sp_list <- find_wd(census_sp_list = tmpHeight_sp_list)
tmpHeight_sp_list <- tmpHeight_sp_list %>% mutate(species=paste(genus,species))

tmpHeight <- left_join(tmpHeight, 
                       tmpHeight_sp_list %>% select(species, wd), by="species")

tmpHeight <- tmpHeight %>% group_by(tree_tag, treeid, family, genus, species) %>% 
  summarize(height=max(height),dbh=max(dbh), wd=mean(wd))

log_fit <- lm(height~log(dbh), data=tmpHeight %>% filter(is.na(wd)==F));
lm_fit <- lm(height~dbh, data=tmpHeight %>% filter(is.na(wd)==F))
nl_fit <- nls(height~ b0*((1/wd)^b1)*dbh^b2, data=tmpHeight %>% filter(is.na(wd)==F),
              start = list(b0=2.643, b1=-0.3, b2=0.5))
mod_list <- list(log_fit, lm_fit, nl_fit)
best_mod <- AIC(log_fit, lm_fit, nl_fit)[,2] %>% which.min()
print("Best mod is: ***********************************************************************")
print(mod_list[[best_mod]])
curve(predict(mod_list[[1]], newdata=data.frame(dbh=x)), 100,1000, ylim=c(0,50))
curve(predict(mod_list[[2]], newdata=data.frame(dbh=x)), 100,1000, add=T, col="blue")
curve(predict(mod_list[[3]], newdata=data.frame(dbh=x, wd=0.65)), 100,1000, add=T, col="darkgreen")
points(height~dbh, data=tmpHeight, pch=20)
legend("topleft",col=c("black","blue","darkgreen"),legend=c("linear","log","non-linear"),lwd=c(1,1,1))

try(census$height_m <- predict(mod_list[[best_mod]],
                        newdata=data.frame(wd=census$wd, dbh=census$dbh)), silent=T)
census <- census[census$f1!="0",] # REMOVE THE DEAD TREES
census <- census %>% filter(is.na(dbh)==F)
}
mod_height <- mod_list[[best_mod]]
#---END CENSUS PREPROCESS -------------------------------------------------------------------
#--------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------
#--- FIX PROBLEMS IN CENSUS AND DENDROMETER DATA 

# predict.lm missing heights and dbh
tags_missing_height <- census[is.na(census$height_m)==T,]$tree_tag
tmp_missing_height <- census %>% filter(tree_tag %in% tags_missing_height) 
if(dim(tmp_missing_height)[1] != 0){
for(i in 1:length(tags_missing_height)){
  tmp_fit <- tmp_missing_height %>% select(dbh,wd,census_date,tree_tag,f1,f2) %>% 
  filter(tree_tag==tags_missing_height[i]) %>% 
  lm(dbh~census_date, data=.)

pred_dbh <- predict(tmp_fit, newdata=data.frame(census_date=
tmp_missing_height %>% select(dbh,wd,height_m, census_date,tree_tag,f1,f2) %>% 
  filter(tree_tag==tags_missing_height[i] & is.na(height_m)==T) %>% 
  pull(census_date)))

census[census$tree_tag==tags_missing_height[i] & is.na(census$height_m)==T,]$dbh <- pred_dbh
census[census$tree_tag==tags_missing_height[i] & is.na(census$height_m)==T,]$height_m <- 
  predict(nl_fit, newdata=data.frame(wd=census[census$tree_tag==tags_missing_height[i] & is.na(census$height_m)==T,]$wd,
                                     dbh=census[census$tree_tag==tags_missing_height[i] & is.na(census$height_m)==T,]$dbh))
}
}
census$height_m %>% is.na %>% table
census$dbh %>% is.na %>% table
census <- census %>% filter(dbh>0) #get rid of dead trees 

#--- END FIX PROBLEMS IN CENSUS AND DENDROMETER DATA 

#--- ESTIMATE CANOPY POSITION ------------------------------------------------
vec_census_dates <- census$census_date %>% unique() %>% sort()
tmp <- census %>% filter(near(census_date,vec_census_dates[1])) %>% ppa()
for(i in 2:length(vec_census_dates)){
  tmp_2 <- census %>% filter(near(census_date,vec_census_dates[i])) %>% ppa()
  tmp <- bind_rows(tmp, tmp_2)
}
census <- tmp

#--- END ESTIMATE CANOPY POSITION --------------------------------------------

#--------------------------------------------------------------------------------------------

# Allometric options: [2-Chave2005,dry eqI.3][3-Chave2005,moist eqI.6][4-Chave2014,wet-eq4]
# Height correction options: [1-estimate height in function if there are 50+ trees]
#                            [2-Feldpauch][3-all trees have heights from pre-processing step]

### --- SUBSET TEST -------------------------------------------------------------------------
# vec_tags <- dendrometer$tree_tag %>% unique()
# vec_tags_ss <- vec_tags[sample.int(n = length(vec_tags), 100)]
# dendrometer_ss <- dendrometer %>% filter(tree_tag %in% vec_tags_ss)
# 
# dend_out_ss <- dendro_calc_agC_Mg_dev(census = census,
#                                    mod_height = mod_height,
#                                    dendrometer = dendrometer_ss,
#                                    plot_code = this_plot_code,
#                                    allometric_option="5",
#                                    ret="nppacw.permonth.perha",
#                                    census_year = census_year,
#                                    are_dend_units_correct = are_dend_units_correct)


dend_out <- dendro_calc_agC_Mg_dev(census = census, 
                            mod_height = mod_height,
                            dendrometer = dendrometer, 
                            plot_code = this_plot_code, 
                            allometric_option="5", 
                            ret="nppacw.permonth.perha", 
                            census_year = census_year, 
                            are_dend_units_correct = are_dend_units_correct)
f_out_path <- paste0("outputs/stem_npp_month/stem_npp_",this_plot_code,"_",Sys.Date(),".csv")
write_csv(dend_out, path = f_out_path)


# diagnostic plot
dend_out %>% 
  ggplot(data=., aes(date, npp_u_agC_month))+
  geom_line(color="black")+
  geom_point()+ # raw mean
  geom_point(aes(date, npp_wu_agC_month), lwd=3, col='orange')+ # weighted mean
  geom_line(aes(date, npp_hybrid_bestMod_agC_month), col="blue", lwd=2)+ # hybrid
  geom_line(aes(date, npp_bestEst_agC_month), col="purple", lwd=2, alpha=0.9)+ # best estimate
  geom_line(aes(date, npp_hybrid_gam_agC_month), col="darkgreen")+ # hybrid with GAM
  geom_line(aes(date, npp_lm_agC_month), col="red")+ # don't use , linear model 
  # geom_point(aes(date,dev),col="red",alpha=0.7)+ # deviance explained
  geom_abline(aes(intercept=0,slope=0),lwd=0.5,alpha=0.5)+
  geom_abline(aes(intercept=0.2,slope=0),lwd=0.5,alpha=0.2)

