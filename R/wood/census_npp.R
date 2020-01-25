library(tidyverse); library(lubridate)
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("collapse", "dplyr")
conflicted::conflict_prefer("lag", "dplyr")

source("R/wood/DendroFunctions.R")

# Chave et al, 2014:
Chave2014_eq4 <- function(diam_cm, density, height_m, wd) {
  density=wd; 
  AGB_est <- 0.0673*(density*((diam_cm)^2)*height_m)^0.976 
  return(AGB_est)
}

nxv_cens <- read_csv("data/wood/census_nxv.csv")
names(nxv_cens) <- gsub(" ","_",names(nxv_cens))

nxv_cens <- nxv_cens %>% 
  rename(tree_tag = tag_number)

nxv_cens <- nxv_cens %>% 
  # select(plot_code, census_date, tag, treeid, d0,d1,d2,d3,d4) %>% 
  distinct() %>% 
  mutate(d0 = if_else(d0>=13000, d0*0.1, d0))

nxv_cens <- nxv_cens %>% 
  filter(str_detect(species, 'Morto')==F) %>%  # Filter out dead trees
  filter(d0 != 0) %>% # filter out dead trees
  # filter(is.na(genus)==T) %>% 
  ungroup() %>% 
  rowwise() %>% 
  mutate(genus = if_else(is.na(genus)==T, 
                        str_split(species, pattern = " ",simplify = T)[1], 
                        genus)) %>% 
  ungroup() 

vec_nxv1_subplots <- nxv_cens %>% filter(plot_code=='NXV-01') %>% 
  filter(census_date==min(census_date)) %>% 
  pull(sub_plot_t1) %>% 
  unique() %>% 
  sort()

# plot basal area
nxv_cens %>% 
  mutate(d0 = if_else(d0>=13000, d0*0.1, d0)) %>% 
  filter(is.na(d0)==F) %>% 
  select(plot_code, census_date, d0, tree_tag) %>% 
  distinct() %>% 
  group_by(plot_code, census_date) %>% 
  summarize(basal_area_m2 = sum(3.141593*(d0/2)**2)/(1000**2)) %>% 
  ungroup() %>% 
  ggplot(data=., aes(census_date, basal_area_m2,color=plot_code))+
  geom_point()+
  geom_line()




# NXV-01 ------------------------------------------------------------------
# live trees only 
nxv1_cens <- nxv_cens %>% 
  filter(plot_code=='NXV-01') %>% 
  filter(d0 < 1300) %>% # Was there really a tree with 1.3 m in diameter at NXV-01? I am counting this as an outlier
  filter(d0>=100)


## Extrapolate Wood Density ***********************************************
# Get wood wd from Zanne global wood wd database
nxv1_cens$wd <- NA
nxv1_cens_sp_list <- nxv1_cens %>% distinct(family, genus, species, continent, tree_tag)

# recast nxv1_cens species list
nxv1_cens_sp_list <- nxv1_cens %>% distinct(family, genus, species, continent)
nxv1_cens_sp_list$wd <- NA
nxv1_cens_sp_list <- find_wd(census_sp_list = nxv1_cens_sp_list)
nxv1_cens_sp_list <- nxv1_cens_sp_list %>% mutate(species=paste(genus,species))

nxv1_cens <- nxv1_cens %>% select(-wd)
nxv1_cens_sp_list <- nxv1_cens_sp_list %>% 
  select(species, wd) %>% 
  group_by(species) %>% 
  summarize(wd=mean(wd,na.rm=T)) %>% 
  ungroup()
nxv1_cens <- left_join(nxv1_cens,nxv1_cens_sp_list, by="species")


# Extrapolate height *******************************************************
log_fit <- lm(height~log(d0), data=nxv1_cens %>% filter(is.na(wd)==F));
lm_fit <- lm(height~d0, data=nxv1_cens %>% filter(is.na(wd)==F))
nl_fit <- nls(height~ b0*((1/wd)^b1)*d0^b2, data=nxv1_cens %>% filter(is.na(wd)==F),
              start = list(b0=2.643, b1=-0.3, b2=0.5))
mod_list <- list(log_fit, lm_fit, nl_fit)
best_mod <- AIC(log_fit, lm_fit, nl_fit)[,2] %>% which.min()
best_mod <- mod_list[[best_mod]]
nxv1_cens <- nxv1_cens %>% 
  mutate(pred_height = predict(best_mod, newdata=.)) %>% 
  mutate(height = if_else(is.na(height)==T, pred_height, height)) 
  

# Estimate Biomass ********************************************************
nxv1_biomass <- nxv1_cens %>% 
  filter(is.na(d0)==F) %>% 
  filter(sub_plot_t1 %in% vec_nxv1_subplots) %>% # first 50 subplots [total is 0.5 ha]
  mutate(biomass = Chave2014_eq4(diam_cm = d0/10, density=wd, height_m = height, wd=wd)) %>% 
  group_by(plot_code, census_date) %>% 
  summarize(tot_biomass_Mg_ha = (sum(biomass)*2)/1000) %>% # multiplied by 2 because plot is 0.5 ha
  ungroup()
nxv1_biomass %>% 
  ggplot(data=., aes(census_date, tot_biomass_Mg_ha,color=plot_code))+
  geom_point()+
  geom_line()


# Estimate NXV-01 NPP -----------------------------------------------------
nxv1_npp <- nxv1_biomass %>% 
  arrange(census_date) %>% 
  mutate(year_diff = census_date - lag(census_date, n=1, order_by = census_date), 
         biomass_diff = tot_biomass_Mg_ha - lag(tot_biomass_Mg_ha, n = 1, order_by = census_date)) %>% 
  mutate(npp_Mg_ha_yr = biomass_diff/year_diff)


# Estimate SEM of NXV-01 NPP ----------------------------------------------
nxv1_sem <- nxv1_npp %>% 
  filter(is.na(npp_Mg_ha_yr)==F) %>% 
  group_by(plot_code) %>% 
  summarize(npp_u = mean(npp_Mg_ha_yr),
            npp_sd = sd(npp_Mg_ha_yr), 
            nobs = n(), 
            npp_sem = npp_sd/sqrt(nobs))
nxv1_sem %>% 
  kable()


nxv1_npp %>% 
  ggplot(data=., aes(census_date, npp_Mg_ha_yr))+
  geom_point()+
  geom_line()+
  labs(y='NPP Mg Biomass ha-1 yr-1')+
  scale_y_continuous(limits = c(0, 5))











# NXV-02 ------------------------------------------------------------------
# live trees only 
nxv2_cens <- nxv_cens %>% 
  filter(plot_code=='NXV-02') %>% 
  filter(d0 < 1300) %>% 
  filter(d0>=100)


## Extrapolate Wood Density ***********************************************
# Get wood wd from Zanne global wood wd database
nxv2_cens$wd <- NA
nxv2_cens_sp_list <- nxv2_cens %>% distinct(family, genus, species, continent, tree_tag)

# recast nxv2_cens species list

nxv2_cens_sp_list <- nxv2_cens %>% 
  distinct(family, genus, species, continent)
nxv2_cens_sp_list$wd <- NA
nxv2_cens_sp_list <- find_wd(census_sp_list = nxv2_cens_sp_list)
nxv2_cens_sp_list <- nxv2_cens_sp_list %>% mutate(species=paste(genus,species))

nxv2_cens <- nxv2_cens %>% select(-wd)
nxv2_cens_sp_list <- nxv2_cens_sp_list %>% 
  select(species, wd) %>% 
  group_by(species) %>% 
  summarize(wd=mean(wd,na.rm=T)) %>% 
  ungroup()
nxv2_cens <- left_join(nxv2_cens,nxv2_cens_sp_list, by="species")

# Extrapolate height *******************************************************
log_fit <- lm(height~log(d0), data=nxv2_cens %>% filter(is.na(wd)==F));
lm_fit <- lm(height~d0, data=nxv2_cens %>% filter(is.na(wd)==F))
nl_fit <- nls(height~ b0*((1/wd)^b1)*d0^b2, data=nxv2_cens %>% filter(is.na(wd)==F),
              start = list(b0=2.643, b1=-0.3, b2=0.5))
mod_list <- list(log_fit, lm_fit, nl_fit)
best_mod <- AIC(log_fit, lm_fit, nl_fit)[,2] %>% which.min()
best_mod <- mod_list[[best_mod]]
nxv2_cens <- nxv2_cens %>% 
  mutate(pred_height = predict(best_mod, newdata=.)) %>% 
  mutate(height = if_else(is.na(height)==T, pred_height, height)) 


# Estimate Biomass ********************************************************
nxv2_biomass <- nxv2_cens %>% 
  filter(is.na(d0)==F) %>% 
  # filter(sub_plot_t1 %in% vec_nxv2_subplots) %>% # first 50 subplots [total is 0.5 ha]
  mutate(biomass = Chave2014_eq4(diam_cm = d0/10, density=wd, height_m = height, wd=wd)) %>% 
  group_by(plot_code, census_date) %>% 
  summarize(tot_biomass_Mg_ha = (sum(biomass))/1000) %>% # multiplied by 2 because plot is 0.5 ha
  ungroup()
nxv2_cens %>% group_by(census_date) %>% summarize(nobs = n()) %>% 
  ggplot(data=.,aes(census_date, nobs))+geom_line()+geom_point()
nxv2_biomass %>% 
  ggplot(data=., aes(census_date, tot_biomass_Mg_ha,color=plot_code))+
  geom_point()+
  geom_line()


# Estimate NXV-02 NPP -----------------------------------------------------
nxv2_npp <- nxv2_biomass %>% 
  arrange(census_date) %>% 
  mutate(year_diff = census_date - lag(census_date, n=1, order_by = census_date), 
         biomass_diff = tot_biomass_Mg_ha - lag(tot_biomass_Mg_ha, n = 1, order_by = census_date)) %>% 
  mutate(npp_Mg_ha_yr = biomass_diff/year_diff)


# Estimate SEM of NXV-02 NPP ----------------------------------------------
nxv2_sem <- nxv2_npp %>% 
  filter(is.na(npp_Mg_ha_yr)==F) %>% 
  group_by(plot_code) %>% 
  summarize(npp_u = mean(npp_Mg_ha_yr),
            npp_sd = sd(npp_Mg_ha_yr), 
            nobs = n(), 
            npp_sem = npp_sd/sqrt(nobs)); 
nxv2_sem

nxv2_npp %>% 
  ggplot(data=., aes(census_date, npp_Mg_ha_yr))+
  geom_point()+
  geom_line()+
  labs(y='NPP Mg Biomass ha-1 yr-1')+
  scale_y_continuous(limits = c(0, 5))



# Print table of Census Woody NPP -----------------------------------------
bind_rows(nxv1_sem, nxv2_sem) %>% 
  kable()









# # Scratch -----------------------------------------------------------------
# 
#   
# # NXV-01 Dead trees ------------------------------------------------------
# nxv_cens %>%
#   filter(plot_code=='NXV-01') %>% 
#   filter(is.na(d0)==F) %>% 
#   filter(sub_plot_t1 %in% vec_nxv1_subplots) %>% # first 50 subplots [total is 0.5 ha]
#   filter(d1 ==0) %>% # means tree is dead according to forest plots manual 
#   group_by(census_date) %>% 
#   summarize(n_dead = sum(d1 == 0)) %>% 
#   ungroup() %>% 
#   ggplot(data=., aes(census_date, n_dead))+
#   geom_point()
#   # mutate(biomass = Chave2014_eq4(diam_cm = d0/10, density=wd, height_m = height, wd=wd)) %>% 
#   # group_by(plot_code, census_date) %>% 
#   # summarize(tot_biomass_Mg_ha = (sum(biomass)*2)/1000) %>% # multiplied by 2 because plot is 0.5 ha
#   # ungroup()
# 
# 
# 
# # Plot basal area  ********************************************************
# nxv1_cens %>% 
#   filter(is.na(d0)==F) %>% 
#   filter(sub_plot_t1 %in% vec_nxv1_subplots) %>% 
#   # select(plot_code, census_date, d0, tag) %>% 
#   # distinct() %>% 
#   group_by(plot_code, census_date) %>% 
#   summarize(basal_area_m2 = sum(3.141593*(d0/2)**2)/(1000**2)) %>% 
#   ungroup() %>% 
#   ggplot(data=., aes(census_date, basal_area_m2,color=plot_code))+
#   geom_point()+
#   geom_line()
# 
# 
# 
# # are the dendrometer bands in the plot? 
# dend_all %>% 
#   filter(plot_code=='NXV-01') %>% 
#   pull(tree_tag) %>% unique %in% 
#   unique(nxv1_cens$tree_tag)
# 
# 
# nxv_cens %>% 
#   filter(is.na(d0)==F) %>% 
#   filter(sub_plot_t1 %in% vec_nxv1_subplots) %>% 
#   # select(plot_code, census_date, d0, tag) %>% 
#   # distinct() %>% 
#   group_by(plot_code, census_date) %>% 
#   summarize(basal_area_m2 = sum(3.141593*(d0/2)**2)/(1000**2)) %>% 
#   ungroup() %>% 
#   ggplot(data=., aes(census_date, basal_area_m2,color=plot_code))+
#   geom_point()+
#   geom_line()
# 
# 
# 
# # number of unique species per census date
# nxv1_cens %>% 
#   group_by(census_date) %>% 
#   summarize(n_spp = length(unique(species))) %>% 
#   ungroup() %>% 
#   ggplot(data=., aes(census_date, n_spp))+
#   geom_point()
# 
# 
# 
# curve(Chave2014_eq4(x, 0.65, height_m = 10, wd = 0.65), 0, 10)
# 
# # plot basal area
# nxv_cens %>% 
#   filter(is.na(d0)==F) %>% 
#   filter(sub_plot_t1 %in% vec_nxv1_subplots) %>% 
#   # select(plot_code, census_date, d0, tag) %>% 
#   # distinct() %>% 
#   group_by(plot_code, census_date) %>% 
#   summarize(basal_area_m2 = sum(3.141593*(d0/2)**2)/(1000**2)) %>% 
#   ungroup() %>% 
#   ggplot(data=., aes(census_date, basal_area_m2,color=plot_code))+
#   geom_point()+
#   geom_line()
# 
# 
# nxv_cens %>% 
#   select(d0,d1,d2,d3,d4) %>% arrange(desc(d0))
#   
# nxv_cens %>% 
#   group_by(plot_code) %>% 
#   select(tag) %>% 
#   distinct()
# 
# 
# nxv_cens$d0 %>% max(na.rm=T)
# 
# 
# nxv_cens %>% 
#   filter(is.na(d0)==F) %>% 
#   filter(d0 > 0) %>% 
#   mutate(d0 = if_else(d0>=13000, d0*0.1, d0)) %>% 
#   group_by(plot_code, census_date) %>% 
#   summarize(val = min(d0, na.rm=T), 
#             nobs =n())
# 
# 
# # plot number of observations (live trees)
# nxv_cens %>% 
#   mutate(d0 = if_else(d0>=13000, d0*0.1, d0)) %>% 
#   filter(is.na(d0)==F) %>% 
#   filter(d0 > 0) %>% 
#   select(plot_code, census_date, d0, treeid) %>% 
#   distinct() %>% 
#   group_by(plot_code, census_date) %>% 
#   summarize(nobs = n()) %>% 
#   ungroup() %>% 
#   ggplot(data=., aes(census_date, nobs,color=plot_code))+
#   geom_point()+
#   geom_line()
# 
# nxv_cens %>% 
#   select(plot_code, d0, d1, d2, d3, d4) %>% 
#   apply(., 2, FUN=function(x) sum(is.na(x)))
# 
# nxv_cens %>% filter(is.na(d0)==T) %>% 
#   select(plot_code, d0, d1, d2, d3, d4) %>% 
# apply(., 2, FUN=function(x) sum(is.na(x)))
# 
