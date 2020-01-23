library(tidyverse); library(lubridate); library(mgcv); 
source('R/litterfall/proc_litterfall_utils.R')
# --- Import and prep clim ----------------------------------------------------------------------------------
lf_o <- read_csv("data/litterfall/NXV_litterfall.csv") %>% 
  mutate(site=substr(plot_code,1,3))
lf_o <- lf_o %>% mutate(date=ymd(paste(year,month,day),tz='UTC'), 
                        doy = yday(date)) %>% 
  mutate(frac=leavesflf_MgC_ha_day/totalflf)

# Nova-Xavantina -------------------------------------------------------------------
# oddities NXV-02 2013/02 2013/07
dat_nxv <- lf_o %>% filter(site=='NXV') %>%
  filter(!(site=='NXV' & date=='2013-02-01' & plot_code=='NXV-01')) %>%  # Bad colletion
  mutate(leavesflf_MgC_ha_day = 
           ifelse(plot_code=='NXV-02' &  
                    date %in% as.POSIXct(c("2013-07-15","2013-07-02","2013-02-01"),tz='UTC'), 
                  leavesflf_MgC_ha_day*0.1, leavesflf_MgC_ha_day))


vec <- rgamma(10000, shape = 3, rate = 2)
hist(vec)
glm(vec~1, family=Gamma(link='log')) %>% summary
1/0.338
o <- rstanarm::stan_glm(vec~1, family=Gamma(link='log'),algorithm='meanfield')
summary(o, probs = c(0.1,0.5,0.9))


dat_nxv %>% 
  group_by(plot_code, date) %>% 
  summarize(lf_mean = mean(leavesflf_MgC_ha_day, na.rm=T), 
            lf_sd = sd(leavesflf_MgC_ha_day, na.rm=T), 
            nobs = n()) %>% 
  ungroup() %>% 
  mutate(year=year(date), month=month(date)) %>% 
  filter(year >= 2013 & year<= 2017) %>% 
  group_by(plot_code, year) %>% 
  summarize(lf_annual = mean(lf_mean, na.rm=T)*365, 
            lf_sd = sd(lf_mean, na.rm=T), 
            lf_p05 = quantile(lf_mean, 0.05, na.rm=T)*365,
            lf_p95 = quantile(lf_mean, 0.95, na.rm=T)*365,
            nobs = n()) %>% 
  ungroup() %>%   
  ggplot(data=., aes(year, lf_annual,color=plot_code))+
  # geom_ribbon(aes(x=year, 
  #                 ymax=lf_annual+2*lf_sd, 
  #                 ymin=lf_annual-2*lf_sd, 
  #                 fill=plot_code), 
  #             alpha=0.3, lty=0)+
  geom_ribbon(aes(x=year, 
                  ymax=lf_p95, 
                  ymin=lf_p05, 
                  fill=plot_code), 
              alpha=0.3, lty=0)+
  geom_point()+
  geom_line()
  

dat_nxv %>% 
  group_by(plot_code, date) %>% 
  summarize(lf_mean = mean(leavesflf_MgC_ha_day, na.rm=T), 
            lf_sd = sd(leavesflf_MgC_ha_day, na.rm=T), 
            nobs = n()) %>% 
  ungroup() %>% 
  ggplot(data=., aes(date, lf_mean, color=plot_code))+
  geom_ribbon(aes(x=date, ymax=lf_mean+2*lf_sd, ymin=lf_mean-2*lf_sd, fill=plot_code), 
              alpha=0.3, lty=0)+
  geom_line()

dat_nxv %>% names
dat_nxv %>% ggplot(data=., aes(date, leavesflf_MgC_ha_day,color=plot_code))+geom_smooth()

dat_nxv <- est_leaf_frac(dat=dat_nxv, site='NXV', max_val=0.05)
o_nxv01 <- rollsmooth_mean_lf_mo(dat_nxv, gem_plot_code = 'NXV-01')
o_nxv02 <- rollsmooth_mean_lf_mo(dat_nxv, gem_plot_code = 'NXV-02')

lf_rollsmooth <- bind_rows(o_nxv01, o_nxv02)
write_csv(lf_rollsmooth, paste0("data/litterfall/flf_rollsmooth_MgC_ha_mo_",Sys.Date(),'.csv'))

