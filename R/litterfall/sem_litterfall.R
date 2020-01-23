library(tidyverse); library(lubridate); library(mgcv); 
source('R/litterfall/proc_litterfall_utils.R')
# --- Import and prep clim ----------------------------------------------------------------------------------
lf_o <- read_csv("data/litterfall/NXV_litterfall.csv") %>% 
  mutate(site=substr(plot_code,1,3))
lf_o <- lf_o %>% mutate(date=ymd(paste(year,month,day),tz='UTC'), 
                        doy = yday(date)) %>% 
  mutate(frac=leavesflf_MgC_ha_day/totalflf)

# Nova-Xavantina -------------------------------------------------------------------
lf_o %>% 
  filter(site=='NXV') %>% 
  filter(meas_interval_days > 7) %>% 
  group_by(plot_code, date) %>% 
  summarize(tot = sum(totalflf)) %>% 
  ungroup() %>% 
  ggplot(data=., aes(date, tot,color=plot_code))+
  geom_point()


# distribution approach
vec_lf_nxv1 <- lf_o %>% filter(plot_code=='NXV-01' & 
                                 meas_interval_days>7) %>% 
  mutate(totalflf = if_else(totalflf==0 ,0.001, totalflf)) %>% 
  filter(date >= ymd(paste(2011,02,01)) & 
           date < ymd(paste(2016,02,01))) %>% 
  mutate(flf_MgC_ha_mo = totalflf*30.4) %>% 
  pull(flf_MgC_ha_mo)
hist(vec_lf_nxv1)
fitdistrplus::descdist(vec_lf_nxv1)
f_g <- fitdistrplus::fitdist(vec_lf_nxv1, distr='gamma',method='mle')
f_ln <- fitdistrplus::fitdist(vec_lf_nxv1, distr='lnorm',method='mle')
f_n <- fitdistrplus::fitdist(vec_lf_nxv1, distr='norm',method='mle')
f_w <- fitdistrplus::fitdist(vec_lf_nxv1, distr='weibull',method='mle')
f_e <- fitdistrplus::fitdist(vec_lf_nxv1, distr='exp',method='mle')
f_g$aic
f_ln$aic
f_n$aic
f_w$aic
f_e$aic

summary(f_w)
plot(f_w)

summary(vec_lf_nxv1)
mean(vec_lf_nxv1)*12
f_w$estimate[2]*gamma((1+1/f_w$estimate[1]))
b_w <- fitdistrplus::bootdist(f_w, niter=1000)
summary(b_w)
quantile(f_w, probs = 0.5)
quantile(b_w, probs = 0.5)

lf_o %>% 
  filter(site=='NXV') %>% 
  filter(meas_interval_days > 7) %>% 
  mutate(days_per_month = days_in_month(date)) %>% 
  group_by(plot_code, year,month,date) %>% 
  summarize(flf = mean(totalflf*days_per_month), 
            flf_sd = sd(totalflf), 
            nobs = n()) %>% 
  ungroup() %>% 
  mutate(flf_sem = flf_sd/sqrt(nobs)) %>% 
  group_by(plot_code, year, month) %>% 
  summarize(flf_MgC_ha_mo = mean(flf, na.rm=T)) %>% 
  ungroup() %>% 
  mutate(date=ymd(paste(year,month,1),tz='UTC')) %>% 
  ggplot(data=., aes(date, flf_MgC_ha_mo,color=plot_code))+
  geom_point()

lf_nxv <- lf_o %>% 
  filter(site=='NXV') %>% 
  filter(meas_interval_days > 7) %>% 
  mutate(days_per_month = days_in_month(date)) %>% 
  group_by(plot_code, year,month,date) %>% 
  summarize(flf = mean(totalflf*days_per_month)) %>% 
  ungroup()
lf_nxv %>% group_by(plot_code) %>% summarize(val = mean(flf)*12)
lf_nxv %>% 
  group_by(plot_code, year, month) %>% 
  summarize(flf_MgC_ha_mo = mean(flf, na.rm=T)) %>% 
  ungroup() %>% 
  mutate(date=ymd(paste(year,month,1),tz='UTC')) %>% 
  ggplot(data=., aes(date, flf_MgC_ha_mo,color=plot_code))+
  geom_point()


dat_nxv %>% 
  group_by(plot_code, year, month) %>% 
  summarize(val = mean(leavesflf_MgC_ha_day*days_in_month(date), na.rm=T), 
            val_sd = sd(leavesflf_MgC_ha_day*days_in_month(date), na.rm=T), 
            nobs = n(), 
            val_sem = val_sd/sqrt(nobs)) %>% 
  ungroup() %>% 
  group_by(year) %>% 
  summarize(tot = sum(val, na.rm=T), 
            tot_sd = sqrt(sum(val_sd**2)),
            nobs = n()) %>% 
  ungroup() %>% 
  mutate(tot_sem = tot_sd/sqrt(nobs)) %>% 
  ggplot(data=., aes(year, tot))+
  geom_ribbon(aes(year, ymax=tot+1*tot_sem, ymin=tot-1*tot_sem),alpha=0.25,lty=0)+
  geom_point()


dat_nxv %>% 
  group_by(plot_code, year, month) %>% 
  summarize(val = mean(leavesflf_MgC_ha_day*days_in_month(date), na.rm=T), 
            val_sd = sd(leavesflf_MgC_ha_day*days_in_month(date), na.rm=T)) %>% 
  ungroup() %>% 
  group_by(plot_code, year) %>% 
  summarize(tot = sum(val, na.rm=T), 
            tot_sd = sqrt(sum(val_sd**2)),
            nobs = n()) %>% 
  ungroup() %>% 
  mutate(tot_sem = tot_sd/sqrt(nobs)) %>% 
  group_by(plot_code) %>% 
  summarize(u = mean(tot),
            u_sd = sd(tot),
            nobs=n(),
            u_sem = u_sd/sqrt(nobs)) %>% 
  ungroup()


# simplest way to estimate overall annual mean of total fine litterfall
lf_o %>% 
  filter(site=='NXV' & meas_interval_days>7) %>% 
  group_by(plot_code, year, month) %>% 
  summarize(val = mean(totalflf*days_in_month(date), na.rm=T), 
            val_sd = sd(totalflf*days_in_month(date), na.rm=T)) %>% 
  ungroup() %>% 
  group_by(plot_code, year) %>% 
  summarize(tot = sum(val, na.rm=T), 
            tot_sd = sqrt(sum(val_sd**2)),
            nobs = n()) %>% 
  ungroup() %>% 
  mutate(tot_sem = tot_sd/sqrt(nobs)) %>% 
  group_by(plot_code) %>% 
  summarize(u = mean(tot),
            u_sd = sd(tot),
            nobs=n(),
            u_sem = u_sd/sqrt(nobs)) %>% 
  ungroup()
