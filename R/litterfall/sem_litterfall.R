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
  filter(meas_interval_days > 7) %>% # This removes some outliers! I don't know if this is correct, but I suspect not
  group_by(plot_code, date) %>% 
  summarize(tot = sum(totalflf)) %>% 
  ungroup() %>% 
  ggplot(data=., aes(date, tot,color=plot_code))+
  geom_point()



#--- arithmetic way to estimate mean of total fine litterfall across months, then sum to year---
# (1) calc total litterfall across each trap for each year
# (2) calc plot level mean litterfall and SEM
library(lubridate); library(kable)
lf_o %>% 
  filter(site=='NXV' & meas_interval_days>7) %>%  # remove outliers
  filter(date >= ymd(paste(2011,02,01)) &         # subset to a period of 5 full years
           date < ymd(paste(2016,02,01))) %>% 
  mutate(date = ymd(paste(year,month,1))) %>%     # recreate date column
  mutate(fake_date = date-months(1)) %>%      # in order to create year column corresponding to full 12 months, for grouping
  mutate(year=year(fake_date)) %>%  # recreate year col
  group_by(plot_code, litter_trap_num, year) %>% 
  summarize(val = mean(totalflf,na.rm=T)*365) %>% 
  ungroup() %>% 
  group_by(plot_code) %>% 
  summarize(nobs = n(), 
            val_u = mean(val), 
            val_sd = sd(val)) %>%
  ungroup() %>% 
  mutate(val_sem = val_sd/sqrt(nobs)) %>% 
  kable()


#--- alt arithmetic way to estimate overall annual mean of total fine litterfall ---
# (1) calc mean & sd across traps for each year & month
# (2) calc annual values of total littefall and sd across years (~5 years)
# (3) calc standard error across years
lf_o %>% 
  filter(site=='NXV' & meas_interval_days>7) %>% # remove outliers
  filter(date >= ymd(paste(2011,02,01)) & 
           date < ymd(paste(2016,02,01))) %>% 
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
  summarize(nobs=n(), 
            val_u = mean(tot),
            val_sd = sd(tot)) %>% 
  ungroup() %>% 
  mutate(val_sem = val_sd/sqrt(nobs)) %>% 
  kable()


# Distribution estimation approach -----------------------------------------
vec_lf_nxv1 <- lf_o %>% 
  filter(plot_code=='NXV-01') %>% 
  filter(site=='NXV' & meas_interval_days>7) %>%  # remove outliers
  filter(date >= ymd(paste(2011,02,01)) &         # subset to a period of 5 full years
           date < ymd(paste(2016,02,01))) %>% 
  mutate(date = ymd(paste(year,month,1))) %>%     # recreate date column
  mutate(fake_date = date-months(1)) %>%      # in order to create year column corresponding to full 12 months, for grouping
  mutate(year=year(fake_date)) %>%  # recreate year col
  group_by(plot_code, litter_trap_num, year) %>% 
  summarize(val = mean(totalflf,na.rm=T)*365) %>% 
  ungroup() %>% 
  pull(val)

vec_lf_nxv2 <- lf_o %>% 
  filter(plot_code=='NXV-02') %>% 
  filter(site=='NXV' & meas_interval_days>7) %>%  # remove outliers
  filter(date >= ymd(paste(2011,02,01)) &         # subset to a period of 5 full years
           date < ymd(paste(2016,02,01))) %>% 
  mutate(date = ymd(paste(year,month,1))) %>%     # recreate date column
  mutate(fake_date = date-months(1)) %>%      # in order to create year column corresponding to full 12 months, for grouping
  mutate(year=year(fake_date)) %>%  # recreate year col
  group_by(plot_code, litter_trap_num, year) %>% 
  summarize(val = mean(totalflf,na.rm=T)*365) %>% 
  ungroup() %>% 
  pull(val)

# compare distributions
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

hist(vec_lf_nxv2)
fitdistrplus::descdist(vec_lf_nxv2)
f2_g <- fitdistrplus::fitdist(vec_lf_nxv2, distr='gamma',method='mle')
f2_ln <- fitdistrplus::fitdist(vec_lf_nxv2, distr='lnorm',method='mle')
f2_n <- fitdistrplus::fitdist(vec_lf_nxv2, distr='norm',method='mle')
f2_w <- fitdistrplus::fitdist(vec_lf_nxv2, distr='weibull',method='mle')
f2_e <- fitdistrplus::fitdist(vec_lf_nxv2, distr='exp',method='mle')
f2_g$aic
f2_ln$aic
f2_n$aic
f2_w$aic
f2_e$aic


# Log-normal is best fit for both NXV-01 and NXV-02
b1 <- bootdist(f_ln, niter = 1000)
summary(b1)
o1 <- quantile(b1, probs = 0.5, CI.type = 'two.sided')
o1$quantiles %>% as.numeric
o1$quantCI

b2 <- bootdist(f2_ln, niter = 1000)
summary(b2)
o2 <- quantile(b2, probs = 0.5, CI.type = 'two.sided')
o2$quantiles %>% as.numeric
o2$quantCI$`p=0.5`

tibble(plot_code = c("NXV-01", "NXV-02"), 
       median_est = c(as.numeric(o1$quantiles), 
                      as.numeric(o2$quantiles)), 
       CI_2.5 = c(o1$quantCI$`p=0.5`[1] ,o2$quantCI$`p=0.5`[1] ), 
       CI_97.5 = c(o1$quantCI$`p=0.5`[2] ,o2$quantCI$`p=0.5`[2] )
) %>% 
  kable()


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








# Distribution approach ---------------------------------------------------

fitdistrplus::fitdist()

t_lf_distribution <- lf_o %>% 
  filter(site=='NXV' & meas_interval_days>7) %>% 
  filter(date >= ymd(paste(2011,02,01)) & 
           date < ymd(paste(2016,02,01))) %>% 
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
