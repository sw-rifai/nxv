# Use the mean or the median of litterfall traps? ------------------------
# I guess use the mean
lf_o %>% 
  filter(plot_code=='NXV-01') %>% 
  ggplot(data=., aes(leavesflf_MgC_ha_day))+
  geom_histogram()+
  geom_vline(aes(xintercept=median(leavesflf_MgC_ha_day)), col='red')+
  geom_vline(aes(xintercept=mean(leavesflf_MgC_ha_day)), col='blue')+
  facet_wrap(~month)


dat <- dat_nxv;  gem_plot_code <- 'NXV-01'; min_nobs <- 15

rollsmooth_mean_lf_mo <- function(dat, gem_plot_code, min_nobs=15){
  library(dplyr); library(lubridate); library(RcppRoll)
  junk <- dat %>%
    filter(plot_code==gem_plot_code) %>%
    group_by(site, plot_code, date) %>% 
    summarize(nobs=n(), 
              tot=mean(est_leaf_lf, na.rm=T)) %>% 
    ungroup() %>% 
    filter(nobs>=min_nobs) 
  
  tmp_df <- dat %>%
    filter(plot_code==gem_plot_code) %>% 
    group_by(site,plot_code, date) %>%
    summarize(u=mean(est_leaf_lf, na.rm=T),
              lf_sd = sd(est_leaf_lf, na.rm=T),
              nobs=n()) %>% ungroup() %>%
    filter(nobs>=min_nobs) %>%
    right_join(.,tibble(date=seq(min(dat$date,na.rm=T) %m-% months(1), max(dat$date, na.rm=T) %m+% months(1), by='1 day')), by='date') %>%
    arrange(date) %>%
    fill(u) %>%
    fill(lf_sd) %>% 
    mutate(u_med = roll_mean(u, n=31, align='center', fill = NA)) %>%
    mutate(u_sum = roll_sum(u_med, n=31, align = 'left', fill=NA)) %>%
    mutate(u_sd = roll_mean(lf_sd, n=31, align = 'left', fill=NA)) %>% 
    mutate(month=month(date), year=year(date), doy=yday(date)) %>%
    left_join(., junk, by=c('site','plot_code','date')) %>%
    ungroup()
  
  out_df <- tmp_df %>% 
    group_by(site, plot_code, year,month) %>%
    summarize(n_obs = sum(is.na(tot)==F),
              lf_mo = mean(u_sum, na.rm=T),
              lf_sd_mo = mean(u_sd, na.rm=T), 
              lf_raw_u_mo = mean(tot, na.rm=T)*30.4) %>%
    ungroup() %>%
    filter(n_obs>=1) %>%
    mutate(date=ymd(paste(year,month,15),tz='UTC'))
  return(out_df)    
}

out_df %>% 
  ggplot(data=., aes(date, lf_mo))+
  geom_line()+
  geom_point(aes(date, lf_sd_mo))


dat_nxv %>% names
dat_nxv %>% 
  ggplot(data=., aes(pred_leaf_lf, leavesflf_MgC_ha_day))+
  geom_point()+
  geom_smooth(method='lm')

dat_nxv %>% dim
dat_nxv$leavesflf_MgC_ha_day %>% is.na %>% table
