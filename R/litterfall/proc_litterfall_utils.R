library(tidyverse); library(lubridate); library(mgcv)

# Processing functions specific for litterfall in MgC/day -----------------
est_leaf_frac <- function(dat, sites, start_date, end_date,max_val=0.07, plot_fig=F){
  tmp_dat <- dat %>% dplyr::filter(site %in% sites) 
  if(missing(start_date)){start_date <- min(tmp_dat$date, na.rm=T)};
  if(missing(end_date)){end_date <- max(tmp_dat$date, na.rm=T)};
  
  tmp_dat %>% 
    dplyr::filter(date >= start_date & date < end_date) %>%
    dplyr::mutate(totalflf = ifelse(totalflf>max_val, totalflf*0.1, totalflf)) %>% 
    dplyr::mutate(leavesflf_MgC_ha_day = 
                  ifelse(leavesflf_MgC_ha_day>(max_val*0.9), leavesflf_MgC_ha_day*0.1, leavesflf_MgC_ha_day)) %>% 
    dplyr::filter(totalflf <= max_val) %>% 
    dplyr::filter(twigsflf <= 0.1) 
  tmp_fit <- mgcv::gam(frac~te(doy, totalflf), 
                       data=tmp_dat, 
                       method='REML',
                       family=gaussian(link='identity'),
  )
  if(plot_fig==T){
    visreg::visreg(tmp_fit, 'doy', by='totalflf')
    }
  # visreg::visreg(tmp_fit, 'totalflf', by='doy')
  tmp_df <- dat %>% 
    dplyr::filter(site %in% sites) %>%
    # dplyr::filter(totalflf <= max_val) %>% 
    dplyr::mutate(totalflf = ifelse(totalflf>max_val, totalflf*0.1, totalflf)) %>% 
    dplyr::mutate(leavesflf_MgC_ha_day = 
                  ifelse(leavesflf_MgC_ha_day>(max_val*0.9), 
                         leavesflf_MgC_ha_day*0.1, leavesflf_MgC_ha_day)) %>% 
    dplyr::filter(is.na(twigsflf)==T || twigsflf<0.01) %>% 
    dplyr::mutate(pred_frac = predict(tmp_fit, newdata=., type='response')) %>% 
    dplyr::mutate(pred_frac = ifelse(pred_frac<0,0,pred_frac), 
                  pred_frac = ifelse(pred_frac>1,1,pred_frac)) %>% 
    dplyr::mutate(pred_leaf_lf = pred_frac*totalflf) %>% 
    dplyr::mutate(pred_leaf_lf = ifelse(pred_leaf_lf<0,0,pred_leaf_lf)) %>% 
    dplyr::mutate(pred_leaf_lf = ifelse(pred_leaf_lf <=0, 0.01, pred_leaf_lf), 
                  plot_code=as.factor(plot_code), 
                  date_num =as.numeric(date)) %>% 
    dplyr::mutate(est_leaf_lf = ifelse(is.na(leavesflf_MgC_ha_day)==T, pred_leaf_lf, leavesflf_MgC_ha_day))
  return(tmp_df)
}

est_leaf_frac_simple <- function(dat, sites, start_date, end_date,max_val=0.07){
  tmp_dat <- dat %>% dplyr::filter(site %in% sites) 
  if(missing(start_date)){start_date <- min(tmp_dat$date, na.rm=T)};
  if(missing(end_date)){end_date <- max(tmp_dat$date, na.rm=T)};
  
  tmp_dat %>% 
    dplyr::filter(date >= start_date & date < end_date) %>%
    dplyr::filter(totalflf <= max_val) %>% 
    dplyr::filter(twigsflf <= 0.1) 
  tmp_dat %>% mutate(frac =)
  tmp_fit <- mgcv::gam(frac~s(doy), 
                       data=tmp_dat, 
                       family=gaussian(link='identity')
  )
  visreg::visreg(tmp_fit, 'doy', by='totalflf')
  # visreg::visreg(tmp_fit, 'totalflf', by='doy')
  tmp_df <- dat %>% 
    dplyr::filter(site %in% sites) %>%
    filter(is.na(twigsflf)==T || twigsflf<0.01) %>% 
    dplyr::mutate(pred_frac = predict(tmp_fit, newdata=., type='response')) %>% 
    dplyr::mutate(pred_frac = ifelse(pred_frac<0,0,pred_frac), 
                  pred_frac = ifelse(pred_frac>1,1,pred_frac)) %>% 
    dplyr::mutate(pred_leaf_lf = pred_frac*totalflf) %>% 
    dplyr::mutate(pred_leaf_lf = ifelse(pred_leaf_lf<0,0,pred_leaf_lf)) %>% 
    dplyr::mutate(pred_leaf_lf = ifelse(pred_leaf_lf <=0, 0.01, pred_leaf_lf), 
                  plot_code=as.factor(plot_code), 
                  date_num =as.numeric(date)) %>% 
    dplyr::mutate(est_leaf_lf = ifelse(is.na(leavesflf_MgC_ha_day)==T, pred_leaf_lf, leavesflf_MgC_ha_day))
  return(tmp_df)
}

gamsmooth_lf_mo <- function(dat, max_val){
  tmp_dat <- dat %>% group_by(plot_code,date) %>% 
    summarize(u=median(est_leaf_lf,na.rm=T),
              meas_interval_days = median(meas_interval_days, na.rm=T),
              nobs_date=n()) %>% 
    ungroup() %>% 
    filter(nobs_date >= 15) %>% 
    mutate(doy=lubridate::yday(date), 
           month=lubridate::month(date),
           year=lubridate::year(date),
           date_num=as.numeric(date), plot_code=as.factor(plot_code))
  nobs_df <- dat %>% group_by(year,month,plot_code) %>%
    summarize(nobs_traps=n(),
              lf_obs_med=median(est_leaf_lf,na.rm=T)*30.4) %>%
    ungroup() %>%
    mutate(date=lubridate::ymd(paste(year,month,15),tz='UTC'))
  tmp_mod <- mgcv::gam(u~ 
                         s(meas_interval_days)+
                         s(date_num,by=plot_code)+s(month,by=plot_code),     
                       # te(scale(year),scale(month),scale(doy), by=plot_code),
                       # s(month,bs='cc')+ # doy because collection date is often at beginning/end of month
                       # s(scale(date_num), by=plot_code, bs='ad'), # ad & GCV.Cp overfit smooth, which is ideal for this purpose
                       # method='GCV.Cp',
                       data=tmp_dat %>% mutate(u=u+0.001), 
                       gamma = 0.1,
                       family=Gamma(link='log'))
  summary(tmp_mod) %>% print()
  # plot(tmp_mod)
  # gam.check(tmp_mod);
  # abline(0,1,col='red')
  tmp_dat <- tmp_dat %>%
    expand(plot_code,
           meas_interval_days=27,
           date=seq(min(tmp_dat$date)-months(1), max(tmp_dat$date)+months(1), by='1 day')) %>% 
    mutate(date_num = as.numeric(date), 
           month=lubridate::month(date), 
           doy=lubridate::yday(date)) %>% 
    mutate(pred_lf=predict(tmp_mod,newdata=.,type='response')-0.001, 
           pred_lf_se=predict(tmp_mod,newdata=.,type='response',se.fit = T)$se.fit) 
  tmp_dat <- tmp_dat %>% mutate(year=lubridate::year(date), month=lubridate::month(date)) %>% 
    group_by(plot_code,year,month) %>% 
    summarize(pred_lf = sum(pred_lf, na.rm=T), 
              pred_lf_se = sum(pred_lf_se, na.rm=T)) %>% 
    ungroup() %>% 
    mutate(date=lubridate::ymd(paste(year,month,15),tz='UTC'))
  tmp_dat <- left_join(tmp_dat, nobs_df, by=c('date','year','month','plot_code'))
  return(tmp_dat)
}

rollsmooth_lf_mo <- function(dat, gem_plot_code, min_nobs=15){
  library(dplyr); library(lubridate); library(RcppRoll)
  junk <- dat %>%
    filter(plot_code==gem_plot_code) %>%
    group_by(site, plot_code, date) %>% 
    summarize(nobs=n(), 
              tot=median(est_leaf_lf, na.rm=T)) %>% 
    ungroup() %>% 
    filter(nobs>=min_nobs) 
  
  tmp_df <- dat %>%
    filter(plot_code==gem_plot_code) %>% 
    group_by(site,plot_code, date) %>%
    summarize(u=mean(est_leaf_lf, na.rm=T),
              nobs=n()) %>% ungroup() %>%
    filter(nobs>=min_nobs) %>%
    right_join(.,tibble(date=seq(min(dat$date,na.rm=T) %m-% months(1), max(dat$date, na.rm=T) %m+% months(1), by='1 day')), by='date') %>%
    arrange(date) %>%
    fill(u) %>%
    mutate(u_med = roll_mean(u, n=31, align='center', fill = NA)) %>%
    mutate(u_sum = roll_sum(u_med, n=31, align = 'left', fill=NA)) %>%
    mutate(month=month(date), year=year(date), doy=yday(date)) %>%
    left_join(., junk, by=c('site','plot_code','date')) %>%
    ungroup()
  out_df <- tmp_df %>% 
    group_by(site, plot_code, year,month) %>%
    summarize(n_obs = sum(is.na(tot)==F),
              lf_mo = mean(u_sum, na.rm=T),
              lf_raw_u_mo = mean(tot, na.rm=T)*30.4) %>%
    ungroup() %>%
    filter(n_obs>=1) %>%
    mutate(date=ymd(paste(year,month,15),tz='UTC'))
  return(out_df)    
}
