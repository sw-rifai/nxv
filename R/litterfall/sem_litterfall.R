
dat_nxv %>% 
  group_by(plot_code, year, month) %>% 
  summarize(val = mean(leavesflf_MgC_ha_day*days_in_month(date), na.rm=T), 
            val_sd = sd(leavesflf_MgC_ha_day*days_in_month(date), na.rm=T)) %>% 
  ungroup() %>% 
  group_by(year) %>% 
  summarize(tot = sum(val, na.rm=T), 
            tot_sd = sqrt(sum(val_sd**2)),
            nobs = n()) %>% 
  ungroup() %>% 
  mutate(tot_sem = tot_sd/sqrt(nobs)) %>% 
  ggplot(data=., aes(year, tot))+
  geom_ribbon(aes(year, ymax=tot+2*tot_sd, ymin=tot-2*tot_sd),alpha=0.25,lty=0)+
  geom_ribbon(aes(year, ymax=tot+2*tot_sem, ymin=tot-2*tot_sem),alpha=0.25,lty=0)+
  geom_point()


dat_nxv %>% 
  group_by(plot_code, year, month) %>% 
  summarize(val = mean(leavesflf_MgC_ha_day*days_in_month(date), na.rm=T), 
            val_sd = sd(leavesflf_MgC_ha_day*days_in_month(date), na.rm=T)) %>% 
  ungroup() %>% 
  group_by(year) %>% 
  summarize(tot = sum(val, na.rm=T), 
            tot_sd = sqrt(sum(val_sd**2)),
            nobs = n()) %>% 
  ungroup() %>% 
  mutate(tot_sem = tot_sd/sqrt(nobs)) %>% 
  summarize(u = mean(tot),
            u_sd = sd(tot),
            nobs=n(),
            u_sem = u_sd/sqrt(nobs))
