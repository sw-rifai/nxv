# AggregateDendro & calibrate dates
library(tidyverse); library(lubridate); library(RcppRoll)

# DON'T SUBTRACT DAYS FROM CAXIUANA DATES!!! 
# I don't have the raw dendroband data so I have to work with the date intervals 
# that the processed dendro data came in. 
# THIS ALL NEEDS A REWRITE IF WE EVER GET THE ORIGINAL DBAND DATA FOR CAX
out_list <- list.files("outputs/stem_npp_month/", pattern = "stem_npp", full.names = T)
first_fname <- out_list[str_which(out_list, regex("CAX", ignore_case=T))]
rest_fname <- out_list[-str_which(out_list, regex("CAX", ignore_case=T))]
out_list <- c(first_fname, rest_fname)

atmp <- read_csv(out_list[1]) %>% 
  select(-npp_lm_agC_month) #%>% 
  # rename(nobs=nobs, date_col=date) %>% 
  # mutate(year=year(date_col), month=month(date_col)) %>% 
  # mutate(date = date_col-days(floor(date_diff/2))) %>% 
  # mutate(year = year(date), month=month(date)) 
 
tmp <- expand.grid(year = seq(min(atmp$year)-1, max(atmp$year)),
                   month=1:12) %>% as.tibble()
atmp <- full_join(atmp, tmp, by=c("year","month"))
atmp <- atmp %>% mutate(date = parse_date_time(paste(year,month,15), "ymd"))
atmp <- atmp %>% mutate(mu3=npp_bestEst_agC_month, mu3_hybrid=npp_bestEst_agC_month)

#
# atmp <- atmp %>% arrange(date) %>%
#   mutate(mu3=roll_meanr(npp_bestEst_agC_month,n = 3, fill = NA, na.rm=T)) %>% 
#   mutate(mu3_hybrid = ifelse(date_diff<45, mu3, npp_bestEst_agC_month)) 
  

for(i in 2:length(out_list)){
  btmp <- read_csv(out_list[i]) %>% #%>% select(date,site,plot_code, npp_bestEst_agC_month, nobs.x) %>% rename(nobs=nobs.x)
    select(-npp_lm_agC_month) %>% 
    rename(nobs=nobs.x, date_col=date) %>% 
    mutate(year=year(date_col), month=month(date_col)) %>% 
    mutate(date = date_col-days(floor(date_diff/2))) %>% 
    mutate(year = year(date), month=month(date)) %>% 
    mutate(newdate=parse_date_time(paste(year,month,15),"ymd")) %>% 
    select(-date) %>% rename(date=newdate)

  tmp <- expand.grid(year = seq(min(atmp$year)-1, max(atmp$year)),
                     month=1:12) %>% as.tibble()
  btmp <- full_join(btmp, tmp, by=c("year","month"))
  btmp <- btmp %>% mutate(date = parse_date_time(paste(year,month,15), "ymd"))
  btmp <- btmp %>% arrange(date) %>%
    mutate(mu3=roll_meanr(npp_bestEst_agC_month,n = 3, fill = NA, na.rm=T)) %>% 
    mutate(mu3_hybrid = ifelse(date_diff<45, mu3, npp_bestEst_agC_month))
    
  atmp <- bind_rows(atmp, btmp)
}

atmp <- atmp %>% filter(is.na(mu3_hybrid)==F)

atmp <- atmp %>% 
  mutate(mu3_hybrid = ifelse(plot_code=="BOB-06" & date==as.POSIXct("2014-08-15",tz="UTC"), NA, mu3_hybrid)) %>% 
  mutate(mu3_hybrid = ifelse(plot_code=="BOB-04" & date==as.POSIXct("2014-07-15",tz="UTC"), NA, mu3_hybrid)) %>%
  mutate(mu3_hybrid = ifelse(plot_code=="BOB-05" & date==as.POSIXct("2014-07-15",tz="UTC"), NA, mu3_hybrid)) %>%
  mutate(mu3_hybrid = ifelse(plot_code=="BOB-01" & date==as.POSIXct("2014-07-15",tz="UTC"), NA, mu3_hybrid)) %>%
  mutate(mu3_hybrid = ifelse(plot_code=="BOB-02" & date==as.POSIXct("2014-07-15",tz="UTC"), NA, mu3_hybrid)) %>%
  mutate(mu3_hybrid = ifelse(plot_code=="BOB-03" & date==as.POSIXct("2014-07-15",tz="UTC"), NA, mu3_hybrid)) %>%
  mutate(mu3_hybrid = ifelse(plot_code=="KOG-02" & date==as.POSIXct("2014-10-15",tz="UTC"), NA, mu3_hybrid)) %>%
  mutate(mu3_hybrid = ifelse(plot_code=="KOG-03" & date==as.POSIXct("2014-10-15",tz="UTC"), NA, mu3_hybrid)) %>%
  mutate(mu3_hybrid = ifelse(plot_code=="KOG-04" & date==as.POSIXct("2014-10-15",tz="UTC"), NA, mu3_hybrid)) %>%
  mutate(mu3_hybrid = ifelse(plot_code=="KOG-05" & date==as.POSIXct("2014-10-15",tz="UTC"), NA, mu3_hybrid)) %>%
  mutate(mu3_hybrid = ifelse(plot_code=="KOG-06" & date==as.POSIXct("2014-10-15",tz="UTC"), NA, mu3_hybrid))
  

fname <- paste0("DataTable_NPPliveWood_",gsub("-","",Sys.Date()),".csv")
write_csv(atmp, 
          path = paste0("../El_Niño_pantropical_forest_climate_anomalies/data/new_in_progress/DataTablesforModel20171102/",fname))
write_csv(atmp, 
          path = paste0("../El_Nino_StemNPP/data/gem_stem_npp/",fname))
write_csv(atmp, 
          path = paste0("../El_Nino_ws/data/gem_stem_npp/",fname))

# atmp %>% arrange(date) %>% 
#   mutate(mu3_hybrid = ifelse(date_diff<45, mu3, npp_bestEst_agC_month)) %>% 
#   ggplot(data=., aes(date, npp_bestEst_agC_month))+geom_point(cex=3)+geom_line()+
#   geom_point(aes(date, mu3), col="red", cex=1)+
#   geom_point(aes(date, mu3_hybrid), col="blue")

# library(tidyverse)
# out_list <- list.files("outputs/stem_npp_month/", pattern = "stem_npp", full.names = T)
#  
# atmp <- read_csv(out_list[1]) %>% select(date,site,plot_code, np p_bestEst_agC_month, nobs.x) %>% rename(nobs=nobs.x)
# 
# 
# for(i in 2:length(out_list)){
#   btmp <- read_csv(out_list[i]) %>% select(date,site,plot_code, npp_bestEst_agC_month, nobs.x) %>% rename(nobs=nobs.x)
#   atmp <- bind_rows(atmp, btmp)
# }
# fname <- paste0("DataTable_NPPliveWood_",gsub("-","",Sys.Date()),".csv")
# write_csv(atmp, 
#           path = paste0("../El_Niño_pantropical_forest_climate_anomalies/data/new_in_progress/DataTablesforModel20171102/",fname))
# 
# 
# 
# atmp %>% ggplot(data=., aes(date, npp_bestEst_agC_month))+geom_point()+
#   facet_wrap(~site)
