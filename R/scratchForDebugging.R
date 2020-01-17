census_all$`Plot Code` %>% unique

census_all %>% filter(`Plot Code`=="CAX-01") %>% pull(D4)


(dendrometer$tree_tag %>% unique) %in% census$tree_tag


dend_all %>% filter(plot_code=="KEN-02")
dend_all %>% filter(plot_code=="KEN-01")

census_all %>% filter(`Plot Code`=="KEN-01") %>% View


library(lubridate)
d1 <- parse_date_time(paste(2011,1,1),"ymd")
d2 <- parse_date_time(paste(2012,1,1),"ymd")
d3 <- parse_date_time(paste(2013,6,15),"ymd")

get_time_diffs(c(d1,d2,d3))



dend_all %>% group_by(plot_code) %>% 
  summarize(u_dendrometer_reading_mm = mean(dendrometer_reading_mm,na.rm=T)) %>% View



#######################################################################################
#!--- CHECK IF MEASUREMENTS ARE IN CM OR MM ---!
census_all %>% filter(d4<100) %>% pull(d4) %>% hist # There are a log of 0s, should be NAs 
census_all %>% filter(d4>4000) %>% pull(d4) %>% hist # There are a log of 0s, should be NAs 

dend_all %>% group_by(plot_code) %>% 
  summarize(u_dendrometer_reading_mm = mean(dendrometer_reading_mm,na.rm=T), 
            sd_dendrometer_reading_mm = sd(dendrometer_reading_mm, na.rm=T)) %>% View
  pull(u_dendrometer_reading_mm) %>% hist

dend_all %>% group_by(plot_code,year,month) %>% 
  summarize(n_obs=n()) %>% arrange(plot_code,year,month) %>% View

dend_all %>% filter(plot_code=="TAM-05") %>% filter(year>=2015) %>% View

dend_all %>% filter(plot_code=="BOB-02") %>% 
  filter(dendrometer_reading_mm<=0.1) %>% dim
  # pull(dendrometer_reading_mm) %>% summary

dend_all %>% filter(plot_code=="KOG-02") %>% 
  filter(dendrometer_reading_mm<=0.1) %>% dim
pull(dendrometer_reading_mm) %>% summary

dend_all %>% filter(plot_code=="SAF-05") %>% 
  filter(dendrometer_reading_mm<=0.1) %>% dim
pull(dendrometer_reading_mm) %>% summary


dend_all %>% group_by(plot_code) %>% 
  summarize(n_neg = sum(dendrometer_reading_mm<=0)) %>% View

##################################################################################
# what's going on with TAM census? 
tmp_npp_census_trees
tmp_npp_census

tmp_npp_dend %>% names %>% sort

tmp_npp_dend %>% 
  ggplot(data=., aes(date, npp_weightedAvgTrees_month_dend))+geom_path()



npp_tree %>% names %>% sort

npp_tree %>% ggplot(data=., aes(decDate, agC_Mg_diff))+geom_point()

dend %>% names %>% sort

library(ggjoy)
dend %>% ggplot(data=., aes(x=dendrometer_reading_mm, y=as.factor(date)))+geom_joy()

dend %>% pull(date) %>% unique %>% sort
dend %>% filter(year==2015 & month==7)
dend %>% filter(date==parse_date_time("2015-07-19", "ymd"))


npp_tree %>% group_by(plot_code, year, month) %>% 
  summarize(date=mean(decDate), nobs=n()) %>%
  ggplot(data=., aes(date, nobs))+geom_path()+geom_point()

npp_tree %>% group_by(decDate) %>% 
  summarize(u=mean(dateDiff,na.rm=T),
            b_min=min(dateDiff, na.rm=T), 
            b_max=max(dateDiff,na.rm=T), 
            u_w = mean(weight, na.rm=T), 
            cd = mean(closestCensusDate)) %>% 
  arrange(decDate) %>% 
  ggplot(data=., aes(decDate, cd))+geom_path()+geom_point()+
  theme_bw()

npp_tree %>% group_by(dateDiff) %>% 
  summarize(u=mean(agC_Mg_diff,na.rm=T)) %>% 
  ggplot(data=., aes(dateDiff, u))+geom_point()+geom_path()+theme_bw()

npp_tree %>% group_by(plot_code, year, month) %>% 
  filter(decDate>2015.5) %>% 
  summarize(date=mean(decDate), nobs=n()) %>%
  ggplot(data=., aes(date, nobs, color=month))+geom_path()+geom_point()+
  scale_color_viridis()

npp_tree %>% group_by(plot_code, year, month) %>% 
  # filter(decDate>2010) %>% 
  summarize(date=mean(decDate), 
            nobs=n(), 
            npp_weightedAvgTrees_day_dend = weighted.mean(nppacw_Mg_tree_day, w=weight, na.rm=T), 
            npp_avgTrees_day_dend=mean(nppacw_Mg_tree_day, na.rm=T)) %>%
  ggplot(data=., aes(date, npp_weightedAvgTrees_day_dend, size=nobs, color=month))+
  geom_point()+
  geom_abline(aes(intercept=0, slope=0))+
  scale_color_viridis() + theme_bw()

npp_tree %>% group_by(plot_code, year, month) %>% 
  # filter(decDate>2010) %>% 
  summarize(date=mean(decDate), 
            nobs=n(), 
            npp_weightedAvgTrees_day_dend = weighted.mean(nppacw_Mg_tree_day, w=weight, na.rm=T), 
            npp_avgTrees_day_dend=mean(nppacw_Mg_tree_day, na.rm=T)) %>%
  ggplot(data=., aes(date, npp_weightedAvgTrees_day_dend, size=nobs, color=month))+
  geom_point()+
  geom_abline(aes(intercept=0, slope=0))+
  scale_color_viridis() + theme_bw()



# Check tree_tags are the same in dendrometer and census datasets
unique(dend1$tree_tag) %in% unique(census$tree_tag) %>% table # --- REQUIRES ATTENTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!--- #
print(paste("these tree tags from the dendrometers are not in the census:",dend1$tree_tag[!(dend1$tree_tag %in% cen$tree_tag)] %>% unique))


fillDeadDbh <- function(census, min_dend_date){
  #--- This function:
  #--- 1) finds the max diameter of dead trees from previous census measurement dates 
  #--- 2) replaces the 0s in the original forestPlots census object with the prev max dbh
  #--- 3) filters the census file for being within 2 years of the first dendrometer measurment date
  
  tmp <- census 

  # ! There are a LOT of 0s, which I think are meant to be NAs! All these trees are dead.
  tmp <- tmp %>% mutate(d0=replace(d0, d0==0, NA), 
                 d1=replace(d1, d1==0, NA), 
                 d2=replace(d2, d2==0, NA), 
                 d3=replace(d3, d3==0, NA), 
                 d4=replace(d4, d4==0, NA)) 

  dead_tags <- census_obs %>% filter(f1=="0") %>% pull(tree_tag) %>% unique()
  
  tmp_dead <- tmp %>% filter(tree_tag %in% dead_tags)
  
  tmp_dead <- tmp_dead %>%
    rowwise() %>% mutate(dbh_dead=mean(c(d0,d1,d2,d3,d4),na.rm=T)) %>%
    ungroup() %>%  filter(dbh_dead<5000) %>% 
    group_by(date, tree_tag) %>% 
    summarize(dbh_dead=max(dbh_dead))
  
  # generate the dbh from the mean of the measurements
  tmp <- tmp %>%
    rowwise() %>% mutate(#dbh_mm_sd=sd(c(d0,d1,d2,d3,d4),na.rm=T), 
                         dbh=mean(c(d0,d1,d2,d3,d4),na.rm=T)) %>%
    ungroup() %>%  filter(dbh<5000)
  
  tmp <- left_join(tmp,tmp_dead,by=c("tree_tag","date")) %>% 
    rowwise() %>% 
    mutate(dbh=max(c(dbh, dbh_dead),na.rm=T)) %>% select(-dbh_dead)
  
  
  # filter census for dates within one year of min dend date
  tmp <- tmp %>% filter(date >= (min_dend_date-years(2)))
  
  return(tmp)
}
mutate
max(c(0,10,NA), na.rm=T)

# --- DELETE THIS COMMENTED OUT SECTION LATER -------------------
# jesusNppTree <- left_join(npp_tree, 
# census %>% filter(census_date==min(census_date)) %>% select(tree_tag, family, genus, species, tree_tag, density), 
# by="tree_tag")
# write_csv(jesusNppTree, path = "ouputs/BOB02_npp_agC_MgPerDay_2012_2017.csv")
# 
# jesusNppTree %>% mutate(date=paste(year, month, 15), "ymd") %>% 
#   filter(year>=2015) %>% mutate(y=agC_Mg_diff/agC_Mg) %>% 
#   lm(y ~ scale(density), data=. ) %>% summary 
#-----------------------------------------------------------------

npp_tree$agC_Mg_diff/npp_tree$dateDiff

npp_tree %>% group_by(decDate) %>% summarize(u=mean(dateDiff)) %>% 
  ggplot(data=., aes(decDate, u))+geom_point()

npp_tree %>% group_by(decDate) %>% summarize(u=mean(agC_Mg_diff)) %>% 
  ggplot(data=., aes(decDate, u))+geom_point()

npp_tree %>% filter(decDate>2015.03 & decDate<2015.04) %>% View

npp_tree %>% filter(decDate>2015.03 & decDate<2015.04) %>% 
  pull(agC_Mg_diff) %>% hist

dendrometer %>% group_by(date) %>% 
  summarize(u=mean(dendrometer_reading_mm,na.rm=T)) %>% 
  ggplot(data=., aes(date, u))+geom_path()+
  labs(y="Mean of dendrometer_reading_mm")+
  ggtitle("TAM-05")

dendrometer %>% group_by(date) %>% 
  summarize(u=mean(dendrometer_reading_replaced_mm, na.rm=T)) %>% 
  ggplot(data=., aes(date, u))+geom_path()+
  labs(y="Mean of dendrometer_reading_replaced_mm")+
  ggtitle("TAM-05")

dendrometer %>% group_by(date) %>% 
  summarize(u=sum(is.na(dendrometer_reading_replaced_mm)==F)) %>% 
  ggplot(data=., aes(date, u))+geom_path()+
  labs(y="Number of dendrometer_reading_replaced_mm")+
  ggtitle("TAM-05")

dendrometer$dbh_nodirection_mm %>% hist
dendrometer$dbh_northsouth_mm #NA
dendrometer$dbh_westeast_mm #NA


dendrometer %>% 
  ggplot(data=., aes(dendrometer_reading_mm, dendrometer_reading_replaced_mm))+geom_point()+
  geom_abline(aes(intercept=0,slope=1), col="red")


dendrometer %>% 
  select(comments, dendrometer_reading_mm, dendrometer_reading_replaced_mm) %>% 
  View

library(ggjoy)
dend_all %>%
  ggplot(data=., aes(x=dendrometer_reading_mm, y=plot_code)) + 
  geom_joy()

dend_all %>% group_by(plot_code) %>% summarize(min(dendrometer_reading_mm)) %>% View

jnk <- dend1$tree_tag[!(dend1$tree_tag %in% cen$tree_tag)] %>% unique
dend1 %>% filter(tree_tag %in% jnk) %>% group_by(tree_tag) %>% 
  summarize(nobs=n()) %>% as.data.frame

census %>% filter(tree_tag=="114")


census %>% filter(tree_tag=="282") %>% arrange(date) %>% pull(dbh)
tmp %>% filter(tree_tag=="282") %>% pull(dbh) #[1] 632 632 632 632 632 632 632 632 632
tmp %>% rowwise() %>% 
  mutate(dbh_new = est_dbh_from_dendro(dbh = dbh, dendrometer_reading_mm = dendrometer_reading_mm), 
         dbh_replaced = est_dbh_from_dendro(dbh=dbh, dendrometer_reading_mm = dendrometer_reading_replaced_mm)) %>% 
  filter(tree_tag=="282") %>% 
  ggplot(data=., aes(date, dbh_new))+geom_path()+
  geom_path(aes(date, dbh_replaced, color="blue"))


dend_all %>%
  ggplot(data=., aes(dendrometer_reading_mm, dendrometer_reading_replaced_mm, color=plot_code))+
  geom_point()+geom_smooth(method="lm", se=F)

dend_all %>% group_by(plot_code) %>% 
  summarize(nobs = sum(is.na(dendrometer_reading_replaced_mm)==F)) %>% View

dend %>% filter(tree_tag=="100") %>% View
census %>% filter(tree_tag=="100") %>% select(date,dbh)
# tmp_new_dbands <- 
tmp %>% group_by(date, tree_tag) %>% 
  summarize(u=mean(thisdbh_mm)) %>% 
  filter(tree_tag=="100") %>% 
  ggplot(data=., aes(date,u))+geom_path()

# tree_tag 12.1, 163.1 has multiple dendro replacements

dend %>% filter(tree_tag == "12.1") %>% arrange(date) %>% 
  ggplot(data=., aes(date, current_dbh))+geom_path()+
  geom_vline(aes(xintercept=date_decimal(dbh_replace_date)))

dend %>% filter(tree_tag == "296") %>% arrange(date) %>%
  select(date, dbh_replace_date,
         current_dbh, dendrometer_reading_mm, 
         dendrometer_reading_replaced_mm) %>%
  View

dend %>% filter(tree_tag == "180.1") %>% arrange(date) %>% 
  ggplot(data=., aes(date, current_dbh))+geom_path()+
  geom_vline(aes(xintercept=date_decimal(dbh_replace_date)))
library(viridis)

dend %>% filter(tree_tag=="296") %>%
  arrange(date) %>% 
  ggplot(data=., aes(date, c(NA,diff(dendrometer_reading_mm)), color=tree_tag))+
  geom_path()+theme(legend.position="none")+
  scale_color_viridis(discrete = T)

dend %>% #filter(tree_tag=="180.1") %>%
  arrange(date) %>% 
  ggplot(data=., aes(date, c(NA,diff(dendrometer_reading_replaced_mm)), color=tree_tag))+
  geom_path()+theme(legend.position="none")+
  scale_color_viridis(discrete = T)

dend %>% names %>% sort

dend$dendrometer_reading_mm %>% is.na %>% table

dend %>% filter(tree_tag=="180.1") %>% 
  ggplot(data=., aes(date, delta1))+geom_path()


# 
# [1] "100"   "100.1" "102"   "103"   "104"   "105"   "107"   "11"    "110"  
# [10] "111"   "113"   "116"   "12.1"  "121"   "122"   "124"   "124.1" "125"  
# [19] "127"   "13"    "13.1"  "135.1" "139.1" "140"   "142"   "143"   "147"  
# [28] "147.1" "148"   "149"   "150"   "151"   "151.1" "152"   "154"   "157"  
# [37] "157.1" "162.1" "163.1" "164"   "165"   "169"   "170"   "170.1" "173"  
# [46] "173.1" "174"   "180"   "180.1" "184.1" "184.2" "186"   "188"   "191"  
# [55] "192"   "199"   "202"   "203"   "206"   "207"   "210"   "212"   "217"  
# [64] "218"   "223"   "228"   "230"   "232.2" "232.4" "233"   "235"   "237"  
# [73] "238"   "242"   "243"   "244"   "245"   "247"   "248"   "249"   "250"  
# [82] "250.1" "252"   "26"    "264"   "27"    "27.1"  "271"   "276"   "278"  
# [91] "279"   "281"   "282"   "284"   "285"   "287.1" "291"   "292"   "295"  
# [100] "296"   "298"   "298.1" "299"   "300"   "300.1" "304"   "306"   "307"  
# [109] "311"   "313"   "316"   "327"   "328"   "330"   "332"   "337"   "338"  
# [118] "339"   "342"   "343"   "353"   "358"   "359"   "36"    "360"   "361"  
# [127] "362"   "366"   "370"   "384"   "385"   "388"   "390"   "395"   "40"   
# [136] "401"   "413"   "414"   "415.1" "416"   "424"   "429"   "431"   "434"  
# [145] "435"   "438"   "441"   "442"   "445.1" "446"   "449"   "45.1"  "452"  
# [154] "454"   "456"   "463"   "465"   "477"   "482"   "485.1" "490"   "491"  
# [163] "494"   "495"   "496"   "497"   "51"    "510"   "512"   "517"   "520"  
# [172] "521"   "524"   "526"   "529"   "530"   "532"   "535"   "536.1" "544"  
# [181] "544.2" "545"   "546"   "550"   "554"   "555"   "559"   "56"    "560"  
# [190] "568"   "570"   "570.1" "570.2" "572.1" "578"   "582"   "583"   "584"  
# [199] "585"   "586"   "587"   "6"     "61"    "66"    "68"    "69"    "7"    
# [208] "70"    "72"    "77"    "79"    "8"     "83"    "85"    "88.1"  "9"    
# [217] "92"    "96.1"  "96.2"  "97"    "99"   

library(viridis)
dend %>% filter(tree_tag=="572.1") %>%
  arrange(date) %>% 
  ggplot(data=., aes(date, current_dbh))+
  geom_path()+theme(legend.position="none")+
  scale_color_viridis(discrete = T)

dend %>% filter(tree_tag=="296") %>%
  arrange(date) %>% 
  ggplot(data=., aes(date, current_dbh))+
  geom_path()+theme(legend.position="none")+
  scale_color_viridis(discrete = T)


dend %>% filter(date>parse_date_time(2015, "y")) %>% 
  pull(dendrometer_reading_mm) %>% hist
dend %>% filter(date>parse_date_time(2016, "y")) %>% 
  pull(dendrometer_reading_mm) %>% hist


dend %>% group_by(tree_tag) %>% arrange(date) %>% 
  mutate(baseline_dbh = ifelse(is.na(baseline_dbh)==T, lead(baseline_dbh), baseline_dbh)) %>%
  mutate(baseline_dbh = cummax(baseline_dbh)) %>% 
  filter(tree_tag=="296") %>% pull(baseline_dbh)

  
  mutate(baseline_dbh = cummax(baseline_dbh)) %>% 
  filter(tree_tag=="296") %>% pull(baseline_dbh)

dend %>% filter(tree_tag=="296") %>% pull(baseline_dbh)


dend %>% filter(tree_tag == "296") %>% arrange(date) %>%
  select(date, 
         current_dbh, dendrometer_reading_mm, 
         dendrometer_reading_replaced_mm, delta1, baseline_dbh) %>%
  View

dend

ave(c(1:10), FUN)


dend %>% arrange(date) %>% group_by(tree_tag) %>% rowwise() %>% 
  mutate(baseline_dbh = cummax(baseline_dbh)) %>% 
  filter(tree_tag=="296") %>%   
  select(date, 
         current_dbh, dendrometer_reading_mm, 
         dendrometer_reading_replaced_mm, delta1, baseline_dbh) %>% 
  View

dend %>% arrange(date) %>% group_by(tree_tag) %>% 
  mutate(baseline_dbh = ifelse(is.na(baseline_dbh)==T, lead(baseline_dbh), baseline_dbh)) %>% 
  filter(tree_tag=="296") %>%   
  select(date, 
         current_dbh, dendrometer_reading_mm, 
         dendrometer_reading_replaced_mm, delta1, baseline_dbh) %>% 
  View

ave(dend$baseline_dbh, dend$tree_tag, FUN=cummax, na.rm=T) %>% summary


dend$current_dbh <- est_dbh_from_dendro(dbh=dend$baseline_dbh, dendrometer_reading_mm = dend$dendrometer_reading_mm)

# add first dbh measurement (cm) to dendrometer measurement (cm) = thisdbh
# dend$dendrometer_reading_mm <- as.numeric(dend$dendrometer_reading_mm) # 
# dend$thisdbh_mm <- ((dend$dbh*pi) + (dend$dendrometer_reading_mm))/pi 
# dend$thisdbh_cm <- dend$thisdbh_mm/10


# #--- ADJOIN THE REPLACED DENDROMETER READINGS TO THE FILE ------------------ #
# #--- Not all GEM plots have dendrometer_reading_replaced_mm values. 
# #--- This needs to be checked in the future in case the criteria for baseline_dbh changes
# #--- The following seems inefficient and easy to break. A well thought out rewrite is warrented.
# suppressWarnings(tmpJunk <- dend %>% rowwise() %>% 
#   mutate(new_dbh_date_req = ifelse(is.na(dendrometer_reading_replaced_mm)==F, paste(year,month,day), NA)) %>% 
#   ungroup() %>% mutate(new_dbh_date_req=decimal_date(parse_date_time(new_dbh_date_req, "ymd"))) %>% 
#   group_by(tree_tag) %>% 
#   summarize(dbh_replace_date=min(new_dbh_date_req, na.rm=T)) %>% 
#   mutate(dbh_replace_date=ifelse(is.infinite(dbh_replace_date), NA, dbh_replace_date)))
# tmpJunk <- left_join(tmpJunk,
#           (dend %>% mutate(dbh_replace_date=decimal_date(date)) %>%
#            select(tree_tag, dbh_replace_date, current_dbh)), 
#           by=c("tree_tag","dbh_replace_date")) %>% rename(baseline_dbh2=current_dbh)
# dend <- left_join(dend, tmpJunk, by="tree_tag")
# dend <- dend %>% rowwise() %>% 
#   mutate(current_dbh = ifelse(decimal_date(date)>=dbh_replace_date, 
#                               est_dbh_from_dendro(dbh = baseline_dbh2,
#                                   dendrometer_reading_mm = dendrometer_reading_replaced_mm), 
#                               current_dbh)) %>% ungroup()
#  # --------------------------------------------------------------------------- #  

tmp_npp_dend %>% names %>% sort

census_all %>% filter(plot_code=="KEN-01") %>% pull(tree_tag) %>% unique %>% sort

npp_tree %>% pull(tree_tag) %>% unique %>% class

npp_tree %>% filter(tree_tag=="101") %>% 
  ggplot(data=., aes(decDate, agC_Mg_diff))+geom_path()

npp_tree %>% filter(agC_Mg>0.25) %>%
  ggplot(data=., aes(decDate, agC_Mg_diff, color=tree_tag))+
  geom_path() +
  # geom_smooth(method="lm",se=F)+
  theme_bw()+
  viridis::scale_color_viridis(discrete=T)+
  theme(legend.position = "none")
 
dend %>% group_by(tree_tag) %>% summarize(u=min(current_dbh)) %>% arrange((u))
# tree_tag        u
# <chr>    <dbl>
#   1    263.2 95.50394
# 2      172 98.47211
# 3    156.1 98.56761
# 4     28.1 98.69493
# 5    172.2 98.72676
# 6      425 98.79042
# 7    390.1 99.29972
# 8    163.1 99.52254
# 9    356.1 99.52254
# 10    583.2 99.58620

dend %>% filter(tree_tag=="263.2") %>% 
  ggplot(data=., aes(date, dendrometer_reading_mm))+geom_path()+
  geom_vline(aes(xintercept=as.POSIXct("2015-01-01")),col="red")

dend %>% filter(tree_tag=="172") %>% 
  ggplot(data=., aes(date, dendrometer_reading_mm))+geom_path()+geom_point()+
  geom_vline(aes(xintercept=as.POSIXct("2015-01-01")),col="red")

dend %>% filter(tree_tag=="583.2") %>% 
  ggplot(data=., aes(date, dendrometer_reading_mm))+geom_path()+geom_point()+
  geom_vline(aes(xintercept=as.POSIXct("2015-01-01")),col="red")

dend %>% filter(tree_tag=="356.1") %>% 
  ggplot(data=., aes(date, dendrometer_reading_mm))+geom_path()+geom_point()+
  geom_vline(aes(xintercept=as.POSIXct("2015-01-01")),col="red")

dend %>% filter(tree_tag=="163.1") %>% #good example of problem
  ggplot(data=., aes(date, dendrometer_reading_mm))+geom_path()+geom_point()+
  geom_vline(aes(xintercept=as.POSIXct("2015-01-01")),col="red")

dend %>% filter(tree_tag=="390.1") %>% 
  ggplot(data=., aes(date, dendrometer_reading_mm))+geom_path()+geom_point()+
  geom_vline(aes(xintercept=as.POSIXct("2015-01-01")),col="red")

dend %>% filter(date>as.POSIXct("2014-01-01")) %>% 
  ggplot(data=., aes(date, dendrometer_reading_mm))+geom_point()+
  geom_violin(aes(x = as.factor(), y=dendrometer_reading_mm))



dend %>% filter(tree_tag=="163.1") %>% #good example of problem
  ggplot(data=., aes(date, delta1))+geom_path()+geom_point()+
  geom_vline(aes(xintercept=as.POSIXct("2015-01-01")),col="red")

dev.new()
dend %>% filter(tree_tag=="163.1") %>% #good example of problem
  ggplot(data=., aes(date, current_dbh))+geom_path()+geom_point()+
  geom_vline(aes(xintercept=as.POSIXct("2015-01-01")),col="red")

dend %>% filter(tree_tag=="163.1") %>% #good example of problem
  ggplot(data=., aes(date, baseline_dbh))+geom_path()+geom_point()+
  geom_vline(aes(xintercept=as.POSIXct("2015-01-01")),col="red")

dend %>% ggplot(data=., aes(x = delta1))+geom_histogram()+
  facet_wrap(~date)
dend %>% ggplot(data=., aes(x = current_dbh))+geom_histogram()+
  facet_wrap(~date)

library(viridis)
out_df %>% ggplot(data=., aes(date, weighted_npp_agC_month, color=nobs, size=nobs))+
  geom_point()+scale_color_viridis()
out_df %>% ggplot(data=., aes(date, npp_avgtrees_day_dend_sd, color=nobs, size=nobs))+
  geom_point()+scale_color_viridis()


npp_tree %>% filter(near(decDate,2015.038, tol=0.01)) %>% 
  pull(nppacw_Mg_tree_day) %>% hist

npp_tree %>% group_by(decDate) %>% 
  summarize(npp_max = max(nppacw_Mg_tree_day)) %>% 
  ggplot(data=., aes(decDate, npp_max))+geom_path()


npp_tree %>% filter(near(decDate,2015.038, tol=0.01)) %>% 
 filter(nppacw_Mg_tree_day==max(nppacw_Mg_tree_day))

dend %>% filter(tree_tag=="210") %>% ggplot(data=., aes(date, current_dbh))+geom_path()

dend %>% filter(current_dbh>500) %>% 
  ggplot(data=., aes(date, current_dbh, color=tree_tag))+geom_path()+
  theme_bw()+scale_color_viridis(discrete = T)



install.packages("microbenchmark")
library(microbenchmark)
x <- runif(100)
microbenchmark(
  sqrt(x),
  x ^ 0.5
)

x <- runif(100)
y <- runif(100)
microbenchmark(
  dend %>% arrange(date),
  append(x,y)
)

microbenchmark(
  arrange(dend$date),
  get_time_diffs(tmpDates)
)

microbenchmark(
thistree  <- which(dend$tree_tag == uid[ii]), #dend$tree_tag[ii])
if(length(thistree)==1){next},
agC_Mg       <- dend$agC_Mg[thistree],
tag       <- dend$tree_tag[thistree],
#agCdiff   <- dend$agCdiff[thistree]
agC_Mg_diff   <- ave(dend$agC_Mg[thistree], FUN = function(x) c(NA, diff(x))),
tmpyear      <- dend$year[thistree],
tmpmonth     <- dend$month[thistree],
tmpplot_code <- dend$plot_code[thistree],
tmpDates <- dend$date[thistree],
tmpDayDiffs <- get_time_diffs(tmpDates), # IN THEORY THIS SHOULD GIVE THE NUMBER OF DAYS BETWEEN OBSERVATIONS
tmpDayDiffs <- c(NA,tmpDayDiffs), # need to tack an NA to the front of the dayDiff vector
dbhs <- dend$dbh[thistree],
tree_tags <- dend$tree_tag[thistree],

print(paste("uid[",ii,"]",length(aa)-length(gg))),

aa            <- c(aa, tmpplot_code),
bb            <- c(bb, tag),
cc            <- c(cc, tmpyear),
dd            <- c(dd, tmpmonth),
ee            <- c(ee, agC_Mg),
ff            <- c(ff, agC_Mg_diff),
gg            <- c(gg, tmpDayDiffs),
hh            <- c(hh, dbhs),
tag_vec       <- c(tag_vec, tree_tags)
)

microbenchmark(
diff(tmpDates), 
get_time_diffs(tmpDates)
)

npp_tree$dateDiff %>% summary

diff(tmpDates, units="days")
diff(tmpDates, units="seconds")

diff(as.numeric(tmpDates))/(24*60*60)

diff(as.numeric(c(tmpDates[1],tmpDates)))/(24*60*60)
diff((c(tmpDates[1],tmpDates)))

out_df %>% ggplot(data=., aes(date, npp_avgtrees_day_dend_sd))+
  geom_point()

dend_all %>% filter(plot_code=="YAY-03") %>% pull(date) %>% unique

dend %>% filter(current_dbh<200) %>% 
  ggplot(data=., aes(date, current_dbh))+geom_path()
  
dend %>% filter(current_dbh<200) %>% 
  ggplot(data=., aes(date, dendrometer_reading_mm))+geom_path()

dend %>% pull(date) %>% unique %>% length
dend %>% pull(tree_tag) %>% unique() %>% length()
104*14
dim(dend)

dend %>% arrange(date) %>% group_by(date) %>% 
  summarize(nobs=n()) %>% 
  ggplot(data=., aes(date, nobs))+geom_path()

dend %>% mutate(dendDate=parse_date_time(paste(year,month,15),"ymd")) %>% 
  group_by(dendDate) %>% 
  summarize(nobs=n()) %>% 
  ggplot(data=., aes(dendDate, nobs))+geom_path()

dend %>% filter(current_dbh==max(current_dbh)) %>% pull(tree_tag) %>% unique

dend %>% filter(tree_tag=="451") %>% #good example of problem
  ggplot(data=., aes(date, delta1))+geom_path()+geom_point()

dend %>% filter(tree_tag=="451") %>% #good example of problem
  ggplot(data=., aes(date, current_dbh))+geom_path()+geom_point()

dend %>% filter(tree_tag=="451") %>% #good example of problem
  ggplot(data=., aes(date, baseline_dbh))+geom_path()+geom_point()

dend %>% ungroup %>% filter(current_dbh>1500)


dend1 %>% filter(tree_tag=="537") %>% #good example of problem
  ggplot(data=., aes(date, dendrometer_reading_mm))+geom_path()+geom_point()

dend1 %>% group_by(tree_tag) %>% arrange(tree_tag, date) %>% 
  ggplot(data=., aes(date, dendrometer_reading_mm,color=tree_tag))+geom_path()+geom_point()

dend %>% group_by(tree_tag) %>% arrange(tree_tag, date) %>% 
  ggplot(data=., aes(date, current_dbh,color=tree_tag))+geom_path()+geom_point()


diff(c(0,NA,10,NA,20,30,50))

est_dbh_from_dendro(100, NA)

dend %>% pull(dendrometer_reading_mm) %>% summary

dend$dendrometer_reading_mm %>% max

census$date %>% unique
census %>% arrange(date) %>% group_by(date) %>% summarize(n_trees=n())
census$plot_code


npp_tree$pos_meas_error %>% table
npp_tree %>% filter(nppacw_Mg_tree_day==max(nppacw_Mg_tree_day))

npp_tree %>% filter(tree_tag=="451") %>% 
  ggplot(data=., aes(decDate, nppacw_Mg_tree_day))+geom_path()

dend %>% 
  ggplot(data=., aes(date, current_dbh, color=tree_tag))+geom_path()+geom_point()

dend %>% filter(tree_tag=="451") %>% lm(current_dbh~date, data=.) %>% residuals() %>% abs() %>% max()

dend %>% ggplot(data=., aes(date, delta1,color=tree_tag))+geom_path()


dend$dendrometer_reading_mm %>% is.na %>% table
dend$dendrometer_reading_replaced_mm %>% is.na %>% table


library(broom)

vec_gte3 <- dend %>% group_by(tree_tag) %>% summarize(nobs=n()) %>% 
  filter(nobs>=3) %>% pull(tree_tag)

tmp_lm_df <- dend %>% filter(tree_tag %in% vec_gte3) %>% ungroup()
tmp_lm_fit <- tmp_lm_df %>% group_by(tree_tag) %>% 
  do(growthFit = lm(current_dbh~date, data=.)) %>% ungroup()
tidy(tmp_lm_df, tmp_lm_fit)

tmp_lm_fit %>% names
tmp_lm_fit$growthFit %>% residuals()


# extract max residuals from linear trend, match with tree_tag [ok], and date of max resid[to do]
library(purrr)
vec_max_resids <- dend %>% filter(tree_tag %in% vec_gte3) %>%
  split(.$tree_tag) %>% 
  map(~lm(current_dbh~date, data=.)) %>% 
  map(residuals) %>% 
  map_dbl(max)

df_resids <- tibble(tree_tag=vec_max_resids %>% names() %>% as.character(), 
                        max_resid=vec_max_resids)

vec_max_resids %>% length
vec_gte3 %>% length

get_date_resid <- function(date,current_dbh){
  fit <- lm(current_dbh~date); 
  vec_resids <- residuals(fit); 
  max_resid <- max(vec_resids); 
  pos <- which(vec_resids==max_resid)#position of max residual
  out <- list((x[pos]))
  return(out)
}

gt<-function(x,b){x>b}
gt_sum<-function(x,b){sum(x>b)}


dend %>% filter(tree_tag %in% vec_gte3) %>%
  split(.$tree_tag) %>% 
  map(~lm(current_dbh~date, data=.)) %>% 
  map(residuals) %>% 
  map_dbl(gt_sum,20)


dend %>% filter(tree_tag %in% vec_gte3) %>%
  split(.$tree_tag) %>% 
  map(.$date, .$current_dbh, get_date_resid)  


dend %>% filter(tree_tag %in% vec_gte3) %>%
  split(.$tree_tag) %>% 
  map(~lm(current_dbh~date, data=.)) %>% 
  map(residuals) %>% 
  map(detect_index(gte_10))

dend %>% filter(tree_tag %in% vec_gte3) %>%
  split(.$tree_tag) %>% 
  map(~lm(current_dbh~date, data=.)) %>% 
  map(residuals) %>% 
  map_lgl(., .gt, b=10)

dend %>% filter(tree_tag %in% vec_gte3) %>%
  split(.$tree_tag) %>% 
  map(~lm(current_dbh~date, data=.)) %>% 
  map(residuals) %>%
  as.tibble(default=NA)

dend %>% filter(max_resid >=25) %>% 
  ggplot(data=., aes(date, current_dbh, color=tree_tag))+geom_path()+
  geom_smooth(method='lm', se=F, lty=3)

dend %>% filter(max_resid >=25) %>% 
  ggplot(data=., aes(date, current_dbh, color=tree_tag))+geom_path()+
  geom_path(aes(date, baseline_dbh), lty=3)

dend %>% 
  split(.$tree_tag) %>% 
  map(diff(current_dbh, data=.))

dend %>% filter(max_resid >=20) %>% 
  ggplot(data=., aes(date, current_dbh, color=tree_tag))+geom_path()+
  geom_path(aes(date, baseline_dbh), lty=3)

dend %>% filter(max_resid >=20) %>% 
  ggplot(data=., aes(date, dendrometer_reading_mm, color=tree_tag))+geom_path()+
  geom_path(aes(date, delta1), lty=3)

match(census_all$plot_code=="NXV-02")

jnk <- census_all$plot_code =="NXV-02"

census_all[match(jnk==T),]

census_all %>% filter(plot_code=="NXV-02")
census_all %>% pull(plot_code) %>% `==` "NXV-02"

census_all %>% filter(plot_code=="SAF-01") %>% View

d_pn <- dend_all %>% pull(plot_code) %>% unique() %>% sort 
c_pn <- census_all %>% pull(plot_code) %>% unique() %>% sort
com_pn <- full_join(tibble(plot_code=d_pn, d="dband"), 
                    tibble(plot_code=c_pn, d="census"), by="plot_code")


dend_all$plot_code %>% unique
dend_all %>% filter(plot_code=="STO-03") %>% pull(tree_tag)
census_all %>% pull(plot_code) %>% unique
census_all %>% filter(plot_code=="STO-03") %>% View

dend_all <- readxl::read_excel("data/dendro_20171103.xlsx",sheet=1)
dend_all %>% filter(plot_code=="STO-03") %>% pull(tree_tag)

dend_all <- readxl::read_excel("data/dendro_STM.xlsx",sheet=1)
dend_all %>% filter(plot_code=="STO-03") %>% pull(tree_tag)

dend_all %>% filter(plot_code=="STO-03") %>% pull(tree_tag) %>% unique
dend_all %>% filter(plot_code=="STO-03") %>% pull(dendrometer_reading_mm)



dend_all <- read_csv("data/dendro_STM.csv")
dend_all %>% mutate(dendrometer_reading_mm=as.numeric(dendrometer_reading_mm), 
                    dendrometer_reading_replaced_mm=as.numeric(dendrometer_reading_replaced_mm))

census_all %>% filter(plot_code=="STO-03") %>% names %>% sort
census_all %>% filter(plot_code=="STO-03") %>% 
  select(d0,d1,d2,d3,d4)


dend_all %>% filter(plot_code=="SAF-05") %>% names

tree_tag, date, wood_density, species, baseline_dbh, dates_dendro_replaced


tmp_a <- read.csv(file = "data/Census_Santarem_2014_2016.csv", skip = 1)

tmp_b <- read_csv(file = "data/Census_Santarem_2014_2016.csv", skip = 1)
tmp_b <- tmp_b %>% rename(tree_tag=B357_T9_N1, dbh=`23.8`) %>% 
  select(tree_tag, dbh)

abline(0,1,col="red")

dend_tmp <- read_csv("data/dendro_STM.csv") %>% mutate(dendrometer_reading_mm=as.numeric(dendrometer_reading_mm), 
                    dendrometer_reading_replaced_mm=as.numeric(dendrometer_reading_replaced_mm))
dend_tmp$tree_tag

junk <- left_join(dend_tmp, tmp_b, by="tree_tag")

junk %>% filter(is.na(dendrometer_reading_mm)==F) %>% 
  filter(plot_code=="STO-03") #798 obs
junk %>% filter(is.na(dendrometer_reading_mm)==F) %>% 
  filter(plot_code=="STO-07") #1003 obs
junk %>% filter(is.na(dendrometer_reading_mm)==F) %>% 
  filter(plot_code=="STO-06") #996 obs

junk %>% filter(is.na(dendrometer_reading_mm)==F) %>% 
  filter(plot_code=="STO-03") %>% distinct(year,month,day)
junk %>% filter(is.na(dendrometer_reading_mm)==F) %>% 
  filter(plot_code=="STO-06") %>% distinct(year,month,day)
junk %>% filter(is.na(dendrometer_reading_mm)==F) %>% 
  filter(plot_code=="STO-07") %>% distinct(year,month,day)


Bacia	Transecto	Trans code	Land-use	Data 2014	Data 2016	Parcela	N�	Tree code	Coletada 2014	Fam�lia	G�nero	Esp�cie	Subsp�cie	G�nero/SP/SPP	Wood Density	POM 2014	POM 2016	Forma	DAP 2014 (average)	DAP 2016 (average)	Altura 2014 CORRIGIDA	Altura2016	Deco2014	Deco2016	Observa��es_campo2014	Observa��es_ID2014	Observa��es_campo2016	Observa��es_ID2016	Excluir de an�lises?	Refer�ncia 2014


census_all$plot_code %>% table

# CAX-01 
# 1715   2327   2224   2500   2527   3919   3163   1479   1546   1608   2674 
# CAX-02 CAX-06 CAX-08
# 

dend_all$plot_code %>% table

dend_all %>% filter(plot_code=="CAX-03") %>% pull(tree_tag) %>% table()
dend_all %>% filter(plot_code=="CAX-03") %>% pull(tree_tag) %>% table()

dend_all %>% filter(plot_code=="CAX-03") %>% View


#--- NXV-01 ----------------------------------------------------------------------
dend %>% filter(tree_tag=="1017") %>% 
  ggplot(data=., aes(date, current_dbh))+geom_path()
dend$current_dbh %>% max
npp_tree %>% ggplot(data=., aes(dbh,nppacw_Mg_tree_day))+geom_point()
npp_tree %>% group_by(date) %>% 
  summarize(wu=weighted.mean(nppacw_Mg_tree_day, w=weight))

# (n_trees_first_census*tmp_npp_dend$npp_weightedAvgTrees_day_dend
plot((550*tmp_npp_dend$npp_weightedAvgTrees_day_dend*365) ~ 
         tmp_npp_dend$date, col="darkgreen",lwd=3, type='l', 
       ylab="Annual agC Mg change rate", xlab="Date", ylim=c(0,10), 
       main=census$plot_code %>% unique)
  
census$date %>% unique %>% sort
census %>% filter(date==min(date)) %>% pull(tree_tag) %>% unique %>% length
census %>% filter(date==min(date)) %>% filter(dbh>=100) %>% 
  pull(tree_tag) %>% unique %>% length
census %>% filter(date==min(date)) %>% #filter(dbh>=100) %>% 
  pull(tree_tag) %>% unique %>% length

census %>% filter(date==min(date)) %>% filter(dbh>=100) %>% 
  pull(tree_tag) %>% unique %>% length

jnk <- read_csv("data/NXV01_FP.csv")
jnk %>% filter(`Census Date`==min(`Census Date`)) %>% 
  filter(D0>=100)
jnk %>% filter(`Census Date`==max(`Census Date`)) %>% 
  filter(D0>=100)

jnk %>% filter(D0>=100) %>% 
  rowwise() %>% 
  mutate(ba_m2 = (pi*(D0/2)**2)/(100**2)) %>% 
  ungroup() %>% 
  group_by(`Census Date`) %>%
  summarize(sum_ba_m2 = sum(ba_m2,na.rm=T))

census_all %>% 
  filter(d0>=100) %>% 
  rowwise() %>% 
  mutate(ba_m2 = (3.141593*(d0/2)^2)) %>% 
  ungroup() %>% 
  group_by(plot_code, census_date) %>%
  summarize(sum_ba_m2 = sum(ba_m2,na.rm=T)*1e-4) %>% View


npp_tree %>% group_by(plot_code, year, month) %>% 
  summarize(npp_avgtrees_day_dend=mean(nppacw_Mg_tree_day,na.rm=T), 
            npp_avgtrees_day_dend_sd=sd(nppacw_Mg_tree_day, na.rm=T), 
            npp_weightedAvgTrees_day_dend = weighted.mean(nppacw_Mg_tree_day, w=weight, na.rm=T),
            npp_medianAvgTrees_day_dend=median(nppacw_Mg_tree_day,na.rm=T),
            nobs=n()) %>% 
  mutate(date=parse_date_time(paste(year,month,15), "ymd"), 
         npp_weightedAvgTrees_month_dend = npp_weightedAvgTrees_day_dend*30.4, 
         npp_weightedAvgTrees_yr_dend=npp_weightedAvgTrees_month_dend*12, 
         npp_medianAvgTrees_month_dend = npp_medianAvgTrees_day_dend*30.4, 
         npp_medianAvgTrees_yr_dend=npp_medianAvgTrees_month_dend*12) %>%
  mutate(mgC_ha_month=npp_weightedAvgTrees_month_dend*nobs) %>% 
  ggplot(data=., aes(date, mgC_ha_month))+geom_path()

(census_all %>% filter(plot_code=="NXV-02") %>% pull(tree_tag)) %in%
(dend_all %>% filter(plot_code=="NXV-02") %>% pull(tree_tag) %>% as.character) %>% table

plot(census$d0~census$d3)

census %>% pull(tree_tag) %>% unique
census %>% pull(date) %>% unique





out_df %>% ggplot(data=., aes(date, npp_lm_agC_month))+geom_line()
npp_tree %>% filter(tag %in% c("108","118","18", "436","515","61")) %>% 
  ggplot(data=., aes(date, current_dbh, color=tree_tag))+geom_line()+geom_point()
npp_tree %>% filter(tag %in% c("108")) %>% 
  ggplot(data=., aes(date, current_dbh, color=tree_tag))+geom_line()+geom_point()
npp_tree %>% filter(tag %in% c("515")) %>% 
  ggplot(data=., aes(date, current_dbh, color=tree_tag))+geom_line()+geom_point()

npp_tree %>% filter(year>=2014) %>% 
  ggplot(data=., aes(date, current_dbh, color=tree_tag))+geom_line()+
  theme(legend.position = "none")

npp_tree %>% filter(year>=2014 & current_dbh>=500) %>% 
  ggplot(data=., aes(date, current_dbh, color=tree_tag))+geom_line()


npp_tree %>% filter(tag %in% c("159")) %>% 
  ggplot(data=., aes(date, current_dbh, color=tree_tag))+geom_line()+geom_point()
npp_tree %>% filter(tag %in% c("1")) %>% 
  ggplot(data=., aes(date, dateDiff, color=tree_tag))+geom_line()+geom_point()

dend %>% names()
dend %>% group_by(tree_tag) %>%
  arrange(date) %>% 
  mutate(prior_date = lag(date, order_by = date)) %>%
  filter(is.na(prior_date)==F) %>% 
  mutate(date_diff = (as.numeric(date)-as.numeric(prior_date))) %>% pull(date_diff) %>% summary

dend %>% group_by(tree_tag) %>%
  arrange(date) %>% 
  mutate(prior_date = lag(date, order_by = date)) %>%
  filter(is.na(prior_date)==F) %>% 
  mutate(date_diff = (as.numeric(date)-as.numeric(prior_date))) %>% 
  filter(date_diff==0)

dend_all %>% mutate(site=substr(plot_code, 1,3)) %>% 
  filter(site %in% c("ANK","BOB","KOG")) %>% 
  write_csv(., path = "outputs/dendro_Ghana_20180208.csv")

dend_out %>% filter(nobs>=10) %>% 
  ggplot(data=., aes(date, npp_lm_agC_month))+geom_line()



npp_tree %>% filter(tag %in% c("159")) %>% 
  ggplot(data=., aes(date, current_dbh, color=tree_tag))+geom_line()+geom_point()
dend %>% filter(tree_tag %in% c("159")) %>% 
  ggplot(data=., aes(date, delta1_std, color=tree_tag))+geom_line()+geom_point()
dend %>% filter(tree_tag %in% c("159")) %>% pull(delta1) %>% plot
dend %>% filter(tree_tag %in% c("159")) %>% pull(baseline_dbh) %>% points(col="red")
dend %>% filter(tree_tag %in% c("159")) %>% pull(band_reset) %>% plot(col="blue")



dend %>% #filter(tree_tag %in% c("159")) %>% 
  ggplot(data=., aes(date, delta1_std, color=tree_tag))+geom_line()+geom_point()+
  theme(legend.position = 'none')

dend %>% 
  ggplot(data=., aes(date, dendrometer_reading_mm,color=tree_tag))+geom_line()+
  theme(legend.position = "none")

dend %>% 
  ggplot(data=., aes(date, dendrometer_reading_replaced_mm,color=tree_tag))+geom_line()+
  theme(legend.position = "none")

out_df %>% 
  ggplot(data=., aes(date, npp_lm_agC_month))+geom_line()+geom_point()

dend %>% filter(tree_tag == "159") %>% 
  ggplot(data=., aes(date, height_m))+geom_line()



dend3 %>% 
  # filter(current_dbh>=250 & current_dbh<=400) %>% 
  filter(current_dbh<500 &current_dbh>=400) %>% 
  # filter(year>=2015) %>% 
  # filter(tree_tag=="159") %>% 
  ggplot(data=., aes(date, dendrometer_reading_replaced_mm, color=tree_tag))+geom_line()+geom_point()+
  theme(legend.position = "none")

dend3 %>% filter(is.na(dendrometer_reading_replaced_mm)==F & is.na(dendrometer_reading_mm)==T)
dend3 %>% filter(is.na(dendrometer_reading_replaced_mm)==F & is.na(dendrometer_reading_mm)==F)
dend3 %>% filter(is.na(dendrometer_reading_replaced_mm)==T & is.na(dendrometer_reading_mm)==F)




dend2 %>% ungroup() %>%  group_by(tree_tag) %>% arrange(date) %>% 
  mutate(baseline_max = cummax(baseline_dbh)) %>% 
  filter(tree_tag=="101") %>% 
  # ggplot(data=., aes(date, dbh_max))+geom_line()
  ggplot(data=., aes(date, baseline_dbh2))+geom_line()
# ggplot(data=., aes(date, delta1_std))+geom_line()


mutate(band_reset_1 = ifelse(delta1 < -3, T, F)) %>%   
  mutate(baseline_dbh = ifelse(delta1 < -3, lag(current_dbh,n=1), baseline_dbh)) %>% # judgement call: how many mm before call it a new baseline?
  mutate(band_reset_2 = ifelse(is.na(baseline_dbh)==T, T, F)) %>%   
  mutate(baseline_dbh = ifelse(is.na(baseline_dbh)==T, lead(baseline_dbh), baseline_dbh)) %>%
  mutate(band_reset_3 = ifelse(delta1_std > 0.1, T,F)) %>%   
  mutate(baseline_dbh = ifelse(delta1_std > 0.1, 
                               est_dbh_from_dendro(dbh=baseline_dbh, dendrometer_reading_mm = as.numeric(dendrometer_reading_mm)), 
                               baseline_dbh)) %>% 
  mutate(baseline_dbh = RcppRoll::roll_max(baseline_dbh)) %>% 
  mutate(current_dbh = est_dbh_from_dendro(dbh=baseline_dbh, dendrometer_reading_mm = dendrometer_reading_mm)) %>% 
  mutate(current_dbh = ifelse(current_dbh<(lag(current_dbh,1,order_by = date)+5),    #added
                              lag(current_dbh,n = 1, order_by=date),                 #
                              current_dbh)) %>%                                      #
  # mutate(pos_meas_err = ifelse((current_dbh-lag(current_dbh,n=1)) > 40, T, F)) %>%  # This removes readings of more than 100mm change
  # mutate(pos_meas_err = ifelse((current_dbh-lag(current_dbh,n=1)) < -5, T, F)) %>% 
  mutate(delta1 = current_dbh-lag(current_dbh,order_by = date)) %>% # added - recompute current_dbh
  mutate(delta1_std = delta1/current_dbh)  # added



#-------------------------------------------------------------------------------------------------------------------
# NPP of the plot from the censuses --------------------------------------------------------------------------------
# The monthly npp_agC from the dendrometer extrapolation NEEDS TO BE REALITY CHECKED AGAINST THE CENSUS VALS !!!
# Steps: filter for only live trees, calc biomass, convert to agC_Mg
tmp_npp_census_trees <- census %>% filter(dbh>=90) %>% 
  # filter(f2=="1") %>%  # Comment this out if you want to include dead trees
  rowwise() %>% 
  mutate(biomass_tree_kg=sum(Chave2014_eq4(diam_cm=dbh/10, density=density, height_m=height_m), na.rm=T), 
         basal_area_tree_cm=sum(pi*(0.1*dbh/2)^2, na.rm=T)) %>% 
  mutate(agC_Mg = (1/(2.1097*1000))*biomass_tree_kg) %>% ungroup()  

tmp_npp_census <- tmp_npp_census_trees %>% group_by(year,month) %>% 
  summarize(n_trees=n(), 
            basal_area_cm2 = sum(basal_area_tree_cm,na.rm=T), 
            biomass_kg = sum(biomass_tree_kg, na.rm=T)) %>% 
  mutate(agC_Mg = (1/(2.1097*1000))*biomass_kg, 
         date=parse_date_time(paste(year,month,15), "ymd"))

tmp_npp_census$day_interval <- c(NA,get_time_diffs(tmp_npp_census$date))
tmp_npp_census$agC_Mg_diff <- c(NA, diff(tmp_npp_census$agC_Mg))
tmp_npp_census$agC_Mg_day <- tmp_npp_census$agC_Mg_diff/tmp_npp_census$day_interval

# calculate the daily accumulation rate of Mg C between the first and last census dates
n_days_first_last_census <- get_time_diffs(c(min(tmp_npp_census$date), max(tmp_npp_census$date)))
first_census_agC_Mg <- tmp_npp_census %>% filter(date==min(tmp_npp_census$date)) %>% pull(agC_Mg)
last_census_agC_Mg <- tmp_npp_census %>% filter(date==max(tmp_npp_census$date)) %>% pull(agC_Mg)
overall_census_agC_Mg_daily <- (last_census_agC_Mg - first_census_agC_Mg)/(n_days_first_last_census)


#--- DIAGNOSTIC CHECKS FOR MATCHING BIOMASS CHANGE IN DENDROMETERS WITH THAT FROM CENSUS ------------------# 
n_trees_first_census <- tmp_npp_census %>% filter(date==min(tmp_npp_census$date)) %>% pull(n_trees)
n_trees_last_census <- tmp_npp_census %>% filter(date==max(tmp_npp_census$date)) %>% pull(n_trees)
overall_dendro_naive_agC_Mg_daily <- mean(c(n_trees_first_census, n_trees_last_census),na.rm=T)*
  (tmp_npp_dend %>% pull(npp_avgtrees_day_dend) %>% mean(na.rm=T))
print(paste("Overall Census agC Mg annual accumulation rate: ", overall_census_agC_Mg_daily*365))

# --------------------------------------------------------------------------------------------------------------# 


census_all %>% names %>% sort
census_all %>% filter(substr(plotname,1,3)=="KEN") %>% pull(height)

census_all %>% filter(plot_code=="KEN-01") %>% pull(height) %>% is.na %>% table
census_all %>% filter(plot_code=="KEN-02") %>% pull(height) %>% is.na %>% table


dend_all %>% filter(is.na(date)==T) %>% pull(year) %>% is.na() %>% table()
dend_all %>% filter(is.na(date)==T) %>% pull(month) %>% is.na() %>% table()
dend_all %>% filter(is.na(date)==T) %>% pull(day) %>% is.na() %>% table()
dend_all %>% filter(is.na(date)==T) %>%
  mutate(date=parse_date_time(paste(year,month,day), "ymd"))

census$dbh %>% is.na %>% table
census$dbh %>% is.infinite() %>% table
census$dbh %>% is.infinite() %>% table
tmp_missing_height$dbh
tmp_missing_height$d0
tmp_missing_height$d1
tmp_missing_height$d2
tmp_missing_height$d3

" Problems: 
- dbh in census is coming back as -Inf
- negative growth predictions in SAF
"

dend %>% ggplot(data=., aes(date, current_dbh,color=tree_tag))+geom_line()+
  geom_point()+theme(legend.position = "none")
dend %>% ggplot(data=., aes(date, dendrometer_reading_mm,color=dbh))+geom_point()+
  viridis::scale_color_viridis()#+
  theme(legend.position = "none")

out_df %>% ggplot(data=., aes(date, npp_lm_agC_month))+geom_line() 
out_df %>% ggplot(data=., aes(date, npp_u_agC_month))+geom_line() 

curve(predict(tmp_fit_lm, type="response", newdata=data.frame(dbh=x, wd=0.65, height_pred=16)),100,1000); 
abline(h=0)

npp_tree %>% filter(decDate==vec_dend_date[idx]) %>% filter(dbh>=90) %>% 
  ggplot(data=., aes(dbh, nppacw_Mg_tree_day))+geom_point()

tmp_nl_fit <- nls(nppacw_Mg_tree_day~ b0*((1/wd)^b1)*dbh^b2, 
                  data=npp_tree %>% filter(decDate==vec_dend_date[idx]) %>% filter(dbh>=90),
                     start = list(b0=0, b1=0.2, b2=0.1))

tmp_nl_fit <- nls(nppacw_Mg_tree_day~ b0*((1/wd)^b1)*dbh^b2, 
                  data=npp_tree %>% filter(decDate==vec_dend_date[idx]) %>% filter(dbh>=90),
                  start = list(b0=0, b1=0.2, b2=0.1))

curve(5.104e-08*exp(x**1.531/1500), 100,1000, add=T)


#------------------------------------------------------------------------------------------------------------
library(mgcv)
tmp_dat <- npp_tree %>% filter(dbh>=90) %>% 
  mutate(rat_growth = nppacw_Mg_tree_day/agC_Mg) %>% filter(rat_growth>=-0.0015)

idx <- 3
tmp_dat <- npp_tree %>% filter(decDate==vec_dend_date[idx]) %>% filter(dbh>=90) %>% 
  mutate(rat_growth = nppacw_Mg_tree_day/agC_Mg) %>% filter(rat_growth>=-0.0015)
tmp_lm_fit <- lm(nppacw_Mg_tree_day~dbh*wd, data=tmp_dat)
tmp_nl_fit <- nls(nppacw_Mg_tree_day~ b2+b0*exp(((dbh)^b1)/1500), 
    data=tmp_dat,
    start = list(b0=0.00003, b1=1.15, b2=0),
    control=list(maxiter = 100000, warnOnly=T),
    lower=c(0,0,-1e-04), upper=c(0.001, 3,5e-04), algorithm = "port")
# tmp_gam_fit <- gam(nppacw_Mg_tree_day~s(dbh, bs="gp",k=3)+dbh*wd+s(height_pred,bs="gp",k=3),select=T,data=tmp_dat,gamma=1)
tmp_gam_fit <- gam(nppacw_Mg_tree_day~dbh+wd+s(height_pred),select=T,gamma=2,data=tmp_dat, method="ML")
# tmp_gam_fit <- gam(nppacw_Mg_tree_day~s(scale(dbh),scale(wd))+s(scale(height_pred)),select=T,gamma=1,data=npp_tree)

tmp_dat %>% plot(predict(tmp_gam_fit)~nppacw_Mg_tree_day, data=., col="red"); abline(0,1,col="purple")
tmp_dat %>% points(predict(tmp_lm_fit)~nppacw_Mg_tree_day, data=.)
tmp_dat %>% points(predict(tmp_nl_fit)~nppacw_Mg_tree_day, data=., col="blue")

tmp_dat %>% ggplot(data=., aes(nppacw_Mg_tree_day, rat_growth,col=dbh))+geom_point()+
  viridis::scale_color_viridis()

tmp_dat %>% ggplot(data=., aes(date, nppacw_Mg_tree_day,col=dbh,group=as.factor(year)))+geom_point()+
  viridis::scale_color_viridis()+
  scale_x_datetime(date_breaks = "6 months", date_labels="%m")

tmp_dat %>% ggplot(data=., aes(month, nppacw_Mg_tree_day,col=as.factor(year)))+geom_point()+geom_smooth(se=F)+
  viridis::scale_color_viridis(discrete = T)
  # scale_x_datetime(date_breaks = "6 months", date_labels="%m")

plot(nppacw_Mg_tree_day~dbh, data = tmp_dat)
curve(predict(tmp_lm_fit, newdata=data.frame(wd=0.65, dbh=x)), 100,1300, add=T, col="black")
curve(predict(tmp_nl_fit, newdata=data.frame(wd=0.65, dbh=x)), 100,1300, add=T, col="blue")
curve(predict(tmp_gam_fit, newdata=data.frame(wd=0.65, height_pred=15, dbh=x)), 100,1300, add=T, col="red")

predict(tmp_lm_fit, newdata=census_max_basal_area) %>% sum
predict(tmp_nl_fit, newdata=census_max_basal_area) %>% sum
predict(tmp_gam_fit, type="response", newdata=census_max_basal_area) %>% sum


npp_tree %>% ggplot(data=., aes(dbh, nppacw_Mg_tree_day/agC_Mg, group=interaction(year,month)))+geom_point()+
  viridis::scale_color_viridis()


bbmle::AICctab(tmp_lm_fit, tmp_nl_fit, tmp_gam_fit)
summary(tmp_nl_fit)
summary(tmp_fit_lm)

(predict(tmp_nl_fit) %>% sum)/(predict(tmp_fit_lm) %>% sum)

dendrometer %>% filter(year>=2014) %>% 
  mutate(dendrometer_reading_mm = ifelse(year>=2015, dendrometer_reading_mm*10, dendrometer_reading_mm)) %>% 
  ggplot(data=., aes(y=dendrometer_reading_mm, x=as.factor(date)))+geom_violin()

dendrometer %>% filter(year>=2013) %>% 
  mutate(dendrometer_reading_mm = ifelse(year>=2015, dendrometer_reading_mm*10, dendrometer_reading_mm)) %>% 
  ggplot(data=., aes(y=dendrometer_reading_mm, x=as.factor(date)))+geom_boxplot()


if(exists("tmp_fit_nl")==F){
  tmp_fit_nl <- nls(nppacw_Mg_tree_day~ b2*height_pred+b0*exp(((dbh)^b1)/1500), 
      data=tmp_dat,
      start = list(b0=0.00003, b1=1.15, b2=0),
      control=list(maxiter = 100000, warnOnly=T),
      lower=c(0,0,-1e-04), upper=c(0.001, 3,5e-04), algorithm = "port")
}




tmp_nl_fit0 <- nls(nppacw_Mg_tree_day~ b2+b0*exp(((dbh)^b1)/1500), 
                  data=tmp_dat,
                  start = list(b0=0.00003, b1=1.15, b2=0),
                  control=list(maxiter = 100000, warnOnly=T),
                  lower=c(0,0,-1e-04), upper=c(0.001, 3,5e-04), algorithm = "port")
tmp_nl_fit1 <- nls(nppacw_Mg_tree_day~ b2*height_pred+b0*exp(((dbh)^b1)/1500), 
                  data=tmp_dat,
                  start = list(b0=0.00003, b1=1.15, b2=0),
                  control=list(maxiter = 100000, warnOnly=T),
                  lower=c(0,0,-1e-04), upper=c(0.001, 3,5e-04), algorithm = "port")
bbmle::AICctab(tmp_nl_fit0, tmp_nl_fit1)

plot(nppacw_Mg_tree_day~dbh, data = tmp_dat)
curve(predict(tmp_nl_fit0, newdata=data.frame(wd=0.65, dbh=x)), 100,1300, add=T, col="black")
curve(predict(tmp_nl_fit1, newdata=data.frame(wd=0.65, dbh=x, height_pred=30)), 100,1300, add=T, col="blue")
curve(predict(tmp_nl_fit1, newdata=data.frame(wd=0.65, dbh=x, height_pred=50)), 100,1300, add=T, col="blue")


tmp_dat <- npp_tree %>% mutate(rat_growth = nppacw_Mg_tree_day/agC_Mg) %>% filter(rat_growth>=-0.0015)
fit_a <- lm(nppacw_Mg_tree_day~ dbh, data=tmp_dat)
fit_b <- lm(nppacw_Mg_tree_day~ dbh+wd, data=tmp_dat)
fit_c <- lm(nppacw_Mg_tree_day~ dbh*wd, data=tmp_dat)
fit_d <- lm(nppacw_Mg_tree_day~ height_pred+dbh+wd, data=tmp_dat)
fit_e <- lm(nppacw_Mg_tree_day~ height_pred+dbh*wd, data=tmp_dat)
fit_f <- lm(nppacw_Mg_tree_day~ height_pred+dbh+wd+I(dbh**2), data=tmp_dat)
fit_g <- lm(nppacw_Mg_tree_day~ height_pred+I(pi*(dbh/2)^2)+wd+I(dbh**2), data=tmp_dat)
fit_h <- lm(nppacw_Mg_tree_day~ scale(height_pred)+scale(dbh)+scale(wd)+scale(I(dbh**2)), data=tmp_dat)
fit_i <- lm(nppacw_Mg_tree_day~ scale(height_pred)+scale(dbh)+scale(wd)+scale(I(dbh**2)), data=tmp_dat)
fit_j <- lm(nppacw_Mg_tree_day~ scale(height_pred)+scale(dbh)*scale(wd)+scale(I(dbh**2)), data=tmp_dat)
bbmle::AICctab(fit_a, fit_b, fit_c, fit_d, fit_e, fit_f, fit_g, fit_h,fit_i, fit_j)

fit_gam <- gam(nppacw_Mg_tree_day~s(scale(dbh))+wd+s(height_pred)+s(scale(I(dbh**2))),select=T,gamma=2,data=tmp_dat)
summary(fit_gam)
summary(fit_h)


predict(tmp_gam_fit, type="response", newdata=census_max_basal_area) %>% sum
predict(fit_h, type="response", newdata=census_max_basal_area) %>% sum
dim(census_max_basal_area)[1] * (tmp_dat$nppacw_Mg_tree_day %>% mean())

summary(fit_gam) %>% str
summary(fit_gam)$dev.expl

out_df %>% #filter(npp_nl_agC_month<=2.5) %>% 
  ggplot(data=., aes(date, npp_lm_agC_month))+  geom_line()+
  geom_line(aes(date, npp_gam_agC_month),col="darkgreen")+
  geom_line(aes(date, npp_nl_agC_month),col="blue")+
  geom_line(aes(date, npp_u_agC_month), col="red")

out_df$npp_gam_agC_month %>% summary
out_df$npp_u_agC_month %>% summary
out_df$npp_lm_agC_month %>% summary
out_df$npp_nl_agC_month %>% summary

npp_tree %>% names() %>% sort()

npp_tree %>% filter(tree_tag=="159") %>% 
  ggplot(data=., aes(date, current_dbh))+geom_line()
dend %>% filter(tree_tag=="159") %>% 
  ggplot(data=., aes(date, baseline_dbh_max))+geom_line()
dend %>% filter(tree_tag=="159") %>% 
  ggplot(data=., aes(date, dendrometer_reading_mm))+geom_line()
dend %>% filter(tree_tag=="159") %>% 
  ggplot(data=., aes(date, delta1_std))+geom_line()
dend %>% filter(tree_tag=="159") %>% pull(delta1_std) %>% min(na.rm=T)
  


dend %>% filter(year>=2010 & dbh>750) %>%  
  ggplot(data=., aes(date, current_dbh, color=tree_tag))+geom_line()+
  geom_line(aes(date, baseline_dbh), lty=3)

dend2 %>% ggplot(data=., aes(baseline_dbh, baseline_dbh2))+geom_point()+
  geom_abline(aes(intercept=0,slope=1),color="red")




#--------------------------------------------------------------------
vec_dats <- npp_tree$date %>% unique() %>% sort()
tmp_dat <- npp_tree %>% filter(date==as.POSIXct("2013-01-27 00:00:00 GMT")) %>% 
 mutate(rat_growth = nppacw_Mg_tree_day/agC_Mg) %>% filter(rat_growth>=-0.0015)

library(caret)
fit_caret <- train(nppacw_Mg_tree_day~ ., 
                   data=tmp_dat %>% select(nppacw_Mg_tree_day, wd, dbh, height_pred), 
                   method="glmnet")

fit_lm <- lm(nppacw_Mg_tree_day~wd+dbh+height_pred, data=tmp_dat)       
summary(fit_lm)

fit_lm <- lm(nppacw_Mg_tree_day~wd+I(exp(dbh))+height_pred, data=tmp_dat)       
summary(fit_lm)
predict(fit_lm, )

fit_gam <- gam(nppacw_Mg_tree_day~s(wd)+s(dbh)+s(height_pred)+s(I(dbh*wd)), 
               data=tmp_dat, 
               select = T)       
summary(fit_gam)
plot(fit_gam)
predict(fit_gam) %>% sum()

fit_gam2 <- gam(nppacw_Mg_tree_day~s(wd,dbh)+s(dbh, height_pred, k=3), 
               data=tmp_dat, 
               select = T, 
               method="ML")       
summary(fit_gam2)
plot(fit_gam2)
predict(fit_gam2) %>% hist
predict(fit_gam2) %>% sum

fit_gam3 <- gam(nppacw_Mg_tree_day~ti(wd,dbh,height_pred, k=3), 
                data=tmp_dat, 
                select = T, 
                method="ML")       
summary(fit_gam3)
vis.gam(fit_gam3,view = c("dbh","wd"), theta=-30)
vis.gam(fit_gam3,view = c("height_pred","wd"), theta=120)
vis.gam(fit_gam3,view = c("height_pred","dbh"), theta=120)

predict(fit_gam, type="response", newdata=census_max_basal_area) %>% sum
predict(fit_gam2, type="response", newdata=census_max_basal_area) %>% sum
predict(fit_gam3, type="response", newdata=census_max_basal_area) %>% sum
predict(fit_gam4, type="response", newdata=census_max_basal_area) %>% sum


fit_gam4 <- gam(nppacw_Mg_tree_day~ti(wd,dbh, by=month)+s(height_pred, k=3, by=month), 
                data=npp_tree %>% mutate(month=as.factor(month)), 
                select = T,
                method="ML")       
summary(fit_gam4)
vis.gam(fit_gam4,view = c("dbh","wd"), theta=-30)
vis.gam(fit_gam4,view = c("height_pred","wd"), theta=120)
vis.gam(fit_gam4,view = c("height_pred","dbh"), theta=120)

out_df %>% #filter(np_agC_month<=2.5) %>% 
  ggplot(data=., aes(date, npp_lm_agC_month))+  geom_line()+
  geom_line(aes(date, npp_gam_agC_month),col="darkgreen")+
  geom_line(aes(date, npp_hybrid_bestMod_agC_month),col="blue")+
  geom_line(aes(date, npp_u_agC_month), col="red")+
  geom_abline(aes(intercept=0,slope=0),col="grey",lwd=2)



dend_out %>% filter(npp_nl_agC_month<3) %>% 
  ggplot(data=., aes(date, npp_gam_agC_month))+geom_path(col="darkgreen")+
  geom_path(aes(date, npp_u_agC_month),col="red")+
  geom_path(aes(date, npp_nl_agC_month),col="blue")+
  geom_abline(aes(intercept=0,slope=0,col="red"))

dendrometer$tree_tag %>% unique() %>% length()

vec_tags_ss

dend %>% ggplot(data=., aes(date, current_dbh, group=tree_tag))+geom_line()

dend1 %>% ggplot(data=., aes(date, dbh_max, group=tree_tag))+geom_line()

dend %>% ggplot(data=., aes(date, current_dbh, color=tree_tag))+geom_line()+
  theme(legend.position="none")

# --- THIS FORM SEEMS TO WORK --- 
gam(nppacw_Mg_tree_day~ti(wd,dbh, k=3)+s(height_pred,k=2),
    data=tmp_dat,
    select = T,
    method="ML")


[1:40]
plot(dend_out_ss$npp_gam_agC_month[-41]~dend_out$npp_u_agC_month[-41],cex=1.3,col="darkgreen",pch=20); abline(0,1,col="red")
points(dend_out_ss$npp_lm_agC_month[-41]~dend_out$npp_u_agC_month[-41], col="blue",pch=20);

lm(dend_out_ss$npp_gam_agC_month[-41]~dend_out$npp_u_agC_month[-41]) %>% summary
lm(dend_out_ss$npp_lm_agC_month[-41]~dend_out$npp_u_agC_month[-41]) %>% summary

# ti: 0.8232



data.frame(plot_code=aa, tag=bb, year=cc, month=dd, 
           agC_Mg=ee, agC_Mg_diff=ff, dateDiff=gg, current_dbh=vec_current_dbh,
           tree_tag=vec_tag, 
           # pos_meas_error=vec_pmer, 
           canopyLayer_pred = vec_canopyLayer,
           height_pred = vec_height,
           wd = vec_wd, 
           date=vec_dates)
length(aa)
length(gg)
length(vec_current_dbh)
length(vec_dates)


dend_all$plot_code %>% unique %>% sort
census_all$plot_code %>% unique %>% sort

tmp_a <- gam(nppacw_Mg_tree_day~wd+s(dbh)+s(height_pred,k=2),
                    data=tmp_dat,
                    select = T,
                    method="ML")
tmp_b <- gam(I(nppacw_Mg_tree_day+1)~wd*dbh+s(height_pred,k=2),
             data=tmp_dat,
             family=Gamma(link="inverse"),
             select = T,
             method="ML")

plot((predict(tmp_b, type="response")-1)~predict(tmp_a, type="response")); abline(0,1,col="red")

abs(sum(residuals(tmp_a, type="response")**2))
abs(sum(residuals(tmp_b, type="response")**2))

census %>% names
census %>% group_by(census_date) %>% summarize(nobs=n())

census_all %>% mutate(site=substr(plot_code,1,3)) %>% 
  filter(site%in%c("ANK","BOB","CAX","KEN","KOG","LPG","MLA","MNG","SAF","TAM")) %>%
  filter(d0 >= 100) %>% 
  group_by(plot_code, census_date) %>% summarize(nobs=n()) %>% 
  ggplot(data=., aes(census_date, nobs, color=plot_code))+geom_line()

# problem with the stem count on one of the NXV plots 

census_all %>% filter(plot_code==this_plot_code) %>% dim
census_all %>% filter(plot_code==this_plot_code) %>%
  distinct() %>% dim

dendrometer
census

"CAX-01" "CAX-02" "CAX-04" "CAX-06" "CAX-08"
census_all %>% filter(plot_code=="CAX-04")


############################################################################
# --- CAX census height ------
caxh <- read_csv("data/CAX_heights_fromELDSscavenge.csv")
caxh <- caxh %>% filter(family_abbrev!="morta")
fam_names <- census_all$family %>% unique() %>% sort() %>% tolower()
cax_fams <- caxh$family_abbrev %>% tolower()
caxh$family <- fam_names[pmatch(cax_fams, fam_names, duplicates.ok = T)]
caxh$genus <- stringr::str_split_fixed(caxh$species, " ", n=3)[,1]
caxh$continent <- "South America"
caxh_wd <- find_wd(caxh %>% select(family, genus, species, continent))
cax_heights <- caxh_wd
cax_heights$dbh <- caxh$dbh %>% as.numeric()
cax_heights$height <- caxh$alt_t %>% as.numeric()
cax_heights$dbh <- cax_heights$dbh*10
cax_heights %>% ggplot(data=., aes(dbh, height))+geom_point()+geom_smooth()
write_csv(cax_heights, "data/CAX_heights_cleaned.csv")
############################################################################

unique(dend$tree_tag)[!(unique(dend$tree_tag) %in% unique(census$tree_tag))]

unique(dend$tree_tag)[!(unique(dend$tree_tag) %in% unique(census_all %>% filter(plot_code=="CAX-01") %>% pull(tree_tag)))]
unique(dend$tree_tag)[!(unique(dend$tree_tag) %in% unique(census_all %>% filter(plot_code=="CAX-02") %>% pull(tree_tag)))]
unique(dend$tree_tag)[!(unique(dend$tree_tag) %in% unique(census_all %>% filter(plot_code=="CAX-03") %>% pull(tree_tag)))]
unique(dend$tree_tag)[!(unique(dend$tree_tag) %in% unique(census_all %>% filter(plot_code=="CAX-04") %>% pull(tree_tag)))]
unique(dend$tree_tag)[!(unique(dend$tree_tag) %in% unique(census_all %>% filter(plot_code=="CAX-06") %>% pull(tree_tag)))]
unique(dend$tree_tag)[!(unique(dend$tree_tag) %in% unique(census_all %>% filter(plot_code=="CAX-08") %>% pull(tree_tag)))]

tmp_vec_tags <- dend %>% group_by(tree_tag) %>%
  summarize(dbh=max(current_dbh)) %>%
  arrange(desc(dbh)) %>% pull(tree_tag) %>% unique()
dend %>% filter(tree_tag %in% tmp_vec_tags[4]) %>% 
  ggplot(data=., aes(date, current_dbh,color=tree_tag))+
  geom_point()+geom_line()+
  theme(legend.position = " none")

dend %>% filter(tree_tag %in% tmp_vec_tags[4]) %>% 
  ggplot(data=., aes(date, dendrometer_reading_mm,color=tree_tag))+
  geom_point()+geom_line()+
  theme(legend.position = " none")

dend %>% filter(tree_tag=="60") %>% 
 ggplot(data=., aes(date, dendrometer_reading_mm))+geom_line()

npp_tree %>% filter(tag=="60") %>% 
  ggplot(data=., aes(date, agC_Mg, group=tag))+geom_line()

dend_all %>% filter(plot_code=="CAX-04") %>% pull(dendrometer_reading_mm) %>% summary

dendrometer %>% group_by(tree_tag) %>% 
  summarize(sum_diam=sum(dendrometer_reading_mm, na.rm=T)*100) %>% 
  pull(sum_diam) %>% hist(50)


#--- decrypting the super super super super super messed up caixuana data formatting ---
dat <- read_csv("../../../../Downloads/dend_cax_24Oct17_3.csv")
dat$tree_tag <- as.character(dat$tree_tag)
vec_tags <- dat %>% group_by(tree_tag) %>% summarize(u=mean(dap_jan_2013, na.rm=T)) %>% 
  arrange(desc(u)) %>% pull(tree_tag)

dat %>% 
  mutate(date=as.POSIXct(date, tz="UTC")) %>%
  filter(is.na(value)==F) %>% 
  # filter(tree_tag %in% vec_tags[100:130]) %>%
  group_by(tree_tag) %>%
  arrange(date) %>% 
  mutate(value=value*900) %>% 
  mutate(dbh_increment = cumsum(value)) %>% 
  ungroup() %>% 
  # select(date, dbh_increment, value) %>% 
  ggplot(data=., aes(date, dbh_increment, color=tree_tag))+
  geom_point()+
  geom_line()+
  # geom_line(aes(date, value*900),col="red")+
  # geom_smooth(method='lm', se=F)+
  scale_x_datetime(date_breaks="2 years")+
  theme(legend.position = "none")

dat %>% 
  mutate(date=as.POSIXct(date, tz="UTC")) %>%
  filter(is.na(value)==F) %>% 
  # filter(tree_tag %in% vec_tags[100:130]) %>%
  group_by(tree_tag) %>%
  arrange(date) %>% 
  mutate(value=value*900) %>% 
  mutate(dbh_increment = cumsum(value)) %>% 
  ungroup() %>% 
  group_by(date) %>% 
  summarize(tot=sum(dbh_increment, na.rm=T)) %>% 
  ggplot(data=., aes(date, tot))+
  geom_point()+
  geom_line()+
  # geom_line(aes(date, value*900),col="red")+
  # geom_smooth(method='lm', se=F)+
  scale_x_datetime(date_breaks="2 years")+
  theme(legend.position = "none")

dat2 <- dat %>% 
  mutate(date=as.POSIXct(date, tz="UTC")) %>%
  filter(is.na(value)==F) %>% 
  group_by(tree_tag) %>%
  arrange(date) %>% 
  mutate(value=value*900) %>% 
  mutate(dbh_increment = cumsum(value)) %>% 
  ungroup() 

dat2$wd %>% mean(na.rm=T)

tmp_dat <- npp_tree %>% filter(date==vec_dend_date[idx]) %>%
  mutate(rat_growth = nppacw_Mg_tree_day/agC_Mg) %>% 
  filter(nppacw_Mg_tree_day!=0)
samp_idx <- sample.int(dim(tmp_dat)[1], dim(tmp_dat)[1]/5)

fit <- caret::train(nppacw_Mg_tree_day~ height_pred+dbh*wd+I(dbh**2),
             method="gbm", data=tmp_dat[samp_idx, ],
             trControl=trainControl(method="cv"))
predict(fit, newdata=tmp_census) %>% sum
predict(tmp_fit_gam1, newdata=tmp_census) %>% sum

tmp_dat$nppacw_Mg_tree_day %>% sum()


fit <- caret::train(nppacw_Mg_tree_day~.,
                    method="earth", 
                    data=tmp_dat%>% 
                      select(nppacw_Mg_tree_day,wd,dbh,height_pred),
                    trControl=trainControl(method="repeatedcv", repeats = 10))
summary(fit$finalModel)
predict(fit, newdata=tmp_census[1:323,]) %>% sum
mean(tmp_dat[samp_idx,]$nppacw_Mg_tree_day, na.rm=T)*323
mean(tmp_dat$nppacw_Mg_tree_day, na.rm=T)*323





fit <- gam(nppacw_Mg_tree_day~te(I(dbh*wd),height_pred),
           select=T, 
           data=tmp_dat[samp_idx,])
predict(fit, newdata=tmp_census[1:441,]) %>% sum

tmp_vec <- npp_tree %>% filter(nppacw_Mg_tree_day==0) %>% pull(tree_tag) %>% unique()
  group_by(date) %>% 
  summarize(nobs=n()) %>% 
  ggplot(data=., aes(date, nobs))+geom_line()

npp_tree %>% 
  filter(tree_tag %in% tmp_vec[300:350]) %>% 
  ggplot(data=., aes(date, nppacw_Mg_tree_day)) + geom_line()+
  theme(legend.position = "none")+
  facet_wrap(~tree_tag, scales = "free")

#possibly dead: 112,173, 288, 281, 342, 454, 469

npp_tree %>% group_by(date) %>% summarize(u=mean(dateDiff)) %>% 
  ggplot(data=., aes(date, u))+geom_line()

# dates with 180 day diff
2007-10-01
2013-04-01
2015-04-01
2015-10-01
2016-10-01

census %>% group_by(date) %>% summarize(nobs=n())

census %>% filter(plot_code=="NXV-02")
census %>% pull(d0) %>% hist

out_df %>% ggplot(data=., aes(date, npp_u_agC_month))+geom_line(lwd=3)+
  geom_line(aes(date, npp_hybrid_bestMod_agC_month), col="blue")+
  geom_line(aes(date, npp_bestEst_agC_month), col="purple")+
 geom_line(aes(date, npp_hybrid_gam_agC_month), col="darkgreen")+
 geom_line(aes(date, npp_lm_agC_month), col="red")

census_all %>% filter(plot_code=="CAX-04") %>% 
  filter(d0>=100) %>% 
  group_by(census_date) %>% 
  summarize(nobs=n()) %>% 
  ggplot(data=., aes(census_date, nobs))+geom_line()

census_all %>% filter(plot_code=="NXV-02") %>% 
  pull(census_date) %>% unique

census_all %>% filter(plot_code=="NXV-02") %>% View
  filter(census_date==2011.255) %>% 
  pull(d1) 
  

dend_all %>% filter(plot_code=="NXV-02") %>% 
  group_by(date) %>% summarize(nobs=n())


dat2 %>% ggplot(data=., aes(date, dbh_increment, color=tree_tag))+geom_point()

dat %>% ggplot(data=., aes(date, dendrometer_reading_mm_cum))+
  geom_point()+geom_line()+
  facet_wrap(~tree_tag, scales = "free_y")

dend %>% ggplot(data=., aes(date, current_dbh, color=tree_tag))+
  geom_point()+geom_line()+
  facet_wrap(~tree_tag, scales = "free_y")+theme(legend.position = "none")

dend$height_pred %>% hist

dat$dendrometer_reading_mm %>% hist(50)
dat$dendrometer_reading_mm %>% hist(50)

out_df %>% 
  ggplot(data=., aes(date, npp_wu_agC_month))+geom_line(color="orange")+
  geom_line(aes(date, npp_u_agC_month), lwd=3)+
  geom_line(aes(date, npp_hybrid_bestMod_agC_month), col="blue", lwd=2)+
  geom_line(aes(date, npp_bestEst_agC_month), col="purple")+
  geom_line(aes(date, npp_hybrid_gam_agC_month), col="darkgreen")+
  geom_line(aes(date, npp_lm_agC_month), col="red")



dend_out %>% 
  ggplot(data=., aes(date, npp_wu_agC_month))+geom_line(color="orange")+
  geom_line(aes(date, npp_u_agC_month), data=dend_out, lwd=3)+
  geom_line(aes(date, npp_hybrid_bestMod_agC_month), col="blue", lwd=2)+
  geom_line(aes(date, npp_bestEst_agC_month), col="purple")+
  geom_line(aes(date, npp_hybrid_gam_agC_month), col="darkgreen")+
  geom_line(aes(date, npp_lm_agC_month), col="red")

census_all %>% filter(plot_code=="NXV-02") %>% pull(census_date) %>% unique() %>% sort()

census_all  <- read_csv("data/forestplots20180320.csv", na=c("NA", "NaN", "")) # Census_Santarem_2014_2016.csv
names(census_all) <- tolower(names(census_all))
names(census_all) <- gsub(pattern=" ",replacement="_",names(census_all))
census_all <- census_all %>% rename(tree_tag = tag_number) %>% mutate(tree_tag = as.character(tree_tag))

census_all %>% filter(plot_code=="NXV-02") %>% 
  filter(near(census_date, 2013.001, tol=0.01)) %>% 
  filter(d0>=10)
  pull(d0) 


census_all %>% 
  mutate(site=substr(plot_code,1,3)) %>% 
  distinct() %>% 
  group_by(census_date, plot_code,site) %>% 
  filter(d0>=100) %>% 
  summarize(nobs=n()) %>% 
  ggplot(data=., aes(census_date, nobs, color=plot_code))+geom_point()+geom_line()+
  facet_wrap(~site, scales="free") 



census_all %>% dim
census_all %>% distinct() %>% dim

dendrometer$dendrometer_reading_mm %>% hist

dendrometer %>% ggplot(data=., aes(date, dendrometer_reading_mm))+geom_point()+
  geom_smooth(method="lm")+
  theme(legend.position = "none")

census %>% filter(date==first(date)) %>% dim()
dendrometer$tree_tag %>% unique() %>% length()

census$dbh %>% median
census %>% filter(tree_tag %in% unique(dendrometer$tree_tag)) %>% pull(dbh) %>% median()

dend_out %>% 
  ggplot(data=., aes(date, npp_wu_agC_month))+geom_line(color="orange")+geom_point()+
  geom_line(aes(date, npp_u_agC_month), data=dend_out, lwd=3)+
  geom_line(aes(date, npp_hybrid_bestMod_agC_month), col="blue", lwd=2)+
  geom_line(aes(date, npp_bestEst_agC_month), col="purple")+
  geom_line(aes(date, npp_hybrid_gam_agC_month), col="darkgreen")+
  geom_line(aes(date, npp_lm_agC_month), col="red")


plot(residuals.gam(tmp_fit_gam4, type="response")~tmp_dat$dbh)
abline(lm(residuals.gam(tmp_fit_gam4, type="response")~tmp_dat$dbh))

tmp_wfun <- approxfun(density(census_max_basal_area$dbh))
growth_wu_mgc <- weighted.mean(tmp_dat$nppacw_Mg_tree_day,
                               w = tmp_wfun(tmp_dat$dbh),na.rm=T)

plot(tmp_wfun(tmp_dat$dbh)~tmp_dat$dbh)
gam(nppacw_Mg_tree_day~wd+s(dbh)+s(height_pred),
                    data=tmp_dat,
                    select = T,weights = tmp_wfun(tmp_dat$dbh),
                    method="ML") %>% predict(census_max_basal_area) %>% sum()




################################################################################################
# --- use this to debug TAM-06 --- figuring out how to update the dendrobands baseline_dbh
################################################################################################

dend1 <-dend %>% group_by(tree_tag) %>% arrange(date) %>% 
  mutate(delta1 = current_dbh-lag(current_dbh,order_by = date)) %>% 
  mutate(delta1_std = delta1/current_dbh) %>% 
  mutate(dbh_max = cummax(current_dbh)) %>% 
  mutate(flag_new_dband = ifelse(is.na(dendrometer_reading_replaced_mm)==F & lag(is.na(dendrometer_reading_replaced_mm)==T), 1, 0)) %>% 
  mutate(baseline_dbh2 = ifelse(flag_new_dband==1, current_dbh, baseline_dbh))

dend2 <- dend1 %>% ungroup() %>% #rowwise() %>%
  mutate(baseline_dbh = as.double(baseline_dbh)) %>% 
  mutate(baseline_dbh2 = if_else(delta1_std < -0.0075, dbh_max, baseline_dbh)) %>% 
  mutate(baseline_dbh3 = if_else(is.na(dendrometer_reading_replaced_mm)==F, current_dbh, baseline_dbh))

dend3 <- dend2 %>% group_by(tree_tag) %>% arrange(date) %>% 
  mutate(baseline_dbh2 = if_else(is.na(baseline_dbh2)==T, lead(baseline_dbh2), baseline_dbh2)) %>% 
  mutate(baseline_dbh_max = cummax(baseline_dbh2))

dend3$current_dbh <- est_dbh_from_dendro(dbh=dend3$baseline_dbh_max, dendrometer_reading_mm = as.numeric(dend3$dendrometer_reading_mm))

dend <- dend3
dend$thisdbh_cm <- dend$current_dbh/10
dend$agC_Mg <- NA
dend %>% names()
dend$height_pred <- predict(mod_height, newdata=dend, type="response")
rm(dend1, dend2, dend3)
# END FCKNG GOLD!!!----------------------------------------------------------------------------------------------
################################################################################################
# --- end new shit 

pass1 <- dend %>% 
  filter(tree_tag=="51") %>%
  # filter(current_dbh>450 & current_dbh<480) %>% 
  group_by(tree_tag) %>% 
  arrange(date) %>% 
  mutate(delta1 = current_dbh-lag(current_dbh,order_by = date)) %>% 
  mutate(delta1_std = as.double(delta1/current_dbh)) %>%
  mutate(baseline_dbh = if_else(delta1_std < -0.0075, as.double(current_dbh), as.double(baseline_dbh))) %>% 
  mutate(baseline_dbh = if_else(is.na(baseline_dbh)==T, lead(baseline_dbh), baseline_dbh)) %>% 
  mutate(baseline_dbh1 = cummax(baseline_dbh)) %>%
  mutate(current_dbh = est_dbh_from_dendro(dbh=baseline_dbh1, dendrometer_reading_mm = dendrometer_reading_mm)) 
  
pass2 <- pass1 %>% 
  # filter(tree_tag=="512") %>%
  # filter(current_dbh>450 & current_dbh<480) %>% 
  group_by(tree_tag) %>% 
  arrange(date) %>% 
  mutate(delta1 = current_dbh-lag(current_dbh,order_by = date)) %>% 
  mutate(delta1_std = as.double(delta1/current_dbh)) %>%
  mutate(baseline_dbh = if_else(delta1_std < -0.0075, as.double(current_dbh), as.double(baseline_dbh))) %>% 
  mutate(baseline_dbh = if_else(is.na(baseline_dbh)==T, lead(baseline_dbh), baseline_dbh)) %>% 
  mutate(baseline_dbh1 = cummax(baseline_dbh)) %>%
  mutate(current_dbh = est_dbh_from_dendro(dbh=baseline_dbh1, dendrometer_reading_mm = dendrometer_reading_mm)) 

pass2 %>% 
  filter(current_dbh>470 & current_dbh<480) %>% 
  filter(tree_tag=="51") %>% 
  ggplot(data=., aes(date, dendrometer_reading_mm, color=tree_tag))+geom_line()+geom_point()

dend %>% 
  ggplot(data=., aes(date, current_dbh))+geom_point()+geom_smooth()

dend %>% group_by(date) %>% summarize(nobs=n()) %>% 
  ggplot(data=., aes(date, nobs))+geom_line()

dend %>% filter(tree_tag=="114") %>% select(dendrometer_reading_mm, dendrometer_reading_replaced_mm) %>% View
  
dend %>% filter(tree_tag=="114") %>% 
  ggplot(data=., aes(date, dendrometer_reading_mm))+geom_line()+
  geom_point(aes(date, dendrometer_reading_replaced_mm),color="blue")

dend %>% filter(tree_tag=="114") %>% 
  ggplot(data=., aes(date, current_dbh, color=tree_tag))+
  geom_point()+geom_line()+
  geom_line(aes(date, baseline_dbh))+
  geom_line(aes(date, dbh_max),col="black")+
  geom_line(aes(date, baseline_dbh3),col="purple")+
  theme(legend.position = "none")

dend3 %>% filter(tree_tag=="114") %>% 
  ggplot(data=., aes(date, current_dbh, color=tree_tag))+
  geom_point()+geom_line()

dend_out %>% filter(date>="2012-01-01") %>%
  filter(date_diff < 140) %>% 
  ggplot(data=., aes(date, npp_u_agC_month))+geom_line(color="orange")+geom_point()+
  geom_line(aes(date, npp_u_agC_month), lwd=3)+
  geom_line(aes(date, npp_hybrid_bestMod_agC_month), col="blue", lwd=2)+
  geom_line(aes(date, npp_bestEst_agC_month), col="purple")+
  geom_line(aes(date, npp_hybrid_gam_agC_month), col="darkgreen")+
  geom_line(aes(date, npp_lm_agC_month), col="red")


dend_out %>% ggplot(data=., aes(npp_u_agC_month, dev))+geom_point()+geom_smooth(method="lm")

dend %>% filter(date>="2015-01-01") %>% 
  mutate(category=cut(current_dbh, breaks=c(-Inf,150, 200, 300, 500, 750, Inf))) %>% #, labels=c("small","middle","large","massive"))) %>% 
  filter(is.na(category)==F) %>% 
  ggplot(data=., aes(date, current_dbh, color=tree_tag))+
  geom_point()+geom_line()+
  facet_wrap(~category, scales="free", ncol = 2)+
  theme(legend.position = "none")

npp_tree %>% names

npp_tree %>% #filter(date>="2010-01-01") %>% 
  mutate(agC_Mg_diff_std = as.double(agC_Mg_diff/agC_Mg)) %>% 
  filter(agC_Mg_diff_std > -0.05) %>% 
  mutate(category=cut(current_dbh, breaks=c(-Inf,150, 200, 300, 500, 750, Inf))) %>% #, labels=c("small","middle","large","massive"))) %>% 
  filter(is.na(category)==F) %>% 
  filter(dateDiff<=130) %>% 
  ggplot(data=., aes(date, agC_Mg_diff_std, color=tree_tag))+
  geom_point()+geom_line()+
  facet_wrap(~category, scales="free", ncol = 2)+
  theme(legend.position = "none")

npp_tree %>% filter(date>="2015-01-01") %>%
  filter(current_dbh>250 & current_dbh<272) %>%
  arrange(tree_tag) %>% 
  filter(tree_tag=="93") %>% 
  filter(dateDiff<=130) %>% 
  mutate(agC_Mg_diff_std = as.double(agC_Mg_diff/agC_Mg)) %>% 
  # filter(agC_Mg_diff_std > -0.05) %>% 
  ggplot(data=., aes(date, current_dbh, color=tree_tag))+
  geom_point()+geom_line()

dend %>% filter(tree_tag=="93") %>% filter(date>="2015-01-01") %>%
  ggplot(data=., aes(date, current_dbh))+geom_point()+geom_line()


dend[thistree,] %>% 
  ggplot(data=., aes(date, current_dbh))+geom_point()


npp_tree %>% 
  mutate(agC_Mg_diff_std = as.double(agC_Mg_diff/agC_Mg)) %>% 
  pull(agC_Mg_diff_std) %>% 
  hist(200)


library(mgcv)
tmp_fit <- gam(npp_bestEst_agC_month~s(year,k=50), 
               data=dend_out %>% 
                 mutate(month=month(date),year=year(date), date=as.numeric(date)), 
               method="ML")
plot(predict(tmp_fit)~dend_out$date, type="l")

tmp_fit_l <- loess(npp_bestEst_agC_month~date, data=dend_out)

library(RcppRoll)
dend_out <- dend_out %>% mutate(year=year(date), month=month(date))
tmp <- expand.grid(year = seq(min(dend_out$year), max(dend_out$year)), 
                   month=1:12) %>% as.tibble()
tmp <- left_join(tmp, dend_out, by=c("year","month"))

tmp %>% mutate(date = parse_date_time(paste(year,month,15), "ymd")) %>% 
  arrange(date) %>% pull(npp_u_agC_month)
  
tmp %>%
  mutate(date = parse_date_time(paste(year,month,15), "ymd")) %>% 
  arrange(date) %>% 
  mutate(mu3=roll_meanr(npp_bestEst_agC_month,n = 3, fill = NA)) %>% 
  # filter(date_diff<120) %>% 
  plot(npp_bestEst_agC_month~date, data=.,type='l')
tmp %>% 
  mutate(date = parse_date_time(paste(year,month,15), "ymd")) %>% 
  arrange(date) %>% 
  mutate(mu3=roll_meanr(npp_bestEst_agC_month,n = 3, fill = NA, na.rm=T)) %>% 
  # filter(date_diff<120) %>% 
  lines(mu3~date, data=., col="darkgreen")


vec_ran <- rnorm(100, 0, 3)
plot(vec_ran, type='l')
lines(roll_meanl(vec_ran,3),type='l',col="red")

roll_meanr(c(NA,NA,NA,1:10), n = 3, fill=NaN, na.rm=T)

roll_meanr(c(NA,NA,1,NA,NA,1,NA,NA,1), n = 3, fill=NaN, na.rm=T)

roll_meanr(c(NA,NA,1,NA,NA,1,NA,NA,1,NA,NA,3,0,0), n = 3, fill=NaN, na.rm=T)


npp_tree %>% #filter(date>="2010-01-01") %>% 
  mutate(agC_Mg_diff_std = as.double(agC_Mg_diff/agC_Mg)) %>% 
  filter(agC_Mg_diff_std > 0) %>% 
  ggplot(data=., aes(log(agC_Mg), log(agC_Mg_diff_std)))+geom_point()+geom_smooth(method='lm')

dend_out %>% ggplot(data = ., aes(dev, npp_bestEst_agC_month))+geom_point()





dend %>% filter(date>="2015-01-01") %>% 
  mutate(category=cut(current_dbh, breaks=c(-Inf,150, 200, 300, 500, 750, Inf))) %>% #, labels=c("small","middle","large","massive"))) %>% 
  filter(is.na(category)==F) %>% 
  ggplot(data=., aes(date, current_dbh, color=tree_tag))+
  geom_point()+geom_line()+
  facet_wrap(~category, scales="free", ncol = 2)+
  theme(legend.position = "none")

census %>% group_by(census_date) %>% summarize(nobs=n())

dend %>% group_by(site) %>% 
  summarize(nobs=n()) %>% 
  ggplot(data=., aes(date, nobs))+geom_point()+
  facet_wrap(~site)

dend_all$plot_code %>% unique() %>% sort()
dend_all %>% mutate(site=substr(plot_code,1,3)) %>% 
                  filter(site %in% c("ANK","BOB","CAX","KOG","LPG","TAM","KEN", 
                                     "STO","NXV","JEN","TAN","SAF","MLA")) %>% 
                  dim()
dend_all %>% mutate(site=substr(plot_code,1,3)) %>% 
  filter(site %in% c("ANK","BOB","CAX","KOG","LPG","TAM","KEN", 
                     "STO","NXV","JEN","TAN","SAF","MLA")) %>% 
  group_by(year,month,plot_code) %>% summarize(nobs=n())

dend_all %>% mutate(site=substr(plot_code,1,3)) %>% 
  filter(site %in% c("ANK","BOB","CAX","KOG","LPG","TAM","KEN", 
                     "STO","NXV","JEN","TAN","SAF","MLA")) %>% 
  mutate(plot_code_tree_tag=paste0(plot_code,"_",tree_tag)) %>% pull(plot_code_tree_tag) %>% unique() %>% length()

dend_all %>% mutate(site=substr(plot_code,1,3)) %>% 
  filter(site %in% c("TAN")) %>% group_by(plot_code, date) %>% 
  summarize(nobs=n(), 
            u=mean(dendrometer_reading_mm,na.rm=T)) %>% 
  ggplot(data=., aes(date, u,color=plot_code))+geom_point()


census_all %>% mutate(site=substr(plot_code,1,3)) %>% filter(site=="TAN") %>% 
ggplot(data=., aes(d0,height))+geom_point()

census$wd



junk1 <- read.csv("data/forestplots20180327.csv") %>% as.tibble()
junk1 %>% filter(Plot.Code=="TAN-05") %>% pull(D0) %>% is.na %>% table
junk1 %>% filter(Plot.Code=="TAN-05") %>% pull(D0) %>% summary()

junk2  <- read_csv("data/forestplots20180327.csv", na=c("NA", "NaN", ""), guess_max=100000) # Census_Santarem_2014_2016.csv
junk2 %>% spec
junk2  <- read_csv2("data/forestplots20180327.csv", ) # Census_Santarem_2014_2016.csv
junk2 %>% filter(Plot code=="TAN-05") %>% pull(D0) %>% is.na() %>% table

dend %>% filter(is.na(dendrometer_reading_mm)==F) %>% 
  distinct() %>% 
  group_by(date) %>% 
  summarize(u=mean(dendrometer_reading_mm), 
            p05=quantile(dendrometer_reading_mm,0.05), 
            p95=quantile(dendrometer_reading_mm,0.95), 
            u_dbh=mean(current_dbh,na.rm=T), 
            nobs=n()) %>% 
  ggplot(data=., aes(date, nobs))+geom_line()
  ggplot(data=., aes(date, u))+geom_line()+
  geom_ribbon(aes(date, ymax=p95,ymin=p05),lty=0,alpha=0.25)





# Subset for nxv repo -----------------------------------------------------


  library(tidyverse); library(lubridate)
  
  d <- read_csv('data/dendro_20180327.csv')
  cen <- read_csv("data/forestplots20180327.csv")
  names(cen) <- tolower(names(cen))
  cen <- cen %>% 
    mutate(site = substr(`plot code`,1,3))
  
  cen_nxv <- cen %>% 
    filter(site=='NXV')
  d <- d %>% 
    mutate(site = substr(plot_code,1,3))
  d <- d %>% 
    filter(site == "NXV")
  
  d %>% write_csv(.,"data/dend_nxv.csv")
  cen_nxv %>% write_csv(.,"data/census_nxv.csv")
  