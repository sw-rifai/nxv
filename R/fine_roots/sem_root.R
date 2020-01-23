# simple arithmetic way to calculate the standard error of the mean for 
# overall annual estimate


library(tidyverse); library(lubridate);

nxv_root <- read_csv("data/fine_roots/nxv_root.csv")
nxv_root %>% 
  group_by(plot_code, year) %>% 
  summarize(val = mean(root)*12, 
            val_sd = sqrt(sum(monthlyNPProot_sd))*12, 
            val_sem = val_sd/n(), 
            nobs=n()) %>% 
  ungroup() %>% 
  group_by(plot_code) %>% 
  summarize(val_mean = mean(val), 
            val_sd = sd(val), 
            val_sem = val_sd/n()) %>% 
  ungroup() %>% 
  ggplot(data=., aes(plot_code, y=val_mean))+
  geom_col()+
  geom_errorbar(aes(
                    ymin=val_mean-val_sem, 
                    ymax = val_mean+val_sem))
  
