lf_rollsmooth %>% 
  ggplot(data=., aes(date, lf_mo, color=plot_code))+
  geom_point()+
  geom_line()
