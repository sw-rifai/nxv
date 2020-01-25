library(fitdistrplus); library(tidyverse); library(knitr)
nxv_root <- read_csv("data/fine_roots/nxv_root.csv")


# Both plots, arithmetic quadrature approach ------------------------------
# very little data to estimate SE!!!
nxv_root %>% 
  filter(root <= 0.6) %>% 
  group_by(plot_code) %>% 
  summarize(nobs = n(), 
            root_u = mean(root)*12, # 12*monthly -> annual
            root_sd = sqrt(sum(monthlyNPProot_sd**2))*12) %>% 
  ungroup() %>% 
  mutate(root_sem = root_sd/sqrt(16*2)) %>%  # dividing by sqrt of 16*2 because there were essentially 2 years of observations with 16 ingrowth cores
  kable()

# NXV-01 ------------------------------------------------------------------
vec_nxv1_root <- nxv_root %>% 
  filter(plot_code=='NXV-01') %>% 
  mutate(date=ymd(paste(year,month,1))) %>% 
  arrange(date) %>% 
  filter(date != min(date)) %>%  # The first measurement was MUCH higher so I think it was a stock and not an NPP measurement
  pull(root)
plot(vec_nxv1_root, type='b') # There is precious little data, so we'll use it all.

# fit distribution for NXV-01 roots
# There is little root data so to determine which distribution is best, I'll combine plots
f_gamma <- fitdistrplus::fitdist(vec_nxv1_root, distr='gamma')
f_w <- fitdistrplus::fitdist(vec_nxv1_root, distr='weibull')
f_ln <- fitdistrplus::fitdist(vec_nxv1_root, distr='lnorm')
f_n <- fitdistrplus::fitdist(vec_nxv1_root, distr='norm')

f_gamma$aic
f_w$aic
f_ln$aic # lowest AIC for the log-normal distribution so we'll use that
f_n$aic


# plot log normal distribution of fine root from NXV-01
exp(f_ln$estimate["meanlog"]+(f_ln$estimate['sdlog']**2)/2) # mean
curve(dlnorm(x, meanlog = f_ln$estimate['meanlog'], 
             sdlog = f_ln$estimate['sdlog']),0,0.5, 
      ylab='probability density')
points(x=vec_nxv1_root,
       y=rep(1, length(vec_nxv1_root)),
       pch=20,col=rgb(0,0,0,alpha=0.3), 
       cex=3)
abline(v=exp(f_ln$estimate["meanlog"]+(f_ln$estimate['sdlog']**2)/2))
abline(v=exp(f_ln$estimate["meanlog"]+1.96*f_ln$sd["meanlog"]+(f_ln$sd['sdlog']**2)/2),lty=3)
abline(v=exp(f_ln$estimate["meanlog"]-1.96*f_ln$sd["meanlog"]+(f_ln$sd['sdlog']**2)/2),lty=3)

out_root_df <- data.frame(plot_code = c("NXV-01","NXV-02"), 
                          mean_fine_root_MgC_ha_yr = NA, 
                          sem_fine_root_MgC_ha_yr = NA)
out_root_df[out_root_df$plot_code=='NXV-01',]$mean_fine_root_MgC_ha_yr <- as.numeric(12*exp(f_ln$estimate["meanlog"]+(f_ln$estimate['sdlog']**2)/2))
out_root_df[out_root_df$plot_code=='NXV-01',]$sem_fine_root_MgC_ha_yr <- as.numeric(12*exp(f_ln$estimate["meanlog"]+(f_ln$estimate['sdlog']**2)/2))

f1 <- fitdistrplus::bootdist(f_ln, niter=1000)
plot(f1,enhance=T)
summary(f1)
quantile(f1,probs = c(0.5))

# NXV-02 ------------------------------------------------------------------
vec_nxv2_root <- nxv_root %>% 
  filter(plot_code=='NXV-02') %>% 
  mutate(date=ymd(paste(year,month,1))) %>% 
  arrange(date) %>% 
  filter(date != min(date)) %>%  # The first measurement was MUCH higher so I think it was a stock and not an NPP measurement
  pull(root)
plot(vec_nxv2_root, type='b') # There is precious little data, so we'll use it all.

# fit distribution for NXV-02 roots
# There is little root data so to determine which distribution is best, I'll combine plots
f_gamma <- fitdistrplus::fitdist(vec_nxv2_root, distr='gamma')
f_w <- fitdistrplus::fitdist(vec_nxv2_root, distr='weibull')
f_ln <- fitdistrplus::fitdist(vec_nxv2_root, distr='lnorm')
f_n <- fitdistrplus::fitdist(vec_nxv2_root, distr='norm')

f_gamma$aic
f_w$aic
f_ln$aic # lowest AIC for the log-normal distribution so we'll use that
f_n$aic


# plot log normal distribution of fine root from NXV-01
exp(f_ln$estimate["meanlog"]+(f_ln$estimate['sdlog']**2)/2) # mean
curve(dlnorm(x, meanlog = f_ln$estimate['meanlog'], 
             sdlog = f_ln$estimate['sdlog']),0,0.75, 
      ylab='probability density')
points(x=vec_nxv2_root,
       y=rep(1, length(vec_nxv1_root)),
       pch=20,col=rgb(0,0,0,alpha=0.3), 
       cex=3)
abline(v=exp(f_ln$estimate["meanlog"]+(f_ln$estimate['sdlog']**2)/2))
abline(v=exp(f_ln$estimate["meanlog"]+1.96*f_ln$sd["meanlog"]+(f_ln$sd['sdlog']**2)/2),lty=3)
abline(v=exp(f_ln$estimate["meanlog"]-1.96*f_ln$sd["meanlog"]+(f_ln$sd['sdlog']**2)/2),lty=3)

out_root_df[out_root_df$plot_code=='NXV-02',]$mean_fine_root_MgC_ha_yr <- as.numeric(12*exp(f_ln$estimate["meanlog"]+(f_ln$estimate['sdlog']**2)/2))
out_root_df[out_root_df$plot_code=='NXV-02',]$sem_fine_root_MgC_ha_yr <- as.numeric(12*exp(f_ln$estimate["meanlog"]+(f_ln$estimate['sdlog']**2)/2))

f1 <- fitdistrplus::bootdist(f_ln, niter=1000)
plot(f1,enhance=T)
summary(f1)
quantile(f1,probs = c(0.5))




tibble(plot_code = c("NXV-01", "NXV-02"), 
       median_est = c(as.numeric(f1$quantiles), 
                      as.numeric(f2$quantiles)), 
       CI_2.5 = c(o1$quantCI$`p=0.5`[1] ,o2$quantCI$`p=0.5`[1] ), 
       CI_97.5 = c(o1$quantCI$`p=0.5`[2] ,o2$quantCI$`p=0.5`[2] )
) %>% 
  kable()



# SCRATCH -----------------------------------------------------------------
# 
# 
# lm(x~1) %>% summary
# hist(x)
# sd(x)/sqrt(length(x))
# 
# x <- rgamma(1000, shape=10, rate=2);
# summary(x)
# glm(x~1, family=Gamma(link='identity')) %>% summary
# 
# 
# p_size <- 1e5 # population
# alpha_p <- rlnorm(p_size, meanlog = -1.96, sdlog=0.184) # full pop
# hist(alpha_p)
# summary(alpha_p)
# sd(alpha_p)/sqrt(length(alpha_p))
# 
# s_size <- 5
# fn_sem <- function(s_size){
#   sd(alpha_p[sample.int(p_size, s_size)])/sqrt(s_size)
# }
# curve(fn_sem(x),2,1000,col='red', ylim=c(0,0.005))
# for(i in 1:1000){
#   curve(fn_sem(x),2,1000, add=T, lwd=0.1)
# }
# 
# fn_sd <- function(s_size){
#   sd(alpha_p[sample.int(p_size, s_size)])
# }
# fn_sd <- Vectorize(fn_sd)
# curve(fn_sd(x),2,1000, ylim=c(0,0.05),lty=0.1)
# for(i in 1:100){
#   curve(fn_sd(x),2,1000, add=T, lwd=0.1)
# }
# abline(h=sqrt((exp(0.184**2)-1)*exp(2*-1.96+0.184**2)),col='red')
# 
# 
# # (1) Fit of a gamma distribution to serving size data
# # using default method (maximum likelihood estimation)
# # followed by parametric bootstrap
# #
# data(groundbeef)
# x1 <- groundbeef$serving
# f1 <- fitdist(x1, "gamma")
# b1 <- bootdist(f1, niter=51)
# print(b1)
# plot(b1)
# plot(b1, enhance=TRUE)
# summary(b1)
# quantile(b1)
# CIcdfplot(b1, CI.output = "quantile")
