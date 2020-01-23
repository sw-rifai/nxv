# Scaling NPP from Nova Xavantina

## Uncertainty propagation 
- Arithmetic approach: 
(1) estimate monthly totals and associated standard error of the mean
(2) estimate annual totals and associated standard error by using the
'quadrature' method.
(3) estimate multi-year estimate of annual NPP, again using the 
'quadrature' method. 
- Attempting to follow method of Malhi et al., (2017) New Phytologist 


## Package requirements
library(dplyr); library(lubridate); library(RcppRoll)

### Fine Litterfall
This part of the NPP data has a large sample size collected between 
2011/Feb through 2016/Feb. I think the arithmetic approach is
sufficiently robust here. 

### Woody NPP from census data
Notes: This is an arithmetic calculation on adult trees greater than 100 mm diameter. I think the census data file I am using did not come from ForestPlots. It has many errors with d0 and species. I have tried to fix these in the code, but there are probably still errors so inspect carefully. s

NXV-01: NPP Mg Biomass ha-1 yr-1
  npp_u npp_sd  nobs npp_sem
  2.79   1.01     5   0.452

NXV-02: NPP Mg Biomass ha-1 yr-1
  npp_u npp_sd  nobs npp_sem
  2.77   1.16     4   0.581


### Woody NPP from dendrometer data
This is the approach I used in the Rifai et al 2018 Phil Trans Roy Soc B paper. There are some stong assumptions here: 
(1) NPP is estimated using dendrometer bands and one fixed census observation. The dendrometer bands are dynamic through time and are used to
fit a model of woody NPP as a function of size and wood density for each measurement date. The model is then upscaled to the plot level with the fixed census observation.


### Fine root NPP 
There is a very limited amount of non-normally distributed data here so I estimate this using an arthmetic approach, and through estimating a best fit parametric distribution. 
Because the sample size is quite small and has a skewed distribution, I suspect the distribution estimate is more robust. 


FOR COMPARISON ------------------------------------------------------
From Marina's manuscript --------------------------------------------
Forest showed higher NPP for canopy (forest: 5.17 and savanna: 2.72 Mg C ha-1 year-1) and wood (forest: 1.35 and savanna: 0.51 Mg C ha-1 year-1) components, but savanna showed higher fine roots NPP (1.78 and 2.33 Mg C ha-1 year-1 for forest and savanna, respectively; Fig. 3a)