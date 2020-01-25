# Scaling NPP from Nova Xavantina

## Uncertainty propagation 
### I. Arithmetic approach: 
(1) estimate monthly totals and associated standard error of the mean
(2) estimate annual totals and associated standard error by using the
'quadrature' method.
(3) estimate multi-year estimate of annual NPP, again using the 
'quadrature' method. 

### II. Distribution & Bootstrap approach

### Package requirements
dplyr\
lubridate\
RcppRoll\
kable\
fitdistrplus\

## NPP Components: 

### Canopy Total Fine Litterfall $(Mg~C~ha^-1~yr^-1)$
This part of the NPP data has a large sample size collected between 
2011/Feb through 2016/Feb, so the more simpole arithmetic approach
might be sufficient here. However the order of calculation makes a difference.\
see: R/litterfall/sem_litterfall.R \


calculation order:\
(1) estimate annual sum of litterfall from each trap\
(2) calculate mean annual litterfall and standard error of the mean\

|plot_code | nobs|    val_u|   val_sd|   val_sem|
|:---------|----:|--------:|--------:|---------:|
|NXV-01    |  154| 2.376703| 1.006774| 0.0811282|
|NXV-02    |  154| 4.778234| 1.444391| 0.1163924|


calculation order:\
(1) calc year/month litterfall from each trap\
(2) calc yearly sum of litterfall\
(3) calc overall annual mean and sem of litterfall\

|plot_code | nobs|    val_u|    val_sd|   val_sem|
|:---------|----:|--------:|---------:|---------:|
|NXV-01    |    6| 1.956995| 0.9362134| 0.3822075|
|NXV-02    |    6| 3.935002| 1.8262495| 0.7455632|

\
\

distribution approach using 1000 bootstraped iterations from the log-normal \

|plot_code | median_est|   CI_2.5|  CI_97.5|
|:---------|----------:|--------:|--------:|
|NXV-01    |   2.182000| 2.042683| 2.317462|
|NXV-02    |   4.579182| 4.372251| 4.795130|

### Woody NPP from census data $(Mg~Biomass~ha^-1~yr^-1)$ 
Notes: This is an arithmetic calculation on adult trees greater than 100 mm diameter. I think the census data file I am using did not come from ForestPlots. It has many errors with d0 and species. I have tried to fix these in the code, but there are probably still errors so inspect carefully. s

|plot_code |    npp_u|   npp_sd| nobs|   npp_sem|
|:---------|--------:|--------:|----:|---------:|
|NXV-01    | 2.791126| 1.009742|    5| 0.4515705|
|NXV-02    | 2.774994| 1.162931|    4| 0.5814653|

\
\

### Woody NPP from dendrometer data $(Mg~C~ha^-1~yr^-1)$ 
**NOT WORKING YET**
This is the approach I used in the Rifai et al 2018 Phil Trans Roy Soc B paper. There are some stong assumptions here:\
(1) NPP is estimated using dendrometer bands and one fixed census observation. The dendrometer bands are dynamic through time and are used to
fit a model of woody NPP as a function of size and wood density for each measurement date. The model is then upscaled to the plot level with the fixed census observation.

\
\

### Fine root NPP $(Mg~C~ha^-1~yr^-1)$ 
There is a very limited amount of non-normally distributed data here so I estimate this using an arthmetic approach, but if I can find the individual ingrowth core data, I can
better estimate this through bootstrapping the best fit parametric distribution. 
The root NPP data generally has a skewed distribution, so I suspect the distribution estimate is more robust. 

arithmetic approach, summing errors in quadrature with an outlier removed from NXV-01:\

|plot_code | nobs|   root_u|  root_sd|  root_sem|
|:---------|----:|--------:|--------:|---------:|
|NXV-01    |    7| 2.221714| 3.180091| 0.5621659|
|NXV-02    |    6| 2.166000| 4.134100| 0.7308126|


arithmetic approach, summing errors in quadrature:\

|plot_code | nobs|   root_u|  root_sd|  root_sem|
|:---------|----:|--------:|--------:|---------:|
|NXV-01    |    7| 2.221714| 3.180091| 0.5621659|
|NXV-02    |    7| 2.916000| 6.110711| 1.0802312|

\
\

### FOR COMPARISON ------------------------------------------------------ 

From Marina's manuscript --------------------------------------------\
canopy (forest: 5.17 and savanna: 2.72 Mg C ha-1 year-1) 
wood (forest: 1.35 and savanna: 0.51 Mg C ha-1 year-1) 
fine roots NPP (1.78 and 2.33 Mg C ha-1 year-1 for forest and savanna