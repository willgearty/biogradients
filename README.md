This repository contains the required R code and data files to produce all of the analyses and plots in ‘Metabolic tradeoffs control biodiversity gradients through geologic time’   
Thomas H. Boag, William Gearty & Richard G. Stockey  
  
Before running these R scripts, we suggest that you download this folder and set it as your R working directory. Required data files should then load when called in each script, and plot files will save within the same folder.  

## Diversity-temperature analyses (main text figures 1 & 2, supplementary figures 1-13):  
Run `diversity_temp.R` to produce all plots and analyses. Note that the analyses can take many hours to run.  
You can load `fossil_results.RData` into your environment to skip running all of the analyses.  

## Metabolic model (main text figure 3, supplementary figures 14-15):  
Figure 3 - run `Metabolic_model.10000.R`, then run `Metabolic_model_plot.10000.R`. Note that to fully reproduce Figure 3 'Temperature_data_extraction.R' will also need to be run (more details below), otherwise the ggplot2 object 'model' can be saved individually.
Figure S14 - run `Metabolic_model.10000.supp.R`, then run `Metabolic_model_plot.10000.supp.R`  
Figure S15 - run `Distribution_summary_plots.R`  
Note that `Metabolic_model.10000.R` and `Metabolic_model.10000.supp.R` have long runtimes. To reduce these runtimes, the number of subsamples taken in the Monte Carlo can be reduced (e.g. 10,000 to 1000 or lower) to produce similar model results in shorter runtime. 

## Climate model equatorial temperature distributions (upper panel main text figure 3):  
To reproduce the equatorial temperature distributions plotted in the upper panel of Figure 3, run 'Temperature_data_extraction.R' after downloading the NetCDF files listed therein from the following webpages:
http://data.ceda.ac.uk/badc/cmip5/data/cmip5/output1/MOHC/HadGEM2-ES/historical/mon/ocean/Omon/r1i1p1/files/thetao_20110916 (Historical)
http://data.ceda.ac.uk/badc/cmip5/data/cmip5/output1/MOHC/HadGEM2-ES/rcp45/mon/ocean/Omon/r1i1p1/v20111206/thetao (RCP/ECP 4.5)
http://data.ceda.ac.uk/badc/cmip5/data/cmip5/output1/MOHC/HadGEM2-ES/rcp85/mon/ocean/Omon/r1i1p1/files/thetao_20111218 (RCP/ECP 8.5)

## To replicate the analyses and plots presented here, the following R packages are required:  
AICcmodvg  
broom  
deeptime  
deSolve  
dispeRse  
dplyr  
geosphere   
ggplot2  
Hmisc  
MASS  
mgcv  
ncdf4
paleoMap  
plyr  
raster  
rgbif  
RNetCDF
rworldmap  
segmented  
viridis  

The deeptime, paleoMap and dispeRse packages are currently only available on GitHub, and will need to be installed from there to reproduce
these analyses and plots. This can be achieved by running the following commands in your R console (ignore the first line 
if you already have devtools installed).  
```r
install.packages("devtools")  
devtools::install_github("willgearty/deeptime")  
devtools::install_github("willgearty/paleoMap")  
devtools::install_github("laurasoul/dispeRse")  
```
  
  
All other packages can be installed from CRAN. These scripts have been tested using R version 4.0.3 - 
Copyright (C) 2020 The R Foundation for Statistical Computing.
