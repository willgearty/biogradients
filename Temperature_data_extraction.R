library(RNetCDF)
library(dplyr)
library(ncdf4)

# This script is used to generate summary distributions of annual mean equatorial sea surface temperature for a range of 
# future climate scenarios and historical models of preindustrial sea surface temperature using the HadGEM2-ES model

###############################################################################
## Historical model
###############################################################################

# open NetCDF
HADGEM.T <- open.nc("thetao_Omon_HadGEM2-ES_historical_r1i1p1_185912-186911.nc")

# print NetCDF summary
print.nc(HADGEM.T)

# extract variables and assign variable names
temp <- var.get.nc(HADGEM.T, "thetao") # units = "K"
depth.T <- var.get.nc(HADGEM.T, "lev") # units = "m"
lat.T <- var.get.nc(HADGEM.T, "lat") # units = "degrees_north" 
lon.T <- var.get.nc(HADGEM.T, "lon") # units = "degrees_east"
time.T <- var.get.nc(HADGEM.T, "time") # units = "days since 1859-12-01", calendar = "360_day"

# How many monthly increments into NetCDF timeseries is the time interval of interest?
months.all.1860_start <- 1

# Generate summary statistics of annual equatorial (-0.34 to 0.34 degrees) mean temperature in surface ocean (top 10m)
temp.eq.1860.K.dist <- quantile(temp[,108:109,1,months.all.1860_start:(months.all.1860_start+11)], c(0.05, .25, .5, .75, .95), na.rm=T)

# Generate dataframe to combine with other model summaries at end of script and rename headers
temp.eq.1860.K.sum <- as.data.frame(cbind(1860, t(temp.eq.1860.K.dist), "Historical"))
names(temp.eq.1860.K.sum) <- c("year", "temp.K.5", "temp.K.25", "temp.K.50", "temp.K.75", "temp.K.95", "scenario")

###############################################################################
## 2100 - RCP 8.5
###############################################################################

# open NetCDF
HADGEM.T <- open.nc("thetao_Omon_HadGEM2-ES_rcp85_r1i1p1_209912-210911.nc")

# print NetCDF summary
print.nc(HADGEM.T)

# extract variables and assign variable names
temp <- var.get.nc(HADGEM.T, "thetao") # units = "K"
depth.T <- var.get.nc(HADGEM.T, "lev") # units = "m"
lat.T <- var.get.nc(HADGEM.T, "lat") # units = "degrees_north" 
lon.T <- var.get.nc(HADGEM.T, "lon") # units = "degrees_east"
time.T <- var.get.nc(HADGEM.T, "time") # units = "days since 1859-12-01", calendar = "360_day"

# How many monthly increments into NetCDF timeseries is the time interval of interest?
months.all.2100_start <- 1

# Generate summary statistics of annual equatorial (-0.34 to 0.34 degrees) mean temperature in surface ocean (top 10m)
temp.eq.2100.K.dist <- quantile(temp[,108:109,1,months.all.2100_start:(months.all.2100_start+11)], c(0.05, .25, .5, .75, .95), na.rm=T)

# Generate dataframe to combine with other model summaries at end of script and rename headers
temp.eq.2100.K.sum.85 <- as.data.frame(cbind(2100, t(temp.eq.2100.K.dist), "RCP8.5"))
names(temp.eq.2100.K.sum.85) <- c("year", "temp.K.5", "temp.K.25", "temp.K.50", "temp.K.75", "temp.K.95", "scenario")

###############################################################################
## 2100 - RCP 4.5
###############################################################################

# open NetCDF
HADGEM.T <- open.nc("thetao_Omon_HadGEM2-ES_rcp45_r1i1p1_201512-202511.nc")

# print NetCDF summary
print.nc(HADGEM.T)

# extract variables and assign variable names
temp <- var.get.nc(HADGEM.T, "thetao") # units = "K"
depth.T <- var.get.nc(HADGEM.T, "lev") # units = "m"
lat.T <- var.get.nc(HADGEM.T, "lat") # units = "degrees_north" 
lon.T <- var.get.nc(HADGEM.T, "lon") # units = "degrees_east"
time.T <- var.get.nc(HADGEM.T, "time") # units = "days since 1859-12-01", calendar = "360_day"

# How many monthly increments into NetCDF timeseries is the time interval of interest?
months.all.2020_start <- 1+12+12+12+12

# Generate summary statistics of annual equatorial (-0.34 to 0.34 degrees) mean temperature in surface ocean (top 10m)
temp.eq.2020.K.dist <- quantile(temp[,108:109,1,months.all.2020_start:(months.all.2020_start+11)], c(0.05, .25, .5, .75, .95), na.rm=T)

# Generate dataframe to combine with other model summaries at end of script and rename headers
temp.eq.2020.K.sum.45 <- as.data.frame(cbind(2020, t(temp.eq.2020.K.dist), "RCP4.5"))
names(temp.eq.2020.K.sum.45) <- c("year", "temp.K.5", "temp.K.25", "temp.K.50", "temp.K.75", "temp.K.95", "scenario")

###############################################################################
## 2299 - RCP 8.5
###############################################################################
# Notably - the year 2299 is spread across two NetCDF files - so two temperature arrays are merged below

# open NetCDF
HADGEM.T.1 <- open.nc("thetao_Omon_HadGEM2-ES_rcp85_r1i1p1_228912-229911.nc")

# print NetCDF summary
print.nc(HADGEM.T)

# extract variables and assign variable names
temp.1 <- var.get.nc(HADGEM.T.1, "thetao") # units = "K"
depth.T <- var.get.nc(HADGEM.T.1, "lev") # units = "m"
lat.T <- var.get.nc(HADGEM.T.1, "lat") # units = "degrees_north" 
lon.T <- var.get.nc(HADGEM.T.1, "lon") # units = "degrees_east"
time.T <- var.get.nc(HADGEM.T.1, "time") # units = "days since 1859-12-01", calendar = "360_day"

HADGEM.T.2 <- open.nc("thetao_Omon_HadGEM2-ES_rcp85_r1i1p1_229912-229912.nc")

print.nc(HADGEM.T.2)

# extract variables and assign variable names
temp.2 <- var.get.nc(HADGEM.T.2, "thetao") # units = "K"

# generate temperature array
temp <- array(dim=c(360,216,40,121))
temp[,,,1:120] <- temp.1
temp[,,,121] <- temp.2

# How many monthly increments into NetCDF timeseries is the time interval of interest?
months.all.2300_start <- 110

# Generate summary statistics of annual equatorial (-0.34 to 0.34 degrees) mean temperature in surface ocean (top 10m)
temp.eq.2300.K.dist <- quantile(temp[,108:109,1,months.all.2300_start:(months.all.2300_start+11)], c(0.05, .25, .5, .75, .95), na.rm=T)

# Generate dataframe to combine with other model summaries at end of script and rename headers
temp.eq.2300.K.sum.85 <- as.data.frame(cbind(2299, t(temp.eq.2300.K.dist), "RCP8.5"))
names(temp.eq.2300.K.sum.85) <- c("year", "temp.K.5", "temp.K.25", "temp.K.50", "temp.K.75", "temp.K.95", "scenario")

###############################################################################
## 2299 - RCP 4.5
###############################################################################

# open NetCDF
HADGEM.T <- open.nc("thetao_Omon_HadGEM2-ES_rcp45_r1i1p1_209912-210911.nc")

# print NetCDF summary
print.nc(HADGEM.T)

# extract variables and assign variable names
temp <- var.get.nc(HADGEM.T, "thetao") # units = "K"
depth.T <- var.get.nc(HADGEM.T, "lev") # units = "m"
lat.T <- var.get.nc(HADGEM.T, "lat") # units = "degrees_north" 
lon.T <- var.get.nc(HADGEM.T, "lon") # units = "degrees_east"
time.T <- var.get.nc(HADGEM.T, "time") # units = "days since 1859-12-01", calendar = "360_day"

# How many monthly increments into NetCDF timeseries is the time interval of interest?
months.all.2100_start <- 1

# Generate summary statistics of annual equatorial (-0.34 to 0.34 degrees) mean temperature in surface ocean (top 10m)
temp.eq.2100.K.dist <- quantile(temp[,108:109,1,months.all.2100_start:(months.all.2100_start+11)], c(0.05, .25, .5, .75, .95), na.rm=T)

# Generate dataframe to combine with other model summaries at end of script and rename headers
temp.eq.2100.K.sum.45 <- as.data.frame(cbind(2100, t(temp.eq.2100.K.dist), "RCP4.5"))
names(temp.eq.2100.K.sum.45) <- c("year", "temp.K.5", "temp.K.25", "temp.K.50", "temp.K.75", "temp.K.95", "scenario")

###############################################################################
## 2299 - RCP 4.5
###############################################################################
# Notably - the year 2299 is spread across two NetCDF files - so two temperature arrays are merged below

# open NetCDF
HADGEM.T.1 <- open.nc("thetao_Omon_HadGEM2-ES_rcp85_r1i1p1_228912-229911.nc")

# print NetCDF summary
print.nc(HADGEM.T)

# extract variables and assign variable names
temp.1 <- var.get.nc(HADGEM.T.1, "thetao") # units = "K"
depth.T <- var.get.nc(HADGEM.T.1, "lev") # units = "m"
lat.T <- var.get.nc(HADGEM.T.1, "lat") # units = "degrees_north" 
lon.T <- var.get.nc(HADGEM.T.1, "lon") # units = "degrees_east"
time.T <- var.get.nc(HADGEM.T.1, "time") # units = "days since 1859-12-01", calendar = "360_day"

# open NetCDF
HADGEM.T.2 <- open.nc("thetao_Omon_HadGEM2-ES_rcp45_r1i1p1_229912-229912.nc")

# print NetCDF summary
print.nc(HADGEM.T.2)

# extract variables and assign variable names
temp.2 <- var.get.nc(HADGEM.T.2, "thetao") # units = "K"

# generate temperature array
temp <- array(dim=c(360,216,40,121))
temp[,,,1:120] <- temp.1
temp[,,,121] <- temp.2

# How many monthly increments into NetCDF timeseries is the time interval of interest?
months.all.2300_start <- 110

# Generate summary statistics of annual equatorial (-0.34 to 0.34 degrees) mean temperature in surface ocean (top 10m)
temp.eq.2300.K.dist <- quantile(temp[,108:109,1,months.all.2300_start:(months.all.2300_start+11)], c(0.05, .25, .5, .75, .95), na.rm=T)

# Generate dataframe to combine with other model summaries at end of script and rename headers
temp.eq.2300.K.sum.45 <- as.data.frame(cbind(2299, t(temp.eq.2300.K.dist), "RCP4.5"))
names(temp.eq.2300.K.sum.45) <- c("year", "temp.K.5", "temp.K.25", "temp.K.50", "temp.K.75", "temp.K.95", "scenario")

###############################################################################
## Generate full summary dataframe and save
###############################################################################
temp.sum <- rbind(temp.eq.1860.K.sum, 
      temp.eq.2100.K.sum.85,
      temp.eq.2300.K.sum.85, 
      temp.eq.2100.K.sum.45,
      temp.eq.2300.K.sum.45)

save(temp.sum, file="Equatorial.temp.sum.RData")