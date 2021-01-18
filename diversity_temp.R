#Packages, etc.####
library(rgbif)
library(plyr)
library(dplyr)
library(raster)
#library(devtools)
#install_github("ropensci/paleobioDB")
#install_github("willgearty/paleoMap")
library(paleoMap)
#install_github("laurasoul/dispeRse")
library(dispeRse)
library(mgcv)

#Shareholder quorom subsampling and rarefaction
source("SQS-3-3.R")

#Load PBDB data
load("PBDBdata.RData")

#Global settings####
#Clean workspace without deleting PBDB data or SQS function
rm(list = ls()[-c(grep("sqs",ls()),grep("moll_raw",ls()))])

#SQS quorum level
#Set this to a single value if you want to do further plotting/analyses
#If you want to aggregate results over many SQS quorum levels, set this to a vector of multiple values
SQS_levels <- c(.25, .5, .75)
#SQS_levels <- .25

#Number of SQS trials to run for each diversity estimate
SQS_trials <- 100

#Number of total latitudinal bins
#Set this to a single value if you want to do further plotting/analyses
#If you want to aggregate results over many bin sizes, set this to a vector of multiple values
lat_bin_nums <- 180/c(1, 2, 4, 7.5, 15)
#lat_bin_nums <- 180/c(.5, .75, 1, 2, 5, 7.5, 10, 15)
#lat_bin_nums <- 24

#Which type of latitudinal bands to use ['equal-width' or 'equal-area']
#Set this to a single value if you want to do further plotting/analyses
#If you want to aggregate results over both band types, set this to a vector of both values
band_types <- c("equal-area", "equal-width")
#band_types <- "equal-area"

#need enough rows to include all band types, numbers of bins, types of diversity metrics,
#whether we are including the modern or not, and whether we are using raw or proportional diversity
fossil_results <- data.frame(matrix(NA, nrow = length(band_types) * length(lat_bin_nums) *
                                      (2 + length(SQS_levels)) * 2 * 2, ncol = 19))
colnames(fossil_results) <- c("band_type", "num_bands", "raw_prop", "div_metric", "inc_modern",
                              "breakpoint", "breakpoint_se",
                              "slope1", "slope1_se", "intercept1", "intercept1_se", "r1_2", "p1",
                              "slope2", "slope2_se", "intercept2", "intercept2_se", "r2_2", "p2")
div_tot_list <- vector(mode="list")
mod_div_tot_list <- vector(mode="list")
T_max_K <- vector(mode="list")
degfs <- vector(mode="list")

pb_main <- tcltk::tkProgressBar(max = nrow(fossil_results))
k <- 1
for(band_type in band_types){
  for(lat_bin_num in lat_bin_nums){
    #Setup latitudinal bins####
    if(band_type == "equal-width"){
      #Equal-width latitudinal bins
      r <- raster(nrows=lat_bin_num, ncols=180, xmn=-180, xmx=180, ymn=-90, ymx=90)
      area <- rowSums(as.matrix(raster::area(r)))
      
      lat_bounds <- seq(-90 , 90, length.out = lat_bin_num + 1)
    } else if(band_type == "equal-area"){
      #Equal-area latitudinal bins
      lat_bounds <- EqualAreaRectangularGrid(N_longitude = 1, N_latitude = lat_bin_num)$latitude_breaks
      
      #set up dataframe for latitudinal temperatures
      area = rep((4*pi*6371^2)/lat_bin_num, lat_bin_num)
    }
    lat_bands <- data.frame(lat_min = lat_bounds[1:(length(lat_bounds) - 1)], lat_max = lat_bounds[2:length(lat_bounds)])
    lat_bands$lat_mid <- (lat_bands$lat_min + lat_bands$lat_max)/2
    lat_bands$lat_bin <- levels(cut(-90:90, lat_bounds, include.lowest = TRUE))
    lat_bands$area <- area
    
    #Modern####
    print("Analyzing Modern...", quote = FALSE)
    #Process modern temperature data
    print("  Analyzing Temperature Data", quote = FALSE)
    temp <- read.csv("woa13_decav_t00an01v2.csv", skip = 1, fill = TRUE, row.names = NULL, header = TRUE, check.names=FALSE)
    colnames(temp)[1:3] <- c("Latitude","Longitude","SST_C")
    temp$SST_K <- temp$SST_C + 273.15
    
    ## gam regression and prediction
    set.seed(17)
    fit <- gam(SST_K ~ s(Latitude), data = temp)
    print(paste("    ", round(sum(fit$edf[fit$smooth[[1]]$first.para:fit$smooth[[1]]$last.para]), digits = 2), "degrees of freedom"), quote=FALSE)
    
    pred <- predict(fit, newdata = data.frame(Latitude = lat_bands$lat_mid), se.fit = TRUE)
    cells <- lat_bands
    cells$SST_K <- pred[["fit"]]
    cells$SST_K_se <- pred[["se.fit"]]
    cells$inv_K <- 1000/cells$SST_K
    
    T_max_K$modern <- max(cells$SST_K, na.rm = TRUE)
    degfs$modern <- round(sum(fit$edf[fit$smooth[[1]]$first.para:fit$smooth[[1]]$last.para]), digits = 2)
    lat_bands$modern_SST_K <- cells$SST_K[match(lat_bands$lat_bin, cells$lat_bin)]
    lat_bands$modern_SST_K_se <- cells$SST_K_se[match(lat_bands$lat_bin, cells$lat_bin)]
    
    #clades <- c("Gastropods","Bivalves","Molluscs")
    #download_ids <- c("0079232-160910150852091", "0080034-160910150852091", "0080041-160910150852091")
    #names(download_ids) <- clades
    
    clades <- c("Molluscs")
    #occurrences up to 200m compiled/downloaded on 12/12/20
    download_ids <- c("Molluscs" = "0135238-200613084148143")
    
    #perform SQS for each of the clades and save diversity as columns in cells df
    for(clade in clades){
      print(paste("  Analyzing Modern", clade), quote = FALSE)
      if(!file.exists(paste(download_ids[clade],".csv", sep=""))){
        print("    Downloading GBIF Data", quote = FALSE)
        occ_download_get(download_ids[clade], overwrite = TRUE)
        unzip(paste(download_ids[clade],".zip", sep=""))
      }
      print("    Analyzing GBIF Data", quote = FALSE)
      specimens <- read.table(paste(download_ids[clade],".csv", sep=""), sep = "\t", header = TRUE, fill = TRUE, quote = NULL, stringsAsFactors = FALSE)
      specimens <- subset(specimens, !is.na(decimalLatitude))
      
      genera <- subset(specimens, taxonRank %in% c("GENUS", "SPECIES", "SUBGENUS", "SUBSPECIES") & genus !="")
      
      genera_loc <- genera[,c("genus","decimalLatitude","decimalLongitude")]
      
      genera_loc$lat_bin <- cut(genera_loc$decimalLatitude, lat_bounds, include.lowest = TRUE)
      counts <- ddply(genera_loc, .(lat_bin), summarize, div = length(unique(genus)), counts = length(genus))
      
      #record specimen counts
      cells[[paste("modern_specimens", clade, sep = "_")]] <- as.numeric(counts$counts[match(cells$lat_bin, counts$lat_bin)])
      
      #Calculate raw diversity
      cells[[paste("modern_raw", clade, sep = "_")]] <- as.numeric(NA)
      mod_div_tot_list[[paste("modern_raw", clade, sep = "_")]] <- mod_div_tot_list[[paste("modern_range", clade, sep = "_")]] <- as.numeric(length(unique(genera_loc$genus)))
      cells[[paste("modern_raw", clade, sep = "_")]] <- as.numeric(counts$div[match(cells$lat_bin, counts$lat_bin)])
      
      #Calculate diversity based on ranges
      ranges <- genera_loc %>%
        group_by(genus) %>%
        dplyr::summarise(
          min_lat = min(decimalLatitude),
          max_lat = max(decimalLatitude)
        )
      
      cells[[paste("modern_range", clade, sep = "_")]] <- as.numeric(NA)
      for(j in 1:length(cells$lat_bin)){
        div <- nrow(subset(ranges, max_lat >= cells$lat_min[j] & min_lat <= cells$lat_max[j]))
        if(div > 0){
          cells[[paste("modern_range", clade, sep = "_")]][j] <- div
        }
      }
      
      #Calculate SQS diversity
      for(SQS_level in SQS_levels){
        print(paste0("      SQS Level: ", SQS_level), quote = FALSE)
        capture.output(div <- sqs(as.numeric(table(as.character(genera_loc$genus))), q = SQS_level, trials = SQS_trials))
        mod_div_tot_list[[paste("modern_SQS", SQS_level, clade, sep = "_")]] <- as.numeric(div["subsampled richness"])
        pb <- txtProgressBar(max = length(cells$lat_bin), style = 3)
        cells[[paste(clade, "SQS", SQS_level, sep = "_")]] <- as.numeric(NA)
        for(n in 1:length(cells$lat_bin)){
          sub_occs <- subset(genera_loc, lat_bin == cells$lat_bin[n])
          if(nrow(sub_occs) > 0){
            tab <- as.numeric(table(as.character(sub_occs$genus)))
            capture.output(div <- sqs(tab, q = SQS_level, trials = SQS_trials))
            cells[[paste(clade, "SQS", SQS_level, sep = "_")]][n] <- as.numeric(div["subsampled richness"])
          }
          setTxtProgressBar(pb, n)
        }
        close(pb)
      }
    }
    
    #Save mollusc data, in case we are looking at more clades
    for(SQS_level in SQS_levels){
      lat_bands[[paste("modern_div_SQS",SQS_level, sep = "_")]] <- cells[[paste("Molluscs_SQS",SQS_level, sep = "_")]][match(lat_bands$lat_bin, cells$lat_bin)]
      div_tot_list[[paste("modern_SQS", SQS_level, sep = "_")]] <- mod_div_tot_list[[paste("modern_SQS", SQS_level, "Molluscs", sep = "_")]]
    }
    div_tot_list[["modern_raw"]] <- mod_div_tot_list[["modern_raw_Molluscs"]]
    div_tot_list[["modern_range"]] <- mod_div_tot_list[["modern_range_Molluscs"]]
    lat_bands$modern_div_raw <- cells$modern_raw_Molluscs[match(lat_bands$lat_bin, cells$lat_bin)]
    lat_bands$modern_specimens <- cells$modern_specimens_Molluscs[match(lat_bands$lat_bin, cells$lat_bin)]
    lat_bands$modern_div_range <- cells$modern_range_Molluscs[match(lat_bands$lat_bin, cells$lat_bin)]
    
    print("  Modern Done!", quote = FALSE)
    
    #Paleo Analyses Setup####
    paleotemps <- read.csv("paleotemp_data.csv", stringsAsFactors = FALSE, na.strings = c("NA", "/", ""))
    paleotemps$SST_K <- paleotemps$Temp_C + 273.15
    # periods <- c("Pliocene", "Miocene", "Early Oligocene", "Middle Eocene", "Late Paleocene-Early Eocene",
    #              "Latest Cretaceous", "Early Cretaceous")
    # periods_abbr <- c("plio", "mio", "ear_olig", "mid_eoc", "late_pal_ear_eoc", "late_cret", "ear_cret")
    # periods_PBDB <- list(c("Pliocene"), c("Late Miocene"), c("Early Oligocene"), c("Middle Eocene"),
    #                      c("Late Paleocene", "Early Eocene"), c("Campanian", "Maastrichtian"),
    #                      c("Early Cretaceous"))
    
    #Aptian-Albian doesn't have enough data
    stages <- c("Late Pliocene", "Late Miocene", "Early Oligocene", "Middle Eocene", "Early Eocene", "Late Paleocene",
                "Maastrichtian", "Campanian", "Berriasian-Barremian")
    stages_abbr <- c("late_plio", "late_mio", "ear_olig", "mid_eoc", "ear_eoc", "late_pal", "maas", "camp", "berr_barr")
    stages_PBDB <- list(c("Late Pliocene"), c("Serravallian", "Tortonian", "Messinian"), c("Early Oligocene"), c("Middle Eocene"),
                        c("Early Eocene"), c("Late Paleocene"), c("Maastrichtian"), c("Campanian"),
                        c("Berriasian", "Valanginian", "Hauterivian", "Barremian"))
    names(stages_abbr) <- names(stages_PBDB) <- stages
    
    #Paleo Analyses####
    for(i in 1:length(stages)){
      print(paste0("Analyzing ", stages[i], "..."), quote = FALSE)
      #Get paleotemperature subset for this period
      print("  Analyzing Paleotemperature Data", quote = FALSE)
      paleotemp <- subset(paleotemps, Stage == stages[i]) %>% mutate(Latitude = Paleolatitude)
      ## gam regression and prediction
      set.seed(17)
      fit <- gam(SST_K ~ s(Latitude, k = 9), data = paleotemp)
      print(paste("    ", round(sum(fit$edf[fit$smooth[[1]]$first.para:fit$smooth[[1]]$last.para]), digits = 2), "degrees of freedom"), quote=FALSE)
      
      pred <- predict(fit, newdata = data.frame(Latitude = lat_bands$lat_mid), se.fit = TRUE)
      lat_bands[[paste(stages_abbr[i],"SST_K", sep = "_")]] <- pred[["fit"]]
      lat_bands[[paste(stages_abbr[i],"SST_K_se", sep = "_")]] <- pred[["se.fit"]]
      T_max_K[[stages_abbr[i]]] <- max(predict(fit, newdata = data.frame(Latitude = seq(-90,90,.1))))
      degfs[[stages_abbr[i]]] <- round(sum(fit$edf[fit$smooth[[1]]$first.para:fit$smooth[[1]]$last.para]), digits = 2)
      
      #Get PBDB data if we don't have it yet
      if(!exists(paste(stages_abbr[i],"moll_raw", sep = "_"))){
        print("  Downloading PBDB Data", quote = FALSE)
        assign(paste(stages_abbr[i],"moll_raw", sep = "_"),
               do.call(rbind, lapply(stages_PBDB[[i]], pm_getdata, base_name = "Mollusca")))
      }
      
      print("  Analyzing PBDB Data", quote = FALSE)
      #Make temporary data.frame so we don't mess up the raw data
      occs <- get(paste(stages_abbr[i],"moll_raw", sep = "_"))
      occs$lat_bin <- cut(occs$paleolat, lat_bounds, include.lowest = TRUE)
      occs <- subset(occs, matched_rank %in% c("genus", "subgenus", "species", "subspecies") & !is.na(paleolat))
      counts <- ddply(occs, .(lat_bin), summarize, div = length(unique(genus)), counts = length(genus))
      
      #record specimen counts
      lat_bands[[paste(stages_abbr[i],"specimens", sep = "_")]] <- as.numeric(counts$counts[match(lat_bands$lat_bin, counts$lat_bin)])
      
      #Calculate raw diversity
      div_tot_list[[paste(stages_abbr[i],"raw", sep = "_")]] <- div_tot_list[[paste(stages_abbr[i],"range", sep = "_")]] <- as.numeric(length(unique(occs$genus)))
      lat_bands[[paste(stages_abbr[i],"div_raw", sep = "_")]] <- as.numeric(counts$div[match(lat_bands$lat_bin, counts$lat_bin)])
      
      #Calculate diversity based on ranges
      ranges <- occs %>%
        group_by(genus) %>%
        dplyr::summarise(
          min_lat = min(paleolat),
          max_lat = max(paleolat)
        )
      
      lat_bands[[paste(stages_abbr[i],"div_range", sep = "_")]] <- as.numeric(NA)
      for(j in 1:length(lat_bands$lat_bin)){
        div <- nrow(subset(ranges, max_lat >= lat_bands$lat_min[j] & min_lat <= lat_bands$lat_max[j]))
        if(div > 0){
          lat_bands[[paste(stages_abbr[i],"div_range", sep = "_")]][j] <- div
        }
      }
      
      #Get SQS diversity
      for(SQS_level in SQS_levels){
        print(paste0("      SQS Level: ", SQS_level), quote = FALSE)
        capture.output(div <- sqs(as.numeric(table(as.character(occs$genus))), q = SQS_level, trials = SQS_trials))
        div_tot_list[[paste(stages_abbr[i],"SQS",SQS_level, sep = "_")]] <- div["subsampled richness"]
        counts[[paste("SQS",SQS_level,"div", sep = "_")]] <- NA
        pb <- txtProgressBar(max = length(counts$lat_bin), style = 3)
        for(n in 1:length(counts$lat_bin)){
          sub_occs <- subset(occs, lat_bin == counts$lat_bin[n])
          if(nrow(sub_occs) > 0){
            tab <- as.numeric(table(as.character(sub_occs$genus)))
            capture.output(div <- sqs(tab, q = SQS_level, trials = SQS_trials))
            counts[[paste("SQS", SQS_level, "div", sep = "_")]][n] <- div["subsampled richness"]
          }
          setTxtProgressBar(pb, n)
        }
        close(pb)
        lat_bands[[paste(stages_abbr[i],"div_SQS", SQS_level, sep = "_")]] <- as.numeric(counts[[paste("SQS",SQS_level,"div", sep = "_")]][match(lat_bands$lat_bin, counts$lat_bin)])
      }
      print(paste(" ", stages[i], "Done!"), quote = FALSE)
    }
    
    #Melt Data####
    print("Melting Data", quote = FALSE)
    library(ggplot2)
    library(reshape2)
    melted_SST <- melt(lat_bands, id.vars = c("lat_bin", "lat_min", "lat_mid", "lat_max", "area"), measure.vars = paste(c("modern", stages_abbr), "SST_K", sep = "_"),
                       variable.name = "time", value.name = "SST_K")
    melted_SST$time <- gsub("_SST_K", "", melted_SST$time)
    melted_div <- melt(lat_bands, id.vars = c("lat_bin"), measure.vars = as.character(outer(paste(c("modern", stages_abbr), "div",sep = "_"), c("raw", "range", paste("SQS", SQS_levels, sep = "_")), FUN = "paste", sep = "_")),
                           variable.name = "time", value.name = "div")
    melted_div$metric <- gsub("^.*_div_", "", melted_div$time)
    melted_div$time <- gsub("_div.*$", "", melted_div$time)
    
    melted_specimens <- melt(lat_bands, id.vars = c("lat_bin"), measure.vars = paste(c("modern", stages_abbr), "specimens", sep = "_"),
                             variable.name = "time", value.name = "specimens")
    melted_specimens$time <- gsub("_specimens", "", melted_specimens$time)
    
    melted_SST_div <- merge(merge(melted_SST, melted_div, by = c("time", "lat_bin")), melted_specimens, by = c("time", "lat_bin"))
    
    degf_df <- data.frame(time = names(degfs), degf = unlist(degfs))
    div_tot_df <- data.frame(time = names(div_tot_list), div_tot = unlist(div_tot_list))
    div_tot_df$metric <- gsub(paste0("(",paste(c("modern", stages_abbr), collapse = "|"),")_"), "", div_tot_df$time)
    div_tot_df$time <- gsub("_(raw|range|SQS).*$", "", div_tot_df$time)
    
    melted_SST_div <- merge(melted_SST_div, div_tot_df, by = c("time", "metric"))
    melted_SST_div <- merge(melted_SST_div, degf_df, by = "time")
    
    melted_SST_div$time <- factor(melted_SST_div$time, levels = c("modern", stages_abbr), ordered = TRUE)
    
    #Calculate Allen stats####
    print("Calculating Stats", quote = FALSE)
    melted_SST_div$Allen.SST <- 1000/melted_SST_div$SST_K
    melted_SST_div$div.prop <- melted_SST_div$div/melted_SST_div$div_tot
    
    #Save melted data for later
    melted_SST_div$band_type <- band_type
    melted_SST_div$num_bands <- lat_bin_num
    if(exists("fossil_data")) {
      fossil_data <- rbind(fossil_data, melted_SST_div)
    } else {
      fossil_data <- melted_SST_div
    }
    
    for(met in levels(factor(melted_SST_div$metric))){
      #Proportional diversity with modern####
      #Find breakpoint
      library(segmented)
      reg <- lm(log(div.prop) ~ Allen.SST, data = subset(melted_SST_div, metric == met))
      reg_seg <- with(subset(melted_SST_div, metric == met), segmented(reg, seg.Z = ~ Allen.SST))
      
      #Find regressions on both sides of breakpoint
      brk <- reg_seg$psi[2:3]
      lm_1 <- lm(log(div.prop) ~ Allen.SST, data = subset(melted_SST_div,Allen.SST < brk[1] & metric == met))
      summ1 <- summary(lm_1)
      lm_2 <- lm(log(div.prop) ~ Allen.SST, data = subset(melted_SST_div,Allen.SST > brk[1] & metric == met))
      summ2 <- summary(lm_2)
      
      #Save allen stats
      fossil_results$band_type[k] <- band_type
      fossil_results$num_bands[k] <- lat_bin_num
      fossil_results$div_metric[k] <- met
      fossil_results$inc_modern[k] <- "yes"
      fossil_results$raw_prop[k] <- "prop"
      fossil_results$breakpoint[k] <- brk[1]
      fossil_results$breakpoint_se[k] <- brk[2]
      fossil_results$slope1[k] <- summ1$coefficients["Allen.SST","Estimate"]
      fossil_results$slope1_se[k] <- summ1$coefficients["Allen.SST","Std. Error"]
      fossil_results$intercept1[k] <- summ1$coefficients["(Intercept)","Estimate"]
      fossil_results$intercept1_se[k] <- summ1$coefficients["(Intercept)","Std. Error"]
      fossil_results$r1_2[k] <- summ1$r.squared
      fossil_results$p1[k] <- summ1$coefficients["Allen.SST","Pr(>|t|)"]
      fossil_results$slope2[k] <- summ2$coefficients["Allen.SST","Estimate"]
      fossil_results$slope2_se[k] <- summ2$coefficients["Allen.SST","Std. Error"]
      fossil_results$intercept2[k] <- summ2$coefficients["(Intercept)","Estimate"]
      fossil_results$intercept2_se[k] <- summ2$coefficients["(Intercept)","Std. Error"]
      fossil_results$r2_2[k] <- summ2$r.squared
      fossil_results$p2[k] <- summ2$coefficients["Allen.SST","Pr(>|t|)"]
      k <- k+ 1
      
      #Proportional diversity without modern####
      #Find breakpoint
      library(segmented)
      reg <- lm(log(div.prop) ~ Allen.SST, data = subset(melted_SST_div, metric == met & time != "modern"))
      reg_seg <- with(subset(melted_SST_div, metric == met & time != "modern"), segmented(reg, seg.Z = ~ Allen.SST))
      
      #Find regressions on both sides of breakpoint
      brk <- reg_seg$psi[2:3]
      lm_1 <- lm(log(div.prop) ~ Allen.SST, data = subset(melted_SST_div,Allen.SST < brk[1] & metric == met & time != "modern"))
      summ1 <- summary(lm_1)
      lm_2 <- lm(log(div.prop) ~ Allen.SST, data = subset(melted_SST_div,Allen.SST > brk[1] & metric == met & time != "modern"))
      summ2 <- summary(lm_2)
      
      #Save allen stats
      fossil_results$band_type[k] <- band_type
      fossil_results$num_bands[k] <- lat_bin_num
      fossil_results$div_metric[k] <- met
      fossil_results$inc_modern[k] <- "no"
      fossil_results$raw_prop[k] <- "prop"
      fossil_results$breakpoint[k] <- brk[1]
      fossil_results$breakpoint_se[k] <- brk[2]
      fossil_results$slope1[k] <- summ1$coefficients["Allen.SST","Estimate"]
      fossil_results$slope1_se[k] <- summ1$coefficients["Allen.SST","Std. Error"]
      fossil_results$intercept1[k] <- summ1$coefficients["(Intercept)","Estimate"]
      fossil_results$intercept1_se[k] <- summ1$coefficients["(Intercept)","Std. Error"]
      fossil_results$r1_2[k] <- summ1$r.squared
      fossil_results$p1[k] <- summ1$coefficients["Allen.SST","Pr(>|t|)"]
      fossil_results$slope2[k] <- summ2$coefficients["Allen.SST","Estimate"]
      fossil_results$slope2_se[k] <- summ2$coefficients["Allen.SST","Std. Error"]
      fossil_results$intercept2[k] <- summ2$coefficients["(Intercept)","Estimate"]
      fossil_results$intercept2_se[k] <- summ2$coefficients["(Intercept)","Std. Error"]
      fossil_results$r2_2[k] <- summ2$r.squared
      fossil_results$p2[k] <- summ2$coefficients["Allen.SST","Pr(>|t|)"]
      k <- k + 1
      
      #Raw diversity with modern####
      #Find breakpoint
      library(segmented)
      reg <- lm(log(div) ~ Allen.SST, data = subset(melted_SST_div, metric == met))
      reg_seg <- with(subset(melted_SST_div, metric == met), segmented(reg, seg.Z = ~ Allen.SST))
      
      #Find regressions on both sides of breakpoint
      brk <- reg_seg$psi[2:3]
      lm_1 <- lm(log(div) ~ Allen.SST, data = subset(melted_SST_div,Allen.SST < brk[1] & metric == met))
      summ1 <- summary(lm_1)
      lm_2 <- lm(log(div) ~ Allen.SST, data = subset(melted_SST_div,Allen.SST > brk[1] & metric == met))
      summ2 <- summary(lm_2)
      
      #Save allen stats
      fossil_results$band_type[k] <- band_type
      fossil_results$num_bands[k] <- lat_bin_num
      fossil_results$div_metric[k] <- met
      fossil_results$inc_modern[k] <- "yes"
      fossil_results$raw_prop[k] <- "raw"
      fossil_results$breakpoint[k] <- brk[1]
      fossil_results$breakpoint_se[k] <- brk[2]
      fossil_results$slope1[k] <- summ1$coefficients["Allen.SST","Estimate"]
      fossil_results$slope1_se[k] <- summ1$coefficients["Allen.SST","Std. Error"]
      fossil_results$intercept1[k] <- summ1$coefficients["(Intercept)","Estimate"]
      fossil_results$intercept1_se[k] <- summ1$coefficients["(Intercept)","Std. Error"]
      fossil_results$r1_2[k] <- summ1$r.squared
      fossil_results$p1[k] <- summ1$coefficients["Allen.SST","Pr(>|t|)"]
      fossil_results$slope2[k] <- summ2$coefficients["Allen.SST","Estimate"]
      fossil_results$slope2_se[k] <- summ2$coefficients["Allen.SST","Std. Error"]
      fossil_results$intercept2[k] <- summ2$coefficients["(Intercept)","Estimate"]
      fossil_results$intercept2_se[k] <- summ2$coefficients["(Intercept)","Std. Error"]
      fossil_results$r2_2[k] <- summ2$r.squared
      fossil_results$p2[k] <- summ2$coefficients["Allen.SST","Pr(>|t|)"]
      k <- k+ 1
      
      #Raw diversity without modern####
      #Find breakpoint
      library(segmented)
      reg <- lm(log(div) ~ Allen.SST, data = subset(melted_SST_div, metric == met & time != "modern"))
      reg_seg <- with(subset(melted_SST_div, metric == met & time != "modern"), segmented(reg, seg.Z = ~ Allen.SST))
      
      #Find regressions on both sides of breakpoint
      brk <- reg_seg$psi[2:3]
      lm_1 <- lm(log(div) ~ Allen.SST, data = subset(melted_SST_div,Allen.SST < brk[1] & metric == met & time != "modern"))
      summ1 <- summary(lm_1)
      lm_2 <- lm(log(div) ~ Allen.SST, data = subset(melted_SST_div,Allen.SST > brk[1] & metric == met & time != "modern"))
      summ2 <- summary(lm_2)
      
      #Save allen stats
      fossil_results$band_type[k] <- band_type
      fossil_results$num_bands[k] <- lat_bin_num
      fossil_results$div_metric[k] <- met
      fossil_results$inc_modern[k] <- "no"
      fossil_results$raw_prop[k] <- "raw"
      fossil_results$breakpoint[k] <- brk[1]
      fossil_results$breakpoint_se[k] <- brk[2]
      fossil_results$slope1[k] <- summ1$coefficients["Allen.SST","Estimate"]
      fossil_results$slope1_se[k] <- summ1$coefficients["Allen.SST","Std. Error"]
      fossil_results$intercept1[k] <- summ1$coefficients["(Intercept)","Estimate"]
      fossil_results$intercept1_se[k] <- summ1$coefficients["(Intercept)","Std. Error"]
      fossil_results$r1_2[k] <- summ1$r.squared
      fossil_results$p1[k] <- summ1$coefficients["Allen.SST","Pr(>|t|)"]
      fossil_results$slope2[k] <- summ2$coefficients["Allen.SST","Estimate"]
      fossil_results$slope2_se[k] <- summ2$coefficients["Allen.SST","Std. Error"]
      fossil_results$intercept2[k] <- summ2$coefficients["(Intercept)","Estimate"]
      fossil_results$intercept2_se[k] <- summ2$coefficients["(Intercept)","Std. Error"]
      fossil_results$r2_2[k] <- summ2$r.squared
      fossil_results$p2[k] <- summ2$coefficients["Allen.SST","Pr(>|t|)"]
      k <- k + 1
      
      tcltk::setTkProgressBar(pb_main, k - 1)
    }
  }
}
close(pb_main)

#Figure S1####
#Plot occurrences
library(deeptime)
library(ggplot2)
library(viridis)
stages_map <- c("Pliocene", "Tortonian", "Rupelian", "Lutetian", "Ypresian", "Thanetian",
                "Maastrichtian", "Campanian", "Valanginian")

colsea = "#00509005"
colland = "#66666660"

occs <- do.call(rbind, lapply(1:length(stages), function(i) get(paste(stages_abbr[i], "moll_raw", sep = "_")) %>%
                                mutate(stage = stages[i]))) %>% mutate(stage = factor(stage, stages))
gg_maps <- do.call(rbind, lapply(1:length(stages), function(i) fortify(get(stages_map[i])) %>%
                                   mutate(stage = stages[i]))) %>% mutate(stage = factor(stage, stages))

time_ages <- data.frame(stages_abbr, start = c(3.333, 13.82, 33.9, 47.8, 56, 59.2, 72.1, 83.6, 145),
                        end = c(2.58, 5.333, 27.82, 37.8, 47.8, 56, 66, 72.1, 125))
stages_with_ages <- setNames(paste0(stages, "\n(", time_ages$start, " - ", time_ages$end, " Ma)"), unique(gg_maps$stage))

ggplot(occs) +
  geom_polygon(data = gg_maps, aes(long, lat, group = group), fill = colland, color = "black") +
  geom_hex(aes(paleolng, paleolat), binwidth = 5) +
  coord_fixed(xlim = c(-180, 180), ylim = c(-90, 90)) +
  scale_x_continuous(expand = c(0,0), breaks = c(seq(-180,180,30)), minor_breaks = NULL) +
  scale_y_continuous(expand = c(0,0), breaks = c(seq(-90,90,30)), minor_breaks = NULL) +
  scale_fill_viridis(name = "# Occs", trans= "log10") +
  facet_wrap(vars(stage), labeller = as_labeller(stages_with_ages)) +
  theme_bw(base_size = 20) +
  theme(legend.position = "right", axis.title = element_blank(), axis.text = element_text(colour = "black")) +
  theme(panel.background = element_rect(fill = colsea), plot.margin = unit(c(.5,.5,.5,.5), "cm"),
        panel.spacing.y = unit(.5, "lines"), panel.spacing.x = unit(3, "lines"),
        strip.background = element_blank(), strip.placement = "outside",
        strip.text = element_text(size = 30, margin = unit(c(.5,.5,.5,.5), "lines")))
ggsave("Fossil Occurences by Stage.pdf", width = 28, height = 15)

#Sampling Sensitivity####
#Figure S3####
ggplot(fossil_data, aes(specimens, div)) +
  geom_vline(xintercept = c(30, 50, 200), size = 2, linetype = "11", color = "grey50") +
  geom_point(shape = 21, alpha = .3) +
  #geom_smooth() +
  scale_x_continuous(name = "# of Specimens in Latitudinal Bin", trans = "log10",
                     breaks = c(1,10,100,1000,10000,100000), labels = c("1","10","100","1000","10000","100000")) +
  scale_y_continuous(name = "Estimated Diversity", trans = "log10") +
  facet_wrap(~metric, ncol = 1, scales = "free_y",
             labeller = as_labeller(c("range" = "Rangethrough", "raw" = "Raw Occurrences", "SQS_0.25" = "SQS (q = .25)", "SQS_0.5" = "SQS (q = .5)", "SQS_0.75" = "SQS (q = .75)"))) +
  theme_classic(base_size = 24, base_family = "") +
  theme(axis.ticks = element_line(color = "black"),
        axis.line = element_blank(), axis.text = element_text(colour = "black"),
        plot.margin = unit(c(1,1,1,1), "lines"), panel.border = element_rect(fill = NA)) +
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(margin = unit(c(1,.5,.5,.5), "lines")))
ggsave("Sampling Sensitivity.pdf", width = 8, height = 12)

#Figure S2####
#Latitude and Temperature Sampling
library(deeptime)
g1 <- ggplot(subset(fossil_data, metric == "raw" & band_type == "equal-area"), aes(SST_K - 273.15, specimens)) +
  geom_point(alpha = .3) +
  scale_x_continuous(name = expression("Mean Sea Surface Temperature of Bin ("*degree*"C)")) +
  scale_y_continuous(name = "# Occurrences", trans = "log10",
                     breaks = c(1,100,10000), labels = c("1","100","10000")) +
  theme_classic(base_size = 24, base_family = "") +
  theme(axis.ticks = element_line(color = "black"),
        axis.line = element_blank(), axis.text = element_text(colour = "black"),
        plot.margin = unit(c(1,1,1,1), "lines"), panel.border = element_rect(fill = NA)) +
  facet_wrap(vars(num_bands), ncol = 1, labeller = as_labeller(function(x) paste(x, "latitude bins")))

g2 <- ggplot(subset(fossil_data, metric == "raw" & band_type == "equal-area"), aes(lat_mid, specimens)) +
  geom_point(alpha = .3) +
  scale_x_continuous(name = "Mean Paleolatitude of Bin") +
  scale_y_continuous(name = "# Occurrences", trans = "log10",
                     breaks = c(1,100,10000), labels = c("1","100","10000")) +
  theme_classic(base_size = 24, base_family = "") +
  theme(axis.ticks = element_line(color = "black"),
        axis.line = element_blank(), axis.text = element_text(colour = "black"),
        plot.margin = unit(c(1,1,1,1), "lines"), panel.border = element_rect(fill = NA)) +
  facet_wrap(vars(num_bands), ncol = 1, labeller = as_labeller(function(x) paste(x, "latitude bins")))
gg <- ggarrange2(g1, g2, nrow = 1, draw = FALSE)
ggsave("Temperature and Latitude Sampling Intensity.pdf", gg, width = 16, height = 14)

#Figure S10####
#model fits before and after cutoffs
test_models <- function(div.prop, SST_K, n_tries = 20) {
  library(segmented)
  library(AICcmodavg)
  df <- data.frame(model = c("linear", "one_break", "two_breaks", "three_breaks", "four_breaks", "sec_ord", "third_ord", "fourth_order", "gam"),
                   r2.adj = NA, AICc = NA, AICc_delta = NA, AICc_weight = NA)
  reg <- lm(log(div.prop) ~ SST_K)
  df[1,2:3] <- c(summary(reg)$adj.r.squared, AICc(reg))
  ind <- 1
  while(ind <= n_tries) {
    tryCatch({
      reg_seg <- segmented(reg, seg.Z = ~ SST_K, npsi = 1)
      ind <- n_tries + 1
    }, warning = function(e){
    }, error = function(e){
    }, finally = {
      if(ind == (n_tries + 1)) df[2,2:3] <- c(summary(reg_seg)$adj.r.squared, AICc(reg_seg))
    })
    ind <- ind + 1
  }
  ind <- 1
  while(ind <= n_tries) {
    tryCatch({
      reg_seg <- segmented(reg, seg.Z = ~ SST_K, npsi = 2)
      ind <- n_tries + 1
    }, warning = function(e){
    }, error = function(e){
    }, finally = {
      if(ind == (n_tries + 1)) df[3,2:3] <- c(summary(reg_seg)$adj.r.squared, AICc(reg_seg))
    })
    ind <- ind + 1
  }
  ind <- 1
  while(ind <= n_tries) {
    tryCatch({
      reg_seg <- segmented(reg, seg.Z = ~ SST_K, npsi = 3)
      ind <- n_tries + 1
    }, warning = function(e){
    }, error = function(e){
    }, finally = {
      if(ind == (n_tries + 1)) df[4,2:3] <- c(summary(reg_seg)$adj.r.squared, AICc(reg_seg))
    })
    ind <- ind + 1
  }
  ind <- 1
  while(ind <= n_tries) {
    tryCatch({
      reg_seg <- segmented(reg, seg.Z = ~ SST_K, npsi = 4)
      ind <- n_tries + 1
    }, warning = function(e){
    }, error = function(e){
    }, finally = {
      if(ind == (n_tries + 1)) df[5,2:3] <- c(summary(reg_seg)$adj.r.squared, AICc(reg_seg))
    })
    ind <- ind + 1
  }
  
  reg <- lm(log(div.prop) ~ poly(SST_K, 2))
  df[6,2:3] <- c(summary(reg)$adj.r.squared, AICc(reg))
  reg <- lm(log(div.prop) ~ poly(SST_K, 3))
  df[7,2:3] <- c(summary(reg)$adj.r.squared, AICc(reg))
  reg <- lm(log(div.prop) ~ poly(SST_K, 4))
  df[8,2:3] <- c(summary(reg)$adj.r.squared, AICc(reg))
  
  gam_fit <- gam(log(div.prop) ~ s(SST_K))
  df[9,2:3] <- c(summary(gam_fit)$r.sq, AICc(gam_fit))
  
  df$AICc_delta[!is.na(df$AICc)] <- df$AICc[!is.na(df$AICc)] - min(df$AICc[!is.na(df$AICc)])
  model_lik <- exp(-0.5*df$AICc_delta[!is.na(df$AICc)])
  df$AICc_weight[!is.na(df$AICc)] <- model_lik/sum(model_lik)
  
  return(df)
}

# run a set of models for each metric with varying sampling cutoffs
set.seed(12345)
model_fits <- do.call(rbind, lapply(1:100, function (i) {
  rbind(fossil_data %>% filter(num_bands == 24, band_type == "equal-area") %>%
                      mutate(cutoff = 0),
                    fossil_data %>% filter(num_bands == 24, band_type == "equal-area", specimens > 30) %>%
                      mutate(cutoff = 30),
                    fossil_data %>% filter(num_bands == 24, band_type == "equal-area", specimens > 50) %>%
                      mutate(cutoff = 50),
                    fossil_data %>% filter(num_bands == 24, band_type == "equal-area", specimens > 200) %>%
                      mutate(cutoff = 200)) %>%
  mutate(cutoff = factor(cutoff)) %>%
  group_by(cutoff, num_bands, metric) %>%
  summarize(test_models(div.prop, SST_K, n_tries = 20)) %>%
  arrange(cutoff, num_bands, metric, AICc) %>%
  mutate(model = factor(model, levels = c("linear", "one_break", "two_breaks", "three_breaks", "four_breaks", "sec_ord", "third_ord", "fourth_order", "gam")),
         rep = i)
}))

library(ggplot2)
library(viridis)
ggplot(model_fits) +
  geom_boxplot(aes(x = cutoff, y = AICc_weight, color = model, fill = model), size = .25, width = .99, position = position_dodge2(padding = .1)) +
  stat_summary(aes(x = cutoff, y = AICc_weight, group = interaction(cutoff, model)), fun = median, geom = "crossbar", size = .1, width = .99, position = position_dodge2(padding = .075)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5)) +
  xlab("Sampling Cutoff (# Specimens)") +
  ylab("AICc Weight") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(lim = c(0,1)) +
  scale_fill_viridis(name = "Model", discrete = TRUE, breaks = c("linear", "one_break", "two_breaks", "three_breaks", "four_breaks", "sec_ord", "third_ord", "fourth_order", "gam"),
                     labels = c("Linear", "Linear with 1 Breakpoint", "Linear with 2 Breakpoints", "Linear with 3 Breakpoints",
                                "Linear with 4 Breakpoints", "Quadratic", "Cubic", "Quartic", "GAM"),
                     guide = guide_legend(nrow = 3)) +
  scale_color_viridis(name = "Model", discrete = TRUE, breaks = c("linear", "one_break", "two_breaks", "three_breaks", "four_breaks", "sec_ord", "third_ord", "fourth_order", "gam"),
                     labels = c("Linear", "Linear with 1 Breakpoint", "Linear with 2 Breakpoints", "Linear with 3 Breakpoints",
                                "Linear with 4 Breakpoints", "Quadratic", "Cubic", "Quartic", "GAM"),
                     guide = guide_legend(nrow = 3)) +
  facet_wrap(vars(metric), ncol = 1,
             labeller = as_labeller(c("range" = "Rangethrough", "raw" = "Raw Occurrences", "SQS_0.25" = "SQS (q = .25)", "SQS_0.5" = "SQS (q = .5)", "SQS_0.75" = "SQS (q = .75)"))) +
  theme_classic(base_size = 20, base_family = "") +
  theme(axis.ticks = element_line(color = "black"),
        axis.line = element_blank(), axis.text = element_text(colour = "black"), legend.position = "top",
        plot.margin = unit(c(1,1,1,1), "lines"), panel.border = element_rect(fill = NA)) +
  theme(panel.spacing = unit(0, "lines"),  legend.margin = margin(0,0,0,0, "lines"),
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(margin = unit(c(1,.5,.5,.5), "lines")))
ggsave("Model Fits by Sampling Cutoff.pdf", width = 10, height = 12)

#Table 1####
#pull out summary stats for main analysis model fitting
model_fits %>% filter(cutoff == 0, metric == "SQS_0.25") %>%
  group_by(model) %>%
  summarize(AICc_d = round(median(AICc_delta), 2), AICc_w = round(median(AICc_weight), 2))

#minimum delta AICc for linear model across all sensitivity permutations
model_fits %>% filter(model == "linear") %>%
  group_by(metric, cutoff) %>%
  summarize(test = min(AICc_delta))

#Figure S12####
#plot allen-like results
library(ggplot2)
library(viridis)
ann_text <- data.frame(inc_modern = "yes", num_bands = 24, raw_prop = "prop",
                       div_metric = "SQS_0.25", band_type = "equal-area", slope2 = 3.5,
                       lab = "This Study")

ggplot(fossil_results, aes(x = as.character(num_bands), y = slope2, shape = band_type)) +
  scale_x_discrete(limits = as.character(sort(lat_bin_nums)), name = "Number of Latitudinal Bands") +
  scale_y_continuous(name = "Estimated Slope", sec.axis = sec_axis(~., labels = NULL)) +
  geom_rect(xmin = 0, xmax = 50, ymin = -7.65, ymax = -6.71, color = "grey80", fill = "grey80") +
  geom_segment(x = 0, xend = 50, y = -7.65, yend = -7.65, aes(linetype = "roy"), color = "grey80") +
  geom_rect(xmin = 0, xmax = 50, ymin = -8.12, ymax = -6.96, linetype = "dashed", color = "black", fill = NA, size = .5) +
  geom_segment(x = 0, xend = 50, y = -6.96, yend = -6.96, aes(linetype = "theory"), color = NA, size = .5) +
  geom_point(aes(color = div_metric), size = 2, position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(color = div_metric, ymin = slope2 - 1.96*slope2_se, ymax = slope2 + 1.96*slope2_se), width = .8, position = position_dodge(width = 0.8)) +
  geom_text(data = ann_text, aes(label = lab), show.legend = FALSE,
            color = "black", hjust = .55, size = 5) +
  geom_segment(data = ann_text, x = 1.96, xend = 1.96, y = 1.5, yend = -3.5, show.legend = FALSE,
               color = "black", arrow = arrow(length = unit(0.03, "npc")), size = .75) +
  theme_classic(base_size = 24, base_family = "") +
  theme(axis.ticks = element_line(color = "black"), axis.ticks.y.right = element_blank(),
        axis.line = element_line(linetype = "solid", colour = "black", lineend = "square"), axis.line.y = element_line(lineend = "butt"),
        axis.text = element_text(colour = "black"),
        plot.margin = unit(c(1,1,1,1), "lines"), panel.border = element_rect(fill = NA),
        legend.key.size = unit(1.75, "lines"), legend.spacing.y = unit(.25, "lines")) +
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(margin = unit(c(.5,.5,.5,.5), "lines"))) +
  scale_color_viridis(discrete=TRUE, name = "Diversity Metric", labels = c("range" = "Rangethrough", "raw" = "Raw Occurrences", "SQS_0.25" = "SQS (q = .25)", "SQS_0.5" = "SQS (q = .5)", "SQS_0.75" = "SQS (q = .75)"), guide = guide_legend(order = 1)) +
  scale_shape(name = "Latitudinal Band Type", labels = c("equal-area" = "Equal Area", "equal-width" = "Equal Width"), guide = guide_legend(order = 3)) +
  scale_linetype_manual(name = NULL, values = c("roy" = "solid", "theory" = "dashed"),
                        labels = c("roy" = "Roy et al 1998\nempirical estimate", "theory" = "Gillooly and Allen 2007\ntheoretical range"),
                        guide = guide_legend(keyheight=unit(2, "cm"), order = 3, override.aes = list(size = c(4,.5), color = c("grey80", "black")))) +
  guides(color = guide_legend(order = 1), shape = guide_legend(order = 2)) +
  facet_grid(inc_modern~raw_prop, labeller = labeller(inc_modern = as_labeller(c('yes' = 'Including~Modern', 'no' = 'Not~Including~Modern'), label_parsed),
                                                      raw_prop = as_labeller(c('prop' = 'Proportional~Diversity', 'raw' = 'Raw~Diversity'), label_parsed)))
ggsave("Allen Fossil Variations Slope.pdf", device = "pdf", width = 16, height = 10.5)

#Figure S11####
ann_text <- data.frame(inc_modern = "yes", num_bands = 24, raw_prop = "prop",
                       div_metric = "SQS_0.25", band_type = "equal-area", breakpoint = 3.54,
                       lab = "This Study")
ann_rect <- fossil_results %>%
  filter(inc_modern == "yes", num_bands == 24, raw_prop == "prop", div_metric == "SQS_0.25", band_type == "equal-area") %>%
  summarize(ymin = breakpoint - 1.96*breakpoint_se, ymax = breakpoint + 1.96*breakpoint_se)

ggplot(fossil_results, aes(x = as.character(num_bands), y = breakpoint, shape = band_type)) +
  scale_x_discrete(limits = as.character(sort(lat_bin_nums)), name = "Number of Latitudinal Bands") +
  scale_y_continuous(name = "Estimated Breakpoint (1000/K)",
                     sec.axis = sec_axis(~(1000/. - 273.15), name = expression("Estimated Breakpoint ("*degree*"C)"), labels = )) +
  geom_rect(xmin = 0, xmax = 50, ymin = ann_rect$ymin, ymax = ann_rect$ymax, color = "grey90", fill = "grey90") +
  geom_segment(x = 0, xend = 50, y = ann_rect$ymin, yend = ann_rect$ymin, aes(linetype = "this_study"), color = "grey90") +
  geom_point(aes(color = div_metric), size = 2, position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(color = div_metric, ymin = breakpoint - 1.96*breakpoint_se, ymax = breakpoint + 1.96*breakpoint_se), width = .8, position = position_dodge(width = 0.8)) +
  geom_text(data = ann_text, aes(label = lab), show.legend = FALSE,
            color = "black", hjust = .55, size = 5) +
  geom_segment(data = ann_text, x = 1.96, xend = 1.96, y = 3.53, yend = 3.49, show.legend = FALSE,
               color = "black", arrow = arrow(length = unit(0.03, "npc")), size = .75) +
  theme_classic(base_size = 24, base_family = "") +
  theme(axis.ticks = element_line(color = "black"),
        axis.line = element_line(linetype = "solid", colour = "black", lineend = "square"), axis.line.y = element_line(lineend = "butt"),
        axis.text = element_text(colour = "black"),
        plot.margin = unit(c(1,1,1,1), "lines"), panel.border = element_rect(fill = NA),
        legend.key.size = unit(1.75, "lines"), legend.spacing.y = unit(.25, "lines")) +
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(), strip.placement = "ouside", strip.text = element_text(margin = unit(c(.5,.5,.5,.5), "lines"))) +
  scale_color_viridis(discrete=TRUE, name = "Diversity Metric", labels = c("range" = "Rangethrough", "raw" = "Raw Occurrences", "SQS_0.25" = "SQS (q = .25)", "SQS_0.5" = "SQS (q = .5)", "SQS_0.75" = "SQS (q = .75)"), guide = guide_legend(order = 1)) +
  scale_shape(name = "Latitudinal Band Type", labels = c("equal-area" = "Equal Area", "equal-width" = "Equal Width"), guide = guide_legend(order = 3)) +
  scale_linetype_manual(name = NULL, values = c("this_study" = "solid"),
                        labels = c("this_study" = "Range for This Study"),
                        guide = guide_legend(keyheight=unit(2, "cm"), order = 3, override.aes = list(size = c(10), color = c("grey90")))) +
  guides(color = guide_legend(order = 1), shape = guide_legend(order = 2)) +
  facet_grid(inc_modern~raw_prop, labeller = labeller(inc_modern = as_labeller(c('yes' = 'Including~Modern', 'no' = 'Not~Including~Modern'), label_parsed),
                                                      raw_prop = as_labeller(c('prop' = 'Proportional~Diversity', 'raw' = 'Raw~Diversity'), label_parsed)))
ggsave("Allen Fossil Variations Breakpoint.pdf", device = "pdf", width = 16, height = 10.5)

#Plotting####
#Black, Orange, Sky Blue, Bluish green, Yellow, Blue, Vermillion, Reddish purple
colors8 <- c(rgb(0,0,0), rgb(230/255, 159/255, 0/255), rgb(86/255, 180/255, 233/255), rgb(0/255, 158/255, 115/255), rgb(240/255,228/255,66/255), rgb(0/255, 114/255, 178/255), rgb(213/255,94/255,0/255), rgb(204/255,121/255,167/255))
times <- c("modern", stages_abbr)
time_ages <- data.frame(stages_abbr, start = c(3.333, 13.82, 33.9, 47.8, 56, 59.2, 72.1, 83.6, 145),
                        end = c(2.58, 5.333, 27.82, 37.8, 47.8, 56, 66, 72.1, 125))

library(viridis)
#need 10 different labels, so split between 4 shapes and 3 colors
time_shapes <- setNames(c(rep_len(c(21, 22, 23, 24), 3), rep_len(c(21, 22, 23, 24), 4), rep_len(c(21, 22, 23, 24), 3)), times)
vir_colors <- viridis_pal()(3)
time_colors <- setNames(c(rep(vir_colors[1], 3), rep(vir_colors[2], 4), rep(vir_colors[3], 3)), times)
time_names <- setNames(c("Modern", stages), times)

library(deeptime)

#We like these settings
fossil_data_sub <- subset(fossil_data, num_bands == 24 & metric == "SQS_0.25" & band_type == "equal-area")

#Proportional diversity
ggplot(data = fossil_data_sub, aes(x = SST_K, y = div/div_tot, color = time)) +
  geom_point() +
  scale_x_continuous(name = "SST (K)") +
  scale_y_continuous(name = paste0("Percent Total Diversity"), trans = "log10") +
  coord_cartesian(ylim = c(.01,1), xlim = c(265,315)) +
  scale_color_manual(values = time_colors, labels = time_names, name = "Time Interval") +
  theme_bw()

#Raw diversity
ggplot(data = fossil_data_sub, aes(x = SST_K, y = div, color = time)) +
  geom_point() +
  scale_x_continuous(name = "SST (K)") +
  scale_y_continuous(name = paste0("Raw Diversity")) +
  coord_cartesian(xlim = c(265,315)) +
  scale_color_manual(values = time_colors, labels = time_names, name = "Time Interval") +
  theme_bw()

#Like Allen
#Proportional diversity
ggplot(data = fossil_data_sub, aes(x = Allen.SST, y = div/div_tot, color = time)) +
  geom_point() +
  scale_x_continuous(name = "SST (1000/K)") +
  scale_y_continuous(name = paste0("Percent Total Diversity")) +
  coord_trans(y = "log") +
  scale_color_manual(values = time_colors, labels = time_names, name = "Time Interval") +
  theme_bw()

#Raw diversity
ggplot(data = fossil_data_sub, aes(x = Allen.SST, y = div, color = time)) +
  geom_point() +
  scale_x_continuous(name = "SST (1000/K)") +
  scale_y_continuous(name = paste0("Raw Diversity")) +
  coord_trans(y = "log") +
  scale_color_manual(values = time_colors, labels = time_names, name = "Time Interval") +
  theme_bw()

#plot latitudinal diversity gradients (with regression lines)
#Proportional Diversity
ggplot(data = fossil_data_sub, aes(x = lat_mid, y = div/div_tot, color = time)) +
  geom_point() +
  geom_smooth(method=lm, formula = y ~ I(x^2) + I(x^4)) +
  scale_x_continuous(name = "Latitude") +
  scale_y_continuous(name = "Percent Total Diversity") +
  scale_color_manual(values = time_colors, labels = time_names, name = "Time Interval") +
  theme_bw()

#Raw diversity
ggplot(data = fossil_data_sub, aes(x = lat_mid, y = div, color = time)) +
  geom_point() +
  geom_smooth(method=lm, formula = y ~ I(x^2) + I(x^4)) +
  scale_x_continuous(name = "Latitude") +
  scale_y_continuous(name = "Raw Diversity") +
  scale_color_manual(values = time_colors, labels = time_names, name = "Time Interval") +
  theme_bw()

#or split using facet_wrap
#Proportional diversity
ggplot(data = fossil_data_sub, aes(x = lat_mid, y = div/div_tot)) +
  geom_point() +
  geom_smooth(method=lm, formula = y ~ I(x^2) + I(x^4)) +
  scale_x_continuous(name = "Latitude") +
  scale_y_continuous(name = "Percent Total Diversity") +
  theme_bw() +
  facet_wrap(~time, ncol = 1)

#Raw diversity
ggplot(data = fossil_data_sub, aes(x = lat_mid, y = div)) +
  geom_point() +
  geom_smooth(method=lm, formula = y ~ I(x^2) + I(x^4)) +
  scale_x_continuous(name = "Latitude") +
  scale_y_continuous(name = "Raw Diversity") +
  theme_bw() +
  facet_wrap(~time, ncol = 1)

#plot latitudinal temperature gradients (with regression lines)
ggplot(data = fossil_data_sub, aes(x = lat_mid, y = SST_K, color = time)) +
  #geom_point() +
  geom_smooth(method=lm, formula = y ~ ns(x, df = min(10, length(x) - 5))) +
  scale_x_continuous(name = "Latitude") +
  scale_y_continuous(name = "Sea Surface Temperature (K)") +
  scale_color_manual(values = time_colors, labels = time_names, name = "Time Interval") +
  theme_bw()

#or split using facet_wrap
ggplot(data = fossil_data_sub, aes(x = lat_mid, y = SST_K)) +
  #geom_point() +
  geom_smooth(method=lm, formula = y ~ ns(x, df = min(10, length(x) - 5))) +
  scale_x_continuous(name = "Latitude") +
  scale_y_continuous(name = "Sea Surface Temperature (K)") +
  theme_bw() +
  facet_wrap(~time, ncol = 1)

#Figure S4####
#Plot overlay of temperature and diversity by time bin
#Proportional diversity
# Adjust data for different y-axis scales
ylim.prim <- c(-13, 47) # Temperature range
ylim.sec <- c(-.10, 1) # Diversity proportion range

b <- diff(ylim.prim)/diff(ylim.sec)
a <- ylim.prim[1] - b*ylim.sec[1]

time_abbrvs <- setNames(names(time_names), time_names)
paleotemps$time <- time_abbrvs[paleotemps$Stage]
temp$time <- "modern"

time_names_with_ages <- setNames(c("Modern", paste0(stages, "\n(", time_ages$start, " - ", time_ages$end, " Ma)")), times)

ggplot(data = fossil_data_sub, aes(x = lat_mid)) +
  geom_smooth(aes(y = a + b*(div.prop), linetype = "dashed", color = "div", fill = "div"), method="gam", formula = y ~ s(x, k = 9), se = T, method.args = list(gamma = .5)) +
  geom_point(aes(y = a + b*(div.prop), shape = "Diversity", color = "div", fill = "div"), size = 3) +
  geom_smooth(data = subset(paleotemps, !is.na(time)), aes(x = Paleolatitude, y = SST_K - 273.15, linetype = "solid", color = "temp", fill = "temp"), method="gam", formula = y ~ s(x, k = 9), se = T, fullrange = T, method.args = list(gamma = 1)) +
  geom_point(data = subset(paleotemps, !is.na(time)), aes(x = Paleolatitude, y = SST_K - 273.15, shape = Marine.Terrestrial, color = "temp"), size = 3, alpha = .5) +
  geom_point(data = temp, aes(x = Latitude, y = SST_K - 273.15, color = "temp", fill = "temp"), alpha = .005, shape = 16) +
  geom_smooth(data = temp, aes(x = Latitude, y = SST_K - 273.15, linetype = "solid"), method="gam", formula = y ~ s(x, k = 9), color = "black", fill = "grey", se = T, fullrange = T) +
  #geom_quantile(aes(y = a + b*(div.prop)), formula = y ~ ns(x, df = 3), quantiles = c(.05, .95)) +
  scale_x_continuous(name = "Latitude", expand = c(0,0), limits = c(-80, 90)) +
  scale_y_continuous(expand = c(0,0), name = expression("Sea Surface Temperature ("*degree*"C)"),
                     sec.axis = sec_axis(~ (. - a)/b, name = "Proportion of Total SQS Genera")) +
  coord_cartesian(ylim = ylim.prim, xlim = c(-90, 90)) +
  scale_shape_manual(name = NULL, values = c("Marine" = 16, "Terrestrial" = 21, "Diversity" = 17),
                     labels = c("Marine" = "Marine\nTemperature Proxy", "Terrestrial" = "Terrestrial\nTemperature Proxy", "Diversity" = "Proportional\nSQS Diversity"),
                     limits = c("Marine", "Terrestrial", "Diversity")) +
  scale_linetype_manual(name = NULL, values = c("solid" = "solid", "dashed" = "dashed"),
                        labels = c("solid" = "Temperature Regression", "dashed" = "Proportional SQS Diversity\nRegression"),
                        limits = c("solid", "dashed")) +
  scale_color_manual(name = NULL, values = c("temp" = time_colors[[1]], "div" = time_colors[[5]]),
                     labels = c("temp" = "Temperature", "div" = "Proportional Diversity"),
                     limits = c("temp", "div")) +
  scale_fill_manual(name = NULL, values = c("temp" = time_colors[[1]], "div" = time_colors[[5]]),
                    labels = c("temp" = "Temperature", "div" = "Proportional SQS Diversity"),
                    limits = c("temp", "div")) +
  facet_wrap(~factor(time, levels = times), ncol = 3, labeller = as_labeller(time_names_with_ages)) +
  theme_classic(base_size = 24) +
  theme(axis.ticks = element_line(color = "black"), axis.line = element_blank(), axis.text = element_text(colour = "black"),
        plot.margin = unit(c(1,1,1,1), "lines"), panel.border = element_rect(fill = NA),
        legend.key.size = unit(1.75, "lines"), legend.spacing.y = unit(.35, "lines"),
        legend.position = c(.7, .125), legend.box.background = element_blank(),
        legend.background = element_blank(), legend.direction = "vertical", legend.box = "horizontal",
        axis.title.y.right = element_text(margin = unit(c(0,0,0,1), "lines")),
        panel.spacing = unit(1.5, "lines"), strip.background = element_blank(),
        strip.placement = "outside", strip.text = element_text(margin = unit(c(0,0,.5,0), "lines"))) +
  guides(shape=guide_legend(order = 1, keyheight = 3.5, default.unit="lines", override.aes=list(alpha=c(.25,.25,1), color = c(time_colors[[1]], time_colors[[1]], time_colors[[5]]))),
         linetype=guide_legend(byrow = TRUE, keyheight = 3, override.aes=list(alpha = .15, fill=c(time_colors[[1]], time_colors[[5]]), color = c(time_colors[[1]], time_colors[[5]]))),
         color = FALSE, fill = FALSE)
ggsave("Temperature and Diversity Overlay (Prop).pdf", width = 15, height = 18)

#Figure 1####
#Modern and Early Cretaceous only
ggplot(data = subset(fossil_data_sub, time %in% c("modern", "berr_barr")), aes(x = lat_mid)) +
  geom_smooth(aes(y = a + b*(div.prop), linetype = "dashed", color = "div", fill = "div"), method="gam", formula = y ~ s(x, k = 9), se = T, method.args = list(gamma = .5)) +
  geom_smooth(data = subset(paleotemps, time %in% c("modern", "berr_barr")), aes(x = Paleolatitude, y = SST_K - 273.15, linetype = "solid", color = "temp", fill = "temp"), method="gam", formula = y ~ s(x, k = 9), se = T, fullrange = T) +
  geom_smooth(data = temp, aes(x = Latitude, y = SST_K - 273.15, linetype = "solid", color = "temp"), method="gam", formula = y ~ s(x), se = T, fullrange = T) +
  #geom_quantile(aes(y = a + b*(div.prop)), formula = y ~ ns(x, df = 3), quantiles = c(.05, .95)) +
  scale_x_continuous(name = "Latitude", expand = c(0,0), limits = c(-80, 90)) +
  scale_y_continuous(expand = c(0,0), name = expression("Sea Surface Temperature ("*degree*"C)"),
                     sec.axis = sec_axis(~ (. - a)/b, name = "Proportion of Total SQS Genera")) +
  coord_cartesian(ylim = ylim.prim, xlim = c(-90, 90)) +
  scale_linetype_manual(name = NULL, values = c("solid" = "solid", "dashed" = "dashed"),
                        labels = c("solid" = "Temperature", "dashed" = "Proportional SQS Diversity"),
                        limits = c("solid", "dashed")) +
  scale_color_manual(name = NULL, values = c("temp" = time_colors[[1]], "div" = time_colors[[5]]),
                     labels = c("temp" = "Temperature", "div" = "Proportional SQS Diversity"),
                     limits = c("temp", "div")) +
  scale_fill_manual(name = NULL, values = c("temp" = time_colors[[1]], "div" = time_colors[[5]]),
                    labels = c("temp" = "Temperature", "div" = "Proportional SQS Diversity"),
                    limits = c("temp", "div")) +
  facet_wrap(~factor(time, levels = times), ncol = 3,
             labeller = as_labeller(c("modern" = "Icehouse Example\n(Modern)",
                                      "berr_barr" = "Greenhouse Example\n(Berriasian-Barremian,\n145 - 125 Ma)"))) +
  theme_classic(base_size = 24) +
  theme(axis.ticks = element_line(color = "black"), axis.line = element_blank(), axis.text = element_text(colour = "black"),
        plot.margin = unit(c(1,1,1,1), "lines"), panel.border = element_rect(fill = NA),
        legend.position = "bottom", legend.box.margin = margin(0,0,0,0), legend.key.width = unit(4, "lines"),
        legend.margin = margin(0,0,0,0), strip.background = element_blank(),
        strip.placement = "outside", strip.text = element_text(size = 28, margin = unit(c(0,0,.5,0), "lines")),
        axis.title.y.right = element_text(margin = unit(c(0,0,0,1), "lines"))) +
  guides(linetype=guide_legend(override.aes=list(alpha = .15, fill=c(time_colors[[1]], time_colors[[5]]), color = c(time_colors[[1]], time_colors[[5]]))))
ggsave("Temperature and Diversity Overlay Comparison.pdf", width = 15, height = 7.5)

#Figure S5####
#Raw diversity
# Adjust data for different y-axis scales
ylim.prim <- c(-13, 47) # Temperature range
ylim.sec <- c(-5, 85) # Diversity raw range

b <- diff(ylim.prim)/diff(ylim.sec)
a <- ylim.prim[1] - b*ylim.sec[1]

ggplot(data = fossil_data_sub, aes(x = lat_mid)) +
  geom_smooth(aes(y = a + b*(div), linetype = "dashed", color = "div", fill = "div"), method="gam", formula = y ~ s(x, k = 9), se = T, method.args = list(gamma = .5)) +
  geom_point(aes(y = a + b*(div), shape = "Diversity", color = "div", fill = "div"), size = 3) +
  geom_smooth(data = subset(paleotemps, !is.na(time)), aes(x = Paleolatitude, y = SST_K - 273.15, linetype = "solid", color = "temp", fill = "temp"), method="gam", formula = y ~ s(x, k = 9), se = T, fullrange = T, method.args = list(gamma = 1)) +
  geom_point(data = subset(paleotemps, !is.na(time)), aes(x = Paleolatitude, y = SST_K - 273.15, shape = Marine.Terrestrial, color = "temp"), size = 3, alpha = .5) +
  geom_point(data = temp, aes(x = Latitude, y = SST_K - 273.15, color = "temp", fill = "temp"), alpha = .005, shape = 16) +
  geom_smooth(data = temp, aes(x = Latitude, y = SST_K - 273.15, linetype = "solid"), method="gam", formula = y ~ s(x), color = "black", fill = "grey", se = T, fullrange = T) +
  #geom_quantile(aes(y = a + b*(div)), formula = y ~ ns(x, df = 3), quantiles = c(.05, .95)) +
  scale_x_continuous(name = "Latitude", expand = c(0,0), limits = c(-80, 90)) +
  scale_y_continuous(expand = c(0,0), name = expression("Sea Surface Temperature ("*degree*"C)"),
                     sec.axis = sec_axis(~ (. - a)/b, name = "# SQS Genera")) +
  coord_cartesian(ylim = ylim.prim, xlim = c(-90, 90)) +
  scale_shape_manual(name = NULL, values = c("Marine" = 16, "Terrestrial" = 21, "Diversity" = 17),
                     labels = c("Marine" = "Marine\nTemperature Proxy", "Terrestrial" = "Terrestrial\nTemperature Proxy", "Diversity" = "SQS Diversity"),
                     limits = c("Marine", "Terrestrial", "Diversity")) +
  scale_linetype_manual(name = NULL, values = c("solid" = "solid", "dashed" = "dashed"),
                        labels = c("solid" = "Temperature Regression", "dashed" = "SQS Diversity Regression"),
                        limits = c("solid", "dashed")) +
  scale_color_manual(name = NULL, values = c("temp" = time_colors[[1]], "div" = time_colors[[5]]),
                     labels = c("temp" = "Temperature", "div" = "SQS Diversity"),
                     limits = c("temp", "div")) +
  scale_fill_manual(name = NULL, values = c("temp" = time_colors[[1]], "div" = time_colors[[5]]),
                    labels = c("temp" = "Temperature", "div" = "SQS Diversity"),
                    limits = c("temp", "div")) +
  facet_wrap(~factor(time, levels = times), ncol = 3, labeller = as_labeller(time_names_with_ages)) +
  theme_classic(base_size = 24) +
  theme(axis.ticks = element_line(color = "black"), axis.line = element_blank(), axis.text = element_text(colour = "black"),
        plot.margin = unit(c(1,1,1,1), "lines"), panel.border = element_rect(fill = NA),
        legend.key.size = unit(1.75, "lines"), legend.spacing.y = unit(.35, "lines"),
        legend.position = c(.7, .125), legend.box.background = element_blank(),
        legend.background = element_blank(), legend.direction = "vertical", legend.box = "horizontal",
        axis.title.y.right = element_text(margin = unit(c(0,0,0,1), "lines")),
        panel.spacing = unit(1.5, "lines"), strip.background = element_blank(),
        strip.placement = "outside", strip.text = element_text(margin = unit(c(0,0,.5,0), "lines"))) +
  guides(shape=guide_legend(order = 1, keyheight = 3.5, default.unit="lines", override.aes=list(alpha=c(.25,.25,1), color = c(time_colors[[1]], time_colors[[1]], time_colors[[5]]))),
         linetype=guide_legend(byrow = TRUE, keyheight = 3, override.aes=list(alpha = .15, fill=c(time_colors[[1]], time_colors[[5]]), color = c(time_colors[[1]], time_colors[[5]]))),
         color = FALSE, fill = FALSE)
ggsave("Temperature and Diversity Overlay (Raw).pdf", width = 15, height = 18)

#Find segments
library(segmented)
reg <- lm(div/div_tot ~ SST_K, data = fossil_data_sub)
reg_seg <- with(fossil_data_sub, segmented(reg, seg.Z = ~ SST_K))

# manual lm's
brk <- reg_seg$psi[2:3]
lm_1 <- lm(div/div_tot ~ SST_K, data = subset(fossil_data_sub,SST_K < brk[1]))
summary(lm_1)
lm_2 <- lm(div/div_tot ~ SST_K, data = subset(fossil_data_sub,SST_K > brk[1]))
summary(lm_2)

#plot regression based on breakpoint
ggplot() +
  geom_rect(aes(xmin=brk[1]-1.96*brk[2], xmax=brk[1]+1.96*brk[2], ymin=-5, ymax=5), fill="black", alpha=0.2) +
  geom_vline(xintercept = brk[1]) +
  geom_point(data = fossil_data_sub, aes(x = SST_K, y = div/div_tot, color = time)) +
  geom_smooth(data = subset(fossil_data_sub, SST_K < brk[1]), aes(x = SST_K, y = div/div_tot), method=lm, formula = y ~ x, se = FALSE) +
  geom_smooth(data = subset(fossil_data_sub, SST_K > brk[1]), aes(x = SST_K, y = div/div_tot), method=lm, formula = y ~ x, se = FALSE) +
  scale_x_continuous(name = "SST (1000/K)") +
  scale_y_continuous(name = "Percent Total Diversity") +
  coord_cartesian(ylim = c(0,1)) +
  theme_bw()

#plot quantile linear regressions based on breakpoint
ggplot(data = fossil_data_sub, aes(x = SST_K, y = div/div_tot)) +
  geom_point(aes(color = time)) +
  geom_quantile(quantiles = c(.05,.95), aes(group = SST_K > brk[1])) +
  scale_x_continuous(name = "SST (K)") +
  scale_y_continuous(name = "Percent Total Diversity") +
  theme_bw()

#Figure 2####
#breakpoint regression
library(segmented)
library(broom)
reg <- lm(log(div.prop) ~ SST_K, data = fossil_data_sub)
glance(reg)
davies.test(reg, seg.Z = ~ SST_K)
reg_seg <- with(fossil_data_sub, segmented(reg, seg.Z = ~ SST_K, npsi = 1))
seg_pred <- as.data.frame(cbind(SST_K = seq(267, 313, .01), predict(reg_seg, newdata = data.frame(SST_K = seq(267, 313, .01)), interval = "confidence", level = .95)))
#slopes on both sides of the breakpoint
slope(reg_seg)
#use broom to get p-value and r-sqaured
glance(reg_seg)

brk <- reg_seg$psi[2:3]
print(paste0("Linear breakpoint at ", round(brk[1], 1), "K (", round(brk[1] - 273.15, 1), "C) +- ", round(brk[2] * 1.96, 1)), quote = FALSE)

#3rd order polynomial regression
reg <- lm(log(div.prop) ~ poly(SST_K, 3), data = fossil_data_sub)
pred <- as.data.frame(cbind(SST_K = seq(267, 313, .01), predict(reg, newdata = data.frame(SST_K = seq(267, 313, .01)), interval = "confidence", level = .95)))
print(paste0("3rd order polynomial peak at ", round(pred$SST_K[which.max(pred$fit)], 1), "K (",
             round(pred$SST_K[which.max(pred$fit)] - 273.15, 1), "C)"), quote = FALSE)
glance(reg)

#gam regression
gam_fit <- gam(log(div.prop) ~ s(SST_K), data = fossil_data_sub)
gam_pred <- cbind(SST_K = seq(267, 313, .01), as.data.frame(predict(gam_fit, newdata = data.frame(SST_K = seq(267, 313, .01)), se.fit = TRUE)))
print(paste("gam peak at", round(gam_pred$SST_K[which.max(gam_pred$fit)], 1), "K"), quote = FALSE)

library(viridis)
ggplot() +
  geom_rect(aes(xmin=(brk[1] + 1.96 * brk[2]), xmax=(brk[1] - 1.96 * brk[2]), ymin=0.001, ymax=1), fill="white", linetype = "dashed", color = "grey40", size = 1, inherit.aes = FALSE) +
  geom_vline(xintercept = brk[1], size = 1) +
  geom_point(data = fossil_data_sub, aes(x = SST_K, y = div.prop, fill = time, shape = time), size = 3) +
  #geom_smooth(data = fossil_data_sub, aes(x = SST_K, y = div.prop), color = "black", size = 1,
  #            method="gam", formula = y ~ s(x, k = 9), se = T, method.args = list(gamma = .5)) +
  geom_ribbon(data = seg_pred, aes(x = SST_K, ymin = exp(lwr), ymax = exp(upr)), alpha = .25) +
  geom_line(data = seg_pred, aes(x = SST_K, y = exp(fit)), size = 1) +
  scale_shape_manual(values = time_shapes, labels = gsub("\n", " ", time_names_with_ages), name = NULL,
                     guide = guide_legend(nrow = 5, byrow = TRUE, label.hjust = 0)) +
  scale_fill_manual(values = time_colors, labels = gsub("\n", " ", time_names_with_ages), name = NULL,
                    guide = guide_legend(nrow = 5, byrow = TRUE, label.hjust = 0)) +
  scale_x_continuous(name = expression("Sea Surface Temperature ("*degree*"C)"), expand = c(0,0),
                     breaks = c(0, 10, 20, 30) + 273.15, labels = c(0, 10, 20, 30), sec.axis = sec_axis(~.)) +
  scale_y_continuous(name = "Percent Total SQS Generic Diversity", expand = c(0,0), sec.axis = sec_axis(~.), trans = "log",
                     breaks = c(0.01, 0.03, 0.1, 0.3, 1), labels = c(0.01, 0.03, 0.1, 0.3, 1) * 100) +
  coord_cartesian(xlim = c(-3, 38) + 273.15, ylim = c(0.0075, 1)) +
  theme_classic(base_size = 24) +
  theme(axis.ticks = element_line(color = "black"), axis.text = element_text(colour = "black"),
        axis.line = element_line(color = "black", lineend = "square"), 
        axis.ticks.x.top = element_blank(), axis.text.x.top = element_blank(),
        axis.ticks.y.right = element_blank(), axis.text.y.right = element_blank(),
        plot.margin = margin(1,1,1,1, "lines"), panel.border = element_rect(fill = NA),
        legend.key.size = unit(.25, "lines"), legend.spacing.y = unit(0.35, "lines"), legend.spacing.x = unit(3, "lines"),
        legend.position = "bottom", legend.margin = margin(0,8,0,0, "lines"))
ggsave("Diversity and Temperature with Breakpoint and Regression.pdf", device = "pdf", width = 12, height = 12)

#Figure S6####
#With Raw diversity
reg <- lm(log(div) ~ SST_K, data = fossil_data_sub)
glance(reg)
davies.test(reg, seg.Z = ~ SST_K)
reg_seg <- with(fossil_data_sub, segmented(reg, seg.Z = ~ SST_K))
seg_pred <- as.data.frame(cbind(SST_K = seq(267, 313, .01), predict(reg_seg, newdata = data.frame(SST_K = seq(267, 313, .01)), interval = "confidence", level = .95)))
brk <- reg_seg$psi[2:3]
print(paste0("Linear breakpoint at ", round(brk[1], 1), "K (", round(brk[1] - 273.15, 1), "C) +- ", round(brk[2] * 1.96, 1)), quote = FALSE)

#3rd order polynomial regression
reg <- lm(log(div) ~ poly(SST_K, 3), data = fossil_data_sub)
pred <- as.data.frame(cbind(SST_K = seq(267, 313, .01), predict(reg, newdata = data.frame(SST_K = seq(267, 313, .01)), interval = "confidence", level = .95)))
print(paste0("3rd order polynomial peak at ", round(pred$SST_K[which.max(pred$fit)], 1), "K (",
             round(pred$SST_K[which.max(pred$fit)] - 273.15, 1), "C)"), quote = FALSE)
glance(reg)

#gam regression
gam_fit <- gam(log(div) ~ s(SST_K), data = fossil_data_sub)
gam_pred <- cbind(SST_K = seq(267, 313, .01), as.data.frame(predict(gam_fit, newdata = data.frame(SST_K = seq(267, 313, .01)), se.fit = TRUE)))
print(paste("gam peak at", round(gam_pred$SST_K[which.max(gam_pred$fit)], 1), "K"), quote = FALSE)

library(viridis)
ggplot() +
  geom_rect(aes(xmin=brk[1]-1.96*brk[2], xmax=brk[1]+1.96*brk[2], ymin=0.01, ymax=120), fill="white", linetype = "dashed", color = "grey40", size = 1, inherit.aes = FALSE) +
  geom_vline(xintercept = brk[1], size = 1) +
  geom_point(data = fossil_data_sub, aes(x = SST_K, y = div, fill = time, shape = time), size = 3) +
  #geom_smooth(data = fossil_data_sub, aes(x = SST_K, y = div), color = "black", size = 1,
  #            method="gam", formula = y ~ s(x, k = 9), se = T, method.args = list(gamma = .5)) +
  geom_ribbon(data = seg_pred, aes(x = SST_K, ymin = exp(lwr), ymax = exp(upr)), alpha = .25) +
  geom_line(data = seg_pred, aes(x = SST_K, y = exp(fit)), size = 1) +
  scale_shape_manual(values = time_shapes, labels = gsub("\n", " ", time_names_with_ages), name = NULL,
                     guide = guide_legend(nrow = 5, byrow = TRUE, label.hjust = 0)) +
  scale_fill_manual(values = time_colors, labels = gsub("\n", " ", time_names_with_ages), name = NULL,
                    guide = guide_legend(nrow = 5, byrow = TRUE, label.hjust = 0)) +
  scale_x_continuous(name = expression("Sea Surface Temperature ("*degree*"C)"), expand = c(0,0),
                     breaks = c(0, 10, 20, 30) + 273.15, labels = c(0, 10, 20, 30), sec.axis = sec_axis(~.)) +
  scale_y_continuous(name = "Raw SQS Generic Diversity (q = .25)", expand = c(0,0), sec.axis = sec_axis(~.),
                     trans = "log", breaks = c(1, 3, 10, 30, 100)) +
  coord_cartesian(xlim = c(-3, 38) + 273.15, ylim = c(.5, 100)) +
  theme_classic(base_size = 24) +
  theme(axis.ticks = element_line(color = "black"), axis.text = element_text(colour = "black"),
        axis.line = element_line(color = "black", lineend = "square"), 
        axis.ticks.x.top = element_blank(), axis.text.x.top = element_blank(),
        axis.ticks.y.right = element_blank(), axis.text.y.right = element_blank(),
        plot.margin = margin(1,1,1,1, "lines"), panel.border = element_rect(fill = NA),
        legend.key.size = unit(.25, "lines"), legend.spacing.y = unit(0.35, "lines"), legend.spacing.x = unit(3, "lines"),
        legend.position = "bottom", legend.margin = margin(0,8,0,0, "lines"))
ggsave("Raw Diversity and Temperature with Breakpoint and Regression.pdf", device = "pdf", width = 12, height = 12)

#proportional diversity without modern
# 3rd order polynomial regression
reg <- lm(div.prop ~ SST_K, data = subset(fossil_data_sub, time != "modern"))
reg_seg <- with(subset(fossil_data_sub, time != "modern"), segmented(reg, seg.Z = ~ SST_K))
brk <- reg_seg$psi[2:3]

reg <- lm(div.prop ~ poly(SST_K, 3), data = subset(fossil_data_sub, time != "modern"))

pred <- cbind(SST_K = seq(265, 315, .1), as.data.frame(predict(reg, newdata = data.frame(SST_K = seq(265, 315, .1)), interval = "confidence")))

ggplot() +
  geom_rect(aes(xmin=brk[1]-1.96*brk[2], xmax=brk[1]+1.96*brk[2], ymin=-5, ymax=5), fill="white", linetype = "dashed", color = "grey40", size = 1, inherit.aes = FALSE) +
  geom_vline(xintercept = brk[1], size = 1) +
  geom_point(data = subset(fossil_data_sub, time != "modern"), aes(x = SST_K, y = div/div_tot, color = time), size = 3) +
  geom_ribbon(data = pred, aes(x = SST_K, ymin = lwr, ymax = upr), fill = "black", alpha = .2) +
  geom_line(data = pred, aes(x = SST_K, y = fit), color = "black", size = 1) +
  scale_color_manual(values = time_colors, labels = time_names, name = NULL,
                     guide = guide_legend(nrow = 3, byrow = TRUE, label.hjust = 0)) +
  scale_x_continuous(name = "SST (K)", expand = c(0,0)) +
  scale_y_continuous(name = "Percent Total Diversity") +
  coord_cartesian(ylim = c(0,.8), xlim = c(265, 315)) +
  theme_classic(base_size = 24) +
  theme(axis.ticks = element_line(color = "black"), axis.text = element_text(colour = "black"),
        axis.line = element_line(color = "black", lineend = "square"), 
        axis.ticks.x.top = element_blank(), axis.text.x.top = element_blank(),
        axis.ticks.y.right = element_blank(), axis.text.y.right = element_blank(),
        plot.margin = margin(1,1,1,1, "lines"), panel.border = element_rect(fill = NA),
        legend.key.size = unit(.25, "lines"), legend.spacing.y = unit(0.25, "lines"),
        legend.position = "top", legend.margin = margin(0,0,0,0, "lines"))
ggsave("Diversity and Temperature with Breakpoint and Regression without Modern.pdf", device = "pdf", width = 12, height = 12)

#Figure S7####
#Breakpoint Bootstrap
#bootstrap observations to determine uncertainty of breakpoint
library(segmented)
nreps <- 1000
breaks <- matrix(NA, nreps, 2)
set.seed(42)

fossil_data_sub_clean <- subset(fossil_data_sub, !is.na(div))
for(i in 1:nreps){
  smpl <- sample(seq(1:nrow(fossil_data_sub_clean)), nrow(fossil_data_sub_clean), replace = TRUE)
  btstrp <- fossil_data_sub_clean[smpl,]
  reg <- lm(log(div.prop) ~ SST_K, data = btstrp)
  tryCatch({
    reg_seg <- with(btstrp, segmented(reg, seg.Z = ~ SST_K))
    breaks[i,] <- reg_seg$psi[2:3]
  }, error = function(e){
  }, finally = {
  })
}

library(Hmisc)
weights <- 1/(breaks[breaks[,1] < 300,2]^2)
mean.wtd <- wtd.mean(breaks[breaks[,1] < 300,1], weights)
stddev <- sqrt(wtd.var(breaks[breaks[,1] < 300,1], weights))

reg <- lm(log(div.prop) ~ SST_K, data = fossil_data_sub)
reg_seg <- with(fossil_data_sub, segmented(reg, seg.Z = ~ SST_K, psi = mean.wtd,
                control = seg.control(it.max = 1, maxit.glm = 0, n.boot=0)))
seg_pred <- as.data.frame(cbind(SST_K = seq(267, 313, .01), predict(reg_seg, newdata = data.frame(SST_K = seq(267, 313, .01)), interval = "confidence", level = .95)))

ggplot() +
  geom_rect(aes(xmin=(mean.wtd + 1.96 * stddev), xmax=(mean.wtd - 1.96 * stddev), ymin=0.001, ymax=1), fill="white", linetype = "dashed", color = "grey40", size = 1, inherit.aes = FALSE) +
  geom_vline(xintercept = mean.wtd, size = 1) +
  geom_point(data = fossil_data_sub, aes(x = SST_K, y = div.prop, fill = time, shape = time), size = 3) +
  #geom_smooth(data = fossil_data_sub, aes(x = SST_K, y = div.prop), color = "black", size = 1,
  #            method="gam", formula = y ~ s(x, k = 9), se = T, method.args = list(gamma = .5)) +
  geom_ribbon(data = seg_pred, aes(x = SST_K, ymin = exp(lwr), ymax = exp(upr)), alpha = .25) +
  geom_line(data = seg_pred, aes(x = SST_K, y = exp(fit)), size = 1) +
  scale_shape_manual(values = time_shapes, labels = gsub("\n", " ", time_names_with_ages), name = NULL,
                     guide = guide_legend(nrow = 5, byrow = TRUE, label.hjust = 0)) +
  scale_fill_manual(values = time_colors, labels = gsub("\n", " ", time_names_with_ages), name = NULL,
                    guide = guide_legend(nrow = 5, byrow = TRUE, label.hjust = 0)) +
  scale_x_continuous(name = expression("Sea Surface Temperature ("*degree*"C)"), expand = c(0,0),
                     breaks = c(0, 10, 20, 30) + 273.15, labels = c(0, 10, 20, 30), sec.axis = sec_axis(~.)) +
  scale_y_continuous(name = "Percent Total SQS Generic Diversity", expand = c(0,0), sec.axis = sec_axis(~.), trans = "log",
                     breaks = c(0.01, 0.03, 0.1, 0.3, 1), labels = c(0.01, 0.03, 0.1, 0.3, 1) * 100) +
  coord_cartesian(xlim = c(-3, 38) + 273.15, ylim = c(0.0075, 1)) +
  theme_classic(base_size = 24) +
  theme(axis.ticks = element_line(color = "black"), axis.text = element_text(colour = "black"),
        axis.line = element_line(color = "black", lineend = "square"), 
        axis.ticks.x.top = element_blank(), axis.text.x.top = element_blank(),
        axis.ticks.y.right = element_blank(), axis.text.y.right = element_blank(),
        plot.margin = margin(1,1,1,1, "lines"), panel.border = element_rect(fill = NA),
        legend.key.size = unit(.25, "lines"), legend.spacing.y = unit(0.35, "lines"), legend.spacing.x = unit(3, "lines"),
        legend.position = "bottom", legend.margin = margin(0,8,0,0, "lines"))
ggsave("Diversity and Temperature with Bootstrap Breakpoint.pdf", device = "pdf", width = 12, height = 12)

#subset analysis####
library(plyr)
n_subsets <- 1000
subset_sizes <- seq(10,40,10)

fossil_data_sub_clean$SST_bin <- cut(fossil_data_sub_clean$SST_K -273.15, breaks = c(-10, 10, 20, 30, 50),
                                     labels = c("< 10", "10 - 20", "20 - 30", "> 30"))

SST_bins <- ddply(fossil_data_sub_clean, .(SST_bin), summarise, N = length(lat_bin))
SST_bins <- subset(SST_bins, N >= 5)
n_bins <- length(SST_bins$SST_bin)

subset_means <- matrix(NA, n_subsets, n_bins)
colnames(subset_means) <- SST_bins$SST_bin

#sampling with replacement
subset_results <- data.frame(SST_bin = rep(SST_bins$SST_bin, length(subset_sizes)), sample = rep(subset_sizes, each = n_bins), mean = NA, std_dev = NA)
set.seed(42)

for(k in 1:length(subset_sizes)){
  subset_size <- subset_sizes[k]
  subset_means[,] <- NA
  for(j in 1:ncol(subset_means)){
    dat <- subset(fossil_data_sub_clean, SST_bin == colnames(subset_means)[j])
    for(i in 1:n_subsets){
      subset_means[i,j] <- mean(sample(log(dat$div.prop), subset_size, replace = TRUE))
    }
  }
  subset_results$mean[((k-1)*n_bins + 1):(k*n_bins)] <- colMeans(subset_means)
  subset_results$std_dev[((k-1)*n_bins + 1):(k*n_bins)] <- unlist(colwise(sd)(as.data.frame(subset_means)))
}

ggplot(data = subset_results, aes(x = SST_bin, y = mean, color = sample, group = sample)) +
  geom_point(position=position_dodge(width = 1)) +
  geom_linerange(aes(ymin = mean - 1.96*std_dev, ymax = mean + 1.96*std_dev), position=position_dodge(width = 1)) +
  scale_x_discrete(name = expression("Sea Surface Temperature Bin ("*degree*"C)")) +
  scale_y_continuous(name = "Percent Total SQS Generic Diversity", expand = c(0,0)) +
  theme_classic(base_size = 24) +
  theme(axis.ticks = element_line(color = "black"), axis.text = element_text(colour = "black"),
        axis.line = element_line(color = "black", lineend = "square"), 
        axis.ticks.x.top = element_blank(), axis.text.x.top = element_blank(),
        axis.ticks.y.right = element_blank(), axis.text.y.right = element_blank(),
        plot.margin = margin(1,1,1,1, "lines"), panel.border = element_rect(fill = NA),
        legend.key.size = unit(.25, "lines"), legend.spacing.y = unit(0.25, "lines"),
        legend.position = "right", legend.margin = margin(0,0,0,0, "lines"))
ggsave("Diversity and Temperature Subsetting with Replacement.pdf", device = "pdf", width = 12, height = 12)

#Figure S8####
#sampling without replacement
subset_sizes <- seq(10,40,10)
subset_results <- data.frame(SST_bin = rep(SST_bins$SST_bin, length(subset_sizes)), sample = rep(subset_sizes, each = n_bins), mean = NA, std_dev = NA)
set.seed(42)

for(k in 1:length(subset_sizes)){
  subset_size <- subset_sizes[k]
  subset_means[,] <- NA
  for(j in 1:ncol(subset_means)){
    dat <- subset(fossil_data_sub_clean, SST_bin == colnames(subset_means)[j])
    if(nrow(dat) < subset_size){
      next
    }
    else{
      for(i in 1:n_subsets){
        subset_means[i,j] <- mean(sample(log(dat$div.prop), subset_size, replace = FALSE))
      }
    }
  }
  subset_results$mean[((k-1)*n_bins + 1):(k*n_bins)] <- colMeans(subset_means)
  subset_results$std_dev[((k-1)*n_bins + 1):(k*n_bins)] <- unlist(colwise(sd)(as.data.frame(subset_means)))
}

ggplot(data = subset_results, aes(x = SST_bin, y = exp(mean), color = factor(sample))) +
  geom_point(position=position_dodge(width = .9)) +
  geom_linerange(aes(ymin = exp(mean - 1.96*std_dev), ymax = exp(mean + 1.96*std_dev)), position=position_dodge(width = .9)) +
  geom_text(data = SST_bins, aes(x = SST_bin, label = paste("n =", N)), size = 6, y = -.8, inherit.aes = FALSE) +
  scale_x_discrete(name = expression("Sea Surface Temperature ("*degree*"C)")) +
  scale_y_continuous(name = "Percent Total SQS Generic Diversity", sec.axis = sec_axis(~.), trans = "log",
                     breaks = c(.03, .1, .3), labels = c(3, 10, 30), lim = c(.025, .45)) +
  scale_color_viridis(name = "Sample Size", discrete = T, breaks = rev(subset_sizes)) +
  theme_classic(base_size = 24) +
  theme(axis.ticks = element_line(color = "black"), axis.text = element_text(colour = "black"),
        axis.line = element_line(color = "black", lineend = "square"), 
        axis.ticks.x.top = element_blank(), axis.text.x.top = element_blank(),
        axis.ticks.y.right = element_blank(), axis.text.y.right = element_blank(),
        plot.margin = margin(1,1,1,1, "lines"), panel.border = element_rect(fill = NA),
        legend.position = "right", legend.margin = margin(0,0,0,0, "lines"),
        legend.key.height = unit(4, "lines"))
ggsave("Diversity and Temperature Subsetting without Replacement.pdf", device = "pdf", width = 12, height = 10)

#fit normal distribution to all diversity data
library(MASS)
normfit <- fitdistr(log(fossil_data_sub_clean$div.prop), "normal")
ks.test(log(fossil_data_sub_clean$div.prop), "pnorm", normfit$estimate)
hist(log(fossil_data_sub_clean$div.prop), breaks = 10, prob = TRUE)
curve(dnorm(x, mean = normfit$estimate[1], sd = normfit$estimate[2]), col = "red", add = TRUE)

#Figure S9####
#bin diversity data
library(MASS)
library(plyr)

bin_width <- 10 #width of the SST bins

fossil_data_sub_clean$SST_bin <- cut(fossil_data_sub_clean$SST_K -273.15, breaks = c(-10, 10, 20, 30, 50),
                                     labels = c("< 10", "10 - 20", "20 - 30", "> 30"))

SST_bins <- ddply(fossil_data_sub_clean, .(SST_bin), summarise, N = length(lat_bin))
SST_bins <- subset(SST_bins, N >= 5)
n_bins <- length(SST_bins$SST_bin)

#plot diversity histograms for each SST bin
ggplot(subset(fossil_data_sub_clean, SST_bin %in% SST_bins$SST_bin)) +
  geom_histogram(aes(x=log(div.prop)), bins = 10) +
  facet_wrap(~SST_bin, ncol = 1)

#simulate data for binned SST using single normal function from entire logged diversity database
normfit <- fitdistr(log(fossil_data_sub_clean$div.prop), "normal")

#how to calculate some statistic that can be summarized over many simulations?
n_sims <- 1000
sim_data <- data.frame(SST_bin=rep(SST_bins$SST_bin, n_sims), diff_mean=rep(NA, n_sims * n_bins), diff_median=rep(NA, n_sims * n_bins), diff_95=rep(NA, n_sims * n_bins), diff_5=rep(NA, n_sims * n_bins), diff_range=rep(NA, n_sims*n_bins))

for(n in 1:n_sims){
  for(i in 1:n_bins){
    sub_data <- subset(fossil_data_sub_clean, SST_bin == SST_bins$SST_bin[i])
    y_obs <- log(sub_data$div.prop)
    y_sim <- rnorm(n = SST_bins$N[i], mean = normfit$estimate[1], sd = normfit$estimate[2])
    sim_data$diff_mean[n_bins*(n-1) + i] <- quantile(y_obs, .5) - quantile(y_sim, .5)
    sim_data$diff_median[n_bins*(n-1) + i] <- median(y_obs) - median(y_sim)
    sim_data$diff_95[n_bins*(n-1) + i] <- quantile(y_obs, .95) - quantile(y_sim, .95)
    sim_data$diff_5[n_bins*(n-1) + i] <- quantile(y_obs, .05) - quantile(y_sim, .05)
    sim_data$diff_range[n_bins*(n-1) + i] <- diff(range(y_obs)) - diff(range(y_sim))
  }
}

tidy_sim_data <- tidyr::gather(sim_data, measure, value, -SST_bin)
tidy_sim_data$measure <- factor(tidy_sim_data$measure, levels = c('diff_mean','diff_median','diff_95','diff_5','diff_range'))

ggplot(data = tidy_sim_data, aes(x = SST_bin, y = value, fill = measure)) +
  geom_hline(yintercept = 0, size = 1, linetype = "dashed", color = "grey") +
  geom_violin() +
  scale_fill_manual(name = NULL, values = c('diff_mean'='black','diff_median'='grey20','diff_95'='grey40','diff_5'='grey60','diff_range'='grey80'),
                    labels = c('diff_mean'='mean','diff_median'='median','diff_95'='95 percentile','diff_5' = '5 percentile', 'diff_range'='range')) +
  scale_y_continuous(name = "Observed - Simulated") +
  scale_x_discrete(name = expression("Sea Surface Temperature Bin ("*degree*"C)")) +
  theme_classic(base_size = 24) +
  theme(axis.ticks = element_line(color = "black"), axis.ticks.y.right = element_blank(),
        axis.line = element_blank(),
        axis.text = element_text(colour = "black"),
        plot.margin = margin(1,1,1,1, "lines"), panel.border = element_rect(fill = NA))
ggsave("Universal Normal Dist vs Binned Raw Data.pdf", device = "pdf", width = 12, height = 12)

#Allen variations####
#Clean up old stuff
rm(list = ls()[-c(grep("sqs",ls()),grep("moll_raw",ls()),
                  grep("fossil_results",ls()),grep("fossil_data",ls()),
                  grep("fossil_results",ls()),grep("fossil_data",ls()))])

#Get modern temperature data
temp <- read.csv("woa13_decav_t00an01v2.csv", skip = 1, fill = TRUE, row.names = NULL, header = TRUE, check.names=FALSE)
colnames(temp)[1:3] <- c("Latitude","Longitude","SST_C")
temp$SST_K <- temp$SST_C + 273.15

## gam regression
set.seed(17)
temp_reg <- gam(SST_K ~ s(Latitude), data = temp)

#get modern Gastropod data
specimens <- read.table("0135238-200613084148143.csv", sep = "\t", header = TRUE, fill = TRUE, quote = NULL, stringsAsFactors = FALSE)
specimens <- subset(specimens, !is.na(decimalLatitude) & class == "Gastropoda")

#Figure out which occurrences would count for Roy et al
#https://www.pnas.org/content/pnas/95/7/3699.full.pdf
#Find the distance between the occurrences and the coast (for things within the Roy range)
library(rworldmap)
library(geosphere)
roypoints <- unique(subset(specimens, decimalLongitude < -25 & decimalLatitude > -5)[,c("decimalLongitude", "decimalLatitude")])
sPDF <- subset(getMap(), continent %in% c("North America", "South America"))
dist2coast <- dist2Line(p = roypoints, line = sPDF)
roypoints$dist2coast <- dist2coast[,"distance"]
#We only want occurrences on the shelf
#This paper says the longest shelf on the East Coast is 420km:
#https://pubs.usgs.gov/pp/0529a/report.pdf
specimens$royrange <- (specimens$decimalLongitude < -25 & specimens$decimalLatitude > -5 &
                         paste(specimens$decimalLongitude, specimens$decimalLatitude) %in% paste(subset(roypoints, dist2coast <= 420000)$decimalLongitude, subset(roypoints, dist2coast <= 420000)$decimalLatitude))

specimens_loc <- specimens[,c("genus", "species","decimalLatitude","decimalLongitude","royrange")]

#Set up band size variations
band_types <- c("equal-width", "equal-area")

#Number of total latitudinal bins
lat_bin_nums <- 180/c(1, 2, 4, 7.5, 15)
#lat_bin_nums <- 180/c(.5, .75, 1, 2, 5, 7.5, 10, 15)

#SQS options
SQS_trials <- 100
SQS_levels <- c(.25, .5, .75)

tax_levels <- c("species","genus")

allen_results <- data.frame(matrix(NA, nrow = length(band_types) * length(lat_bin_nums) *
                               length(tax_levels) * (2 + length(SQS_levels)) * 2, ncol = 11))
colnames(allen_results) <- c("band_type", "num_bands", "tax_level", "div_metric", "royrange",
                       "slope", "slope_se", "intercept", "intercept_se", "r2", "p")

k <- 1
pb <- txtProgressBar(max = nrow(allen_results), style = 3)
for(royrange in c(FALSE, TRUE)){
  for(band_type in band_types){
    for(i in 1:length(lat_bin_nums)){
      #Set up latitude bands
      if(band_type == "equal-width"){
        #Equal-width latitudinal bins
        r <- raster(nrows=lat_bin_nums[i], ncols=180, xmn=-180, xmx=180, ymn=-90, ymx=90)
        area <- rowSums(as.matrix(raster::area(r)))
        
        lat_bounds <- seq(-90 , 90, length.out = lat_bin_nums[i] + 1)
      } else if(band_type == "equal-area"){
        #Equal-area latitudinal bins
        lat_bounds <- EqualAreaRectangularGrid(N_longitude = 1, N_latitude = lat_bin_nums[i])$latitude_breaks
        
        #set up dataframe for latitudinal temperatures
        area <- rep((4*pi*6371^2)/lat_bin_nums[i], lat_bin_nums[i])
      }
      lat_bands <- data.frame(lat_min = lat_bounds[1:(length(lat_bounds) - 1)], lat_max = lat_bounds[2:length(lat_bounds)])
      lat_bands$lat_mid <- (lat_bands$lat_min + lat_bands$lat_max)/2
      lat_bands$lat_bin <- levels(cut(-90:90, lat_bounds, include.lowest = TRUE))
      lat_bands$area <- area
      
      #Get temperature data for bands
      pred <- predict(temp_reg, newdata = data.frame(Latitude = lat_bands$lat_mid), se.fit = TRUE)
      cells <- lat_bands
      cells$SST_K <- pred[["fit"]]
      cells$SST_K_se <- pred[["se.fit"]]
      cells$inv_K <- 1000/cells$SST_K
      
      #Perform analyses for each taxonomic level
      for(tax_level in tax_levels){
        #Bin gastropod data to bands
        dat <- specimens_loc[specimens_loc[,tax_level] != "",]
        if(royrange){
          dat <- subset(dat, royrange)
        }
        
        tot_div <- length(unique(dat[,tax_level]))
        
        dat$lat_bin <- cut(dat$decimalLatitude, lat_bounds, include.lowest = TRUE)
        #Calculate raw diversity
        raw_div <- table(as.character(unique(dat[,c(tax_level,"lat_bin")])$lat_bin))
        cells[[paste0(tax_level,"_raw")]] <- as.numeric(NA)
        cells[[paste0(tax_level,"_raw")]][match(names(raw_div), cells$lat_bin)] <- raw_div
        
        #Calculate regression stuff and save it for raw div
        reg <- lm(log(div) ~ inv_K, data = dplyr::select(cells, inv_K, div = paste0(tax_level,"_raw")))
        summ <- summary(reg)
        allen_results$band_type[k] <- band_type
        allen_results$num_bands[k] <- lat_bin_nums[i]
        allen_results$tax_level[k] <- tax_level
        allen_results$div_metric[k] <- "raw"
        allen_results$royrange[k] <- royrange
        allen_results$slope[k] <- summ$coefficients["inv_K","Estimate"]
        allen_results$slope_se[k] <- summ$coefficients["inv_K","Std. Error"]
        allen_results$intercept[k] <- summ$coefficients["(Intercept)","Estimate"]
        allen_results$intercept_se[k] <- summ$coefficients["(Intercept)","Std. Error"]
        allen_results$r2[k] <- summ$r.squared
        allen_results$p[k] <- summ$coefficients["inv_K","Pr(>|t|)"]
        k <- k+ 1
        setTxtProgressBar(pb, k)
        
        #Do that again but calculating max and min lat for each species, then assuming they exist in all bins between those
        ranges <- dat %>%
          group_by_at(tax_level) %>%
          summarise(
            min_lat = min(decimalLatitude),
            max_lat = max(decimalLatitude),
            .groups = 
          )
        
        cells[[paste0(tax_level,"_range")]] <- as.numeric(NA)
        for(j in 1:length(cells$lat_bin)){
          div <- nrow(subset(ranges, max_lat >= cells$lat_min[j] & min_lat <= cells$lat_max[j]))
          if(div > 0){
            cells[[paste0(tax_level,"_range")]][j] <- div
          }
        }
        
        #Calculate regression stuff and save it for range div
        reg <- lm(log(div) ~ inv_K, data = dplyr::select(cells, inv_K, div = paste0(tax_level,"_range")))
        summ <- summary(reg)
        allen_results$band_type[k] <- band_type
        allen_results$num_bands[k] <- lat_bin_nums[i]
        allen_results$tax_level[k] <- tax_level
        allen_results$div_metric[k] <- "range"
        allen_results$royrange[k] <- royrange
        allen_results$slope[k] <- summ$coefficients["inv_K","Estimate"]
        allen_results$slope_se[k] <- summ$coefficients["inv_K","Std. Error"]
        allen_results$intercept[k] <- summ$coefficients["(Intercept)","Estimate"]
        allen_results$intercept_se[k] <- summ$coefficients["(Intercept)","Std. Error"]
        allen_results$r2[k] <- summ$r.squared
        allen_results$p[k] <- summ$coefficients["inv_K","Pr(>|t|)"]
        k <- k+ 1
        setTxtProgressBar(pb, k)
        
        #SQS diversity
        for(SQS_level in SQS_levels){
          tot_div <- sqs(as.numeric(table(as.character(dat[,tax_level]))), q = .5, trials = SQS_trials)
          tot_div["subsampled richness"]
          
          cells[[paste0(tax_level,"_sqs_",SQS_level)]] <- as.numeric(NA)
          for(j in 1:length(cells$lat_bin)){
            sub_occs <- subset(dat, lat_bin == cells$lat_bin[j])
            if(nrow(sub_occs) > 0){
              tab <- as.numeric(table(as.character(sub_occs[[tax_level]])))
              capture.output(div <- sqs(tab, q = SQS_level, trials = SQS_trials))
              cells[[paste0(tax_level,"_sqs_",SQS_level)]][j] <- as.numeric(div["subsampled richness"])
            }
          }
          #Calculate regression stuff and save it for SQS div
          reg <- lm(log(div) ~ inv_K, data = dplyr::select(cells, inv_K, div = paste0(tax_level,"_sqs_",SQS_level)))
          summ <- summary(reg)
          allen_results$band_type[k] <- band_type
          allen_results$num_bands[k] <- lat_bin_nums[i]
          allen_results$tax_level[k] <- tax_level
          allen_results$div_metric[k] <- paste0("sqs_",SQS_level)
          allen_results$royrange[k] <- royrange
          allen_results$slope[k] <- summ$coefficients["inv_K","Estimate"]
          allen_results$slope_se[k] <- summ$coefficients["inv_K","Std. Error"]
          allen_results$intercept[k] <- summ$coefficients["(Intercept)","Estimate"]
          allen_results$intercept_se[k] <- summ$coefficients["(Intercept)","Std. Error"]
          allen_results$r2[k] <- summ$r.squared
          allen_results$p[k] <- summ$coefficients["inv_K","Pr(>|t|)"]
          k <- k+ 1
          setTxtProgressBar(pb, k)
        }
      }
    }
  }
}
close(pb)

#Figure S13####
library(ggplot2)
library(viridis)

#Plot of slopes
ann_text <- data.frame(royrange = TRUE, num_bands = 180, tax_level = "species",
                       div_metric = "range", band_type = "equal-width", slope = -16,
                       lab = "Roy et al 1998")

ggplot(allen_results, aes(x = as.character(num_bands), y = slope, shape = band_type)) +
  scale_x_discrete(limits = as.character(sort(lat_bin_nums)), name = "Number of Latitudinal Bands") +
  scale_y_continuous(name = "Estimated Slope", sec.axis = sec_axis(~., labels = NULL)) +
  geom_rect(xmin = 0, xmax = 50, ymin = -7.65, ymax = -6.71, color = "grey80", fill = "grey80") +
  geom_segment(x = 0, xend = 50, y = -7.65, yend = -7.65, aes(linetype = "roy"), color = "grey80") +
  geom_rect(xmin = 0, xmax = 50, ymin = -8.12, ymax = -6.96, linetype = "dashed", color = "black", fill = NA, size = .5) +
  geom_segment(x = 0, xend = 50, y = -6.96, yend = -6.96, aes(linetype = "theory"), color = NA, size = .5) +
  geom_point(aes(color = div_metric), size = 2, position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(color = div_metric, ymin = slope - 1.96*slope_se, ymax = slope + 1.96*slope_se), width = .8, position = position_dodge(width = 0.8)) +
  geom_text(data = ann_text, x = 4.725, aes(label = lab), show.legend = FALSE,
            color = "black", hjust = .5, size = 5) +
  geom_segment(data = ann_text, x = 4.725, xend = 4.725, y = -15, yend = -11, show.legend = FALSE,
            color = "black", arrow = arrow(length = unit(0.03, "npc")), size = .75) +
  theme_classic(base_size = 24, base_family = "") +
  theme(axis.ticks = element_line(color = "black"), axis.ticks.y.right = element_blank(),
        axis.line = element_line(linetype = "solid", colour = "black", lineend = "square"), axis.line.y = element_line(lineend = "butt"),
        axis.text = element_text(colour = "black"),
        plot.margin = unit(c(1,1,1,1), "lines"), panel.border = element_rect(fill = NA),
        legend.key.size = unit(1.75, "lines"), legend.spacing.y = unit(.25, "lines")) +
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(margin = unit(c(.5,.5,.5,.5), "lines"))) +
  scale_color_viridis(discrete=TRUE, name = "Diversity Metric", labels = c("range" = "Rangethrough", "raw" = "Raw Occurrences", "sqs_0.25" = "SQS (q = .25)", "sqs_0.5" = "SQS (q = .5)", "sqs_0.75" = "SQS (q = .75)"), guide = guide_legend(order = 1)) +
  scale_shape(name = "Latitudinal Band Type", labels = c("equal-area" = "Equal Area", "equal-width" = "Equal Width"), guide = guide_legend(order = 3)) +
  scale_linetype_manual(name = NULL, values = c("roy" = "solid", "theory" = "dashed"),
                        labels = c("roy" = "Roy et al 1998\nempirical estimate", "theory" = "Gillooly and Allen 2007\ntheoretical range"),
                        guide = guide_legend(keyheight=unit(2, "cm"), order = 3, override.aes = list(size = c(4,.5), color = c("grey80", "black")))) +
  guides(color = guide_legend(order = 1), shape = guide_legend(order = 2)) +
  facet_grid(tax_level~royrange, labeller = labeller(royrange = as_labeller(c(`TRUE` = 'Roy~et~al~1998~Range', `FALSE` = 'Global'), label_parsed)))
ggsave("Allen Variations Slope.pdf", device = "pdf", width = 16, height = 10.5)

#Plot of intercepts
ann_text <- data.frame(royrange = TRUE, num_bands = 180, tax_level = "species",
                       div_metric = "range", band_type = "equal-width", intercept = 57,
                       lab = "Roy et al 1998")

ggplot(allen_results, aes(x = as.character(num_bands), y = intercept, color = div_metric, shape = band_type)) +
  scale_x_discrete(limits = as.character(sort(lat_bin_nums)), name = "Number of Latitudinal Bands") +
  scale_y_continuous(name = "Estimated Intercept", sec.axis = sec_axis(~., labels = NULL)) +
  geom_segment(x = 0, xend = 50, y = 30.51, yend = 30.51, aes(linetype = "roy"), color = "black", size = .5, show.legend = TRUE) +
  geom_point(size = 2, position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = intercept - 1.96*intercept_se, ymax = intercept + 1.96*intercept_se), width = .8, position = position_dodge(width = 0.8)) +
  geom_text(data = ann_text, x = 4.725, aes(label = lab), show.legend = FALSE,
            color = "black", hjust = .5, size = 5) +
  geom_segment(data = ann_text, x = 4.725, xend = 4.725, y = 54, yend = 44, show.legend = FALSE,
               color = "black", arrow = arrow(length = unit(0.03, "npc")), size = .75) +
  theme_classic(base_size = 24, base_family = "") +
  theme(axis.ticks = element_line(color = "black"), axis.ticks.y.right = element_blank(),
        axis.line = element_line(linetype = "solid", colour = "black", lineend = "square"), axis.line.y = element_line(lineend = "butt"),
        axis.text = element_text(colour = "black"),
        plot.margin = unit(c(1,1,1,1), "lines"), panel.border = element_rect(fill = NA),
        legend.key.size = unit(1.75, "lines"), legend.spacing.y = unit(.25, "lines")) +
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(margin = unit(c(.5,.5,.5,.5), "lines"))) +
  scale_color_viridis(discrete=TRUE, name = "Diversity Metric", labels = c("range" = "Rangethrough", "raw" = "Raw Occurrences", "sqs_0.25" = "SQS (q = .25)", "sqs_0.5" = "SQS (q = .5)", "sqs_0.75" = "SQS (q = .75)")) +
  scale_shape(name = "Latitudinal Band Type", labels = c("equal-area" = "Equal Area", "equal-width" = "Equal Width")) +
  scale_linetype_manual(name = NULL, values = c("roy" = "dotted"),
                        labels = c("roy" = "Roy et al 1998\nintercept")) +
  facet_grid(tax_level~royrange, labeller = labeller(royrange = as_labeller(c(`TRUE` = 'Roy~et~al~1998~Range', `FALSE` = 'Global'), label_parsed)))
ggsave("Allen Variations Intercept.pdf", device = "pdf", width = 16, height = 10.5)

#Plot of r-squared
ann_text <- data.frame(royrange = TRUE, num_bands = 180, tax_level = "species",
                       div_metric = "range", band_type = "equal-width", r2 = .8,
                       lab = "Roy et al 1998")

ggplot(allen_results, aes(x = as.character(num_bands), y = r2, color = div_metric, shape = band_type, group = div_metric)) +
  scale_x_discrete(limits = as.character(sort(lat_bin_nums)), name = "Number of Latitudinal Bands") +
  scale_y_continuous(name = expression("R"^2), limits = c(0, 1), sec.axis = sec_axis(~., labels = NULL)) +
  geom_point(size = 2, position = position_dodge(width = 0.8)) +
  geom_text(data = ann_text, aes(label = lab), show.legend = FALSE,
            color = "black", hjust = .65, size = 5) +
  geom_segment(data = ann_text, x = 5.675, xend = 5.675, y = .75, yend = .525, show.legend = FALSE,
               color = "black", arrow = arrow(length = unit(0.03, "npc")), size = .75) +
  theme_classic(base_size = 24, base_family = "") +
  theme(axis.ticks = element_line(color = "black"), axis.ticks.y.right = element_blank(),
        axis.line = element_line(linetype = "solid", colour = "black", lineend = "square"), axis.line.y = element_line(lineend = "butt"),
        axis.text = element_text(colour = "black"),
        plot.margin = unit(c(1,1,1,1), "lines"), panel.border = element_rect(fill = NA),
        legend.key.size = unit(1.75, "lines"), legend.spacing.y = unit(.25, "lines")) +
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(margin = unit(c(.5,.5,.5,.5), "lines"))) +
  scale_color_viridis(discrete=TRUE, name = "Diversity Metric", labels = c("range" = "Rangethrough", "raw" = "Raw Occurrences", "sqs_0.25" = "SQS (q = .25)", "sqs_0.5" = "SQS (q = .5)", "sqs_0.75" = "SQS (q = .75)")) +
  scale_shape(name = "Latitudinal Band Type", labels = c("equal-area" = "Equal Area", "equal-width" = "Equal Width")) +
  facet_grid(tax_level~royrange, labeller = labeller(royrange = as_labeller(c(`TRUE` = 'Roy~et~al~1998~Range', `FALSE` = 'Global'), label_parsed)))
ggsave("Allen Variations R-Squared.pdf", device = "pdf", width = 16, height = 10.5)

#Plot of p-value
ann_text <- data.frame(royrange = TRUE, num_bands = 180, tax_level = "species",
                       div_metric = "range", band_type = "equal-width", p = .4,
                       lab = "Roy et al 1998")

ggplot(allen_results, aes(x = as.character(num_bands), y = p, color = div_metric, shape = band_type, group = div_metric)) +
  scale_x_discrete(limits = as.character(sort(lat_bin_nums)), name = "Number of Latitudinal Bands") +
  scale_y_continuous(name = "p-value", limits = c(0, 1), sec.axis = sec_axis(~., labels = NULL)) +
  geom_hline(yintercept = 0.05, linetype = "dotted") +
  geom_point(size = 2, position = position_dodge(width = 0.8)) +
  geom_text(data = ann_text, aes(label = lab), show.legend = FALSE,
            color = "black", hjust = .65, size = 5) +
  geom_segment(data = ann_text, x = 5.675, xend = 5.675, y = .35, yend = .02, show.legend = FALSE,
               color = "black", arrow = arrow(length = unit(0.03, "npc")), size = .75) +
  theme_classic(base_size = 24, base_family = "") +
  theme(axis.ticks = element_line(color = "black"), axis.ticks.y.right = element_blank(),
        axis.line = element_line(linetype = "solid", colour = "black", lineend = "square"), axis.line.y = element_line(lineend = "butt"),
        axis.text = element_text(colour = "black"),
        plot.margin = unit(c(1,1,1,1), "lines"), panel.border = element_rect(fill = NA),
        legend.key.size = unit(1.75, "lines"), legend.spacing.y = unit(.25, "lines")) +
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(margin = unit(c(.5,.5,.5,.5), "lines"))) +
  scale_color_viridis(discrete=TRUE, name = "Diversity Metric", labels = c("range" = "Rangethrough", "raw" = "Raw Occurrences", "sqs_0.25" = "SQS (q = .25)", "sqs_0.5" = "SQS (q = .5)", "sqs_0.75" = "SQS (q = .75)")) +
  scale_shape(name = "Latitudinal Band Type", labels = c("equal-area" = "Equal Area", "equal-width" = "Equal Width")) +
  facet_grid(tax_level~royrange, labeller = labeller(royrange = as_labeller(c(`TRUE` = 'Roy~et~al~1998~Range', `FALSE` = 'Global'), label_parsed)))
ggsave("Allen Variations P-Value.pdf", device = "pdf", width = 16, height = 10.5)
