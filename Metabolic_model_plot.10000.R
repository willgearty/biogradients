library(ggplot2)
library(viridis)
library(dplyr)
library(zoo)
library(egg)
library(segmented)

# Set working directory as wherever you have metabolic model summary saved
load("Metabolic.model.summary.10000.RData")

# Initiate summary vectors to define plotted envelopes
SD_0.0_rich <- as.numeric()
SD_0.5_rich <- as.numeric()
SD_1.0_rich <- as.numeric()
SD_0.0_med_rich <- as.numeric()

# Generate summary vectors by looping through temperatures 
for (temp in seq(3.685, 3.202, -0.001)){
  
  SD_0.0_rich_temp <- quantile(x = filter(model.summary, Temp.allen <= temp & Temp.allen > temp - 0.001 & Eo.sd < 0.01)$log.richness,
                            probs=c(0.0228, 0.9772, 0.5))
  
  SD_0.0_rich <- as.data.frame(rbind(SD_0.0_rich, cbind(temp, SD_0.0_rich_temp[1], SD_0.0_rich_temp[2])))
  
  SD_0.0_med_rich <- as.data.frame(rbind(SD_0.0_med_rich, cbind(temp, SD_0.0_rich_temp[3])))
  
  SD_0.5_rich_temp <- quantile(x = filter(model.summary, Temp.allen <= temp & Temp.allen > temp - 0.001 & Eo.sd < 0.51 & Eo.sd > 0.49)$log.richness,
                               probs=c(0.0228, 0.9772))
  
  SD_0.5_rich <- as.data.frame(rbind(SD_0.5_rich, cbind(temp, SD_0.5_rich_temp[1], SD_0.5_rich_temp[2])))
  
  SD_1.0_rich_temp <- quantile(x = filter(model.summary, Temp.allen <= temp & Temp.allen > temp - 0.001 & Eo.sd < 1.01 & Eo.sd > 0.99)$log.richness,
                               probs=c(0.0228, 0.9772))
  
  SD_1.0_rich <- as.data.frame(rbind(SD_1.0_rich, cbind(temp, SD_1.0_rich_temp[1], SD_1.0_rich_temp[2])))
  
}

# rename summary vector headers
names(SD_0.0_rich) <- c("Temp.allen", "log_rich_0.05", "log_rich_0.95")

SD_0.0_rich$SD <- rep("0", nrow(SD_0.0_rich))
  
names(SD_0.5_rich) <- c("Temp.allen", "log_rich_0.05", "log_rich_0.95")

SD_0.5_rich$SD <- rep("0.5", nrow(SD_0.5_rich))

names(SD_1.0_rich) <- c("Temp.allen", "log_rich_0.05", "log_rich_0.95")

SD_1.0_rich$SD <- rep("1.0", nrow(SD_1.0_rich))

# combine summary vectors
SD_rich <- rbind(SD_1.0_rich,
                 SD_0.5_rich, 
                 SD_0.0_rich
)

SD_rich$SD <- as.factor(SD_rich$SD)

# relevel factors for plotting
SD_rich$SD <- relevel(SD_rich$SD, "0.5")
SD_rich$SD <- relevel(SD_rich$SD, "1.0")

# rename columns of vector for median model values (used for c.opt analysis below)
names(SD_0.0_med_rich) <- c("Temp.allen", "log_rich_med")

# plot metabolic model overlain by segmented linear regression model from diversity_temp (Figure 2)
# load fossil data (change working directory if necessary)
load("fossil_results.RData")

# subset fossil data
fossil_data_sub <- subset(fossil_data, num_bands == 24 & metric == "SQS_0.25" & band_type == "equal-area")

# segmented linear regression
reg <- lm(log(div.prop) ~ SST_K, data = fossil_data_sub)
reg_seg <- with(fossil_data_sub, segmented(reg, seg.Z = ~ SST_K, npsi = 1))
seg_pred <- as.data.frame(cbind(SST_K = 1000/seq(3.685, 3.202, -0.001), predict(reg_seg, newdata = data.frame(SST_K = 1000/seq(3.685, 3.202, -0.001)), interval = "confidence", level = .95)))

# name columns of regression dataframe
names(seg_pred) <- c("SST_K", "log_div", "lwr_log_div", "upr_log_div")

# test plot to see model fit without intercept adjustment
ggplot()+geom_point(aes(x=seg_pred$log_div, y=SD_0.0_med_rich$log_rich_med)) +
  theme_classic(base_size = 24) +
  theme(axis.ticks = element_line(color = "black"), axis.text = element_text(colour = "black"),
        axis.line = element_line(color = "black", lineend = "square"), 
        axis.ticks.x.top = element_blank(), axis.text.x.top = element_blank(),
        axis.ticks.y.right = element_blank(), axis.text.y.right = element_blank(),
        plot.margin = margin(0,1,1,1, "lines"), panel.border = element_rect(fill = NA),
        legend.key.size = unit(.5, "lines"), legend.spacing.y = unit(0.75, "lines"),
        legend.justification=c(1,1), legend.position=c(.45,1.0),
        legend.background = element_rect(fill=NA),    
        legend.title = element_text(color = "black", size = 18),
        legend.text = element_text(color = "black", size = 16))

# establish 'best' intercept for metabolic model - i.e. the one that allows best correlation with fossil analyses based upon RSS
intercept.test.sum <- as.numeric()
for(c in seq(-3,3,0.001)){
  res <- sum((seg_pred$log_div - SD_0.0_med_rich$log_rich_med + c)^2, na.rm=T)
  intercept.test.sum <- rbind(intercept.test.sum, cbind(c, res))
}

intercept.test.sum <- as.data.frame(intercept.test.sum)
# name 'best' intercept c.opt
c.opt <- intercept.test.sum$c[which.min(intercept.test.sum$res)]

# Plot model and fossil data - this can be independently saved if not interested in replotting climate data
model <- ggplot()+
  geom_ribbon(data=SD_rich, aes(x=(1000/Temp.allen)-273.15, ymin=exp(log_rich_0.05-c.opt), ymax=exp(log_rich_0.95-c.opt), alpha=SD), fill=viridis_pal()(4)[1])+
  scale_alpha_discrete(range=c(0.3, 0.7), name= "Temperature sensitivity\n(SD from mean)")+
  geom_ribbon(data = seg_pred, aes(x = SST_K-273.15, ymin = exp(lwr_log_div), ymax = exp(upr_log_div)), alpha = .25, fill = viridis_pal()(4)[3], alpha=.3) +
  geom_line(data = seg_pred, aes(x = SST_K-273.15, y = exp(log_div)), size = 1, color = viridis_pal()(4)[3]) +
  scale_x_continuous(name = expression("Sea Surface Temperature ("*degree*"C)"), expand = c(0,0),
                     breaks = c(0, 10, 20, 30), labels = c(0, 10, 20, 30), sec.axis = sec_axis(~.)) +
  scale_y_continuous(name = "Percent Total Generic Diversity", expand = c(0,0), sec.axis = sec_axis(~.), trans = "log",
                     breaks = c(0.01, 0.03, 0.1, 0.3, 1), labels = c(0.01, 0.03, 0.1, 0.3, 1) * 100) +
  coord_cartesian(xlim = c(-3, 38), ylim = c(0.0075, 1)) +
  theme_classic(base_size = 24) +
  theme(axis.ticks = element_line(color = "black"), axis.text = element_text(colour = "black"),
        axis.line = element_line(color = "black", lineend = "square"), 
        axis.ticks.x.top = element_blank(), axis.text.x.top = element_blank(),
        axis.ticks.y.right = element_blank(), axis.text.y.right = element_blank(),
        plot.margin = margin(0,1,1,1, "lines"), panel.border = element_rect(fill = NA),
        legend.key.size = unit(.5, "lines"), legend.spacing.y = unit(0.75, "lines"),
        legend.justification=c(1,1), legend.position=c(.45,1.0),
        legend.background = element_rect(fill=NA),    
        legend.title = element_text(color = "black", size = 18),
        legend.text = element_text(color = "black", size = 16))
model

# Import climate model data
load("Equatorial.temp.sum.RData")

# make all dataframe columns numeric
temp.sum$temp.K.5 <- as.numeric(paste(temp.sum$temp.K.5))
temp.sum$temp.K.25 <- as.numeric(paste(temp.sum$temp.K.25))
temp.sum$temp.K.50 <- as.numeric(paste(temp.sum$temp.K.50))
temp.sum$temp.K.75 <- as.numeric(paste(temp.sum$temp.K.75))
temp.sum$temp.K.95 <- as.numeric(paste(temp.sum$temp.K.95))

# relevel for plotting
levels(temp.sum$scenario)[levels(temp.sum$scenario) == "RCP4.5"] <- "RCP/ECP 4.5"
levels(temp.sum$scenario)[levels(temp.sum$scenario) == "RCP8.5"] <- "RCP/ECP 8.5"
temp.sum <- filter(temp.sum, year == "2299" | year == "2100" |  year == "1860" )
  
# plot summary distributions of climate model temperatures
temp <- ggplot(temp.sum, aes(xmiddle=temp.K.50-273.15, y=year))+
  geom_boxplot(stat="identity", width=.9, color= NA, alpha =.4, position=position_dodge2(padding = 0.3, preserve="single"), aes(fill = scenario, xlower=temp.K.5-273.15, xupper=temp.K.95-273.15, xmin= temp.K.5-273.15, xmax= temp.K.95-273.15))+
  geom_boxplot(stat="identity", width=.9, color= NA, alpha =.7, position=position_dodge2(padding = 0.3, preserve="single"), aes(fill = scenario,  xlower=temp.K.25-273.15, xupper=temp.K.75-273.15, xmin= temp.K.25-273.15, xmax= temp.K.75-273.15))+
  geom_boxplot(stat="identity", width=.9, color= "white", alpha =.9, position=position_dodge2(padding = 0.3, preserve="single"), aes(fill = scenario, xlower=temp.K.50-273.15, xupper=temp.K.50-273.15, xmin= temp.K.50-273.15, xmax= temp.K.50-273.15))+
  scale_fill_manual(values=c("grey20", viridis_pal()(20)[9], "goldenrod2"), name ="Model Scenario")+
  scale_x_continuous(name = expression("Sea Surface Temperature ("*degree*"C)"), expand = c(0,0),
                     breaks = c(0, 10, 20, 30), labels = c(0, 10, 20, 30), sec.axis = sec_axis(~.)) +
  coord_cartesian(xlim = c(-3, 38)) +
  theme_classic(base_size = 24) + 
  ylab("Year")+
  geom_text(aes(y=year, label = year, x=21.0), size=5)+
  theme(axis.ticks = element_line(color = "black"), axis.text = element_text(colour = "black"),
      axis.line = element_line(color = NA, lineend = "square"), 
      axis.ticks.x = element_blank(), axis.text.x = element_text(size=0, color=NA),
      axis.ticks.y = element_blank(), 
      axis.text.y = element_blank(),
      axis.title.y = element_text(vjust=-51), 
      axis.title.x = element_text(size=1, color=NA),
      panel.background = element_rect(fill=NA),
      plot.margin = margin(1,1,0,1, "lines"), panel.border = element_rect(fill = NA, color=NA),
      legend.key.size = unit(1, "lines"), legend.spacing.y = unit(.5, "lines"),
      legend.justification=c(1,1), legend.position=c(.35, 1.1), 
      legend.background = element_rect(fill=NA), 
      legend.title = element_text(color = "black", size = 18),
      legend.text = element_text(color = "black", size = 16))

temp

# combine 'model' and 'temp' to generate final figure 3
full.figure <-  ggarrange(temp, model, ncol = 1, heights =c(0.2, 1))

# save as pdf and png
ggsave("Metabolic model 10000 with temp.pdf", full.figure, width = 8, height = 9)
ggsave("Metabolic model 10000 with temp.png", full.figure, width = 8, height = 9)

