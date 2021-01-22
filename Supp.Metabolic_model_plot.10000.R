library(ggplot2)
library(viridis)
library(dplyr)
library(zoo)
library(deeptime)
library(ggpubr)
library(segmented)

# Set working directory as wherever you have metabolic model summary saved
load("Metabolic.model.summary.10000.supp.RData")

# Initiate summary vectors to define plotted envelopes for each taxon
Crassostrea.gigas_rich <- as.numeric()
Gibberulus.gibbosus_rich <- as.numeric()
Mytilus.species_rich <- as.numeric()
Pecten.maximus_rich <- as.numeric()
Dosidicus.gigas_rich <- as.numeric()
Nautilus.pompillius_rich <- as.numeric()
Sepia.officinalis_rich <- as.numeric()

# List Eo values for all molluscs in Penn et al. (2018)
Crassostrea.gigas <- 0.21 # pacific oyster
Gibberulus.gibbosus <- 0.41 # humpbacked conch (sea snail)
Mytilus.species <- 0.16 # bivalve
Pecten.maximus <- 0.33 # great scallops
Dosidicus.gigas <- 0.62 # Humboldt squid
Nautilus.pompillius <- 0.59 # chambered nautilus
Sepia.officinalis <- 0.72 # common cuttlefish

# Generate summary vectors by looping through temperatures 
for (temp in seq(3.685, 3.195, -0.001)){
  
  Crassostrea.gigas_rich_temp <- quantile(x = filter(model.summary, Temp.allen <= temp & Temp.allen > temp - 0.001 & Eo.val == Crassostrea.gigas)$log.richness,
                            probs=c(0.0228, 0.9772))
  
  Crassostrea.gigas_rich <- as.data.frame(rbind(Crassostrea.gigas_rich, cbind(temp, Crassostrea.gigas_rich_temp[1], Crassostrea.gigas_rich_temp[2])))
  
  Gibberulus.gibbosus_rich_temp <- quantile(x = filter(model.summary, Temp.allen <= temp & Temp.allen > temp - 0.001 & Eo.val == Gibberulus.gibbosus)$log.richness,
                                          probs=c(0.0228, 0.9772))
  
  Gibberulus.gibbosus_rich <- as.data.frame(rbind(Gibberulus.gibbosus_rich, cbind(temp, Gibberulus.gibbosus_rich_temp[1], Gibberulus.gibbosus_rich_temp[2])))
  
  Mytilus.species_rich_temp <- quantile(x = filter(model.summary, Temp.allen <= temp & Temp.allen > temp - 0.001 & Eo.val == Mytilus.species)$log.richness,
                                          probs=c(0.0228, 0.9772))
  
  Mytilus.species_rich <- as.data.frame(rbind(Mytilus.species_rich, cbind(temp, Mytilus.species_rich_temp[1], Mytilus.species_rich_temp[2])))
  
  Pecten.maximus_rich_temp <- quantile(x = filter(model.summary, Temp.allen <= temp & Temp.allen > temp - 0.001 & Eo.val == Pecten.maximus)$log.richness,
                                        probs=c(0.0228, 0.9772))
  
  Pecten.maximus_rich <- as.data.frame(rbind(Pecten.maximus_rich, cbind(temp, Pecten.maximus_rich_temp[1], Pecten.maximus_rich_temp[2])))
  
  Dosidicus.gigas_rich_temp <- quantile(x = filter(model.summary, Temp.allen <= temp & Temp.allen > temp - 0.001 & Eo.val == Dosidicus.gigas)$log.richness,
                                            probs=c(0.0228, 0.9772))
  
  Dosidicus.gigas_rich <- as.data.frame(rbind(Dosidicus.gigas_rich, cbind(temp, Dosidicus.gigas_rich_temp[1], Dosidicus.gigas_rich_temp[2])))
  
  Nautilus.pompillius_rich_temp <- quantile(x = filter(model.summary, Temp.allen <= temp & Temp.allen > temp - 0.001 & Eo.val == Nautilus.pompillius)$log.richness,
                                        probs=c(0.0228, 0.9772))
  
  Nautilus.pompillius_rich <- as.data.frame(rbind(Nautilus.pompillius_rich, cbind(temp, Nautilus.pompillius_rich_temp[1], Nautilus.pompillius_rich_temp[2])))
  
  Sepia.officinalis_rich_temp <- quantile(x = filter(model.summary, Temp.allen <= temp & Temp.allen > temp - 0.001 & Eo.val == Sepia.officinalis)$log.richness,
                                       probs=c(0.0228, 0.9772))
  
  Sepia.officinalis_rich <- as.data.frame(rbind(Sepia.officinalis_rich, cbind(temp, Sepia.officinalis_rich_temp[1], Sepia.officinalis_rich_temp[2])))
  
  
}

# rename summary vector headers
names(Crassostrea.gigas_rich) <- c("Temp.allen", "log_rich_0.05", "log_rich_0.95")
names(Gibberulus.gibbosus_rich) <- c("Temp.allen", "log_rich_0.05", "log_rich_0.95")
names(Mytilus.species_rich) <- c("Temp.allen", "log_rich_0.05", "log_rich_0.95")
names(Pecten.maximus_rich) <- c("Temp.allen", "log_rich_0.05", "log_rich_0.95")
names(Dosidicus.gigas_rich) <- c("Temp.allen", "log_rich_0.05", "log_rich_0.95")
names(Nautilus.pompillius_rich) <- c("Temp.allen", "log_rich_0.05", "log_rich_0.95")
names(Sepia.officinalis_rich) <- c("Temp.allen", "log_rich_0.05", "log_rich_0.95")

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

# set c.opt based upon main text analyses -> if editing the analyses in any way refer to c.opt value established in Figure 3 plotting script
c.opt <- -0.62

# Generate individual plot for each molluscan taxon
Crassostrea.gigas_plot <- ggplot(Crassostrea.gigas_rich, aes(x=(1000/Temp.allen)-273.15, ymin=exp(log_rich_0.05-c.opt), ymax=exp(log_rich_0.95-c.opt)))+
  geom_ribbon(alpha=.8, fill=viridis_pal()(4)[1])+
  geom_ribbon(data = seg_pred, aes(x = SST_K-273.15, ymin = exp(lwr_log_div), ymax = exp(upr_log_div)), alpha = .25, fill = viridis_pal()(4)[3], alpha=.3) +
  annotate(geom="line",  x = seg_pred$SST_K-273.15, y = exp(seg_pred$log_div), size = 1, color = viridis_pal()(4)[3]) +
  scale_x_continuous(name = expression("Sea Surface Temperature ("*degree*"C)"), expand = c(0,0),
                     breaks = c(0, 10, 20, 30), labels = c(0, 10, 20, 30), sec.axis = sec_axis(~.)) +
  scale_y_continuous(name = "Percent Total Generic Diversity", expand = c(0,0), sec.axis = sec_axis(~.), trans = "log",
                     breaks = c(0.01, 0.03, 0.1, 0.3, 1), labels = c(0.01, 0.03, 0.1, 0.3, 1) * 100) +
  coord_cartesian(xlim = c(-3, 38), ylim = c(0.0075, 1)) +
  theme_classic(base_size = 24) +
  theme(axis.ticks = element_line(color = "black"), axis.text = element_blank(),
        axis.title = element_blank(),
      axis.line = element_line(color = "black", lineend = "square"), 
      axis.ticks.x.top = element_blank(), axis.text.x.top = element_blank(),
      axis.ticks.y.right = element_blank(), axis.text.y.right = element_blank(),
      plot.margin = margin(.2,.2,.2,.2, "lines"), panel.border = element_rect(fill = NA),
      legend.key.size = unit(.25, "lines"), legend.spacing.y = unit(0.25, "lines"),
      legend.justification=c(1,1), legend.position=c(.4,.98),
      legend.background = element_rect(fill=alpha('white', 0.6)))+
  annotate(geom="text", y=80/100, x=9.5, label="Crassostrea gigas", fontface=3, size=6)
  
Gibberulus.gibbosus_plot <- ggplot(Gibberulus.gibbosus_rich, aes(x=(1000/Temp.allen)-273.15, ymin=exp(log_rich_0.05-c.opt), ymax=exp(log_rich_0.95-c.opt)))+
  geom_ribbon(alpha=.8, fill=viridis_pal()(4)[1])+
  geom_ribbon(data = seg_pred, aes(x = SST_K-273.15, ymin = exp(lwr_log_div), ymax = exp(upr_log_div)), alpha = .25, fill = viridis_pal()(4)[3], alpha=.3) +
  annotate(geom="line",  x = seg_pred$SST_K-273.15, y = exp(seg_pred$log_div), size = 1, color = viridis_pal()(4)[3]) +
  scale_x_continuous(name = expression("Sea Surface Temperature ("*degree*"C)"), expand = c(0,0),
                     breaks = c(0, 10, 20, 30), labels = c(0, 10, 20, 30), sec.axis = sec_axis(~.)) +
  scale_y_continuous(name = "Percent Total Generic Diversity", expand = c(0,0), sec.axis = sec_axis(~.), trans = "log",
                     breaks = c(0.01, 0.03, 0.1, 0.3, 1), labels = c(0.01, 0.03, 0.1, 0.3, 1) * 100) +
  coord_cartesian(xlim = c(-3, 38), ylim = c(0.0075, 1)) +
  theme_classic(base_size = 24) +
  theme(axis.ticks = element_line(color = "black"), axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.line = element_line(color = "black", lineend = "square"), 
        axis.ticks.x.top = element_blank(), axis.text.x.top = element_blank(),
        axis.ticks.y.right = element_blank(), axis.text.y.right = element_blank(),
        plot.margin = margin(.2,.2,.2,.2, "lines"), panel.border = element_rect(fill = NA),
        legend.key.size = unit(.25, "lines"), legend.spacing.y = unit(0.25, "lines"),
        legend.justification=c(1,1), legend.position=c(.4,.98),
        legend.background = element_rect(fill=alpha('white', 0.6)))+
  annotate(geom="text", y=80/100, x=11, label="Gibberulus gibbosus", fontface=3, size=6)

Mytilus.species_plot <- ggplot(Mytilus.species_rich, aes(x=(1000/Temp.allen)-273.15, ymin=exp(log_rich_0.05-c.opt), ymax=exp(log_rich_0.95-c.opt)))+
  geom_ribbon(alpha=.8, fill=viridis_pal()(4)[1])+
  geom_ribbon(data = seg_pred, aes(x = SST_K-273.15, ymin = exp(lwr_log_div), ymax = exp(upr_log_div)), alpha = .25, fill = viridis_pal()(4)[3], alpha=.3) +
  annotate(geom="line",  x = seg_pred$SST_K-273.15, y = exp(seg_pred$log_div), size = 1, color = viridis_pal()(4)[3]) +
  scale_x_continuous(name = expression("Sea Surface Temperature ("*degree*"C)"), expand = c(0,0),
                     breaks = c(0, 10, 20, 30), labels = c(0, 10, 20, 30), sec.axis = sec_axis(~.)) +
  scale_y_continuous(name = "Percent Total Generic Diversity", expand = c(0,0), sec.axis = sec_axis(~.), trans = "log",
                     breaks = c(0.01, 0.03, 0.1, 0.3, 1), labels = c(0.01, 0.03, 0.1, 0.3, 1) * 100) +
  coord_cartesian(xlim = c(-3, 38), ylim = c(0.0075, 1)) +
  theme_classic(base_size = 24) +
  theme(axis.ticks = element_line(color = "black"), axis.text = element_text(colour = "black"),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.line = element_line(color = "black", lineend = "square"), 
        axis.ticks.x.top = element_blank(), axis.text.x.top = element_blank(),
        axis.ticks.y.right = element_blank(), axis.text.y.right = element_blank(),
        plot.margin = margin(.2,.2,.2,.2, "lines"), panel.border = element_rect(fill = NA),
        legend.key.size = unit(.25, "lines"), legend.spacing.y = unit(0.25, "lines"),
        legend.justification=c(1,1), legend.position=c(.4,.98),
        legend.background = element_rect(fill=alpha('white', 0.6)))+
  annotate(geom="text", y=80/100, x=4.5, label="Mytilus sp.", fontface=3, size=6)

Pecten.maximus_plot <- ggplot(Pecten.maximus_rich, aes(x=(1000/Temp.allen)-273.15, ymin=exp(log_rich_0.05-c.opt), ymax=exp(log_rich_0.95-c.opt)))+
  geom_ribbon(alpha=.8, fill=viridis_pal()(4)[1])+
  geom_ribbon(data = seg_pred, aes(x = SST_K-273.15, ymin = exp(lwr_log_div), ymax = exp(upr_log_div)), alpha = .25, fill = viridis_pal()(4)[3], alpha=.3) +
  annotate(geom="line",  x = seg_pred$SST_K-273.15, y = exp(seg_pred$log_div), size = 1, color = viridis_pal()(4)[3]) +
  scale_x_continuous(name = expression("Sea Surface Temperature ("*degree*"C)"), expand = c(0,0),
                     breaks = c(0, 10, 20, 30), labels = c(0, 10, 20, 30), sec.axis = sec_axis(~.)) +
  scale_y_continuous(name = "Percent Total Generic Diversity", expand = c(0,0), sec.axis = sec_axis(~.), trans = "log",
                     breaks = c(0.01, 0.03, 0.1, 0.3, 1), labels = c(0.01, 0.03, 0.1, 0.3, 1) * 100) +
  coord_cartesian(xlim = c(-3, 38), ylim = c(0.0075, 1)) +
  theme_classic(base_size = 24) +
  theme(axis.ticks = element_line(color = "black"), axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line = element_line(color = "black", lineend = "square"), 
        axis.ticks.x.top = element_blank(), axis.text.x.top = element_blank(),
        axis.ticks.y.right = element_blank(), axis.text.y.right = element_blank(),
        plot.margin = margin(.2,.2,.2,.2, "lines"), panel.border = element_rect(fill = NA),
        legend.key.size = unit(.25, "lines"), legend.spacing.y = unit(0.25, "lines"),
        legend.justification=c(1,1), legend.position=c(.4,.98),
        legend.background = element_rect(fill=alpha('white', 0.6)))+
  annotate(geom="text", y=80/100, x=8.5, label="Pecten maximus", fontface=3, size=6)

Dosidicus.gigas_plot <- ggplot(Dosidicus.gigas_rich, aes(x=(1000/Temp.allen)-273.15, ymin=exp(log_rich_0.05-c.opt), ymax=exp(log_rich_0.95-c.opt)))+
  geom_ribbon(alpha=.8, fill=viridis_pal()(4)[1])+
  geom_ribbon(data = seg_pred, aes(x = SST_K-273.15, ymin = exp(lwr_log_div), ymax = exp(upr_log_div)), alpha = .25, fill = viridis_pal()(4)[3], alpha=.3) +
  annotate(geom="line",  x = seg_pred$SST_K-273.15, y = exp(seg_pred$log_div), size = 1, color = viridis_pal()(4)[3]) +
  scale_x_continuous(name = expression("Sea Surface Temperature ("*degree*"C)"), expand = c(0,0),
                     breaks = c(0, 10, 20, 30), labels = c(0, 10, 20, 30), sec.axis = sec_axis(~.)) +
  scale_y_continuous(name = "Percent Total Generic Diversity", expand = c(0,0), sec.axis = sec_axis(~.), trans = "log",
                     breaks = c(0.01, 0.03, 0.1, 0.3, 1), labels = c(0.01, 0.03, 0.1, 0.3, 1) * 100) +
  coord_cartesian(xlim = c(-3, 38), ylim = c(0.0075, 1)) +
  theme_classic(base_size = 24) +
  theme(axis.ticks = element_line(color = "black"), axis.text = element_text(colour = "black"),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.line = element_line(color = "black", lineend = "square"), 
        axis.ticks.x.top = element_blank(), axis.text.x.top = element_blank(),
        axis.ticks.y.right = element_blank(), axis.text.y.right = element_blank(),
        plot.margin = margin(.2,.2,.2,.2, "lines"), panel.border = element_rect(fill = NA),
        legend.key.size = unit(.25, "lines"), legend.spacing.y = unit(0.25, "lines"),
        legend.justification=c(1,1), legend.position=c(.4,.98),
        legend.background = element_rect(fill=alpha('white', 0.6)))+
  annotate(geom="text", y=80/100, x=8, label="Dosidicus gigas", fontface=3, size=6)

Nautilus.pompillius_plot <- ggplot(Nautilus.pompillius_rich, aes(x=(1000/Temp.allen)-273.15, ymin=exp(log_rich_0.05-c.opt), ymax=exp(log_rich_0.95-c.opt)))+
  geom_ribbon(alpha=.8, fill=viridis_pal()(4)[1])+
  geom_ribbon(data = seg_pred, aes(x = SST_K-273.15, ymin = exp(lwr_log_div), ymax = exp(upr_log_div)), alpha = .25, fill = viridis_pal()(4)[3], alpha=.3) +
  annotate(geom="line",  x = seg_pred$SST_K-273.15, y = exp(seg_pred$log_div), size = 1, color = viridis_pal()(4)[3]) +
  scale_x_continuous(name = expression("Sea Surface Temperature ("*degree*"C)"), expand = c(0,0),
                     breaks = c(0, 10, 20, 30), labels = c(0, 10, 20, 30), sec.axis = sec_axis(~.)) +
  scale_y_continuous(name = "Percent Total Generic Diversity", expand = c(0,0), sec.axis = sec_axis(~.), trans = "log",
                     breaks = c(0.01, 0.03, 0.1, 0.3, 1), labels = c(0.01, 0.03, 0.1, 0.3, 1) * 100) +
  coord_cartesian(xlim = c(-3, 38), ylim = c(0.0075, 1)) +
  theme_classic(base_size = 24) +
  theme(axis.ticks = element_line(color = "black"), axis.text = element_text(colour = "black"),
        axis.title = element_blank(),
        axis.line = element_line(color = "black", lineend = "square"), 
        axis.ticks.x.top = element_blank(), axis.text.x.top = element_blank(),
        axis.ticks.y.right = element_blank(), axis.text.y.right = element_blank(),
        plot.margin = margin(.2,.2,.2,.2, "lines"), panel.border = element_rect(fill = NA),
        legend.key.size = unit(.25, "lines"), legend.spacing.y = unit(0.25, "lines"),
        legend.justification=c(1,1), legend.position=c(.4,.98),
        legend.background = element_rect(fill=alpha('white', 0.6)))+
  annotate(geom="text", y=80/100, x=10, label="Nautilus pompillius", fontface=3, size=6)

Sepia.officinalis_plot <- ggplot(Sepia.officinalis_rich, aes(x=(1000/Temp.allen)-273.15, ymin=exp(log_rich_0.05-c.opt), ymax=exp(log_rich_0.95-c.opt)))+
  geom_ribbon(alpha=.8, fill=viridis_pal()(4)[1])+
  geom_ribbon(data = seg_pred, aes(x = SST_K-273.15, ymin = exp(lwr_log_div), ymax = exp(upr_log_div)), alpha = .25, fill = viridis_pal()(4)[3], alpha=.3) +
  annotate(geom="line",  x = seg_pred$SST_K-273.15, y = exp(seg_pred$log_div), size = 1, color = viridis_pal()(4)[3]) +
  scale_x_continuous(name = expression("Sea Surface Temperature ("*degree*"C)"), expand = c(0,0),
                     breaks = c(0, 10, 20, 30), labels = c(0, 10, 20, 30), sec.axis = sec_axis(~.)) +
  scale_y_continuous(name = "Percent Total Generic Diversity", expand = c(0,0), sec.axis = sec_axis(~.), trans = "log",
                     breaks = c(0.01, 0.03, 0.1, 0.3, 1), labels = c(0.01, 0.03, 0.1, 0.3, 1) * 100) +
  coord_cartesian(xlim = c(-3, 38), ylim = c(0.0075, 1)) +
  theme_classic(base_size = 24) +
  theme(axis.ticks = element_line(color = "black"), axis.text = element_text(colour = "black"),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.line = element_line(color = "black", lineend = "square"), 
        axis.ticks.x.top = element_blank(), axis.text.x.top = element_blank(),
        axis.ticks.y.right = element_blank(), axis.text.y.right = element_blank(),
        plot.margin = margin(.2,.2,.2,.2, "lines"), panel.border = element_rect(fill = NA),
        legend.key.size = unit(.25, "lines"), legend.spacing.y = unit(0.25, "lines"),
        legend.justification=c(1,1), legend.position=c(.4,.98),
        legend.background = element_rect(fill=alpha('white', 0.6)))+
  annotate(geom="text", y=80/100, x=8, label="Sepia officinalis", fontface=3, size=6)

# Generate summary panel
summary.fig <-ggarrange2(Mytilus.species_plot, Crassostrea.gigas_plot, Pecten.maximus_plot, Gibberulus.gibbosus_plot, 
                         Nautilus.pompillius_plot, Dosidicus.gigas_plot, Sepia.officinalis_plot, ncol=4)

# Add axis labels
summary.fig.labels <- annotate_figure(summary.fig,
                                      bottom = text_grob(expression("Sea Surface Temperature ("*degree*"C)"), size = 30), 
                                      left = text_grob("Percent Total Generic Diversity", size = 30, rot=90))
                
ggsave("Supp metabolic model 10000 (molluscs only).pdf", summary.fig.labels, width = 16, height = 8)
ggsave("Supp metabolic model 10000 (molluscs only).png", summary.fig.labels, width = 16, height = 8)

