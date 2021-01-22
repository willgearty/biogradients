library(ggplot2)
library(viridis)
library(deeptime)

## Define number of repetitions for random number generation
reps <- 5000
increments <- 0.001

# set seed for random number generation
set.seed(1993)

# input Penn et al. 2018 parameter values
u_E <- 0.41 # eV
sig_E <- 0.29 # eV
c_E <- 1.40 

u_A <- 3.01
sig_A <- 0.49
c_A <- 0.05

u_C <- 1.10
sig_C <- 0.42
c_C <- 0.35


# sample Ao values from PDF
A0_vec <- seq(0.1,80,0.1)
A0.1000 <- dlnorm(A0_vec, meanlog = u_A, sdlog = sig_A)/c_A
A0.sum <- data.frame(A0_vec, A0.1000)

# sample phi crit values from PDF
phi_crit_vec <- seq(0.01,12,.01)
phi_crit.1000 <- dlnorm(phi_crit_vec, meanlog = u_C, sdlog = sig_C)/c_C
phi_crit.sum <- data.frame(phi_crit_vec, phi_crit.1000)

E0_vec <- seq(-0.6,1.5,0.001)
E0.1000 <- dnorm(E0_vec, mean = u_E, sd = sig_E)/c_E
E0.sum <- data.frame(E0_vec, E0.1000)

# Make individual plots
A0.plot <- ggplot(A0.sum, aes(x=A0_vec, ymax=A0.1000, ymin=-Inf))+
geom_ribbon(alpha=.9, fill=viridis_pal()(4)[1])+                    
  theme_classic(base_size = 24) + xlab(expression("A"[o]))+ ylab("Frequency")+
  scale_x_continuous(sec.axis = sec_axis(~.)) +
  scale_y_continuous(sec.axis = sec_axis(~.)) +
  coord_cartesian(ylim=c(0,1), expand=c(0,0))+
  theme(axis.ticks = element_line(color = "black"), axis.text = element_text(colour = "black"),
        axis.line = element_line(color = "black", lineend = "square"), 
        axis.ticks.x.top = element_blank(), axis.text.x.top = element_blank(),
        axis.ticks.y.right = element_blank(), axis.text.y.right = element_blank(),
        plot.margin = margin(1,1,1,1, "lines"), panel.border = element_rect(fill = NA),
        legend.key.size = unit(.25, "lines"), legend.spacing.y = unit(0.25, "lines"),
        legend.justification=c(1,1), legend.position=c(.4,.98),
        legend.background = element_rect(fill=alpha('white', 0.6)))
  
  E0.plot <- ggplot(E0.sum, aes(x=E0_vec, ymax=E0.1000, ymin=-Inf))+
    geom_ribbon(alpha=.9, fill=viridis_pal()(4)[1])+                    
    theme_classic(base_size = 24) + xlab(expression("E"[o]))+ ylab("Frequency")+
    scale_x_continuous(sec.axis = sec_axis(~.)) +
    scale_y_continuous(sec.axis = sec_axis(~.)) +
    coord_cartesian(ylim=c(0,1), expand=c(0,0))+
  theme(axis.ticks = element_line(color = "black"), axis.text = element_text(colour = "black"),
        axis.line = element_line(color = "black", lineend = "square"), 
        axis.ticks.x.top = element_blank(), axis.text.x.top = element_blank(),
        axis.ticks.y.right = element_blank(), axis.text.y.right = element_blank(),
        plot.margin = margin(1,1,1,1, "lines"), panel.border = element_rect(fill = NA),
        legend.key.size = unit(.25, "lines"), legend.spacing.y = unit(0.25, "lines"),
        legend.justification=c(1,1), legend.position=c(.4,.98),
        legend.background = element_rect(fill=alpha('white', 0.6)),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())
  
  phi.crit.plot <- ggplot(phi_crit.sum, aes(x=phi_crit_vec, ymax=phi_crit.1000, ymin=-Inf))+
    geom_ribbon(alpha=.9, fill=viridis_pal()(4)[1])+                    
    theme_classic(base_size = 24) + xlab(expression(phi[crit]))+ ylab("Frequency")+
    scale_x_continuous(sec.axis = sec_axis(~.)) +
    scale_y_continuous(sec.axis = sec_axis(~.)) +
    coord_cartesian(ylim=c(0,1), expand=c(0,0))+
  theme(axis.ticks = element_line(color = "black"), axis.text = element_text(colour = "black"),
        axis.line = element_line(color = "black", lineend = "square"), 
        axis.ticks.x.top = element_blank(), axis.text.x.top = element_blank(),
        axis.ticks.y.right = element_blank(), axis.text.y.right = element_blank(),
        plot.margin = margin(1,1,1,1, "lines"), panel.border = element_rect(fill = NA),
        legend.key.size = unit(.25, "lines"), legend.spacing.y = unit(0.25, "lines"),
        legend.justification=c(1,1), legend.position=c(.4,.98),
        legend.background = element_rect(fill=alpha('white', 0.6)),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())

# Combine plots
freq.sums <- ggarrange2(A0.plot, E0.plot, phi.crit.plot, ncol=3)

# Save plots
ggsave("Metabolic index distribution summary.pdf", freq.sums, width = 20, height = 7)
  
    