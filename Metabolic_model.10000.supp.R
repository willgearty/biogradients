library(ggplot2)
library(dplyr)
library(viridis)
library(deSolve)
library(dplyr)
library(segmented)

## Define number of repetitions for random number generation
reps <- 10000
increments <- 0.001

# set seed for random number generation
set.seed(1993)

# input Penn et al. 2018 parameter values
u_E <- 0.41 # eV
sig_E <- 0.29 # eV

u_A = 3.01
sig_A = 0.49

u_C <- 1.10
sig_C <- 0.42

# sample Ao values from PDF
A0.10000 <- rlnorm(reps, meanlog = u_A, sdlog = sig_A)

# sample phi crit values from PDF
phi_crit.10000 <- rlnorm(reps, meanlog = u_C, sdlog = sig_C)

# Set Penn et al. 2018 constants
kB = 8.6173303*10^(-5)
Tref <- 15+273.15

# randomly generate pO2 values
pO2.10000 <- runif(reps, 0.1957, 0.2146) ## 5th and 95th percentiles of Boag et al 2018 surface data (and uniform distribution...)
## Units: atm

## Want to monte carlo initial value from paleo data
## Load paleo data
## Set appropriate working directory to load "fossil_results.RData"
load("fossil_results.RData")
fossil_data_sub <- subset(fossil_data, num_bands == 24 & metric == "SQS_0.25" & band_type == "equal-area")

# segmented linear regression
reg <- lm(log(div.prop) ~ SST_K, data = fossil_data_sub)
reg_seg <- with(fossil_data_sub, segmented(reg, seg.Z = ~ SST_K, npsi = 1))
seg_pred <- as.data.frame(cbind(SST_K = seq(267, 313, .01), predict(reg_seg, newdata = data.frame(SST_K = seq(267, 313, .01)), interval = "confidence", level = .95)))

# Calculate initial value and mean error across entire segmented linear model
mean_fossil_seed <- seg_pred$fit[1]
SD_fossil_seed <- mean(seg_pred$fit-seg_pred$lwr) # using 95% conf interval as SD of distribution to give approximately equal error quantification on both models

# Make vector of initial diversity values
S.vec <- rnorm(100, mean = mean_fossil_seed, sd = SD_fossil_seed)

# Make vector of Allen E values for uncertainty
E.vec <- runif(n = 100, min = 0.6, max = 0.7)

# Calculate lowest temperature in fossil data compilation
Min.temp <- summary(fossil_data_sub$SST_K)[1]
Min.temp.Allen <- 1000/Min.temp
model.summary <- data.frame(Temp.allen = double(), Allen.E=double(), log.richness=double(), Eo.cum.freq = double())

# Loop through Eo values for all molluscs in Penn et al. dataset
for(E0.val in c(0.21, 0.41, 0.16, 0.33, 0.62, 0.59, 0.72)){ 
  # Initiate temperature-specific summary data frame
  phi.temp.sum <- data.frame(A0.10000 = double(), E0.10000=double(), phi_crit.10000=double(), temp=double(), pO2.10000=double(), phi.10000=double(), iteration=double())

  # E0 values now prescribed in loop itself for this supplementary dataset - so two lines below (from main text script) commented out.
  # set Eo value based on place in cumulative frequency loop
  # E0.val <- qnorm(E0.quantile, mean = u_E, sd = sig_E)

  # How many SDs from the mean is the current Eo value? (for plotting purposes later)
  E0.val.sd <- abs(u_E-E0.val)/sig_E

  # Throughout this script, "Allen Temp" is used as shorthand to refer to 1000/T, with temperature in Kelvin
  # Allen temps of 3.685 TO 3.2 give a celsius temp span of -1.78 TO 39.35
  # Loop through these temperature values
  for(Allen.Temp in seq(3.685, 3.2, -increments)){
    # convert Allen temp to degrees K
    Temp <- (1000/Allen.Temp)

    # Calculate 1000 phi values, based on Ao sampling, for each temp looped through and current Eo value
    phi.10000 <- A0.10000 * (pO2.10000/(exp((-E0.val/kB)*((1/Temp)-(1/Tref)))))

    # Summarise this temperature
    phi.temp <- as.data.frame(cbind(A0.10000, rep(E0.val, reps), phi_crit.10000, rep(Allen.Temp, reps), rep(Temp, reps), pO2.10000, phi.10000, seq(1:reps)))
    # Combine with summary of all temperatures
    phi.temp.sum <- rbind(phi.temp.sum, phi.temp)
  }
  # Remove scenarios for which phi < phi crit (again, with Penn distributions)
  phi.crit.xxx.sum <- filter(phi.temp.sum, phi.10000>=phi_crit.10000)

  # Generate summary table - how many of the 10,000 ecotypes are still here
  survival.by.temp <- as.data.frame(table(phi.crit.xxx.sum$V4))

  # Make temperature a numerical vector
  survival.by.temp$Var1 <- as.numeric(paste(survival.by.temp$Var1))

  # Calculate the slope in log(proportional diversity) per temperature increment
  slope.survival <- diff(log(survival.by.temp$Freq/reps))

  # Solve Eo-specific differential equations
  # Loop through uncertainty in E values and initial values from fossil data
  for (Allen.E.no in 1:100){ # Take into account uncertainty on E parameter in Allen et al 2002 using Gillooly & Allen 2007
    for (init_val in 1:100){
    Allen.E  <- E.vec[Allen.E.no]
      pars <- c(kinet.slope = Allen.E/(1000*kB)*increments) # set kinetic slope based on Allen et al. 2002 eqn
      # because we are 'counting' downwards in Allen space, we remove the negative sign from the equation in
      # Allen et al. 2002, effectively reversing the axes in their study.

      solve.model <- function(pars){
        derivs <- function(t, y, pars){
          with(as.list(c(y,pars)), {

            aerob.slope <- -1 * slope.survival[match(round(Temp.allen, digits=3),
                                                round(survival.by.temp$Var1, digits=3))] # sign of slope reversed as we are 'counting' downwards in Allen space (as discussed above)

            # define model equations
            dTemp <- -increments # Cannot have "time" moving backwards in deSolve, so set temp here
            dS <- aerob.slope + kinet.slope

            # return values
            return(list(c(Temp.allen=dTemp, S = dS),
                        aerob.slope = aerob.slope , kinet.slope = kinet.slope, temp=t))
          })
        }

        # Define initial conditions based on compilation and starting temperature
        initial.conditions <- c(Temp.allen = (3.685-increments), S = S.vec[init_val])

        time.step <-  as.numeric(seq(0, (length(seq(3.2, 3.685, increments))-1),1))

        return(ode(y = initial.conditions, t = time.step, func = derivs, parms=pars, d="vode", maxsteps=50000))
      }
      x <- as.data.frame(solve.model(pars))

      model.summary <- rbind(model.summary,
                             cbind(x$Temp.allen, Allen.E, x$S, E0.val, E0.val.sd))
    }}
  print(E0.val)
}

names(model.summary) <- c("Temp.allen", "Allen.E", "log.richness", "Eo.val", "Eo.sd")

save(model.summary, file = "Metabolic.model.summary.10000.supp.RData")
