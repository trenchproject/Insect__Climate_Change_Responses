#################################################################################
#### This R script simulates population dynamics of hemipteran insects in their native habitats
#################################################################################



#################################################
#### Load packages and set working directory ####
#################################################
library(tidyr)
library(ggplot2)
library(deSolve)
library(cowplot)
library(dplyr)
library(tidyverse)
#library(Cairo)
#library(grid)
#library(gridExtra)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

######################################
#### Simulate population dynamics ####
######################################

# STEP 1: Assign parameter values
# read data and select insect
data <- read_csv("Temperature response data.csv")

# select an insect by removing # in front of name and placing # in front of other species
#sp.data <- subset(data, Species == "Clavigralla shadabi")
sp.data <- subset(data, Species == "Clavigralla tomentosicollis")
# define model parameters
params <- c(sp.data[2], sp.data[3], sp.data[4], sp.data[5],sp.data[6], sp.data[7], sp.data[8],
            sp.data[9], sp.data[10], sp.data[11], sp.data[12], sp.data[13], sp.data[14], sp.data[15],
            sp.data[16], sp.data[17], sp.data[18], sp.data[19], sp.data[20], sp.data[21], sp.data[22],
            sp.data[23], sp.data[24])
# time
yrs <- 10; # Number of years over which to run the model
timestep <- 1; # Length of time-steps (days)
times <- seq(0, yrs*365, by = timestep)
# temperature parameters
meanT.incr = 0 # Increase in mean temperature over next 50 years
amplT.incr = 0 # Increase in amplitude of temperature fluctuations over next 50 years
dMeanT = meanT.incr/(yrs*365) # simulate x degree increase in mean temperature over 50 years
dAmplT = amplT.incr/(yrs*365) # simulate x degree increase in the amplitude of temperature fluctuations over 50 years


# STEP 2: Specify model
model = function(t, State, Pars){
  with(as.list(c(State, Pars)),
       {
         # Temperature regime
         #meanT.incr = 0 # Increase in mean temperature over next 50 years
         #amplT.incr = 0 # Increase in amplitude of temperature fluctuations over next 50 years
         #dMeanT = meanT.incr/(yrs*365) # simulate x degree increase in mean temperature over 50 years
         #dAmplT = amplT.incr/(yrs*365) # simulate x degree increase in the amplitude of temperature fluctuations over 50 years
         # temperature function
         temp = (meanT + dMeanT*times) + (amplT + dAmplT*times)*sin((2*pi*times + deltaT)/365)
         # time-series of temperature values
         signal = as.data.frame(list(times = times, temp = rep(0, length(times))))
         signal$temp = temp
         # create interpolating function
         input = approxfun(signal, rule = 2 ) # rule = 2 sets any extrapolated point to the closest data extreme
         # import interpolated temperature function
         T =  input(t)
         
         # Temperature responses
         b = bTopt*exp(-(T-Toptb)^2/(2*sb^2)) # Birth rate
         q = qTR*exp(Aq*(1/TR-1/Tmax))*exp(-(T-Toptq)^2/(2*sq^2)) # Density-dependence
         m = mTR*(T/TR)*exp(A*(1/TR-1/T))/(1+skew*(exp(AL*(1/TL-1/T))+exp(AH*(1/TH-1/T)))) # Maturation rate
         dJ = dJTR*exp(AdJ*(1/TR-1/T)) # Juvenile mortality
         dA = dATR*exp(AdA*(1/TR-1/T)) # Adult mortality
         
         # Juvenile and adult densities
         dJdt = b*(A >= 0.2)*A*exp(-q*A) - m*J - dJ*J # Allee threshold of 10% of equilibrium adult density
         dAdt = m*J - dA*A
         
         return(list(c(dJdt, dAdt), signal = T))
       })
}


# STEP 3: Model simulation and plotting population dynamics
init = c(J = 1, A = 1)
# run model and input into data frame
densities = as.data.frame(ode(func = model, y = init, parms = params, times = times, method = "ode45"))
# model variables (time, life stage, density)
output = densities[seq(0, dim(densities)[1], by=timestep), ]
output = gather(output, key=Variable, value=Output, -time)

# Plot population dynamics
pop.dyn = ggplot(output[output$Variable %in% c("J", "A"), ], aes(x=time, y=Output, color=Variable)) + 
  geom_line(size=1.2) +
  scale_color_manual(values=c("J"="#d1495b", "A"="#30638e")) + 
  labs(x="Time", y="Density") +
  # scale_x_continuous(expand=c(0, 5)) +
  scale_y_log10(limits=c(0.1, 100)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none", 
        axis.text = element_text(size=13),
        axis.title = element_text(size=20)) 

# Plot temperature regime
temp.regime  = ggplot(output[output$Variable %in% c("signal"), ], aes(x=time, y=Output, color=Variable)) + 
  geom_line(size=1.2) +
  scale_color_manual(values=c("signal"="#d1495b")) + 
  labs(x="Time", y="Temperature (K)") +
  # scale_x_continuous(expand=c(0, 5)) +
  scale_y_continuous(limits=c(298 - amplT.incr, 302 + meanT.incr + amplT.incr)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none", 
        axis.text = element_text(size=13),
        axis.title = element_text(size=20))

Plots = ggdraw() +
  draw_plot(pop.dyn, x = 0, y = 0.5, width = 1, height = 0.5) +
  draw_plot(temp.regime, x = 0, y = 0, width = 1, height = 0.5) #+
#draw_plot_label(label = c("(c)", ""), size = 18, x = c(0, 0), y = c(1, 0.5))
Plots


####################################################################
#### Plot empirical time-series for Clavigralla tomentosicollis ####
####################################################################

# Read data
data.plotA <- read_csv("Egwuatu_1977a.csv")

# Convert to data frame
data.A <- as.data.frame(data.plotA) %>%
  gather("stage","popdens",4:5)

# Plot time-series data
plotA = ggplot(data.A, aes(x=time, y=popdens, color=stage)) + 
  geom_line(size=1.2) +
  scale_color_manual(values=c("J"="#d1495b", "A"="#30638e")) + 
  labs(x="Time", y="Density") +
  # scale_x_continuous(expand=c(0, 5)) +
  scale_y_log10(limits=c(0.2, 10)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none", 
        axis.text = element_text(size=13),
        axis.title = element_text(size=20)) 
plotA