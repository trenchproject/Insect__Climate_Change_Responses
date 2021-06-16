####################################################################################################
#### This R script simulates population dynamics of hemipteran insects in their native habitats ####
####################################################################################################


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
data <- read_csv("Temperature response parameters.csv")

# select an insect by removing # in front of name and placing # in front of other species
#sp.data <- subset(data, Species == "Clavigralla shadabi")
#sp.data <- subset(data, Species == "Clavigralla tomentosicollis Benin")
sp.data <- subset(data, Species == "Clavigralla tomentosicollis Nigeria B")
#sp.data <- subset(data, Species == "Clavigralla tomentosicollis Burkina Faso A")

# Read time-series data (if applicable)
data.density <- read_csv("Egwuatu_1977.csv")
# Set the plot (A, B, or C)
data.plot <- subset(data.density, Plot=="B")

# define model parameters
params <- c(sp.data[2], sp.data[3], sp.data[4], sp.data[5],sp.data[6], sp.data[7], sp.data[8],
            sp.data[9], sp.data[10], sp.data[11], sp.data[12], sp.data[13], sp.data[14], sp.data[15],
            sp.data[16], sp.data[17], sp.data[18], sp.data[19], sp.data[20], sp.data[21], sp.data[22],
            sp.data[23], sp.data[24], sp.data[25], sp.data[26], sp.data[27], sp.data[28])
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
         # temperature function
         temp = (meanT + dMeanT*times) + (amplT + dAmplT*times)*sin(2*pi*(times + shiftT)/365)
         # resource variation parameters
         if(Rend > Rstart) {
           Rlength = (Rend - Rstart)/365 # proportion of year when resource is available
         } else { Rlength = (365 - Rstart + Rend)/365 }
         Rshift = cos(pi*Rlength) # shift used to model resource availability via sine wave
         # resource availability function
         dum = 1/(1-Rshift)*(-Rshift + sin(2*pi*(times-Rstart)/365 + asin(Rshift))) # sine wave describing resource availability
         if(Res == 1) { res = 0.5*(abs(dum) + dum) } else { res = 1 }
         # time-series of temperature and resource values
         signalT = as.data.frame(list(times = times, temp = rep(0, length(times))))
         signalT$temp = temp
         signalR = as.data.frame(list(times = times, res = rep(0, length(times))))
         signalR$res = res
         # create interpolating functions
         inputT = approxfun(signalT, rule = 2 ) # rule = 2 sets any extrapolated point to the closest data extreme
         inputR = approxfun(signalR, rule = 2 ) # rule = 2 sets any extrapolated point to the closest data extreme
         # import interpolated functions
         T =  inputT(t)
         R =  inputR(t)
         
         # Temperature responses
         b = R*bTopt*exp(-(T-Toptb)^2/(2*sb^2)) # Birth rate
         #q = 1 # Density-dependence (temperature-independent for simplicity)
         q = qTR*exp(Aq*(1/TR-1/Tmax))*exp(-(T-Toptq)^2/(2*sq^2)) # Density-dependence (temperature-dependent)
         m = mTR*(T/TR)*exp(A*(1/TR-1/T))/(1+skew*(exp(AL*(1/TL-1/T))+exp(AH*(1/TH-1/T)))) # Maturation rate
         dJ = dJTR*exp(AdJ*(1/TR-1/T)) # Juvenile mortality
         dA = dATR*exp(AdA*(1/TR-1/T)) # Adult mortality
         
         # Juvenile and adult densities
         dJdt = b*(A >= 0.1)*A*exp(-q*A) - m*J - dJ*J # Allee threshold of 0.1
         dAdt = m*J - dA*A
         
         return(list(c(dJdt, dAdt), signalT = T))
       })
}


# STEP 3: Model simulation and plotting population dynamics
init = c(J = 1, A = 1)
# run model and input into data frame
densities = as.data.frame(ode(func = model, y = init, parms = params, times = times, method = "ode45"))
# gather model variables (time, temp, life stage, density)
output = densities[seq(0, dim(densities)[1], by=timestep), ]
output = gather(output, key=Variable, value=Output, -time)

# Plot population dynamics
pop.dyn <- ggplot(output[output$Variable %in% c("J", "A"), ], aes(x=time, y=Output, color=Variable)) + 
  geom_line(size=1.2) +
  scale_color_manual(values=c("J"="#d1495b", "A"="#30638e")) + 
  labs(x="Time", y="Density") +
  scale_y_continuous(limits=c(0, 100)) +
  #scale_y_log10(limits=c(1, 10)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none", 
        axis.text = element_text(size=13),
        axis.title = element_text(size=20)) 

# Plot temperature regime
temp.regime <- ggplot(output[output$Variable %in% c("signalT"), ], aes(x=time, y=Output, color=Variable)) + 
  geom_line(size=1.2) +
  scale_color_manual(values=c("signalT"="#d1495b")) + 
  labs(x="Time", y="Temperature (K)") +
  scale_y_continuous(limits=c(params$meanT - params$amplT - amplT.incr - 1, params$meanT + params$amplT + meanT.incr + amplT.incr + 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none", 
        axis.text = element_text(size=13),
        axis.title = element_text(size=20))

ggdraw() +
  draw_plot(pop.dyn, x = 0, y = 0.5, width = 1, height = 0.5) +
  draw_plot(temp.regime, x = 0, y = 0, width = 1, height = 0.5)


####################################################################
#### Plot empirical time-series for Clavigralla tomentosicollis ####
####################################################################
# Plot options
xmin <- 200 # days
xmax <- 750 # days
ymin <- 0
ymax <- 10

# Plot time-series data
plot.J <- ggplot(data.plot, aes(x=time, y=J, ymin=J-J_SE, ymax=J+J_SE)) + 
  geom_pointrange(size=1.2, color="#d1495b") +
  #geom_line(size=0.8, linetype="longdash", color="#d1495b") +
  labs(x="Time", y="Density") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(ymin, ymax)) +
  #scale_y_log10(limits=c(1, 10)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20))
#plot.J

plot.A <- ggplot(data.plot, aes(x=time, y=A, ymin=A-A_SE, ymax=A+A_SE)) + 
  geom_pointrange(size=1.2, color="#30638e") +
  #geom_line(size=0.8, linetype="longdash", color="#30638e") +
  labs(x="Time", y="Density") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(ymin, ymax)) +
  #scale_y_log10(limits=c(1, 10)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20)) 
#plot.A

# Plot model dynamics (starting 8 years in to avoid transients)
plot.model <- ggplot(subset(output[output$Variable %in% c("J", "A"), ], time>=8*365+xmin & time<=8*365+xmax),
                          aes(x=time, y=Output, color=Variable)) + 
  geom_line(size=1.2) +
  scale_color_manual(values=c("J"="#d1495b", "A"="#30638e")) + 
  labs(x="Time", y="Density") +
  scale_y_continuous(limits=c(ymin, ymax)) +
  #scale_y_log10(limits=c(1, 10)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20),
        axis.text.x = element_text(colour = "white"), axis.ticks.x = element_blank())
#plot.model

# Plot temperature regime
plot.temp <- ggplot(subset(output[output$Variable %in% c("signalT"), ], time>=xmin & time<=xmax),
                      aes(x=time, y=Output, color=Variable)) + 
  geom_line(size=1.2) +
  scale_color_manual(values=c("signalT"="#d1495b")) + 
  labs(x="Time", y="T (K)") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(params$meanT - params$amplT - amplT.incr - 1, params$meanT + params$amplT + meanT.incr + amplT.incr + 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20)) 
#plot.temp

ggdraw()  +
    draw_plot(plot.temp, x = 0, y = 0, width = 1, height = 0.3) +
    draw_plot(plot.model, x = 0, y = 0.3, width = 1, height = 0.7) +
    draw_plot(plot.J, x = 0, y = 0.3, width = 1, height = 0.7) +
    draw_plot(plot.A, x = 0, y = 0.3, width = 1, height = 0.7)
