########################################################################
#### This R script fits model outputs to empirical time-series data ####
########################################################################


# Load packages and set working directory
library(tidyr)
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyverse)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# READ IN TEMPERATURE RESPONSE PARAMETERS AND TIME-SERIES DATA
# Temperature response parameters
data <- as.data.frame(read_csv("Temperature response parameters.csv"))

# Time-series data
# Select an file by removing # in front of name and placing # in front of other files
data.density <- read_csv("Population data Nigeria.csv")
#data.density <- read_csv("Population data China.csv")

# Select time-series data
data.TS <- subset(data.density, Plot=="B") # select plot A, B, or C
#data.TS <- subset(data.density, location=="Dafeng" & species=="Apolygus_lucorum") # select location and species


# READ IN MODEL OUTPUT
# From DDE population dynamics.py 
# Select an insect by removing # in front of name and placing # in front of other species
#data.model <- as.data.frame(read_csv("Time Series Clavigralla shadabi.csv"))
#sp.data <- subset(data, Species == "Clavigralla shadabi")
#data.model <- as.data.frame(read_csv("Time Series Clavigralla tomentosicollis Benin.csv"))
#sp.data <- subset(data, Species == "Clavigralla tomentosicollis Benin")
data.model <- as.data.frame(read_csv("Time Series Clavigralla tomentosicollis Nigeria B.csv"))
sp.data <- subset(data, Species == "Clavigralla tomentosicollis Nigeria B")
#data.model <- as.data.frame(read_csv("Time Series Clavigralla tomentosicollis Burkina Faso.csv"))
#sp.data <- subset(data, Species == "Clavigralla tomentosicollis Burkina Faso")
#data.model <- as.data.frame(read_csv("Time Series Apolygus lucorum China Dafeng.csv"))
#sp.data <- subset(data, Species == "Apolygus lucorum China Dafeng")

# Population dynamics with climate change
clim.data <- data.model


# SET PLOT OPTIONS
# Default: plot last 2 year of model
# For historical time period
xmin <- 0 #200
xmax <- 750 #720
ymin <- 0
ymax <- 10 #60
# For climate change time period
xmin.CC <- 0 #200
xmax.CC <- 750 #750
ymin.CC <- 0
ymax.CC <- 10 #60
yr <- 360 # days in a year (using 360 for simplicity)
init_yrs <- 10 # number of years to initiate the model (from Python DDE model)
TS.length <- xmax - xmin # length of time-series data
end <- nrow(data.model)


# FORMAT MODEL OUTPUT TO ALIGN WITH TIME-SERIES DATA
# Remove all rows before time-series data starts

# for Clavigralla tomentosicollis
data.model <- data.model[c(-1:-((1972-1961+init_yrs)*yr + xmin)), ]
# for Apolygus lucorum
#data.model <- data.model[c(-1:-((2005-1961+init_yrs)*yr + xmin)), ]
# climate change period (remove all but last 2 years of data)
clim.data <- clim.data[c(-1:-(end-2*yr + xmin.CC)), ]

# Remove all rows after xmax days
data.model <- data.model[c(-(xmax+1):-end), ]

# Re-scale time to start at xmin
data.model <- sweep(data.model, 2, c(data.model[[1,1]]-xmin,0,0,0,0))
# with climate change
CC.steps <- clim.data[[1,1]] + xmin.CC # model time-steps before climate change period
clim.data <- sweep(clim.data, 2, c(CC.steps,0,0,0,0))

# Log transform density if needed
data.model$J <- log(data.model$J, 10)
data.model$A <- log(data.model$A, 10)
clim.data$J <- log(clim.data$J, 10)
clim.data$A <- log(clim.data$A, 10)



# PLOT TIME-SERIES DATA
# Juvenile density
plot.J = ggplot(data.TS, aes(x=time, y=J, ymin=J-J_SE, ymax=J+J_SE)) + 
  geom_pointrange(size=0.8, color="#d1495b") +
  #geom_line(size=0.8, linetype="longdash", color="#d1495b") +
  labs(x="Time", y="Density") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(ymin, ymax)) +
  #scale_y_log10(limits=c(ymin, ymax)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20))

# Adult density
plot.A = ggplot(data.TS, aes(x=time, y=A, ymin=A-A_SE, ymax=A+A_SE)) + 
  geom_pointrange(size=0.8, color="#30638e") +
  #geom_line(size=0.8, linetype="longdash", color="#30638e") +
  labs(x="Time", y="Density") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(ymin, ymax)) +
  #scale_y_log10(limits=c(ymin, ymax)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20)) 



# PLOT MODEL DATA
# Historical time period
# Juvenile density
model.J = ggplot(data.model, aes(x=Time, y=J)) + 
  geom_line(size=0.8, color="#d1495b") +
  labs(x="Time", y="Density") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(ymin, ymax)) +
  #scale_y_log10(limits=c(ymin, ymax)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20))

# Adult density
model.A = ggplot(data.model, aes(x=Time, y=A)) + 
  geom_line(size=0.8, color="#30638e") +
  labs(x="Time", y="Density") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(ymin, ymax)) +
  #scale_y_log10(limits=c(ymin, ymax)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20)) 

# Insect density (juveniles + adults)
model.I = ggplot(data.model, aes(x=Time, y=J+A)) + 
  geom_line(size=0.8, color="black") +
  labs(x="Time", y="Density") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(ymin, ymax)) +
  #scale_y_log10(limits=c(ymin, ymax)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20)) 


# Climate change time period
# Juvenile density
model.J.CC = ggplot(clim.data, aes(x=Time, y=J)) + 
  geom_line(size=0.8, color="#d1495b") +
  labs(x="Time", y="Density") +
  scale_x_continuous(limits=c(xmin.CC, xmax.CC)) +
  scale_y_continuous(limits=c(ymin.CC, ymax.CC)) +
  #scale_y_log10(limits=c(ymin.CC, ymax.CC)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20))

# Adult density
model.A.CC = ggplot(clim.data, aes(x=Time, y=A)) + 
  geom_line(size=0.8, color="#30638e") +
  labs(x="Time", y="Density") +
  scale_x_continuous(limits=c(xmin.CC, xmax.CC)) +
  scale_y_continuous(limits=c(ymin.CC, ymax.CC)) +
  #scale_y_log10(limits=c(ymin.CC, ymax.CC)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20)) 

# Insect density (juveniles + adults)
model.I.CC = ggplot(clim.data, aes(x=Time, y=J+A)) + 
  geom_line(size=0.8, color="black") +
  labs(x="Time", y="Density") +
  scale_x_continuous(limits=c(xmin.CC, xmax.CC)) +
  scale_y_continuous(limits=c(ymin.CC, ymax.CC)) +
  #scale_y_log10(limits=c(ymin.CC, ymax.CC)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20))



# PLOT HABITAT TEMPERATURE FUNCTION
# Historical time period
plot.temp <- ggplot() +
  geom_function(fun = function(t) (sp.data$meanT + sp.data$delta_mean*t)  + (sp.data$amplT + sp.data$delta_ampl*t) * sin(2*pi*(t + sp.data$shiftT)/yr),
                size=0.8, linetype="longdash", color="#d1495b") +
  labs(x="Time", y="T (K)") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(sp.data$meanT - abs(sp.data$amplT) - 1, sp.data$meanT + abs(sp.data$amplT) + 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20))

# climate change time period
plot.temp.CC <- ggplot() +
  geom_function(fun = function(t) (sp.data$meanT + sp.data$delta_mean*(t+CC.steps))  + (sp.data$amplT + sp.data$delta_ampl*(t+CC.steps)) * sin(2*pi*((t+CC.steps) + sp.data$shiftT)/yr),
                size=0.8, linetype="longdash", color="#d1495b") +
  labs(x="Time", y="T (K)") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(sp.data$meanT - abs(sp.data$amplT) + sp.data$delta_mean*CC.steps, sp.data$meanT + abs(sp.data$amplT) + + sp.data$delta_mean*CC.steps)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20))


# DRAW FINAL PLOTS
# Historical time period
# Juveniles and adults
plot <- ggdraw()  +
  draw_plot(plot.temp, x = 0, y = 0, width = 1, height = 0.3) +
  draw_plot(plot.J, x = 0, y = 0.3, width = 1, height = 0.7) +
  draw_plot(plot.A, x = 0, y = 0.3, width = 1, height = 0.7) +
  draw_plot(model.J, x = 0, y = 0.3, width = 1, height = 0.7) +
  draw_plot(model.A, x = 0, y = 0.3, width = 1, height = 0.7)
plot

# Total insects
#plot <- ggdraw()  +
  #draw_plot(plot.temp, x = 0, y = 0, width = 1, height = 0.3) +
  #draw_plot(plot.A, x = 0, y = 0.3, width = 1, height = 0.7) +
  #draw_plot(model.J, x = 0, y = 0.3, width = 1, height = 0.7) +
  #draw_plot(model.A, x = 0, y = 0.3, width = 1, height = 0.7) +
  #draw_plot(model.I, x = 0, y = 0.3, width = 1, height = 0.7)
#plot
#plot

# Climate change time period
# Juveniles and adults
plot.CC <- ggdraw()  +
draw_plot(plot.temp.CC, x = 0, y = 0, width = 1, height = 0.3) +
draw_plot(model.J.CC, x = 0, y = 0.3, width = 1, height = 0.7) +
draw_plot(model.A.CC, x = 0, y = 0.3, width = 1, height = 0.7)
plot.CC

# Total insects
#plot.CC <- ggdraw()  +
#  draw_plot(plot.temp.CC, x = 0, y = 0, width = 1, height = 0.3) +
#  draw_plot(model.J.CC, x = 0, y = 0.3, width = 1, height = 0.7) +
#  draw_plot(model.A.CC, x = 0, y = 0.3, width = 1, height = 0.7) +
#  draw_plot(model.I.CC, x = 0, y = 0.3, width = 1, height = 0.7)
#plot.CC
