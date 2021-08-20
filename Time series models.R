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


# READ IN TIME-SERIES DATA AND MODEL OUTPUT
# Select an file by removing # in front of name and placing # in front of other files
#data.density <- read_csv("Population data Nigeria.csv")
data <- as.data.frame(read_csv("Temperature response parameters.csv"))
data.density <- read_csv("Population data China.csv")

# Select time-series data
#data.TS <- subset(data.density, Plot=="B") # select plot A, B, or C
data.TS <- subset(data.density, location=="Dafeng" & species=="Apolygus_lucorum") # select location and species


# Read in model output
# Select an insect by removing # in front of name and placing # in front of other species
#data.model <- as.data.frame(read_csv("Time Series Clavigralla shadabi.csv"))
#sp.data <- subset(data, Species == "Clavigralla shadabi")
#data.model <- as.data.frame(read_csv("Time Series Clavigralla tomentosicollis Benin.csv"))
#sp.data <- subset(data, Species == "Clavigralla tomentosicollis Benin")
data.model <- as.data.frame(read_csv("Time Series Clavigralla tomentosicollis Nigeria B.csv"))
sp.data <- subset(data, Species == "Clavigralla tomentosicollis Nigeria B")
#data.model <- as.data.frame(read_csv("Time Series Clavigralla tomentosicollis Burkina Faso.csv"))
#sp.data <- subset(data, Species == "Clavigralla tomentosicollis Burkina Faso")
data.model <- as.data.frame(read_csv("Time Series Apolygus lucorum China Dafeng.csv"))
sp.data <- subset(data, Species == "Apolygus lucorum China Dafeng")


# SET PLOT OPTIONS
# Default: plot last 2 year of model
xmin <- 0 #200
xmax <- 720 #750
ymin <- 0
ymax <- 100 #10
yr <- 360
init_yrs <- 10 # number of years to initiate the model (from Python DDE model)
TS.length <- xmax - xmin # length of time-series data
end <- nrow(data.model)


# FORMAT MODEL OUTPUT TO ALIGN WITH TIME-SERIES DATA
# Remove all rows before time-series data starts
# for Clavigralla tomentosicollis
#endCut <- (1972-1961+init_yrs)*yr + xmin
#data.model <- data.model[c(-1:-endCut), ]
# for Apolygus lucorum
endCut <-(2005-1961+init_yrs)*yr + xmin
data.model <- data.model[c(-1:-endCut), ]

# Remove all rows after xmax days
startCut <- xmax+1
data.model <- data.model[c(-startCut:-end), ]

# Re-scale time so that last x years are plotted as the first x years
#data.model <- sweep(data.model, 2, c((1972-1961+init_yrs)*yr+xmin,0,0,0,0)) # for Clavigralla tomentosicollis
data.model <- sweep(data.model, 2, c((2005-1961+init_yrs)*yr+xmin,0,0,0,0)) # for Apolygus lucorum


# PLOT TIME-SERIES DATA
# Juvenile density
plot.J = ggplot(data.TS, aes(x=time, y=J, ymin=J-J_SE, ymax=J+J_SE)) + 
  geom_pointrange(size=0.8, color="#d1495b") +
  #geom_line(size=0.8, linetype="longdash", color="#d1495b") +
  labs(x="Time", y="Log(Density)") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(ymin, ymax)) +
  #scale_y_log10(limits=c(ymin, ymax)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20))
#plot.J

# Adult density
plot.A = ggplot(data.TS, aes(x=time, y=A, ymin=A-A_SE, ymax=A+A_SE)) + 
  geom_pointrange(size=0.8, color="#30638e") +
  #geom_line(size=0.8, linetype="longdash", color="#30638e") +
  labs(x="Time", y="Log(Density)") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(ymin, ymax)) +
  #scale_y_log10(limits=c(ymin, ymax)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20)) 
#plot.A


# PLOT MODEL DATA
# Juvenile density
model.J = ggplot(data.model, aes(x=Time, y=J)) + 
  geom_line(size=0.8, color="#d1495b") +
  labs(x="Time", y="Log(Density)") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(ymin, ymax)) +
  #scale_y_log10(limits=c(ymin, ymax)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20))
#model.J

# Adult density
model.A = ggplot(data.model, aes(x=Time, y=A)) + 
  geom_line(size=0.8, color="#30638e") +
  labs(x="Time", y="Log(Density)") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(ymin, ymax)) +
  #scale_y_log10(limits=c(ymin, ymax)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20)) 
#model.A

# Insect density (juveniles + adults)
model.I = ggplot(data.model, aes(x=Time, y=J+A)) + 
  geom_line(size=0.8, color="#30638e") +
  labs(x="Time", y="Log(Density)") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(ymin, ymax)) +
  #scale_y_log10(limits=c(ymin, ymax)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20)) 
#model.I


# PLOT HABITAT TEMPERATURE FUNCTION
plot.temp <- ggplot() +
  geom_function(fun = function(t) sp.data$meanT + sp.data$amplT*sin(2*pi*(t + sp.data$shiftT)/yr),
                size=0.8, linetype="longdash", color="#d1495b") +
  labs(x="Time", y="T (K)") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(sp.data$meanT - sp.data$amplT - 1, sp.data$meanT + sp.data$amplT + 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20))
#plot.temp


# DRAW FINAL PLOTS
ggdraw()  +
  draw_plot(plot.temp, x = 0, y = 0, width = 1, height = 0.3) +
  #draw_plot(plot.J, x = 0, y = 0.3, width = 1, height = 0.7) +
  draw_plot(plot.A, x = 0, y = 0.3, width = 1, height = 0.7) +
  draw_plot(model.J, x = 0, y = 0.3, width = 1, height = 0.7) +
  draw_plot(model.A, x = 0, y = 0.3, width = 1, height = 0.7)
