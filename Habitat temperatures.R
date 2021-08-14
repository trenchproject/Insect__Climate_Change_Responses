#############################################################################
#### This R script fits sinusoidal functions to habitat temperature data ####
#############################################################################


# Load packages and set working directory
library(tidyr)
library(ggplot2)
library(dplyr)
library(tidyverse)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


#### Fit temperature functions to data from Climate Wizard using Fourier analysis ####
# Select a location by removing # in front of name and placing # in front of other locations
data <- as.data.frame(read_csv("Nigeria climate data.csv"))

# Fit sinusoidal function to habitat temperature data
fit <- nls(T_K ~ (meanT + dMean * days) + (amplT + dAmpl * days) * sin(2*pi*(days + shiftT)/360), data = data,
           start = list(meanT = 300, dMean = 0.1, amplT = 5, dAmpl = 0.1, shiftT = 30))
summary(fit)

# Plot model fit to TEMPERATURE data
xmin <- 0
xmax <- 360*(1991-1961)
ggplot(data, aes(x=days, y=T_K)) + 
  geom_point(size=3, color="black") +
  geom_line(size=0.8) +
  #geom_function(fun = function(t){(coef(fit)[1] + coef(fit)[2] * t) + (coef(fit)[3] + coef(fit)[4] * t) * sin(2*pi*(t + coef(fit)[5])/360)},
  #              size=1, color="#d1495b") +
  labs(x="Time (days)", y="Mean Temperature (K)") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(coef(fit)[1] - coef(fit)[2] - 6, coef(fit)[1] + coef(fit)[2] + 6)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20))

# Plot fit to TEMPERATURE data pooled across years
ggplot(data, aes(x=month, y=T_K)) + 
  geom_point(size=3, color="black") +
  scale_x_continuous(limits=c(0, 12)) +
  scale_y_continuous(limits=c(coef(fit)[1] - coef(fit)[2] - 6, coef(fit)[1] + coef(fit)[2] + 6)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20))


# Plot model fit to PRECIPITATION data
ymin <- 0
ymax <- 400
ggplot(data, aes(x=days, y=Prec_mm)) + 
  geom_point(size=3, color="blue") +
  geom_line(size=0.8, color="blue") +
  labs(x="Time (days)", y="Precipitation (mm)") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(ymin, ymax)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20))

# Plot fit to PRECIPITATION data pooled across years
ggplot(data, aes(x=month, y=Prec_mm)) + 
  geom_point(size=3, color="blue") +
  scale_x_continuous(limits=c(0, 12)) +
  scale_y_continuous(limits=c(ymin, ymax)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20))


########## Fit temperature functions to data from Habitat temperatures.csv ##########
# Read data
temp.data <- as.data.frame(read_csv("Habitat temperatures.csv"))

# Select an insect by removing # in front of name and placing # in front of other species
#sp.data <- subset(temp.data, Species == "Clavigralla shadabi")
#sp.data <- subset(temp.data, Species == "Clavigralla tomentosicollis Benin")
#sp.data <- subset(temp.data, Species == "Clavigralla tomentosicollis Nigeria A")
#sp.data <- subset(temp.data, Species == "Clavigralla tomentosicollis Nigeria B")
#sp.data <- subset(temp.data, Species == "Clavigralla tomentosicollis Nigeria C")
sp.data <- subset(temp.data, Species == "Clavigralla tomentosicollis Burkina Faso")

# Remove columns that do not contain temperature data
sp.data <- sp.data[-1:-7]
rownames(sp.data) <- c("T_K")

# Make data table with time and mean monthly temperature
data <- matrix(nrow = length(sp.data), ncol = 1)
for(i in 1:length(data)) {data[i] <- 30*(i-1)+15} # mid-point of month
colnames(data) <- c("day")
data <- cbind(data, t(sp.data)) # mean monthly temperature data

# Remove blank cells
data <- data[rowSums(is.na(data)) == 0,]
data <- as.data.frame(data)

# Fit sinusoidal function to habitat temperature data
fit <- nls(T_K ~ meanT + amplT*sin(2*pi*(day + shiftT)/365), data = data,
           start = list(meanT = 293, amplT = 2, shiftT = 30))
summary(fit)

# Plot model data
ggplot(data, aes(x=day, y=T_K)) + 
  geom_point(size=5, color="#d1495b") +
  geom_function(fun = function(t) coef(fit)[1] + coef(fit)[2]*sin(2*pi*(t + coef(fit)[3])/365),
                size=0.8, linetype="longdash", color="#d1495b") +
  labs(x="Time (days)", y="Mean Temperature (K)") +
  scale_x_continuous(limits=c(0, 720)) +
  scale_y_continuous(limits=c(coef(fit)[1] - coef(fit)[2] - 2, coef(fit)[1] + coef(fit)[2] + 3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20))

