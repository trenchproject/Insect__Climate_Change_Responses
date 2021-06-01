########################################################################
#### This R script fits model outputs to empirical time-series data ####
########################################################################


# Load packages and set working directory
library(tidyr)
library(ggplot2)
library(dplyr)
library(tidyverse)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Select an insect by removing # in front of name and placing # in front of other species
#read.data <- read_csv("Time Series Clavigralla shadabi.csv")
#read.data <- read_csv("Time Series Clavigralla tomentosicollis Benin.csv")
read.data <- read_csv("Time Series Clavigralla tomentosicollis Nigeria A.csv")
data <- as.data.frame(read.data)


# Plot model data
ggplot(data, aes(x=Time, y=J)) + 
  geom_point(size=5, color="#d1495b") +
  geom_function(fun = function(t) coef(fit)[1] + coef(fit)[2]*sin((2*pi*t + coef(fit)[3])/365),
                size=0.8, linetype="longdash", color="#d1495b") +
  labs(x="Time", y="Density") +
  scale_x_continuous(limits=c(0, 720)) +
  scale_y_log10(limits=c(0.1, 10)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20))

# Plot time-series data
plotA.J = ggplot(data.plotA, aes(x=time, y=J, ymin=J-J_SE, ymax=J+J_SE)) + 
  geom_pointrange(size=1.2, color="#d1495b") +
  #geom_line(size=0.8, linetype="longdash", color="#d1495b") +
  labs(x="Time", y="Density") +
  scale_x_continuous(limits=c(200, 500)) +
  scale_y_log10(limits=c(0.3, 10)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20))
