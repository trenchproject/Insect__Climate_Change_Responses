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


# Read in time-series data
data.density <- read_csv("Egwuatu_1977.csv")
# Set the plot (A, B, or C)
data.TS <- subset(data.density, Plot=="A")


# Select an insect by removing # in front of name and placing # in front of other species
#read.data <- read_csv("Time Series Clavigralla shadabi.csv")
#read.data <- read_csv("Time Series Clavigralla tomentosicollis Benin.csv")
read.data <- read_csv("Time Series Clavigralla tomentosicollis Nigeria A.csv")
data <- as.data.frame(read.data)


# Set plot options (default: plot last 2 year of model)
#xmin <- nrow(data) - 2*365
#xmax <- nrow(data)
xmin <- 0
xmax <-730
ymin <- -50
ymax <- 50


# Plot time-series data
plot.J = ggplot(data.TS, aes(x=time, y=J, ymin=J-J_SE, ymax=J+J_SE)) + 
  geom_pointrange(size=0.8, color="#d1495b") +
  #geom_line(size=0.8, linetype="longdash", color="#d1495b") +
  labs(x="Time", y="Density") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(ymin, ymax)) +
  #scale_y_log10(limits=c(0.3, 10)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20))
#plot.J

plot.A = ggplot(data.TS, aes(x=time, y=A, ymin=A-A_SE, ymax=A+A_SE)) + 
  geom_pointrange(size=0.8, color="#30638e") +
  #geom_line(size=0.8, linetype="longdash", color="#30638e") +
  labs(x="Time", y="Density") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(ymin, ymax)) +
  #scale_y_log10(limits=c(0.3, 10)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20)) 
#plot.A

# Plot model data
model.J = ggplot(data, aes(x=Time, y=J)) + 
  geom_line(size=0.8, color="#d1495b") +
  labs(x="Time", y="Density") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(ymin, ymax)) +
  #scale_y_log10(limits=c(0.3, 10)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20))
#model.J

model.A = ggplot(data, aes(x=Time, y=A)) + 
  geom_line(size=0.8, color="#30638e") +
  labs(x="Time", y="Density") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(ymin, ymax)) +
  #scale_y_log10(limits=c(0.3, 10)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20)) 
#model.A

# Draw final plot
ggdraw()  +
  #draw_plot(plot.temp, x = 0, y = 0, width = 1, height = 0.3) +
  draw_plot(plot.J, x = 0, y = 0.3, width = 1, height = 0.7) +
  draw_plot(plot.A, x = 0, y = 0.3, width = 1, height = 0.7) +
  draw_plot(model.J, x = 0, y = 0.3, width = 1, height = 0.7) +
  draw_plot(model.A, x = 0, y = 0.3, width = 1, height = 0.7)