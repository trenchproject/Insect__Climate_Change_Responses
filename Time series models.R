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


# Read in time-series data and temperature response data
data.density <- read_csv("Egwuatu_1977.csv")
data <- as.data.frame(read_csv("Temperature response parameters.csv"))


# Select time-series data (for Nigeria data, select plot A, B, or C)
data.TS <- subset(data.density, Plot=="C")


# Read in model output
# Select an insect by removing # in front of name and placing # in front of other species
#read.data <- as.data.frame(read_csv("Time Series Clavigralla shadabi.csv"))
#sp.data <- subset(sp.data, Species == "Clavigralla shadabi")
#read.data <- as.data.frame(read_csv("Time Series Clavigralla tomentosicollis Benin.csv"))
#sp.data <- subset(sp.data, Species == "Clavigralla tomentosicollis Benin")
#data.model <- as.data.frame(read_csv("Time Series Clavigralla tomentosicollis Nigeria A.csv"))
#sp.data <- subset(data, Species == "Clavigralla tomentosicollis Nigeria A")
data.model <- as.data.frame(read_csv("Time Series Clavigralla tomentosicollis Burkina Faso C.csv"))
sp.data <- subset(data, Species == "Clavigralla tomentosicollis Burkina Faso C")


# Set plot options (default: plot last 2 year of model)
xmin <- 200
xmax <- 750
ymin <- 0
ymax <- 10
TS.length <- xmax - xmin # length of time-series data
end <- nrow(data.model)

# Format model output to align with time-series data
# Remove all rows before last 10 years + xmin days
data.model <- data.model[c(-1:-(end - 10*365 + xmin)), ]

# Remove all rows after xmax days
data.model <- data.model[c(-(xmax - xmin + 2):-end), ]

# Re-scale time so that last x years are plotted as the first x years
data.model <- sweep(data.model, 2, c(end - 10*365,0,0,0,0))


# Plot time-series data
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
#plot.J

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
#plot.A

# Plot model data
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
#model.J

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
#model.A

# Plot habitat temperature function
plot.temp <- ggplot() +
  geom_function(fun = function(t) sp.data$meanT + sp.data$amplT*sin(2*pi*(t + sp.data$shiftT)/365),
                size=0.8, linetype="longdash", color="#d1495b") +
  labs(x="Time", y="T (K)") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(sp.data$meanT - sp.data$amplT - 1, sp.data$meanT + sp.data$amplT + 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20))
#plot.temp


# Draw final plot
ggdraw()  +
  draw_plot(plot.temp, x = 0, y = 0, width = 1, height = 0.3) +
  draw_plot(plot.J, x = 0, y = 0.3, width = 1, height = 0.7) +
  draw_plot(plot.A, x = 0, y = 0.3, width = 1, height = 0.7) +
  draw_plot(model.J, x = 0, y = 0.3, width = 1, height = 0.7) +
  draw_plot(model.A, x = 0, y = 0.3, width = 1, height = 0.7)
