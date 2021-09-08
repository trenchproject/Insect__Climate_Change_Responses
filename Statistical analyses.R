#############################################################################
#### This R script analyzes the thermal performance curve and model data ####
#############################################################################


# Load packages and set working directory
library(tidyr)
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyverse)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# READ IN TEMPERATURE RESPONSE PARAMETERS
# Temperature response parameters
data <- as.data.frame(read_csv("Temperature response parameters.csv"))


# EXTINCTION THRESHOLD (INCREASE IN MEAN T) VERSUS THERMAL SAFETY MARGIN AND WARMING TOLERANCE
# Linear regressions
# thermal safety margin
tsm <- lm(ext_meanT ~ 0 + TSM, data=data)
summary (tsm)
# warming tolerance
wt <- lm(ext_meanT ~ 0 + WT, data=data)
summary (wt)

# Plots
# thermal safety margin
Xmin <- 0
Xmax <- 12
Ymin <- 0
Ymax <- 15
plot(data$TSM, data$ext_meanT, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), main=expression("extinction threshold (meanT) vs. TSM"))
points(seq(Xmin,Xmax,1), coef(tsm)[1]*seq(Xmin,Xmax,1), type="l", col="blue")
# warming tolerance
Xmin <- 0
Xmax <- 25
Ymin <- 0
Ymax <- 15
plot(data$WT, data$ext_meanT, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), main=expression("extinction threshold (meanT) vs. WT"))
points(seq(Xmin,Xmax,1), coef(wt)[1]*seq(Xmin,Xmax,1), type="l", col="blue")


