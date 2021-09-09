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


# EXTINCTION THRESHOLD (increase in meanT that leads to extinction in the model)
# Versus intrinsic growth rate (r)
# all parameters: rTopt significant only in full model
r <- lm(ext_meanT ~ rTopt + Toptr + rTmax, data=data)
summary (r)
r.Topt <- lm(ext_meanT ~ rTopt, data=data)
summary (r.Topt)

# thermal safety margin: non significant
tsm <- lm(ext_meanT ~ TSM, data=data)
summary (tsm)

# warming tolerance: significant
wt <- lm(ext_meanT ~ WT, data=data)
wt <- lm(ext_meanT ~ 0 + WT, data=data) # remove intercept
summary (wt)

# warming tolerance (active period): significant
wt.active <- lm(ext_meanT ~ WT_active, data=data)
wt.active <- lm(ext_meanT ~ 0 + WT_active, data=data) # remove intercept
summary (wt.active)


# Versus net reproductive rate (R0)
# all parameters: sR0 significant
R0 <- lm(ext_meanT ~ R0Topt + ToptR0 + sR0, data=data)
summary (R0)

# R0: Tmax (temperature at which R0 = 1): marginally significant
R0.Tmax <- lm(ext_meanT ~ R0Tmax, data=data)
summary (R0.Tmax)

# R0 thermal safety margin: significant
R0.tsm <- lm(ext_meanT ~ R0TSM, data=data)
R0.tsm <- lm(ext_meanT ~ 0 + R0TSM, data=data) # remove intercept
summary (R0.tsm)

# R0 warming tolerance: significant
R0.wt <- lm(ext_meanT ~ R0WT, data=data)
R0.wt <- lm(ext_meanT ~ 0 + R0WT, data=data) # remove intercept
summary (R0.wt)

# warming tolerance (active period): significant
R0.wt.active <- lm(ext_meanT ~ R0WT_active, data=data)
R0.wt.active <- lm(ext_meanT ~ 0 + R0WT_active, data=data) # remove intercept
summary (R0.wt.active)



# PLOTS
# Intrinsic growth rate (r)
# thermal safety margin: non significant
Xmin <- 0
Xmax <- 12
Ymin <- 0
Ymax <- 15
#plot(data$TSM, data$ext_meanT, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), main=expression("extinction threshold (meanT) vs. TSM"))
#points(seq(Xmin,Xmax,1), coef(tsm)[1]*seq(Xmin,Xmax,1), type="l", col="blue")
#abline(0, 1)

# warming tolerance: significant (slope = 0.81)
Xmin <- 0
Xmax <- 25
Ymin <- 0
Ymax <- 15
plot(data$WT, data$ext_meanT, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), main=expression("extinction threshold (meanT) vs. WT"))
points(seq(Xmin,Xmax,1), coef(wt)[1]*seq(Xmin,Xmax,1), type="l", col="blue")
abline(0, 1)

# warming tolerance (active period): significant (slope = 1.01)
Xmin <- 0
Xmax <- 20
Ymin <- 0
Ymax <- 15
plot(data$WT_active, data$ext_meanT, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), main=expression("extinction threshold (meanT) vs. WT (active)"))
points(seq(Xmin,Xmax,1), coef(wt.active)[1]*seq(Xmin,Xmax,1), type="l", col="blue")
abline(0, 1)



# Net reproductive rate (R0)
# R0 thermal safety margin: significant (slope = 0.46)
Xmin <- 10
Xmax <- 30
Ymin <- 0
Ymax <- 15
plot(data$R0TSM, data$ext_meanT, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), main=expression("extinction threshold (meanT) vs. R0 TSM"))
points(seq(Xmin,Xmax,1), coef(R0.tsm)[1]*seq(Xmin,Xmax,1), type="l", col="blue")
abline(0, 1)

# warming tolerance: significant (slope = 0.38)
Xmin <- 0
Xmax <- 40
Ymin <- 0
Ymax <- 15
plot(data$R0WT, data$ext_meanT, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), main=expression("extinction threshold (meanT) vs. R0 WT"))
points(seq(Xmin,Xmax,1), coef(R0.wt)[1]*seq(Xmin,Xmax,1), type="l", col="blue")
abline(0, 1)

# warming tolerance (active period): significant (slope = 0.5)
Xmin <- 0
Xmax <- 30
Ymin <- 0
Ymax <- 15
plot(data$R0WT_active, data$ext_meanT, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), main=expression("extinction threshold (meanT) vs. R0 WT (active)"))
points(seq(Xmin,Xmax,1), coef(R0.wt.active)[1]*seq(Xmin,Xmax,1), type="l", col="blue")
abline(0, 1)




