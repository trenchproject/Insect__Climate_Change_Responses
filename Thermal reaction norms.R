##############################################################
#### This R script plots species' thermal reaction norms #####
##############################################################


# Load packages and set working directory
library(tidyr)
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyverse)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# READ IN DATA
LH.data <- as.data.frame(read_csv("Temperature response parameters.csv"))
temp.data <- as.data.frame(read_csv("Temperature parameters Tave.csv"))

# EXCLUDE DATA
LH.data <- LH.data[-c(12,13,15,18),]
temp.data <- temp.data[-c(12,13,15,18),]

# SEPARATE BY HABITAT
trop <- LH.data[LH.data$Habitat == "Tropical",]
subtrop <- LH.data[LH.data$Habitat == "Subtropical",]
temp <- LH.data[LH.data$Habitat == "Temperate",]
trop.T <- temp.data[temp.data$Habitat == "Tropical",]
subtrop.T <- temp.data[temp.data$Habitat == "Subtropical",]
temp.T <- temp.data[temp.data$Habitat == "Temperate",]


########################################### PLOTS ############################################
# FITNESS
# Plot options
Tmin <- 275
Tmax <- 315
ymin <- 0
ymax <- 0.5
# Thermal reaction norm
r <- function(T) { ifelse(T <= param$rTopt, param$rMax*exp(-1*((T-param$rTopt)/(2*param$rs))^2),
                            param$rMax*(1 - ((T-param$rTopt)/(param$rTopt-param$rTmax))^2)) }
# Plot
param <- trop[1,]
plot(seq(Tmin,Tmax,0.1), r(seq(Tmin,Tmax,0.1)), type="l", lwd=1, col="red", xlim=c(Tmin,Tmax), ylim=c(ymin,ymax), xlab="T", ylab="m(T)")
for(i in 2:nrow(trop)) {
  param <- trop[i,]
  points(seq(Tmin,Tmax,0.1), r(seq(Tmin,Tmax,0.1)), type="l", lwd=1, col="red", xlim=c(Tmin,Tmax), ylim=c(ymin,ymax), xlab="T", ylab="m(T)")
}
for(i in 1:nrow(subtrop)) {
  param <- subtrop[i,]
  points(seq(Tmin,Tmax,0.1), r(seq(Tmin,Tmax,0.1)), type="l", lwd=1, col="orange", xlim=c(Tmin,Tmax), ylim=c(ymin,ymax), xlab="T", ylab="m(T)")
}
for(i in 1:nrow(temp)) {
  param <- temp[i,]
  points(seq(Tmin,Tmax,0.1), r(seq(Tmin,Tmax,0.1)), type="l", lwd=1, col="blue", xlim=c(Tmin,Tmax), ylim=c(ymin,ymax), xlab="T", ylab="m(T)")
}
# Habitat temperatures
for(i in 1:nrow(trop)) {
  t.param <- trop.T[i,]
  abline(v = t.param$meanT.f + t.param$delta_mean.f*365*75, col="red", lwd=1, lty="longdash")
}
for(i in 1:nrow(subtrop)) {
  t.param <- subtrop.T[i,]
  abline(v = t.param$meanT.f + t.param$delta_mean.f*365*75, col="orange", lwd=1, lty="longdash")
}
for(i in 1:nrow(temp)) {
  t.param <- temp.T[i,]
  abline(v = t.param$meanT.f + t.param$delta_mean.f*365*75, col="blue", lwd=1, lty="longdash")
}


# DEVELOPMENT RATE (LEFT-SKEWED)
# Plot options
Tmin <- 275
Tmax <- 315
ymin <- 0
ymax <- 0.5
# Thermal reaction norm
m <- function (T) { param$mTR*(T/param$TR)*exp(param$AmJ*(1/param$TR-1/T))/(1+exp(param$AL*(1/param$TL-1/T))+exp(param$AH*(1/param$TH-1/T))) }
# Plot
param <- trop[1,]
plot(seq(Tmin,Tmax,0.1), m(seq(Tmin,Tmax,0.1)), type="l", lwd=1, col="red", xlim=c(Tmin,Tmax), ylim=c(ymin,ymax), xlab="T", ylab="m(T)")
for(i in 2:nrow(trop)) {
  param <- trop[i,]
  points(seq(Tmin,Tmax,0.1), m(seq(Tmin,Tmax,0.1)), type="l", lwd=1, col="red", xlim=c(Tmin,Tmax), ylim=c(ymin,ymax), xlab="T", ylab="m(T)")
}
for(i in 1:nrow(subtrop)) {
  param <- subtrop[i,]
  points(seq(Tmin,Tmax,0.1), m(seq(Tmin,Tmax,0.1)), type="l", lwd=1, col="orange", xlim=c(Tmin,Tmax), ylim=c(ymin,ymax), xlab="T", ylab="m(T)")
}
for(i in 1:nrow(temp)) {
  param <- temp[i,]
  points(seq(Tmin,Tmax,0.1), m(seq(Tmin,Tmax,0.1)), type="l", lwd=1, col="blue", xlim=c(Tmin,Tmax), ylim=c(ymin,ymax), xlab="T", ylab="m(T)")
}
# Habitat temperatures
for(i in 1:nrow(trop)) {
  t.param <- trop.T[i,]
  abline(v = t.param$meanT.f + t.param$delta_mean.f*365*75, col="red", lwd=1, lty="longdash")
}
for(i in 1:nrow(subtrop)) {
  t.param <- subtrop.T[i,]
  abline(v = t.param$meanT.f + t.param$delta_mean.f*365*75, col="orange", lwd=1, lty="longdash")
}
for(i in 1:nrow(temp)) {
  t.param <- temp.T[i,]
  abline(v = t.param$meanT.f + t.param$delta_mean.f*365*75, col="blue", lwd=1, lty="longdash")
}


# DEVELOPMENT RATE (STEP FUNCTION)
# Plot options
Tmin <- 275
Tmax <- 315
ymin <- 0
ymax <- 0.5
# Thermal reaction norm
m <- function (T) { ifelse(T <= param$Topt, param$mTR*(T/param$TR)*exp(param$AmJ*(1/param$TR-1/T))/(1+exp(param$AH*(1/param$TH-1/T))),
                           ifelse(T <= param$Tmax, param$mTR*(param$Topt/param$TR)*exp(param$AmJ*(1/param$TR-1/param$Topt))/(1+exp(param$AH*(1/param$TH-1/param$Topt))), 0)) }
# Plot
param <- trop[1,]
plot(seq(Tmin,Tmax,0.1), m(seq(Tmin,Tmax,0.1)), type="l", lwd=1, col="red", xlim=c(Tmin,Tmax), ylim=c(ymin,ymax), xlab="T", ylab="m(T)")
for(i in 2:nrow(trop)) {
  param <- trop[i,]
  points(seq(Tmin,Tmax,0.1), m(seq(Tmin,Tmax,0.1)), type="l", lwd=1, col="red", xlim=c(Tmin,Tmax), ylim=c(ymin,ymax), xlab="T", ylab="m(T)")
}
for(i in 1:nrow(subtrop)) {
  param <- subtrop[i,]
  points(seq(Tmin,Tmax,0.1), m(seq(Tmin,Tmax,0.1)), type="l", lwd=1, col="orange", xlim=c(Tmin,Tmax), ylim=c(ymin,ymax), xlab="T", ylab="m(T)")
}
for(i in 1:nrow(temp)) {
  param <- temp[i,]
  points(seq(Tmin,Tmax,0.1), m(seq(Tmin,Tmax,0.1)), type="l", lwd=1, col="blue", xlim=c(Tmin,Tmax), ylim=c(ymin,ymax), xlab="T", ylab="m(T)")
}
# Habitat temperatures
for(i in 1:nrow(trop)) {
  t.param <- trop.T[i,]
  abline(v = t.param$meanT.f + t.param$delta_mean.f*365*75, col="red", lwd=1, lty="longdash")
}
for(i in 1:nrow(subtrop)) {
  t.param <- subtrop.T[i,]
  abline(v = t.param$meanT.f + t.param$delta_mean.f*365*75, col="orange", lwd=1, lty="longdash")
}
for(i in 1:nrow(temp)) {
  t.param <- temp.T[i,]
  abline(v = t.param$meanT.f + t.param$delta_mean.f*365*75, col="blue", lwd=1, lty="longdash")
}


# THERMAL REACTION NORMS
# r <- function (T) { ifelse(T <= param$rTopt, param$rMax*exp(-1*((T-param$rTopt)/(2*param$rs))^2),
#                            param$rMax*(1 - ((T-param$rTopt)/(param$rTopt-param$rTmax))^2)) }
# b <- function (T) { param$bTopt*exp(-((T-param$Toptb)^2)/(2*param$sb^2)) }
# dJ <- function (T) { param$dJTR*exp(param$AdJ*(1/param$TR-1/T)) }
# dA <- function (T) { param$dATR*exp(param$AdA*(1/param$TR-1/T)) }
#R0 <- function (T) { b/dA * m/(m+dJ) }
#f <- function (T) { b/dA }
#s <- function (T) { exp(-dJ/m) }

