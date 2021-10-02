##################################################################################
#### This R script calculates a species' intrinsic per capita growth rate (r) ####
##################################################################################

# Load packages and set working directory
library(tidyr)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(cubature)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# USER: enter location, time period, and insect species
location <- "Benin"
species <- "Clavigralla shadabi"


################################## TPC: HISTORICAL CLIMATE ###################################
# Read in climate data
temp.h <- as.data.frame(read_csv(paste0("Historical climate data ",location,".csv")))
# Read in temperature response and temperature parameters, and temperature response data for selected insect
param <- subset(as.data.frame(read_csv("Temperature response parameters.csv")), Species == paste(species,location))
# Read in temperature parameters
t.param <- subset(as.data.frame(read_csv("Temperature parameters.csv")), Species == paste(species,location))
# Read in observed r data
r.data <- subset(as.data.frame(read_csv("Temperature response data.csv")), Species == species)[c(10,47)] # col 10 = T_K, col 47 = r

# Brute force method: average r(T) for all climate data T
r.TPC.h <- 0
n.TPC.h <- nrow(temp.h)
for(i in 1:n.TPC.h) {
  r.TPC.h <- r.TPC.h + ifelse(temp.h$T[i] <= param$rTopt, param$rMax*exp(-1*((temp.h$T[i]-param$rTopt)/(2*param$rs))^2),
         param$rMax*(1 - ((temp.h$T[i]-param$rTopt)/(param$rTopt-param$rTmax))^2)) # from Deutsch et al. 2008
}
(r.TPC.h <- r.TPC.h/n.TPC.h)

# Integrate across r(T(t))
T.h <- function(t) { (t.param$meanT.h+t.param$delta_mean.h*t) + (t.param$amplT.h+t.param$delta_ampl.h*t)*sin(2*pi*(t + t.param$shiftT.h)/365) - t.param$amplD.h*cos(2*pi*t) }
r.h <- function(t) {
  ifelse(T.h(t) <= param$rTopt, param$rMax*exp(-1*((T.h(t)-param$rTopt)/(2*param$rs))^2),
                                 param$rMax*(1 - ((T.h(t)-param$rTopt)/(param$rTopt-param$rTmax))^2)) # from Deutsch et al. 2008
}
n.TPC.h <- 3650 # integrate over 10 years
(r.TPC.h <- cubintegrate(r.h, lower = 0, upper = n.TPC.h, method = "pcubature")$integral/n.TPC.h)


# PLOT
Tmin <- 285
Tmax <- 315
ymin <- 0
ymax <- 0.2
hist(temp.h$T, xlim=c(Tmin,Tmax), ylim=c(ymin,ymax), breaks=seq(from=Tmin, to=Tmax, by=1), ylab="r", col="blue", freq=FALSE, main = NULL)
hist(temp.f$T, xlim=c(Tmin,Tmax), ylim=c(ymin,ymax), breaks=seq(from=Tmin, to=Tmax, by=1), ylab="r", col=rgb(255,0,0, max = 255, alpha = 130), freq=FALSE, main = NULL, add=TRUE)
points(seq(Tmin,Tmax,1), ifelse(seq(Tmin,Tmax,1) <= param$rTopt, param$rMax*exp(-1*((seq(Tmin,Tmax,1)-param$rTopt)/(2*param$rs))^2),
                                                           param$rMax*(1 - ((seq(Tmin,Tmax,1)-param$rTopt)/(param$rTopt-param$rTmax))^2)), type="l", lwd=5, col="black")


##################################### TPC: FUTURE CLIMATE ####################################
# Read in climate data
temp.h <- as.data.frame(read_csv(paste0("Future climate data ",location,".csv")))
# Read in temperature response and temperature parameters, and temperature response data for selected insect
param <- subset(as.data.frame(read_csv("Temperature response parameters.csv")), Species == paste(species,location))
# Read in temperature parameters
t.param <- subset(as.data.frame(read_csv("Temperature parameters.csv")), Species == paste(species,location))
# Read in observed r data
r.data <- subset(as.data.frame(read_csv("Temperature response data.csv")), Species == species)[c(10,47)] # col 10 = T_K, col 47 = r


# Brute force method: average r(T) for all climate data T
r.TPC.f <- 0
n.TPC.f <- nrow(temp.f)
for(i in 1:n.TPC.f) {
  r.TPC.f <- r.TPC.f + ifelse(temp.f$T[i] <= param$rTopt, param$rMax*exp(-1*((temp.f$T[i]-param$rTopt)/(2*param$rs))^2),
                  param$rMax*(1 - ((temp.f$T[i]-param$rTopt)/(param$rTopt-param$rTmax))^2)) # from Deutsch et al. 2008
}
r.TPC.f <- r.TPC.f/n.TPC.f
r.TPC.f

# Integrate across r(T(t))
T.f <- function(t) { (t.param$meanT.f+t.param$delta_mean.f*t) + (t.param$amplT.f+t.param$delta_ampl.f*t)*sin(2*pi*(t + t.param$shiftT.f)/365) - t.param$amplD.f*cos(2*pi*t) }
r.f <- function(t) {
  ifelse(T.f(t) <= param$rTopt, param$rMax*exp(-1*((T.f(t)-param$rTopt)/(2*param$rs))^2),
         param$rMax*(1 - ((T.f(t)-param$rTopt)/(param$rTopt-param$rTmax))^2)) # from Deutsch et al. 2008
}
n.TPC.f <- 365*80 # integrate over 80 years
(r.TPC.f <- cubintegrate(r.f, lower = 0, upper = n.TPC.f, method = "pcubature")$integral/n.TPC.f)


################################# MODEL: HISTORICAL CLIMATE ##################################
# Read in climate param and temperature response parameters for selected insect
TS.h <- as.data.frame(read_csv(paste0("Historical time series ",species," ",location,".csv")))

# Integrate across ln(t+1)/ln(t)
start <- 10*365 # skip first 10 years, which are used for model initialization
r.model.h <- 0
n.model.h <- nrow(TS.h)
for(i in start:start+365) {
  r.model.h <- r.model.h + log(TS.h$A[i+1]+1)/log(TS.h$A[i]+1)
}
r.model.h <- r.model.h/365 #n.model.h
r.model.h


################################### MODEL: FUTURE CLIMATE ####################################
# Read in climate param and temperature response parameters for selected insect
TS.f <- as.data.frame(read_csv(paste0("Future time series ",species," ",location,".csv")))

# Integrate across ln(t+1)/ln(t)
start <- 10*365 # skip first 10 years, which are used for model initialization
r.model.f <- 0
n.model.f <- nrow(TS.f)
for(i in start:start+365) {
  r.model.f <- r.model.f + log(TS.f$A[i+1]+1)/log(TS.f$A[i]+1)
}
r.model.f <- r.model.f/365 #n.model.f
r.model.f

r.TPC.h
r.TPC.f
r.model.h
r.model.f

