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


# USER: enter species and location
species <- "Clavigralla tomentosicollis"
location <- "Nigeria"
species <- "Uroleucon ambrosiae"
location <- "Brazil"

################################## TPC: HISTORICAL CLIMATE ###################################
# Read in climate data
temp.h <- as.data.frame(read_csv(paste0("Climate data/Historical climate data ",location,".csv")))
# Read in temperature response and temperature parameters, and temperature response data for selected insect
param <- subset(as.data.frame(read_csv("Temperature response parameters.csv")), Species == paste(species,location))
# Read in temperature parameters
t.param <- subset(as.data.frame(read_csv("Temperature parameters.csv")), Species == paste(species,location))

# # Brute force method: average r(T) for all climate data T
# r.TPC.h <- 0
# n.TPC.h <- nrow(temp.h)
# for(i in 1:n.TPC.h) {
#   r.TPC.h <- r.TPC.h + ifelse(temp.h$T[i] <= param$rTopt, param$rMax*exp(-1*((temp.h$T[i]-param$rTopt)/(2*param$rs))^2),
#          param$rMax*(1 - ((temp.h$T[i]-param$rTopt)/(param$rTopt-param$rTmax))^2)) # from Deutsch et al. 2008
# }
# (r.TPC.h <- r.TPC.h/n.TPC.h)

# Integrate across r(T(t))
T.h <- function(t) { (t.param$meanT.h+t.param$delta_mean.h*t) - (t.param$amplT.h+t.param$delta_ampl.h*t)*cos(2*pi*(t + t.param$shiftT.h)/365) - 0*t.param$amplD.h*cos(2*pi*t) }
r.h <- function(t) {
  ifelse(T.h(t) <= param$rTopt, param$rMax*exp(-1*((T.h(t)-param$rTopt)/(2*param$rs))^2),
                                 param$rMax*(1 - ((T.h(t)-param$rTopt)/(param$rTopt-param$rTmax))^2)) # from Deutsch et al. 2008
}
start <- 0
end <- 365*80
(r.TPC.h <- cubintegrate(r.h, lower = start, upper = end, method = "pcubature")$integral/(end-start))


##################################### TPC: FUTURE CLIMATE ####################################
# Read in climate data
temp.f <- as.data.frame(read_csv(paste0("Climate data/Future climate data ",location,".csv")))
# Read in temperature response and temperature parameters, and temperature response data for selected insect
param <- subset(as.data.frame(read_csv("Temperature response parameters.csv")), Species == paste(species,location))
# Read in temperature parameters
t.param <- subset(as.data.frame(read_csv("Temperature parameters.csv")), Species == paste(species,location))

# # Brute force method: average r(T) for all climate data T
# r.TPC.f <- 0
# n.TPC.f <- nrow(temp.f)
# for(i in 1:n.TPC.f) {
#   r.TPC.f <- r.TPC.f + ifelse(temp.f$T[i] <= param$rTopt, param$rMax*exp(-1*((temp.f$T[i]-param$rTopt)/(2*param$rs))^2),
#                   param$rMax*(1 - ((temp.f$T[i]-param$rTopt)/(param$rTopt-param$rTmax))^2)) # from Deutsch et al. 2008
# }
# r.TPC.f <- r.TPC.f/n.TPC.f
# r.TPC.f

# Integrate across r(T(t))
T.f <- function(t) { (t.param$meanT.f+t.param$delta_mean.f*t) - (t.param$amplT.f+t.param$delta_ampl.f*t)*cos(2*pi*(t + t.param$shiftT.f)/365) - 0*t.param$amplD.f*cos(2*pi*t) }
r.f <- function(t) {
  ifelse(T.f(t) <= param$rTopt, param$rMax*exp(-1*((T.f(t)-param$rTopt)/(2*param$rs))^2),
         param$rMax*(1 - ((T.f(t)-param$rTopt)/(param$rTopt-param$rTmax))^2)) # from Deutsch et al. 2008
}
start <- 365*70 # start 2090
end <- 365*80 # end 2100
(r.TPC.f <- cubintegrate(r.f, lower = start, upper = end, method = "pcubature")$integral/(end-start))


# PLOT
Tmin <- round(min(temp.h$T,temp.f$T),0) - 1
Tmax <- round(max(temp.h$T,temp.f$T),0) + 1
ymin <- 0
ymax <- round(param$rMax,1) + 0.05
hist(temp.h$T, xlim=c(Tmin,Tmax), ylim=c(ymin,ymax), breaks=seq(from=Tmin, to=Tmax, by=1), ylab="r", col=rgb(0,0,255, max = 255, alpha = 80), border=rgb(0,0,255, max = 255, alpha = 80), freq=FALSE, main = NULL)
hist(temp.f$T, xlim=c(Tmin,Tmax), ylim=c(ymin,ymax), breaks=seq(from=Tmin, to=Tmax, by=1), ylab="r", col=rgb(255,0,0, max = 255, alpha = 80), border=rgb(255,0,0, max = 255, alpha = 80), freq=FALSE, main = NULL, add=TRUE)
points(seq(Tmin,Tmax,1), ifelse(seq(Tmin,Tmax,1) <= param$rTopt, param$rMax*exp(-1*((seq(Tmin,Tmax,1)-param$rTopt)/(2*param$rs))^2),
                                param$rMax*(1 - ((seq(Tmin,Tmax,1)-param$rTopt)/(param$rTopt-param$rTmax))^2)), type="l", lwd=4, col="black")
abline(v = mean(temp.h$T), col="blue", lwd=3, lty=2)
abline(v = mean(temp.f$T), col="red", lwd=3, lty=2)
r.TPC.h
r.TPC.f



################################# MODEL: HISTORICAL CLIMATE ##################################
# Read in climate data and temperature response parameters for selected insect
TS.h <- as.data.frame(read_csv(paste0("Time series data/Historical time series ",species," ",location,".csv")))

# Integrate model time-series across ln(t/(t-1))
# r.model.h <- 0
# count.h <- 0
# for(i in 3:nrow(TS.h)) {
#   if(TS.h$A[i] >= 0 && TS.h$A[i-1] >= 0 && is.na(TS.h$A[i]) == FALSE && is.na(TS.h$A[i-1]) == FALSE) {
#     r.model.h <- r.model.h + log(TS.h$A[i]/TS.h$A[i-1])
#     count.h <- count.h + 1
# }}
# (r.model.h <- r.model.h/count.h)

# Integrate low density per capita population growth rate across ln(t/(t-1))
r.model.h <- 0
start <- nrow(TS.h) - 365*10 # integrate over last 10 years of time-series
end <- nrow(TS.h)
for(i in start:end) { r.model.h <- r.model.h + TS.h$r[i] }
(r.model.h <- r.model.h/(end-start))


################################### MODEL: FUTURE CLIMATE ####################################
# Read in climate data and temperature response parameters for selected insect
TS.f <- as.data.frame(read_csv(paste0("Time series data/Future time series ",species," ",location,".csv")))

# Integrate across ln(t+1/t)
# r.model.f <- 0
# count.f <- 0
# for(i in 3:nrow(TS.f)) {
#   if(TS.f$A[i] >= 0 && TS.f$A[i-1] >= 0 && is.na(TS.f$A[i]) == FALSE && is.na(TS.f$A[i-1]) == FALSE) {
#     r.model.f <- r.model.f + log(TS.f$A[i]/TS.f$A[i-1])
#     count.f <- count.f + 1
#   }}
# (r.model.f <- r.model.f/count.f)

# Integrate low density per capita population growth rate across ln(t/(t-1))
r.model.f <- 0
start <- nrow(TS.f) - 365*10 # integrate over last 10 years of time-series
end <- nrow(TS.f)
for(i in start:end) { r.model.f <- r.model.f + TS.f$r[i] }
(r.model.f <- r.model.f/(end-start))


r.TPC.h
r.TPC.f
r.model.h
r.model.f
r.TPC.f/r.TPC.h
r.model.f/r.model.h


