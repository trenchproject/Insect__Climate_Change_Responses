###################################################################################
#### This R script calculates a species' temperature sensitivity to extinction ####
###################################################################################

# Load packages and set working directory
library(tidyr)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(cubature)
library(lamW)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# USER: include overwintering? (i.e., do not integrate over temperatures below Tmin)
overw <- TRUE

# USER: include diurnal variation?
daily <- FALSE

# USER: run TPC analysis?
TPC <- FALSE


# Read in temperature response and temperature parameters, and temperature response data
param.all <- as.data.frame(read_csv("Temperature response parameters.csv"))
# Read in temperature parameters
ifelse(daily == TRUE, t.param.all <- as.data.frame(read_csv("Temperature parameters.csv")),
       t.param.all <- as.data.frame(read_csv("Temperature parameters Tave.csv")))
# Create array for results
results <- data.frame(param.all[,1], 1:nrow(param.all), 1:nrow(param.all))
names(results) <- c("Species","TPC","Model")


# RUN ANALYSIS FOR ALL SPECIES
#for(s in 1:nrow(param.all)) {
for(s in 1:1) {

# Select species
param <- param.all[s,]
t.param <- t.param.all[s,]


############################################ TPC ############################################
if(TPC == TRUE){
# Increase temperature mean and/or amplitude
delta.mean <- 0
delta.ampl <- 0

# Repeat integration of TPC as mean and/or amplitude increase until r < 0
repeat{

# Integrate across r(T(t))
T <- function(t) { (t.param$meanT.h+delta.mean) - (t.param$amplT.h+delta.ampl)*cos(2*pi*(t + t.param$shiftT.h)/365) - t.param$amplD.h*cos(2*pi*t) }
start <- 0
end <- 10*365

if(overw == FALSE) {
  r <- function(t) {
    ifelse(T(t) <= param$rTopt, param$rMax*exp(-1*((T(t)-param$rTopt)/(2*param$rs))^2),
           param$rMax*(1 - ((T(t)-param$rTopt)/(param$rTopt-param$rTmax))^2)) # from Deutsch et al. 2008
  }
  r.TPC <- cubintegrate(r, lower = start, upper = end, method = "pcubature")$integral/(end-start) # pcubature is faster but cannot be used with overwintering
}
if(overw == TRUE) {
  # r during active season
  r <- function(t) {
    ifelse(T(t) <= param$Tmin, 0,
           ifelse(T(t) <= param$rTopt, param$rMax*exp(-1*((T(t)-param$rTopt)/(2*param$rs))^2),
                  param$rMax*(1 - ((T(t)-param$rTopt)/(param$rTopt-param$rTmax))^2))) # from Deutsch et al. 2008
  }
  # integrate across active season
  season <- end # season length
  ifelse(daily == TRUE, length <- 0.5, length <- 1)
  for(t in seq(0,end,length)) { if(T(t) <= param$Tmin) {season <- season - length }} # number of days when T(t) > Tmin
  r.TPC <- cubintegrate(r, lower = start, upper = end, method = "hcubature")$integral/season
}
#print(r.TPC)

# Evaluate whether r < 0 and if so, break repeat
#if(delta.mean > 1){
if(r.TPC <= 0){
  results[s,2] <- delta.mean
  break
}

# Increase temperature mean and/or amplitude
delta.mean <- delta.mean + 0.1
}}


########################################### MODEL ############################################
TS <- as.data.frame(read_csv(paste0("Time series data Ext/Time series ",param[1],".csv")))

# Time of extinction
ext.time <- TS[(TS$J == 0 & TS$A == 0), "Time"][1]
print(ext.time)


}


