##################################################################################
#### This R script calculates a species' intrinsic per capita growth rate (r) ####
##################################################################################

# Load packages and set working directory
library(tidyr)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(ncdf4)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# USER: enter location, time period, and insect species
location <- "Benin"
species <- "Clavigralla shadabi"


################################## TPC: HISTORICAL CLIMATE ###################################
# Read in climate data and temperature response parameters for selected insect
temp.h <- as.data.frame(read_csv(paste0("Historical climate data ",location,".csv")))
data <- subset(as.data.frame(read_csv("Temperature response parameters.csv")),
                  Species == paste(species,location))

# Integrate across r(T(t))
r.TPC.h <- 0
n.TPC.h <- nrow(temp.h)
for(i in 1:n) {
  r.TPC.h <- r.TPC.h + ifelse(temp.h$T[i] <= data$rTopt, data$rMax*exp(-1*((temp.h$T[i]-data$rTopt)/(2*data$rs))^2),
         data$rMax*(1 - ((temp.h$T[i]-data$rTopt)/(data$rTopt-data$rTmax))^2)) # from Deutsch et al. 2008
}
r.TPC.h <- r.TPC.h/n.TPC.h
r.TPC.h


##################################### TPC: FUTURE CLIMATE ####################################
# Read in climate data and temperature response parameters for selected insect
temp.f <- as.data.frame(read_csv(paste0("Future climate data ",location,".csv")))
data <- subset(as.data.frame(read_csv("Temperature response parameters.csv")),
               Species == paste(species,location))

# Integrate across r(T(t))
r.TPC.f <- 0
n.TPC.f <- nrow(temp.f)
for(i in 1:n) {
  r.TPC.f <- r.TPC.f + ifelse(temp.f$T[i] <= data$rTopt, data$rMax*exp(-1*((temp.f$T[i]-data$rTopt)/(2*data$rs))^2),
                  data$rMax*(1 - ((temp.f$T[i]-data$rTopt)/(data$rTopt-data$rTmax))^2)) # from Deutsch et al. 2008
}
r.TPC.f <- r.TPC.f/n.TPC.f
r.TPC.f


################################# MODEL: HISTORICAL CLIMATE ##################################
# Read in climate data and temperature response parameters for selected insect
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
# Read in climate data and temperature response parameters for selected insect
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

