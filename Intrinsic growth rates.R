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


################################## HISTORICAL CLIMATE DATA ###################################
# Read in climate data and temperature response parameters for selected insect
temp.h <- as.data.frame(read_csv(paste0("Historical climate data ",location,".csv")))
data <- subset(as.data.frame(read_csv("Temperature response parameters.csv")),
                  Species == paste(species,location))

# Integrate across r(T(t))
r.h <- 0
n.h <- nrow(temp.h)
for(i in 1:n) {
  r.h <- r.h + ifelse(temp.h$T[i] <= data$rTopt, data$rMax*exp(-1*((temp.h$T[i]-data$rTopt)/(2*data$rs))^2),
         data$rMax*(1 - ((temp.h$T[i]-data$rTopt)/(data$rTopt-data$rTmax))^2)) # from Deutsch et al. 2008
}
r.h <- r.h/n.h
r.h


#################################### FUTURE CLIMATE DATA ####################################
# Read in climate data and temperature response parameters for selected insect
temp.f <- as.data.frame(read_csv(paste0("Future climate data ",location,".csv")))
data <- subset(as.data.frame(read_csv("Temperature response parameters.csv")),
               Species == paste(species,location))

# Integrate across r(T(t))
r.f <- 0
n.f <- nrow(temp.f)
for(i in 1:n) {
  r.f <- r.f + ifelse(temp.f$T[i] <= data$rTopt, data$rMax*exp(-1*((temp.f$T[i]-data$rTopt)/(2*data$rs))^2),
                  data$rMax*(1 - ((temp.f$T[i]-data$rTopt)/(data$rTopt-data$rTmax))^2)) # from Deutsch et al. 2008
}
r.f <- r.f/n.f
r.f


