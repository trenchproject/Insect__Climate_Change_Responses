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
period <- "Historical" # must be capitalized
species <- "Clavigralla shadabi"


# Read in climate data and temperature response parameters for selected insect
temp <- as.data.frame(read_csv(paste0(period," data ",location,".csv")))
sp.data <- subset(as.data.frame(read_csv("Temperature response parameters.csv")),
                  Species == paste(species,location))


# integrate across r(T(t))
n <- length(s)
r <- for(i in 1:n) {
  ifelse(data$T[i] <= Topt, rMax*exp(-1*((data$T[i]-Topt)/(2*sr))^2),
         rMax*(1 - ((data$T[i]-Topt)/(Topt-Tmax))^2))
}