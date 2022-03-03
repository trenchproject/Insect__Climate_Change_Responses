##############################################################################################
### This R script calculates species' mean density, CV, and activity period from DDE model ###
##############################################################################################

# Load packages and set working directory
library(tidyr)
library(ggplot2)
library(dplyr)
library(tidyverse)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# USER: enter species and location or set "all" to TRUE to run analysis for all species
species <- "Macrosiphum euphorbiae"
location <- "Canada"
all <- F

# USER: include diurnal variation?
daily <- FALSE

# USER: output results in csv (only if all == TRUE)?
output <- TRUE


# READ LIFE HISTORY AND TEMPERTURE PARAMETERS
param.all <- as.data.frame(read_csv("Temperature response parameters.csv"))
ifelse(daily == TRUE, t.param.all <- as.data.frame(read_csv("Temperature parameters.csv")),
       t.param.all <- as.data.frame(read_csv("Temperature parameters Tave.csv")))
# Get parameters for selected species
if(all == FALSE) {
  param <- subset(as.data.frame(read_csv("Temperature response parameters.csv")), Species == paste(species,location))
  ifelse(daily == TRUE, t.param <- subset(as.data.frame(read_csv("Temperature parameters.csv")), Species == paste(species,location)),
         t.param <- subset(as.data.frame(read_csv("Temperature parameters Tave.csv")), Species == paste(species,location)))
}


# CREATE ARRAY FOR RESULTS
if(all == TRUE) {
  results <- data.frame(param.all[,1], param.all[,2], param.all[,3], param.all[,4], 1:nrow(param.all), 1:nrow(param.all), 1:nrow(param.all), 1:nrow(param.all),
                        1:nrow(param.all), 1:nrow(param.all), 1:nrow(param.all), 1:nrow(param.all), 1:nrow(param.all))
  names(results) <- c("Species","Latitude","Habitat","Subfamily","mean.h","mean.f","CV.h","CV.f","active.h","active.f","delta.mean","delta.CV","delta.active")
}


# RUN ANALYSES FOR EACH SPECIES
for(s in 1:nrow(param.all)) {
 
  # Get parameters for selected species
  if(all == TRUE) {
    param <- param.all[s,]
    t.param <- t.param.all[s,]
  }
  
  # Read in DDE model data for selected insect
  if(all == FALSE) {
    ifelse(daily == TRUE, TS.h <- as.data.frame(read_csv(paste0("Time series data Diurnal/Historical time series ",species," ",location,".csv"))),
           TS.h <- as.data.frame(read_csv(paste0("Time series data Tave Dev/Historical time series ",species," ",location,".csv"))))
    ifelse(daily == TRUE, TS.f <- as.data.frame(read_csv(paste0("Time series data Diurnal/Future time series ",species," ",location,".csv"))),
           TS.f <- as.data.frame(read_csv(paste0("Time series data Tave Dev/Future time series ",species," ",location,".csv"))))
  }
  if(all == TRUE) {
    ifelse(daily == TRUE, TS.h <- as.data.frame(read_csv(paste0("Time series data Diurnal/Historical time series ",param[1],".csv"))),
           TS.h <- as.data.frame(read_csv(paste0("Time series data Tave Dev/Historical time series ",param[1],".csv"))))
    ifelse(daily == TRUE, TS.f <- as.data.frame(read_csv(paste0("Time series data Diurnal/Future time series ",param[1],".csv"))),
           TS.f <- as.data.frame(read_csv(paste0("Time series data Tave Dev/Future time series ",param[1],".csv"))))
  }
  
  # Define temperature function and the start and end dates for integration
  init_years <- 0 # from Python DDE model
  num_yr <- 5 # number of years over which to integrate time series
  # habitat temperature
  T.h <- function(t) { (t.param$meanT.h + t.param$delta_mean.h*(t+init_years*365)) - (t.param$amplT.h + t.param$delta_ampl.h*(t+init_years*365))*cos(2*pi*(t + t.param$shiftT.h)/365) - t.param$amplD.h*cos(2*pi*t) }
  T.f <- function(t) { (t.param$meanT.f + t.param$delta_mean.f*(t+init_years*365)) - (t.param$amplT.f + t.param$delta_ampl.f*(t+init_years*365))*cos(2*pi*(t + t.param$shiftT.f)/365) - t.param$amplD.f*cos(2*pi*t) }
  # start and end times for integration
  end.h <- nrow(TS.h)
  end.f <- nrow(TS.f)
  start.h <- end.h - 365*num_yr + 1
  start.f <- end.f - 365*num_yr + 1
  # extract time-series data between start and end
  #TS.h <- TS.h[c(start.h:end.h),]
  #TS.f <- TS.f[c(start.f:end.f),]
  
  
  ################################# MODEL: HISTORICAL CLIMATE ##################################
  # Integrate across data from DDE model
  mean.h <- 0
  ss.h <- 0
  active.h <- 0
  # mean density and active period
  for(i in start.h:end.h) {
    if(T.h(i) >= param$Tmin) {
      mean.h <- mean.h + TS.h$A[i]
      active.h <- active.h + 1 }}
  mean.h <- mean.h/active.h # mean daily adult density
  # CV of density
  for(i in start.h:end.h) {
    if(T.h(i) >= param$Tmin) { ss.h <- ss.h + (TS.h$A[i] - mean.h)^2 }} # sum of squares
  cv.h <- sqrt(ss.h/active.h)/mean.h # CV of daily adult density
  active.h <- active.h/num_yr # activity period per year
  # set mean and CV to zero if insect is extinct
  if(TS.h[TS.h$Time == end.h - 1, "A"] == 0) {
    mean.h <- 0
    cv.h <- 0
  }
  
  
  ################################### MODEL: FUTURE CLIMATE ####################################
  # Integrate across data from DDE model
  mean.f <- 0
  ss.f <- 0
  active.f <- 0
  # mean density and active period
  for(i in start.f:end.f) {
    if(T.f(i) >= param$Tmin) {
      mean.f <- mean.f + TS.f$A[i]
      active.f <- active.f + 1 }}
  mean.f <- mean.f/active.f # mean daily adult density
  # CV of density
  for(i in start.f:end.f) {
    if(T.f(i) >= param$Tmin) { ss.f <- ss.f + (TS.f$A[i] - mean.f)^2 }} # sum of squares
  cv.f <- sqrt(ss.f/active.f)/mean.f # CV of daily adult density
  active.f <- active.f/num_yr # activity period per year
  # set mean and CV to zero if insect is extinct
  if(TS.f[TS.f$Time == end.f - 1, "A"] == 0) {
    mean.f <- 0
    cv.f <- 0
  }
  
  
  ################################# RECORD AND DISPLAY RESULTS ###############################
  # INPUT RESUTS INTO ARRAY
  if(all == TRUE) {
      results[s,5] <- mean.h
      results[s,6] <- mean.f
      results[s,7] <- cv.h
      results[s,8] <- cv.f
      results[s,9] <- active.h
      results[s,10] <- active.f
      results[s,11] <- mean.f - mean.h
      results[s,12] <- cv.f - cv.h
      results[s,13] <- (active.f - active.h)/365
    }

  
  # BREAK FOR LOOP IF ANALYSES ARE RUN FOR A SPECIFIED SPECIES
  if(all == FALSE) { break  }
}


# OUTPUT RESULTS IN CSV FILE
if(output == TRUE && all == TRUE) {
  write_csv(results, "Predictions/Predictions population dynamics.csv") }


# SUMMARIZE RESULTS
if(all == FALSE) { 
  print(paste("mean.h:", mean.h))
  print(paste("mean.f:", mean.f))
  print(paste("CV.h:", cv.h))
  print(paste("CV.f:", cv.f))
  print(paste("active.h:", active.h))
  print(paste("active.f:", active.f))
  print(paste("delta.mean:", mean.f - mean.h))
  print(paste("delta.CV:", cv.f - cv.h))
  print(paste("delta.active:", (active.f - active.h)/365))
}
if(all == TRUE) { print(results) }

