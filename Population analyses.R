##############################################################################################
### This R script calculates species' mean density, CV, and activity period from DDE model ###
##############################################################################################

# Load packages and set working directory
library(tidyverse)

# Set working directory (if necessary)
#setwd() # enter working directory of main downloaded file (containing R project file)

# USER: enter species and location or set "all" to TRUE to run analysis for all species
species <- "Macrosiphum euphorbiae"
location <- "Canada"
all <- TRUE

# USER: output results in csv file? (only if all == TRUE)
output <- TRUE


# READ IN TEMPERTURE RESPONSE PARAMETERS AND HABITAT TEMPERATURE PARAMETERS AND CREATE DATA FRAME FOR RESULTS
# Read in files
tr.param <- as.data.frame(read_csv("Model parameters/Temperature response parameters.csv"))
temp.param <- as.data.frame(read_csv("Model parameters/Habitat temperature parameters.csv"))

# Create data frames for results
if(all == TRUE) {
  results <- data.frame(tr.param[,1], tr.param[,2], tr.param[,3], tr.param[,4], tr.param[,5], 1:nrow(tr.param), 1:nrow(tr.param), 1:nrow(tr.param), 1:nrow(tr.param),
                        1:nrow(tr.param), 1:nrow(tr.param), 1:nrow(tr.param), 1:nrow(tr.param), 1:nrow(tr.param))
  names(results) <- c("Species","Location","Latitude","Habitat","Subfamily","mean.h","mean.f","CV.h","CV.f","active.h","active.f","delta.mean","delta.CV","delta.active")
}
extinctions <- data.frame(populations = NA) # data frame to list populations that have gone extinct in the population model


# RUN ANALYSES FOR EACH SPECIES
for(s in 1:nrow(tr.param)) {
 
  # Get parameters for selected species
  if(all == TRUE) {
    param <- tr.param[s,]
    t.param <- temp.param[s,]
  } else {
    param <- tr.param[tr.param$Species == paste(species,location),]
    t.param <- temp.param[temp.param$Species == paste(species,location),]
  }
  
  # Read in DDE model time-series (TS) data for historical (.h) and future (.f) time periods
  TS.h <- as.data.frame(read_csv(paste0("Time series data/Historical time series ",param[1],".csv")))
  TS.f <- as.data.frame(read_csv(paste0("Time series data/Future time series ",param[1],".csv")))
  
  # Determine if the population is extinct and remove any data after extinction
  ext <- FALSE
  if(length(TS.h[TS.f$J == 0 & TS.h$A == 0,"Time"]) > 0) {
    print(paste("PROBLEM:", param[,1], "is predicted to be extinct in its current climate!"))
    break
  }
  if(length(TS.f[TS.f$J == 0 & TS.f$A == 0,"Time"]) > 0) {
    ext <- TRUE
    extinctions[nrow(extinctions) + 1,] <- param[,1] # record name of population that went extinct in the population model
    TS.f <- TS.f[1:min(TS.f[TS.f$J == 0 & TS.f$A == 0, "Time"] + 1),] # remove any data after extinction
  }
  
  # Define temperature function (Eq. 5) for historical (.h) and future (.f) time periods
  T.h <- function(t) { (t.param$meanT.h + t.param$delta_mean.h*(t+start_yr*365)) - (t.param$amplT.h + t.param$delta_ampl.h*(t+start_yr*365))*cos(2*pi*(t + t.param$shiftT.h)/365) }
  T.f <- function(t) { (t.param$meanT.f + t.param$delta_mean.f*(t+start_yr*365)) - (t.param$amplT.f + t.param$delta_ampl.f*(t+start_yr*365))*cos(2*pi*(t + t.param$shiftT.f)/365) }
  
  # start and end times for integration
  start_yr <- 0 # from Python DDE model
  num_yr <- 5 # number of years over which to integrate time series
  end.h <- nrow(TS.h) # end of integration for historical period
  end.f <- nrow(TS.f) # end of integration for future period
  start.h <- max(1, end.h - 365*num_yr + 1) # start of integration for historical period (must not be less than 1)
  start.f <- max(1, end.f - 365*num_yr + 1) # start of integration for future period (must not be less than 1)
  
  
  ################################# MODEL: HISTORICAL CLIMATE ##################################
  # Integrate across data from DDE model
  mean.h <- 0 # mean adult density
  ss.h <- 0 # sum of squares for coefficient of variation (CV)
  active.h <- 0 # active period of insect
  
  # Quantify mean adult density and active period
  for(i in start.h:end.h) {
    if(T.h(i) >= param$Tmin) {
      mean.h <- mean.h + TS.h$A[i]
      active.h <- active.h + 1 }}
  mean.h <- mean.h/active.h # mean daily adult density
  active.h <- active.h/num_yr # mean activity period per year
  
  # Quantify coefficient of variation (CV) of adult density
  for(i in start.h:end.h) {
    if(T.h(i) >= param$Tmin) { ss.h <- ss.h + (TS.h$A[i] - mean.h)^2 }} # sum of squares
  cv.h <- sqrt(ss.h/active.h)/mean.h # CV of daily adult density
  
  # Set mean adult density and CV to zero if insect is extinct
  if(TS.h[TS.h$Time == end.h - 1, "A"] == 0) {
    mean.h <- 0
    cv.h <- 0
  }
  
  
  ################################### MODEL: FUTURE CLIMATE ####################################
  # Integrate across data from DDE model
  mean.f <- 0 # mean adult density
  ss.f <- 0 # sum of squares for coefficient of variation (CV)
  active.f <- 0 # active period of insect
  
  # Quantify mean adult density and active period
  for(i in start.f:end.f) {
    if(T.f(i) >= param$Tmin) {
      mean.f <- mean.f + TS.f$A[i]
      active.f <- active.f + 1 }}
  mean.f <- mean.f/active.f # mean daily adult density
  active.f <- active.f/num_yr # mean activity period per year
  
  # Quantify coefficient of variation (CV) of adult density
  for(i in start.f:end.f) {
    if(T.f(i) >= param$Tmin) { ss.f <- ss.f + (TS.f$A[i] - mean.f)^2 }} # sum of squares
  cv.f <- sqrt(ss.f/active.f)/mean.f # CV of daily adult density
  
  # Set mean adult density and CV to zero if insect is extinct
  if(TS.f[TS.f$Time == end.f - 1, "A"] == 0) {
    mean.f <- 0
    cv.f <- 0
  }
  
  
  ################################# RECORD AND DISPLAY RESULTS ###############################
  # INPUT RESUTS INTO ARRAY
  if(all == TRUE) {
      results[s,6] <- mean.h
      results[s,7] <- mean.f
      results[s,8] <- cv.h
      results[s,9] <- cv.f
      results[s,10] <- active.h
      results[s,11] <- active.f
      results[s,12] <- (mean.f - mean.h)/mean.h
      results[s,13] <- cv.f - cv.h
      results[s,14] <- (active.f - active.h)/365
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
  print(paste("delta.mean:", (mean.f - mean.h)/mean.h))
  print(paste("delta.CV:", cv.f - cv.h))
  print(paste("delta.active:", (active.f - active.h)/365))
}
if(all == TRUE) { print(results) }