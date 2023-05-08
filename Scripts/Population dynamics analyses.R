###################################################################################
## This R script calculates mean density, CV, and activity period from DDE model ##
###################################################################################

# Load packages and set working directory
library(tidyverse)

# Set working directory (if necessary)
#setwd() # enter working directory of main downloaded file (containing R project file)


# USER: enter species and location or set "all" to TRUE to run analysis for all populations
species <- "Clavigralla shadabi"
location <- "Benin"
all <- TRUE

# USER: save results in csv file?
save <- FALSE


# READ IN TEMPERTURE RESPONSE PARAMETERS AND HABITAT TEMPERATURE PARAMETERS AND CREATE DATA FRAME FOR RESULTS
# Read in files
tr.param <- as.data.frame(read_csv("Model parameters/Temperature response parameters.csv"))
temp.param <- as.data.frame(read_csv("Model parameters/Habitat temperature parameters.csv"))

# Create data frames for results and for list of populations that have gone extinct in the population model
results <- data.frame(tr.param[,1], tr.param[,2], tr.param[,3], tr.param[,4], 1:nrow(tr.param), 1:nrow(tr.param), 1:nrow(tr.param), 1:nrow(tr.param),
                      1:nrow(tr.param), 1:nrow(tr.param), 1:nrow(tr.param), 1:nrow(tr.param), 1:nrow(tr.param))
names(results) <- c("Population","Location","Latitude","Habitat","mean.r","mean.f","CV.r","CV.f","active.r","active.f","delta.mean","delta.CV","delta.active")
extinctions <- data.frame(populations = NA) # data frame to list populations that have gone extinct in the population model


# RUN ANALYSES FOR EACH POPULATION
for(s in 1:nrow(tr.param)) {
 
  # Get parameters for selected population
  if(all == TRUE) {
    param <- tr.param[s,]
    t.param <- temp.param[s,]
  } else {
    param <- tr.param[tr.param$Population == paste(species,location),]
    t.param <- temp.param[temp.param$Population == paste(species,location),]
  }
  
  # Read in DDE model time-series data for recent (.r) and future (.f) time periods
  TS.r <- as.data.frame(read_csv(paste0("Time series data/Recent time series ",param[1],".csv")))
  TS.f <- as.data.frame(read_csv(paste0("Time series data/Future time series ",param[1],".csv")))
  
  # Determine if the population is extinct and remove any data after extinction
  ext <- FALSE
  if(length(TS.r[TS.f$J == 0 & TS.r$A == 0,"Time"]) > 0) {
    print(paste("PROBLEM:", param[,1], "is predicted to be extinct in its current climate!"))
    break
  }
  if(length(TS.f[TS.f$J == 0 & TS.f$A == 0,"Time"]) > 0) {
    ext <- TRUE
    extinctions[nrow(extinctions) + 1,] <- param[,1] # record name of population that went extinct in the population model
    TS.f <- TS.f[1:min(TS.f[TS.f$J == 0 & TS.f$A == 0, "Time"] + 1),] # remove any data after extinction
  }
  
  # Define habitat temperature function (Eq. 5) for recent (.r) and future (.f) time periods
  T.r <- function(t) { (t.param$meanT.r + t.param$delta_mean.r*(t+start_yr*365)) - (t.param$amplT.r + t.param$delta_ampl.r*(t+start_yr*365))*cos(2*pi*(t + t.param$shiftT.r)/365) }
  T.f <- function(t) { (t.param$meanT.f + t.param$delta_mean.f*(t+start_yr*365)) - (t.param$amplT.f + t.param$delta_ampl.f*(t+start_yr*365))*cos(2*pi*(t + t.param$shiftT.f)/365) }
  
  # start and end times for integration
  start_yr <- 0 # from Python DDE model
  num_yr <- 5 # number of years over which to integrate time series
  end.r <- nrow(TS.r) # end of integration for recent period
  end.f <- nrow(TS.f) # end of integration for future period
  start.r <- max(1, end.r - 365*num_yr + 1) # start of integration for recent period (must not be less than 1)
  start.f <- max(1, end.f - 365*num_yr + 1) # start of integration for future period (must not be less than 1)
  
  
  ########################### MODEL: RECENT CLIMATE ############################
  # Quantify mean adult density and active period
  mean.r <- 0 # mean adult density
  active.r <- 0 # active period of insect
  for(i in start.r:end.r) {
    if(T.r(i) >= param$Tmin) {
      mean.r <- mean.r + TS.r$A[i]
      active.r <- active.r + 1 }}
  mean.r <- mean.r/active.r # mean daily adult density
  active.r <- active.r/num_yr # mean activity period per year
  
  # Quantify coefficient of variation (CV) of adult density
  ss.r <- 0 # sum of squares for coefficient of variation (CV)
  for(i in start.r:end.r) {
    if(T.r(i) >= param$Tmin) { ss.r <- ss.r + (TS.r$A[i] - mean.r)^2 }} # sum of squares
  cv.r <- sqrt(ss.r/active.r)/mean.r # CV of daily adult density
  
  # Set mean adult density and CV to zero if population has gone extinct
  if(TS.r[TS.r$Time == end.r - 1, "A"] == 0) {
    mean.r <- 0
    cv.r <- 0
  }
  
  
  ########################### MODEL: FUTURE CLIMATE #############################
  # Quantify mean adult density and active period
  mean.f <- 0 # mean adult density
  active.f <- 0 # active period of insect
  for(i in start.f:end.f) {
    if(T.f(i) >= param$Tmin) {
      mean.f <- mean.f + TS.f$A[i]
      active.f <- active.f + 1 }}
  mean.f <- mean.f/active.f # mean daily adult density
  active.f <- active.f/num_yr # mean activity period per year
  
  # Quantify coefficient of variation (CV) of adult density
  ss.f <- 0 # sum of squares for coefficient of variation (CV)
  for(i in start.f:end.f) {
    if(T.f(i) >= param$Tmin) { ss.f <- ss.f + (TS.f$A[i] - mean.f)^2 }} # sum of squares
  cv.f <- sqrt(ss.f/active.f)/mean.f # CV of daily adult density
  
  # Set mean adult density and CV to zero if population has gone extinct
  if(TS.f[TS.f$Time == end.f - 1, "A"] == 0) {
    mean.f <- 0
    cv.f <- 0
  }
  
  
  ############################### RECORD RESULTS ###############################
  # INPUT RESUTS INTO ARRAY
  if(all == TRUE) {
      results[s,5] <- mean.r
      results[s,6] <- mean.f
      results[s,7] <- cv.r
      results[s,8] <- cv.f
      results[s,9] <- active.r
      results[s,10] <- active.f
      results[s,11] <- (mean.f - mean.r)/mean.r
      results[s,12] <- cv.f - cv.r
      results[s,13] <- (active.f - active.r)/365
    }

  
  # BREAK FOR LOOP IF ANALYSES ARE RUN FOR A SPECIFIED POPULATION
  if(all == FALSE) { break  }
}


# SAVE RESULTS IN CSV FILE
if(save == TRUE) { write_csv(results, "Model predictions/Predictions population dynamics.csv") }