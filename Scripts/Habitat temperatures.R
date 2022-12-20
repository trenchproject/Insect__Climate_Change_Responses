#############################################################################
#### This R script fits sinusoidal functions to habitat temperature data ####
#############################################################################

# Load packages
library(tidyverse)

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('..')

# USER: enter location of climate data (see "Climate station data.xlsx") OR set "all" to TRUE to run analysis for all species
location <- "Benin"
all <- FALSE

# USER: Save model fits?
save <- FALSE

# Read model parameters file
params <- as.data.frame(read_csv("Model parameters/Habitat temperature parameters.csv"))

# Run analyses for all species
for(s in 1:nrow(params)) {
  
  # Input data
  if(all == FALSE) {
    data.h <- as.data.frame(read_csv(paste0("Climate data/Historical climate data ",location,".csv")))
    data.f <- as.data.frame(read_csv(paste0("Climate data/Future climate data ",location,".csv")))
  } else {
    data.h <- as.data.frame(read_csv(paste0("Climate data/Historical climate data ",params[s,]$Location,".csv")))
    data.f <- as.data.frame(read_csv(paste0("Climate data/Future climate data ",params[s,]$Location,".csv")))
  }
    
  # Find model fits for climate in selected location (line 13) or move interatively through rows of “Habitat temperature parameters.csv”
  if(all == FALSE) {
    sp.num <- -1
    for(i in 1:nrow(params)) {
      if(params[i,]$Location == location) {
        sp.num <- i
        break }
    }
    if(sp.num == -1) {
      print("Location not found in Habitat temperature parameters.csv, please update line 13")
    }
  } else { sp.num <- s }
  
  # Quantify daily mean temperatures
  # historical
  data.h$day <- floor(data.h$day) # removes the offset between daily min and daily max from "Read climate data.R"
  data.h.min <- data.h[duplicated(data.h$day),] # selects daily min for all days that also have a daily max
  data.h.max <- data.h[duplicated(data.h$day, fromLast=TRUE),] # selects daily max for all days that also have a daily min
  data.h <- data.frame(data.h.min$day, (data.h.min$T + data.h.max$T)/2) # quantifies daily average temperatures
  names(data.h) <- c("day", "T")
  # future
  data.f$day <- floor(data.f$day) # removes the offset between daily min and daily max from "Read climate data.R"
  data.f.min <- data.f[duplicated(data.f$day),] # selects daily min for all days that also have a daily max
  data.f.max <- data.f[duplicated(data.f$day, fromLast=TRUE),] # selects daily max for all days that also have a daily min
  data.f <- data.frame(data.f.min$day, (data.f.min$T + data.f.max$T)/2) # quantifies daily average temperatures
  names(data.f) <- c("day", "T")
  
  # Fit sinusoidal function (Eq. 5 in manuscript) to climate data
  # historical climate (note that delta_mean and delta_ampl are set to zero in the manuscript, so they are not fit)
  fit.h <- summary(nls(T ~ meanT - amplT*cos(2*pi*(day + shiftT)/365), data = data.h,
                        start = list(meanT = params[sp.num,]$meanT.h, amplT = params[sp.num,]$amplT.h, shiftT = params[sp.num,]$shiftT.h)))
  # future climate
  fit.f <- summary(nls(T ~ (meanT + delta_mean*day) - (amplT + delta_ampl*day) * cos(2*pi*(day + shiftT)/365), data = data.f,
                        start = list(meanT = params[sp.num,]$meanT.f, amplT = params[sp.num,]$amplT.f, shiftT = params[sp.num,]$shiftT.f,
                                     delta_mean = params[sp.num,]$delta_mean.f, delta_ampl = params[sp.num,]$delta_ampl.f)))
  
  # Assess whether delta_mean or delta_ampl are significant in the future climate and set to zero if not
  if(fit.f[["coefficients"]][4,4] > 0.05) { fit.f[["coefficients"]][4,1] <- 0 } # delta_mean
  if(fit.f[["coefficients"]][5,4] > 0.05) { fit.f[["coefficients"]][5,1] <- 0 } # delta_ampl
  
  # Assign model parameters
  params[sp.num,]$meanT.h <- round(fit.h[["coefficients"]][1,1], 1)
  params[sp.num,]$amplT.h <- round(fit.h[["coefficients"]][2,1], 2)
  params[sp.num,]$shiftT.h <- round(fit.h[["coefficients"]][3,1], 1) %% 365 # take modulo to quantify day of year for shift
  params[sp.num,]$meanT.f <- round(fit.f[["coefficients"]][1,1], 1)
  params[sp.num,]$amplT.f <- round(fit.f[["coefficients"]][2,1], 2)
  params[sp.num,]$shiftT.f <- round(fit.f[["coefficients"]][3,1], 1) %% 365 # take modulo to quantify day of year for shift
  params[sp.num,]$delta_mean.f <- round(fit.f[["coefficients"]][4,1], 6)
  params[sp.num,]$delta_ampl.f <- round(fit.f[["coefficients"]][5,1], 6)

  # Break for loop (line 23) if analyses are run for a specified species (all <- FALSE in line 14)
  if(all == FALSE) { break }
}

# Save model parameters to "Habitat temperature parameters.csv" (if desired)
if(save == TRUE) { write.csv(params, "Model parameters/Habitat temperature parameters.csv", row.names = FALSE) }

# Print model fits (if analyses are run for a specified species)
if(all == FALSE) {
  print(fit.h)
  print(fit.f) }