############################################################################
#### This R script calculates species' traits using TPCs and DDE model #####
############################################################################

# Load packages
library(tidyverse)
library(cubature) # for integrating over TPC and model predictions

# Set working directory (if necessary)
#setwd() # enter working directory of main downloaded file (containing R project file)

# USER: enter species and location or set "all" to TRUE to run analysis for all species
species <- "Clavigralla shadabi"
location <- "Benin"
all <- TRUE

# USER: save results in csv file? (only if all == TRUE)
save <- TRUE


# READ IN TEMPERTURE RESPONSE PARAMETERS AND HABITAT TEMPERATURE PARAMETERS AND CREATE DATA FRAME FOR RESULTS
# Read in files
tr.param <- as.data.frame(read_csv("Model parameters/Temperature response parameters.csv"))
temp.param <- as.data.frame(read_csv("Model parameters/Habitat temperature parameters.csv"))

# Create data frames for results
if(all == TRUE) {
  results <- data.frame(tr.param[,1], tr.param[,2], tr.param[,3], tr.param[,4], 1:nrow(tr.param), 1:nrow(tr.param), 1:nrow(tr.param), 1:nrow(tr.param), 1:nrow(tr.param), 1:nrow(tr.param))
  names(results) <- c("Species","Location","Latitude","Habitat","TPC.r","TPC.f","Model.r","Model.f","delta.TPC","delta.model")
  r.data <- results; R0.data <- results; birth.data <- results; dev.data <- results; long.data <- results; surv.data <- results # make empty data frames for each trait
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
  
  # Use density-independent model (DI) time-series data for quantifying r_m (TS.r) for recent (.r) and future (.f) time periods
  TS.r.r <- as.data.frame(read_csv(paste0("Time series data DI/Recent time series ",param[1],".csv")))
  TS.r.f <- as.data.frame(read_csv(paste0("Time series data DI/Future time series ",param[1],".csv")))
  # Use density-dependent model time-series (TS) data for quantifying all other traits for recent (.r) and future (.f) time periods
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
    TS.f <- TS.f[1:min(TS.f[TS.f$J == 0 & TS.f$A == 0, "Time"]),] # remove any data after extinction
  }
  if(length(TS.r.f[TS.r.f$J == 0 & TS.r.f$A == 0,"Time"]) > 0) { TS.r.f <- TS.r.f[1:min(TS.r.f[TS.r.f$J == 0 & TS.r.f$A == 0, "Time"]),] } # remove any data after extinction in future density-independent model
  
  # Define temperature function for models (Eq. 5) for recent (.r) and future (.f) time periods
  T.r <- function(t) { (t.param$meanT.r + t.param$delta_mean.r*t) - (t.param$amplT.r + t.param$delta_ampl.r*t)*cos(2*pi*(t + t.param$shiftT.r)/365) }
  T.f <- function(t) { (t.param$meanT.f + t.param$delta_mean.f*t) - (t.param$amplT.f + t.param$delta_ampl.f*t)*cos(2*pi*(t + t.param$shiftT.f)/365) }

  # Extract data
  TS.r <- TS.r[(nrow(TS.r) - 365*5 + 1):nrow(TS.r),] # Extract last 5 years of density-dependent population dynamics
  TS.f <- TS.f[(max(nrow(TS.f) - 365*5 + 1, nrow(TS.f) %% 365 + 1)):nrow(TS.f),] # Extract last 5 years of density-dependent population dynamics (or up to 5 years before extinction)
  if(ext == FALSE) {
    TS.r.r <- TS.r.r[(nrow(TS.r.r) - 365*2 + 1):nrow(TS.r.r),] # Extract last 2 years of density-independent population dynamics
    TS.r.f <- TS.r.f[(max(nrow(TS.r.f) - 365*2 + 1, nrow(TS.r.f) %% 365 + 1)):nrow(TS.r.f),] # Extract last 2 years of density-independent population dynamics
  } else {
    TS.r.r <- TS.r.r[(nrow(TS.r.r) - 365*1 + 1):nrow(TS.r.r),] # Extract last 2 years of density-independent population dynamics
    TS.r.f <- TS.r.f[(max(nrow(TS.r.f) - 365*1 + 1, nrow(TS.r.f) %% 365 + 1)):nrow(TS.r.f),] # Extract up to 1 year of density-independent population dynamics before extinction
  }

  # Start and end times for integration
  start.r <- TS.r[1, "Time"] # start of integration for recent period
  start.f <- TS.f[1, "Time"] # start of integration for future period
  end.r <- TS.r[nrow(TS.r), "Time"] # end of integration for recent period
  end.f <- TS.f[nrow(TS.f), "Time"] # end of integration for future period
  if(ext == FALSE) {
    start.r.r <- 73*365 # population model begins 73 years into simulation (see "start_yr" in "DDE population dynamics.py")
    start.r.f <- 73*365 # population model begins 73 years into simulation (see "start_yr" in "DDE population dynamics.py")
  } else {
    start.r.r <- end.r - 1*365 + 1 # start 1 year before the population goes extinct
    start.r.f <- end.f - 1*365 + 1 # start 1 year before the population goes extinct
  }
  end.r.r <- end.r
  end.r.f <- end.f
  
  # Calculate length of active season (when habitat temperature is above the minimum developmental temperature, Tmin)
  season.r <- 0
  season.f <- 0
  season.r.r <- 0
  season.r.f <- 0
  for(t in seq(start.r,end.r,1)) { if(T.r(t) >= param$Tmin) { season.r <- season.r + 1 }}
  for(t in seq(start.f,end.f,1)) { if(T.f(t) >= param$Tmin) { season.f <- season.f + 1 }}
  for(t in seq(start.r.r,end.r.r,1)) { if(T.r(t) >= param$Tmin) { season.r.r <- season.r.r + 1 }}
  for(t in seq(start.r.f,end.r.f,1)) { if(T.f(t) >= param$Tmin) { season.r.f <- season.r.f + 1 }}

  
################################## DIRECTLY INTEGRATE TPC: recent CLIMATE ###################################
  # Integrate across fitness metrics and components (set to zero when habitat temperature is below the minimum developmental temperature, Tmin, for integration)
  # Fitness metrics
  r.r <- function(t) { ifelse(T.r(t) < param$Tmin, 0, ifelse(T.r(t) <= param$Toptr, param$rMax*exp(-((T.r(t)-param$Toptr)^2)/(2*param$sr^2)),
                  param$rMax*(1 - ((T.r(t)-param$Toptr)/(param$Toptr-param$Tmaxr))^2))) }
  R0.r <- function(t) { ifelse(T.r(t) < param$Tmin, 0, param$R0Topt*exp(-((T.r(t)-param$ToptR0)^2)/(2*param$sR0^2))) }
  # Fitness components
  b.r <- function(t) { ifelse(T.r(t) < param$Tmin, 0, param$bTopt*exp(-((T.r(t)-param$Toptb)^2)/(2*param$sb^2))) }
  g.r <- function(t) { ifelse(T.r(t) < param$Tmin, 0, ifelse(T.r(t) <= param$Toptg, param$gTR*(T.r(t)/param$TR)*exp(param$Ag*(1/param$TR-1/T.r(t)))/(1+exp(param$AL*(1/param$TL-1/T.r(t)))), # if T(t) < Toptg, use monotonic g(T)
                              ifelse(T.r(t) <= param$Tmaxg, param$gMax, 0))) } # If T(t) < Tmaxg, g(T) = gMax; otherwise, g(T) = 0
  dJ.r <- function(t) { ifelse(T.r(t) < param$Tmin, 0, param$dJTR*exp(param$AdJ*(1/param$TR-1/T.r(t)))) }
  dA.r <- function(t) { ifelse(T.r(t) < param$Tmin, 0, param$dATR*exp(param$AdA*(1/param$TR-1/T.r(t)))) }
  s.r <- function(t) { ifelse(T.r(t) < param$Tmin, 0, exp(-dJ.r(t)/g.r(t))) }
  
  # Integration (Note: hcubature must be used with overwintering)
  # Fitness metrics
  r.TPC.r <- cubintegrate(r.r, lower = start.r.r, upper = end.r.r, method = "hcubature")$integral/season.r.r
  R0.TPC.r <- cubintegrate(R0.r, lower = start.r, upper = end.r, method = "hcubature")$integral/season.r
  # Fitness components
  b.TPC.r <- cubintegrate(b.r, lower = start.r, upper = end.r, method = "hcubature")$integral/season.r
  g.TPC.r <- cubintegrate(g.r, lower = start.r, upper = end.r, method = "hcubature")$integral/(end.r-start.r) # different denominator b/c development time does not stop during overwintering
  dA.TPC.r <- cubintegrate(dA.r, lower = start.r, upper = end.r, method = "hcubature")$integral/season.r
  s.TPC.r <- cubintegrate(s.r, lower = start.r, upper = end.r, method = "hcubature")$integral/season.r
  
  
##################################### DIRECTLY INTEGRATE TPC: FUTURE CLIMATE ####################################
  # Integrate across fitness metrics and components (set to zero when habitat temperature is below the minimum developmental temperature, Tmin, for integration)
  # Fitness metrics
  r.f <- function(t) { ifelse(T.f(t) < param$Tmin, 0, ifelse(T.f(t) <= param$Toptr, param$rMax*exp(-((T.f(t)-param$Toptr)^2)/(2*param$sr^2)),
                                                             param$rMax*(1 - ((T.f(t)-param$Toptr)/(param$Toptr-param$Tmaxr))^2))) }
  R0.f <- function(t) { ifelse(T.f(t) < param$Tmin, 0, param$R0Topt*exp(-((T.f(t)-param$ToptR0)^2)/(2*param$sR0^2))) }
  # Fitness components
  b.f <- function(t) { ifelse(T.f(t) < param$Tmin, 0, param$bTopt*exp(-((T.f(t)-param$Toptb)^2)/(2*param$sb^2))) }
  g.f <- function(t) { ifelse(T.f(t) < param$Tmin, 0, ifelse(T.f(t) <= param$Toptg, param$gTR*(T.f(t)/param$TR)*exp(param$Ag*(1/param$TR-1/T.f(t)))/(1+exp(param$AL*(1/param$TL-1/T.f(t)))), # if T(t) < Toptg, use monotonic g(T)
                                                             ifelse(T.f(t) <= param$Tmaxg, param$gMax, 0))) } # If T(t) < Tmaxg, g(T) = gMax; otherwise, g(T) = 0
  dJ.f <- function(t) { ifelse(T.f(t) < param$Tmin, 0, param$dJTR*exp(param$AdJ*(1/param$TR-1/T.f(t)))) }
  dA.f <- function(t) { ifelse(T.f(t) < param$Tmin, 0, param$dATR*exp(param$AdA*(1/param$TR-1/T.f(t)))) }
  s.f <- function(t) { ifelse(T.f(t) < param$Tmin, 0, exp(-dJ.f(t)/g.f(t))) }

  # Integration (Note: hcubature must be used with overwintering)
  # Fitness metrics
  r.TPC.f <- cubintegrate(r.f, lower = start.r.f, upper = end.r.f, method = "hcubature")$integral/season.r.f
  R0.TPC.f <- cubintegrate(R0.f, lower = start.f, upper = end.f, method = "hcubature")$integral/season.f
  # Fitness components
  b.TPC.f <- cubintegrate(b.f, lower = start.f, upper = end.f, method = "hcubature")$integral/season.f
  g.TPC.f <- cubintegrate(g.f, lower = start.f, upper = end.f, method = "hcubature")$integral/(end.f-start.f) # different denominator b/c development time does not stop during overwintering
  dA.TPC.f <- cubintegrate(dA.f, lower = start.f, upper = end.f, method = "hcubature")$integral/season.f
  s.TPC.f <- cubintegrate(s.f, lower = start.f, upper = end.f, method = "hcubature")$integral/season.f


################################# DDE MODEL: recent CLIMATE ##################################
  # Create column for birth rate, adult mortality, and R0 in time-series data
  b.r <- function(t) { param$bTopt*exp(-((T.r(t)-param$Toptb)^2)/(2*param$sb^2)) }
  TS.r$b <- b.r(TS.r$Time)
  dA.r <- function(t) { param$dATR*exp(param$AdA*(1/param$TR-1/T.r(t))) }
  TS.r$dA <- dA.r(TS.r$Time)
  TS.r$R0 <- ifelse(TS.r$dA != 0, TS.r$b/TS.r$dA*TS.r$S, 0) # R0 = 0 during overwintering
  
  # Integrate across data from DDE model
  # initial values
  r.model.r <- 0 # r_m
  r.max.r <- 0 # maximum r_m (for scaling)
  R0.model.r <- 0 # R0
  b.model.r <- 0 # birth rate
  tau.model.r <- 0 # development time
  dA.model.r <- 0 # adult mortality
  s.model.r <- 0 # juvenile survival
  
  # Sum values across each model time-step during active season
  # r_m
  for(i in 2:nrow(TS.r.r)) {
    if(TS.r.r$A[i] > 0 && TS.r.r$A[i-1] > 0 && T.r(start.r.r + TS.r.r$Time[i]) >= param$Tmin) {
      r.max.r <- max(r.max.r, log(TS.r.r$A[i]/TS.r.r$A[i-1]))
      r.model.r <- r.model.r + log(TS.r.r$A[i]/TS.r.r$A[i-1]) } # calculate sum of daily r_m
  }
  # R0 and fitness components   
  for(i in 1:nrow(TS.r)) {
    if(T.r(TS.r$Time[i]) >= param$Tmin) { # i.e., sum traits during active season
      # R0
      R0.model.r <- R0.model.r + TS.r$R0[i]
      # Fitness components (development rate is quantified below)
      b.model.r <- b.model.r + TS.r$b[i]
      dA.model.r <- dA.model.r + TS.r$dA[i]
      s.model.r <- s.model.r + TS.r$S[i]
    }
    # sum development time over all days
    tau.model.r <- tau.model.r + TS.r$tau[i]
  }
  # average values over active season
  r.model.r <- r.model.r/season.r.r
  if(param[1] == "Brevicoryne brassicae US Columbia") { r.max.r <- param$rMax } # NOTE: use rMax for Brevicoryne brassicae b/c r.max.r is biased by one extreme value 
  R0.model.r <- R0.model.r/season.r
  b.model.r <- b.model.r/season.r
  tau.model.r <- tau.model.r/(end.r-start.r)  # different denominator b/c development time does not stop during overwintering
  dA.model.r <- dA.model.r/season.r
  s.model.r <- s.model.r/season.r
  
  
################################### DDE MODEL: FUTURE CLIMATE ####################################
  # Create column for birth rate, adult mortality, and R0 in time-series data
  b.f <- function(t) { param$bTopt*exp(-((T.f(t)-param$Toptb)^2)/(2*param$sb^2)) }
  TS.f$b <- b.f(TS.f$Time)
  dA.f <- function(t) { param$dATR*exp(param$AdA*(1/param$TR-1/T.f(t))) }
  TS.f$dA <- dA.f(TS.f$Time)
  TS.f$R0 <- ifelse(TS.f$dA != 0, TS.f$b/TS.f$dA*TS.f$S, 0) # R0 = 0 during overwintering
  
  # Integrate across data from DDE model
  # initial values
  r.model.f <- 0 # r_m
  r.max.f <- 0 # maximum r_m (for scaling)
  R0.model.f <- 0 # R0
  b.model.f <- 0 # birth rate
  tau.model.f <- 0 # development time
  dA.model.f <- 0 # adult mortality
  s.model.f <- 0 # juvenile survival
  
  # Sum values across each model time-step during active season
  # r_m
  for(i in 2:nrow(TS.r.f)) {
    if(TS.r.f$A[i] > 0 && TS.r.f$A[i-1] > 0 && ((ext == FALSE && T.f(TS.r.f$Time[i] + 73*365) >= param$Tmin) || # model was run 73 years into the future period (see "start_yr" in "DDE population dynamics.py")
       (ext == TRUE && T.f(TS.r.f$Time[i] + 74*365) >= param$Tmin))) { # data was extracted from last year (instead of last 2 years) if the population went extinct
      r.max.f <- max(r.max.f, log(TS.r.f$A[i]/TS.r.f$A[i-1]))
      r.model.f <- r.model.f + log(TS.r.f$A[i]/TS.r.f$A[i-1]) } # calculate sum of daily r_m
  }
  # R0 and fitness components
  for(i in 1:nrow(TS.f)) {
    if(T.f(TS.f$Time[i]) >= param$Tmin) { # i.e., sum traits during active season
      # R0
      R0.model.f <- R0.model.f + TS.f$R0[i]
      # Fitness components (development rate is quantified below)
      b.model.f <- b.model.f + TS.f$b[i]
      dA.model.f <- dA.model.f + TS.f$dA[i]
      s.model.f <- s.model.f + TS.f$S[i]
    }
    # sum development time over all days
    tau.model.f <- tau.model.f + TS.f$tau[i]
  }
  # average values over active season
  r.model.f <- r.model.f/season.r.r
  if(param[1] == "Brevicoryne brassicae US Columbia") { r.max.f <- param$rMax } # NOTE: use rMax for Brevicoryne brassicae b/c r.max.r is biased by one extreme value 
  R0.model.f <- R0.model.f/season.f
  b.model.f <- b.model.f/season.f
  tau.model.f <- tau.model.f/(end.f-start.f)  # different denominator b/c development time does not stop during overwintering
  dA.model.f <- dA.model.f/season.f
  if(ext == FALSE) { s.model.f <- s.model.f/season.f
  } else { s.model.f <- 0 } # if population goes extinct, set survival to zero

  
################################# RECORD AND DISPLAY RESULTS ###############################
  # Calculate optimum development rate, maximum adult longevity, and survival at mean habitat temperature for scaling results
  gTopt <- param$gTR*(param$Toptg/param$TR)*exp(param$Ag*(1/param$TR-1/param$Toptg))/(1+exp(param$AL*(1/param$TL-1/param$Toptg))) #optimum development rate
  Lmax <- 1/min(param$dATR*exp(param$AdA*(1/param$TR-1/T.r(TS.r$Time)))) # maximum adult longevity

  # Input results into arrays
  if(all == TRUE) {
    # r_m
    r.data[s,5] <- round(r.TPC.r/param$rMax,3)
    r.data[s,6] <- round(r.TPC.f/param$rMax,3)
    r.data[s,7] <- round(r.model.r/r.max.r,3)
    r.data[s,8] <- round(r.model.f/r.max.r,3)
    r.data[s,9] <- round((r.TPC.f - r.TPC.r)/param$rMax,3)
    r.data[s,10] <- round((r.model.f - r.model.r)/r.max.r,3)
    
    # R0
    R0.data[s,5] <- round(R0.TPC.r/param$R0Topt,3)
    R0.data[s,6] <- round(R0.TPC.f/param$R0Topt,3)
    R0.data[s,7] <- round(R0.model.r/max(TS.r$R0),3)
    R0.data[s,8] <- round(R0.model.f/max(TS.r$R0),3)
    R0.data[s,9] <- round((R0.TPC.f - R0.TPC.r)/param$R0Topt,3)
    R0.data[s,10] <- round((R0.model.f - R0.model.r)/max(TS.r$R0),3)
  
    # Birth rate
    birth.data[s,5] <- round(b.TPC.r/param$bTopt,3)
    birth.data[s,6] <- round(b.TPC.f/param$bTopt,3)
    birth.data[s,7] <- round(b.model.r/param$bTopt,3)
    birth.data[s,8] <- round(b.model.f/param$bTopt,3)
    birth.data[s,9] <- round((b.TPC.f - b.TPC.r)/param$bTopt,3)
    birth.data[s,10] <- round((b.model.f - b.model.r)/param$bTopt,3)
  
    # Development
    dev.data[s,5] <- round((1/g.TPC.r)/(1/gTopt),3)
    dev.data[s,6] <- round((1/g.TPC.f)/(1/gTopt),3)
    dev.data[s,7] <- round(tau.model.r/min(TS.r$tau),3)
    dev.data[s,8] <- round(tau.model.f/min(TS.r$tau),3)
    dev.data[s,9] <- round((1/g.TPC.f - 1/g.TPC.r)/(1/gTopt),3)
    dev.data[s,10] <- round((tau.model.f - tau.model.r)/min(TS.r$tau),3)
  
    # Longevity
    long.data[s,5] <- round((1/dA.TPC.r)/Lmax,3)
    long.data[s,6] <- round((1/dA.TPC.f)/Lmax,3)
    long.data[s,7] <- round((1/dA.model.r)/(1/min(TS.r$dA)),3)
    long.data[s,8] <- round((1/dA.model.f)/(1/min(TS.r$dA)),3)
    long.data[s,9] <- round((1/dA.TPC.f - 1/dA.TPC.r)/Lmax,3)
    long.data[s,10] <- round((1/dA.model.f - 1/dA.model.r)/(1/min(TS.r$dA)),3)
  
    # Survival
    surv.data[s,5] <- round(s.TPC.r,3)
    surv.data[s,6] <- round(s.TPC.f,3)
    surv.data[s,7] <- round(s.model.r,3)
    surv.data[s,8] <- round(s.model.f,3)
    surv.data[s,9] <- round((s.TPC.f - s.TPC.r)/s.TPC.r,3)
    surv.data[s,10] <- round((s.model.f - s.model.r)/s.TPC.r,3)
  }

  
  # BREAK FOR LOOP IF ANALYSES ARE RUN FOR A SPECIFIED SPECIES
  if(all == FALSE) { break  }
}


# OUTPUT RESULTS IN CSV FILE
if(save == TRUE && all == TRUE) {
  write_csv(r.data, "Predictions/Predictions rm.csv")
  write_csv(R0.data, "Predictions/Predictions R0.csv")
  write_csv(birth.data, "Predictions/Predictions birth.csv")
  write_csv(dev.data, "Predictions/Predictions development.csv")
  write_csv(long.data, "Predictions/Predictions longevity.csv")
  write_csv(surv.data, "Predictions/Predictions survival.csv")
}


# SUMMARIZE RESULTS
if(all == FALSE) {
  # r_m
  print(paste("r.TPC.r:", r.TPC.r/param$rMax))
  print(paste("r.TPC.f:", r.TPC.f/param$rMax))
  print(paste("r.model.r:", r.model.r/r.max.r))
  print(paste("r.model.f:", r.model.f/r.max.r))
  # R0
  print(paste("R0.TPC.r:", R0.TPC.r/param$R0Topt))
  print(paste("R0.TPC.f:", R0.TPC.f/param$R0Topt))
  print(paste("R0.model.r:", R0.model.r/max(TS.r$R0)))
  print(paste("R0.model.f:", R0.model.f/max(TS.r$R0)))
  # Birth rate
  print(paste("b.TPC.r:", b.TPC.r/param$bTopt))
  print(paste("b.TPC.f:", b.TPC.f/param$bTopt))
  print(paste("b.model.r:", b.model.r/param$bTopt))
  print(paste("b.model.f:", b.model.f/param$bTopt))
  # Development rate
  print(paste("tau.TPC.r:", (1/g.TPC.r)/(1/gTopt)))
  print(paste("tau.TPC.f:", (1/g.TPC.f)/(1/gTopt)))
  print(paste("tau.model.r:", tau.model.r/min(TS.r$tau)))
  print(paste("tau.model.f:", tau.model.f/min(TS.r$tau)))
  # Longevity"
  print(paste("1/dA.TPC.r:", (1/dA.TPC.r)/Lmax))
  print(paste("1/dA.TPC.f:", (1/dA.TPC.f)/Lmax))
  print(paste("1/dA.model.r:", (1/dA.model.r)/(1/min(TS.r$dA))))
  print(paste("1/dA.model.f:", (1/dA.model.f)/(1/min(TS.r$dA))))
  # Survival
  print(paste("s.TPC.r:", s.TPC.r))
  print(paste("s.TPC.f:", s.TPC.f))
  print(paste("s.model.r:", s.model.r))
  print(paste("s.model.f:", s.model.f))
} else {
  print(r.data)
  print(R0.data)
  print(birth.data)
  print(dev.data)
  print(long.data)
  print(surv.data)
}