############################################################################
#### This R script calculates species' traits using TPCs and DDE model #####
############################################################################

# Load packages and set working directory
library(tidyverse)
library(cubature) # for integrating over TPC and model predictions

# Set working directory (if necessary)
#setwd() # enter working directory of main downloaded file (containing R project file)

# USER: enter species and location or set "all" to TRUE to run analysis for all species
species <- "Clavigralla shadabi"
location <- "Benin"
all <- TRUE

# USER: save results in csv file? (only if all == TRUE)
save <- FALSE


# READ IN TEMPERTURE RESPONSE PARAMETERS AND HABITAT TEMPERATURE PARAMETERS AND CREATE DATA FRAME FOR RESULTS
# Read in files
tr.param <- as.data.frame(read_csv("Model parameters/Temperature response parameters.csv"))
temp.param <- as.data.frame(read_csv("Model parameters/Habitat temperature parameters.csv"))

# Create data frames for results
if(all == TRUE) {
  results <- data.frame(tr.param[,1], tr.param[,2], tr.param[,3], tr.param[,4], tr.param[,5], 1:nrow(tr.param), 1:nrow(tr.param), 1:nrow(tr.param), 1:nrow(tr.param), 1:nrow(tr.param), 1:nrow(tr.param))
  names(results) <- c("Species","Location","Latitude","Habitat","Subfamily","TPC.h","TPC.f","Model.h","Model.f","delta.TPC","delta.model")
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
  
  # Use density-independent model (DI) time-series data for quantifying r_m (TS.r) for historical (.h) and future (.f) time periods
  TS.r.h <- as.data.frame(read_csv(paste0("Time series data DI/Historical time series ",param[1],".csv")))
  TS.r.f <- as.data.frame(read_csv(paste0("Time series data DI/Future time series ",param[1],".csv")))
  # Use density-dependent model time-series (TS) data for quantifying all other traits for historical (.h) and future (.f) time periods
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
    TS.f <- TS.f[1:min(TS.f[TS.f$J == 0 & TS.f$A == 0, "Time"]),] # remove any data after extinction
  }
  if(length(TS.r.f[TS.r.f$J == 0 & TS.r.f$A == 0,"Time"]) > 0) { TS.r.f <- TS.r.f[1:min(TS.r.f[TS.r.f$J == 0 & TS.r.f$A == 0, "Time"]),] } # remove any data after extinction in future density-independent model
  
  # Define temperature function for models (Eq. 5) for historical (.h) and future (.f) time periods
  T.h <- function(t) { (t.param$meanT.h + t.param$delta_mean.h*t) - (t.param$amplT.h + t.param$delta_ampl.h*t)*cos(2*pi*(t + t.param$shiftT.h)/365) }
  T.f <- function(t) { (t.param$meanT.f + t.param$delta_mean.f*t) - (t.param$amplT.f + t.param$delta_ampl.f*t)*cos(2*pi*(t + t.param$shiftT.f)/365) }

  # Extract data
  TS.h <- TS.h[(nrow(TS.h) - 365*5 + 1):nrow(TS.h),] # Extract last 5 years of density-dependent population dynamics
  TS.f <- TS.f[(max(nrow(TS.f) - 365*5 + 1, nrow(TS.f) %% 365 + 1)):nrow(TS.f),] # Extract last 5 years of density-dependent population dynamics (or up to 5 years before extinction)
  if(ext == FALSE) {
    TS.r.h <- TS.r.h[(nrow(TS.r.h) - 365*2 + 1):nrow(TS.r.h),] # Extract last 2 years of density-independent population dynamics
    TS.r.f <- TS.r.f[(max(nrow(TS.r.f) - 365*2 + 1, nrow(TS.r.f) %% 365 + 1)):nrow(TS.r.f),] # Extract last 2 years of density-independent population dynamics
  } else {
    TS.r.h <- TS.r.h[(nrow(TS.r.h) - 365*1 + 1):nrow(TS.r.h),] # Extract last 2 years of density-independent population dynamics
    TS.r.f <- TS.r.f[(max(nrow(TS.r.f) - 365*1 + 1, nrow(TS.r.f) %% 365 + 1)):nrow(TS.r.f),] # Extract up to 1 year of density-independent population dynamics before extinction
  }

  # Start and end times for integration
  start.h <- TS.h[1, "Time"] # start of integration for historical period
  start.f <- TS.f[1, "Time"] # start of integration for future period
  end.h <- TS.h[nrow(TS.h), "Time"] # end of integration for historical period
  end.f <- TS.f[nrow(TS.f), "Time"] # end of integration for future period
  if(ext == FALSE) {
    start.r.h <- 73*365 # population model begins 73 years into simulation (see "start_yr" in "DDE population dynamics.py")
    start.r.f <- 73*365 # population model begins 73 years into simulation (see "start_yr" in "DDE population dynamics.py")
  } else {
    start.r.h <- end.h - 1*365 + 1 # start 1 year before the population goes extinct
    start.r.f <- end.f - 1*365 + 1 # start 1 year before the population goes extinct
  }
  end.r.h <- end.h
  end.r.f <- end.f
  
  # Calculate length of active season (when habitat temperature is above the minimum developmental temperature, Tmin)
  season.h <- 0
  season.f <- 0
  season.r.h <- 0
  season.r.f <- 0
  for(t in seq(start.h,end.h,1)) { if(T.h(t) >= param$Tmin) { season.h <- season.h + 1 }}
  for(t in seq(start.f,end.f,1)) { if(T.f(t) >= param$Tmin) { season.f <- season.f + 1 }}
  for(t in seq(start.r.h,end.r.h,1)) { if(T.h(t) >= param$Tmin) { season.r.h <- season.r.h + 1 }}
  for(t in seq(start.r.f,end.r.f,1)) { if(T.f(t) >= param$Tmin) { season.r.f <- season.r.f + 1 }}

  
################################## DIRECTLY INTEGRATE TPC: HISTORICAL CLIMATE ###################################
  # Integrate across fitness metrics and components (set to zero when habitat temperature is below the minimum developmental temperature, Tmin, for integration)
  # Fitness metrics
  r.h <- function(t) { ifelse(T.h(t) < param$Tmin, 0, ifelse(T.h(t) <= param$Toptr, param$rMax*exp(-((T.h(t)-param$Toptr)^2)/(2*param$sr^2)),
                  param$rMax*(1 - ((T.h(t)-param$Toptr)/(param$Toptr-param$Tmaxr))^2))) }
  R0.h <- function(t) { ifelse(T.h(t) < param$Tmin, 0, param$R0Topt*exp(-((T.h(t)-param$ToptR0)^2)/(2*param$sR0^2))) }
  # Fitness components
  b.h <- function(t) { ifelse(T.h(t) < param$Tmin, 0, param$bTopt*exp(-((T.h(t)-param$Toptb)^2)/(2*param$sb^2))) }
  g.h <- function(t) { ifelse(T.h(t) < param$Tmin, 0, ifelse(T.h(t) <= param$Toptg, param$gTR*(T.h(t)/param$TR)*exp(param$Ag*(1/param$TR-1/T.h(t)))/(1+exp(param$AL*(1/param$TL-1/T.h(t)))), # if T(t) < Toptg, use monotonic g(T)
                              ifelse(T.h(t) <= param$Tmaxg, param$gMax, 0))) } # If T(t) < Tmaxg, g(T) = gMax; otherwise, g(T) = 0
  dJ.h <- function(t) { ifelse(T.h(t) < param$Tmin, 0, param$dJTR*exp(param$AdJ*(1/param$TR-1/T.h(t)))) }
  dA.h <- function(t) { ifelse(T.h(t) < param$Tmin, 0, param$dATR*exp(param$AdA*(1/param$TR-1/T.h(t)))) }
  s.h <- function(t) { ifelse(T.h(t) < param$Tmin, 0, exp(-dJ.h(t)/g.h(t))) }
  
  # Integration (Note: hcubature must be used with overwintering)
  # Fitness metrics
  r.TPC.h <- cubintegrate(r.h, lower = start.r.h, upper = end.r.h, method = "hcubature")$integral/season.r.h
  R0.TPC.h <- cubintegrate(R0.h, lower = start.h, upper = end.h, method = "hcubature")$integral/season.h
  # Fitness components
  b.TPC.h <- cubintegrate(b.h, lower = start.h, upper = end.h, method = "hcubature")$integral/season.h
  g.TPC.h <- cubintegrate(g.h, lower = start.h, upper = end.h, method = "hcubature")$integral/(end.h-start.h) # different denominator b/c development time does not stop during overwintering
  dA.TPC.h <- cubintegrate(dA.h, lower = start.h, upper = end.h, method = "hcubature")$integral/season.h
  s.TPC.h <- cubintegrate(s.h, lower = start.h, upper = end.h, method = "hcubature")$integral/season.h
  
  
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


################################# DDE MODEL: HISTORICAL CLIMATE ##################################
  # Create column for birth rate, adult mortality, and R0 in time-series data
  b.h <- function(t) { param$bTopt*exp(-((T.h(t)-param$Toptb)^2)/(2*param$sb^2)) }
  TS.h$b <- b.h(TS.h$Time)
  dA.h <- function(t) { param$dATR*exp(param$AdA*(1/param$TR-1/T.h(t))) }
  TS.h$dA <- dA.h(TS.h$Time)
  TS.h$R0 <- ifelse(TS.h$dA != 0, TS.h$b/TS.h$dA*TS.h$S, 0) # R0 = 0 during overwintering
  
  # Integrate across data from DDE model
  # initial values
  r.model.h <- 0 # r_m
  r.max.h <- 0 # maximum r_m (for scaling)
  R0.model.h <- 0 # R0
  b.model.h <- 0 # birth rate
  tau.model.h <- 0 # development time
  dA.model.h <- 0 # adult mortality
  s.model.h <- 0 # juvenile survival
  
  # Sum values across each model time-step during active season
  # r_m
  for(i in 2:nrow(TS.r.h)) {
    if(TS.r.h$A[i] > 0 && TS.r.h$A[i-1] > 0 && T.h(start.r.h + TS.r.h$Time[i]) >= param$Tmin) {
      r.max.h <- max(r.max.h, log(TS.r.h$A[i]/TS.r.h$A[i-1]))
      r.model.h <- r.model.h + log(TS.r.h$A[i]/TS.r.h$A[i-1]) } # calculate sum of daily r_m
  }
  # R0 and fitness components   
  for(i in 1:nrow(TS.h)) {
    if(T.h(TS.h$Time[i]) >= param$Tmin) { # i.e., sum traits during active season
      # R0
      R0.model.h <- R0.model.h + TS.h$R0[i]
      # Fitness components (development rate is quantified below)
      b.model.h <- b.model.h + TS.h$b[i]
      dA.model.h <- dA.model.h + TS.h$dA[i]
      s.model.h <- s.model.h + TS.h$S[i]
    }
    # sum development time over all days
    tau.model.h <- tau.model.h + TS.h$tau[i]
  }
  # average values over active season
  r.model.h <- r.model.h/season.r.h
  if(param[1] == "Brevicoryne brassicae US Columbia") { r.max.h <- param$rMax } # NOTE: use rMax for Brevicoryne brassicae b/c r.max.h is biased by one extreme value 
  R0.model.h <- R0.model.h/season.h
  b.model.h <- b.model.h/season.h
  tau.model.h <- tau.model.h/(end.h-start.h)  # different denominator b/c development time does not stop during overwintering
  dA.model.h <- dA.model.h/season.h
  s.model.h <- s.model.h/season.h
  
  
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
  r.model.f <- r.model.f/season.r.h
  if(param[1] == "Brevicoryne brassicae US Columbia") { r.max.f <- param$rMax } # NOTE: use rMax for Brevicoryne brassicae b/c r.max.h is biased by one extreme value 
  R0.model.f <- R0.model.f/season.f
  b.model.f <- b.model.f/season.f
  tau.model.f <- tau.model.f/(end.f-start.f)  # different denominator b/c development time does not stop during overwintering
  dA.model.f <- dA.model.f/season.f
  if(ext == FALSE) { s.model.f <- s.model.f/season.f
  } else { s.model.f <- 0 } # if population goes extinct, set survival to zero

  
################################# RECORD AND DISPLAY RESULTS ###############################
  # Calculate optimum development rate, maximum adult longevity, and survival at mean habitat temperature for scaling results
  gTopt <- param$gTR*(param$Toptg/param$TR)*exp(param$Ag*(1/param$TR-1/param$Toptg))/(1+exp(param$AL*(1/param$TL-1/param$Toptg))) #optimum development rate
  Lmax <- 1/min(param$dATR*exp(param$AdA*(1/param$TR-1/T.h(TS.h$Time)))) # maximum adult longevity
  # g.h <- function(t) { ifelse(T.h(t) <= param$Toptg, param$gTR*(T.h(t)/param$TR)*exp(param$Ag*(1/param$TR-1/T.h(t)))/(1+exp(param$AL*(1/param$TL-1/T.h(t)))), # if T(t) < Toptg, use monotonic g(T)
  #                                                            ifelse(T.h(t) <= param$Tmaxg, param$gMax, 0)) } # If T(t) < Tmaxg, g(T) = gMax; otherwise, g(T) = 0
  # dJ.h <- function(t) { param$dJTR*exp(param$AdJ*(1/param$TR-1/T.h(t))) }
  # Smean <- exp(-dJ.h(t.param$meanT.h)*g.h(t.param$meanT.h))
  # Shot <- exp(-dJ.h(t.param$meanT.h+t.param$amplT.h)*g.h(t.param$meanT.h+t.param$amplT.h))
  # Smax <- max(exp(-dJ.h(TS.h$Time)*g.h(TS.h$Time)))
  # Smin <- min(exp(-dJ.h(TS.h$Time)*g.h(TS.h$Time)))

  # Input results into arrays
  if(all == TRUE) {
    # r_m
    r.data[s,6] <- round(r.TPC.h/param$rMax,3)
    r.data[s,7] <- round(r.TPC.f/param$rMax,3)
    r.data[s,8] <- round(r.model.h/r.max.h,3)
    r.data[s,9] <- round(r.model.f/r.max.h,3)
    r.data[s,10] <- round((r.TPC.f - r.TPC.h)/param$rMax,3)
    r.data[s,11] <- round((r.model.f - r.model.h)/r.max.h,3)
    
    # R0
    R0.data[s,6] <- round(R0.TPC.h/param$R0Topt,3)
    R0.data[s,7] <- round(R0.TPC.f/param$R0Topt,3)
    R0.data[s,8] <- round(R0.model.h/max(TS.h$R0),3)
    R0.data[s,9] <- round(R0.model.f/max(TS.h$R0),3)
    R0.data[s,10] <- round((R0.TPC.f - R0.TPC.h)/param$R0Topt,3)
    R0.data[s,11] <- round((R0.model.f - R0.model.h)/max(TS.h$R0),3)
  
    # Birth rate
    birth.data[s,6] <- round(b.TPC.h/param$bTopt,3)
    birth.data[s,7] <- round(b.TPC.f/param$bTopt,3)
    birth.data[s,8] <- round(b.model.h/param$bTopt,3)
    birth.data[s,9] <- round(b.model.f/param$bTopt,3)
    birth.data[s,10] <- round((b.TPC.f - b.TPC.h)/param$bTopt,3)
    birth.data[s,11] <- round((b.model.f - b.model.h)/param$bTopt,3)
  
    # Development
    dev.data[s,6] <- round((1/g.TPC.h)/(1/gTopt),3)
    dev.data[s,7] <- round((1/g.TPC.f)/(1/gTopt),3)
    dev.data[s,8] <- round(tau.model.h/min(TS.h$tau),3)
    dev.data[s,9] <- round(tau.model.f/min(TS.h$tau),3)
    dev.data[s,10] <- round((1/g.TPC.f - 1/g.TPC.h)/(1/gTopt),3)
    dev.data[s,11] <- round((tau.model.f - tau.model.h)/min(TS.h$tau),3)
  
    # Longevity
    long.data[s,6] <- round((1/dA.TPC.h)/Lmax,3)
    long.data[s,7] <- round((1/dA.TPC.f)/Lmax,3)
    long.data[s,8] <- round((1/dA.model.h)/(1/min(TS.h$dA)),3)
    long.data[s,9] <- round((1/dA.model.f)/(1/min(TS.h$dA)),3)
    long.data[s,10] <- round((1/dA.TPC.f - 1/dA.TPC.h)/Lmax,3)
    long.data[s,11] <- round((1/dA.model.f - 1/dA.model.h)/(1/min(TS.h$dA)),3)
  
    # Survival
    surv.data[s,6] <- round(s.TPC.h,3)
    surv.data[s,7] <- round(s.TPC.f,3)
    surv.data[s,8] <- round(s.model.h,3)
    surv.data[s,9] <- round(s.model.f,3)
    surv.data[s,10] <- round((s.TPC.f - s.TPC.h)/s.TPC.h,3)
    surv.data[s,11] <- round((s.model.f - s.model.h)/s.TPC.h,3)
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
  print(paste("r.TPC.h:", r.TPC.h/param$rMax))
  print(paste("r.TPC.f:", r.TPC.f/param$rMax))
  print(paste("r.model.h:", r.model.h/r.max.h))
  print(paste("r.model.f:", r.model.f/r.max.h))
  # R0
  print(paste("R0.TPC.h:", R0.TPC.h/param$R0Topt))
  print(paste("R0.TPC.f:", R0.TPC.f/param$R0Topt))
  print(paste("R0.model.h:", R0.model.h/max(TS.h$R0)))
  print(paste("R0.model.f:", R0.model.f/max(TS.h$R0)))
  # Birth rate
  print(paste("b.TPC.h:", b.TPC.h/param$bTopt))
  print(paste("b.TPC.f:", b.TPC.f/param$bTopt))
  print(paste("b.model.h:", b.model.h/param$bTopt))
  print(paste("b.model.f:", b.model.f/param$bTopt))
  # Development rate
  print(paste("tau.TPC.h:", (1/g.TPC.h)/(1/gTopt)))
  print(paste("tau.TPC.f:", (1/g.TPC.f)/(1/gTopt)))
  print(paste("tau.model.h:", tau.model.h/min(TS.h$tau)))
  print(paste("tau.model.f:", tau.model.f/min(TS.h$tau)))
  # Longevity"
  print(paste("1/dA.TPC.h:", (1/dA.TPC.h)/Lmax))
  print(paste("1/dA.TPC.f:", (1/dA.TPC.f)/Lmax))
  print(paste("1/dA.model.h:", (1/dA.model.h)/(1/min(TS.h$dA))))
  print(paste("1/dA.model.f:", (1/dA.model.f)/(1/min(TS.h$dA))))
  # Survival
  print(paste("s.TPC.h:", s.TPC.h))
  print(paste("s.TPC.f:", s.TPC.f))
  print(paste("s.model.h:", s.model.h))
  print(paste("s.model.f:", s.model.f))
} else {
  print(r.data)
  print(R0.data)
  print(birth.data)
  print(dev.data)
  print(long.data)
  print(surv.data)
}