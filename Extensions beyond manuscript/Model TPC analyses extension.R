############################################################################
#### This R script calculates species' traits using TPCs and DDE model #####
############################################################################

# Load packages and set working directory
library(tidyverse)
library(ggplot2)
library(cubature)
library(lamW)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('..')

# USER: choose "Fitness", "R0", "Fecundity", "Survival",
#               "Birth", "Development", "Longevity", or "Recruitment"
trait <- "Fitness" # NOTE: remember to change lines 63, 65, 68, 70, 78 if using "Fitness" or changing from "Fitness"

# USER: enter species and location or set "all" to TRUE to run analysis for all species
species <- "Clavigralla shadabi"
location <- "Benin"
all <- FALSE

# USER: include overwintering? (i.e., do not integrate over temperatures below Tmin)
overw <- TRUE

# USER: include diurnal variation?
daily <- FALSE

# USER: output results in csv (only if all == TRUE)?
output <- TRUE


# READ LIFE HISTORY AND TEMPERATURE PARAMETERS
param.all <- as.data.frame(read_csv("Model parameters/Temperature response parameters manuscript.csv"))
ifelse(daily == TRUE, t.param.all <- as.data.frame(read_csv("Model parameters/Habitat temperature parameters.csv")),
       t.param.all <- as.data.frame(read_csv("Model parameters/Habitat temperature parameters.csv")))
# Get parameters for selected species
if(all == FALSE) {
  param <- subset(as.data.frame(read_csv("Model parameters/Temperature response parameters manuscript.csv")), Species == paste(species,location))
  ifelse(daily == TRUE, t.param <- subset(as.data.frame(read_csv("Model parameters/Habitat temperature parameters.csv")), Species == paste(species,location)),
         t.param <- subset(as.data.frame(read_csv("Model parameters/Habitat temperature parameters.csv")), Species == paste(species,location)))
}


# CREATE ARRAY FOR RESULTS
if(all == TRUE) {
  results <- data.frame(param.all[,1], param.all[,2], param.all[,3], param.all[,4], param.all[,5], 1:nrow(param.all), 1:nrow(param.all), 1:nrow(param.all), 1:nrow(param.all), 1:nrow(param.all), 1:nrow(param.all))
  names(results) <- c("Species","Location","Latitude","Habitat","Subfamily","TPC.h","TPC.f","Model.h","Model.f","delta.TPC","delta.model")
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
           TS.h <- as.data.frame(read_csv(paste0("Time series data DI/Historical time series ",species," ",location,".csv"))))
    ifelse(daily == TRUE, TS.f <- as.data.frame(read_csv(paste0("Time series data Diurnal/Future time series ",species," ",location,".csv"))),
           TS.f <- as.data.frame(read_csv(paste0("Time series data DI/Future time series ",species," ",location,".csv"))))
  } else {
    ifelse(daily == TRUE, TS.h <- as.data.frame(read_csv(paste0("Time series data Diurnal/Historical time series ",param[1],".csv"))),
           TS.h <- as.data.frame(read_csv(paste0("Time series data DI/Historical time series ",param[1],".csv"))))
    ifelse(daily == TRUE, TS.f <- as.data.frame(read_csv(paste0("Time series data Diurnal/Future time series ",param[1],".csv"))),
           TS.f <- as.data.frame(read_csv(paste0("Time series data DI/Future time series ",param[1],".csv"))))
  }
  
  # If needed, remove any rows after juvenile and adult density both become zero simultaneously (i.e., population is extinct)
  if(min(TS.h$A) == 0) { TS.h <- TS.h[1:min(nrow(TS.h), TS.h[TS.h$J == 0 & TS.h$A == 0,"Time"]),] }
  if(min(TS.f$A) == 0) { TS.f <- TS.f[1:min(nrow(TS.f), TS.f[TS.f$J == 0 & TS.f$A == 0,"Time"]),] }
  
  # Define temperature function and the start and end dates for integration
  init_years <- 65 # from Python DDE model (NOTE: REMEMBER TO CHECK THIS VALUE!!!)
  # habitat temperature
  T.h <- function(t) { (t.param$meanT.h + t.param$delta_mean.h*(t+init_years*365)) - (t.param$amplT.h + t.param$delta_ampl.h*(t+init_years*365))*cos(2*pi*(t + init_years*365 + t.param$shiftT.h)/365) } # - t.param$amplD.h*cos(2*pi*t) }
  T.f <- function(t) { (t.param$meanT.f + t.param$delta_mean.f*(t+init_years*365)) - (t.param$amplT.f + t.param$delta_ampl.f*(t+init_years*365))*cos(2*pi*(t + init_years*365 + t.param$shiftT.f)/365) } # - t.param$amplD.f*cos(2*pi*t) }
  # start and end times for integration
  end.h <- TS.h[nrow(TS.h),1]
  end.f <- TS.f[nrow(TS.f),1]
  start.h <- max(end.h - 365*5 + 1, end.h %% 365) # integrate over last 5 years of time-series (max number of years before end date if <5 years in data)
  start.f <- max(end.f - 365*5 + 1, end.f %% 365) # integrate over last 5 years of time-series (max number of years before end date if <5 years in data)
  if(start.h == end.h) { start.h <- 31 } # if < 1 year, set start = 31 to avoid initial transients
  if(start.f == end.f) { start.f <- 31 } # if < 1 year, set start = 31 to avoid initial transients
  #if(trait == "Fitness") { start.h <- start.h + 2; start.f <- start.f + 2; end.h <- end.h + 1; end.f <- end.f + 1 }
  
  
  ################################## TPC: HISTORICAL CLIMATE ###################################
  # Integrate across life history traits
  if(overw == FALSE) {
    # Fitness and R0
    r.h <- function(t) { ifelse(T.h(t) <= param$Toptr, param$rMax*exp(-((T.h(t)-param$Toptr)^2)/(2*param$sr^2)),
                                                               param$rMax*(1 - ((T.h(t)-param$Toptr)/(param$Toptr-param$Tmaxr))^2)) }
    R0.h <- function(t) { param$R0Topt*exp(-((T.h(t)-param$ToptR0)^2)/(2*param$sR0^2)) }
    # Life history traits
    b.h <- function(t) { param$bTopt*exp(-((T.h(t)-param$Toptb)^2)/(2*param$sb^2)) }
    #m.h <- function(t) { param$gTR*(T.h(t)/param$TR)*exp(param$Ag*(1/param$TR-1/T.h(t)))/(1+exp(param$AL*(1/param$TL-1/T.h(t)))+exp(param$AH*(1/param$TH-1/T.h(t)))) }
    m.h <- function(t) { ifelse(T.h(t) <= param$Toptg, param$gTR*(T.h(t)/param$TR)*exp(param$Ag*(1/param$TR-1/T.h(t)))/(1+exp(param$AL*(1/param$TL-1/T.h(t)))), # if T(t) < Toptg, use monotonic mJ(T)
                                                               ifelse(T.h(t) <= param$Tmaxg, param$gMax, 0)) } # If T(t) < Tmaxg, mJ(T) = gMax; otherwise, mJ(T) = 0
    #m.h <- function(t) { param$gTR*(T.h(t)/param$TR)*exp(param$Ag*(1/param$TR-1/T.h(t))) }
    dJ.h <- function(t) { param$dJTR*exp(param$AdJ*(1/param$TR-1/T.h(t))) }
    dA.h <- function(t) { param$dATR*exp(param$AdA*(1/param$TR-1/T.h(t))) }
    # Lifetime fecundity, survivorship, and recruitment
    f.h <- function(t) { b.h(t)/dA.h(t) }
    s.h <- function(t) { exp(-dJ.h(t)/m.h(t)) }
    #R.h <- function(t) { m.h(t) * lambertW0(b.h(t)/m.h(t) * exp((dA.h(t)-dJ.h(t))/m.h(t))) }
    R.h <- function(t) { b.h(t) * s.h(t) }
        
    # Integration (Note: pcubature is faster but cannot be used with overwintering)
    # fitness, R0, lifetime fecundity, survivorship, and recruitment
    r.TPC.h <- cubintegrate(r.h, lower = start.h, upper = end.h, method = "pcubature")$integral/(end.h-start.h)
    R0.TPC.h <- cubintegrate(R0.h, lower = start.h, upper = end.h, method = "pcubature")$integral/(end.h-start.h)
    f.TPC.h <- cubintegrate(f.h, lower = start.h, upper = end.h, method = "pcubature")$integral/(end.h-start.h)
    s.TPC.h <- cubintegrate(s.h, lower = start.h, upper = end.h, method = "pcubature")$integral/(end.h-start.h)
    R.TPC.h <- cubintegrate(R.h, lower = start.h, upper = end.h, method = "pcubature")$integral/(end.h-start.h)
    # life history traits
    b.TPC.h <- cubintegrate(b.h, lower = start.h, upper = end.h, method = "pcubature")$integral/(end.h-start.h)
    m.TPC.h <- cubintegrate(m.h, lower = start.h, upper = end.h, method = "pcubature")$integral/(end.h-start.h)
    dA.TPC.h <- cubintegrate(dA.h, lower = start.h, upper = end.h, method = "pcubature")$integral/(end.h-start.h)
  }
  if(overw == TRUE) {
    # Fitness and R0
    r.h <- function(t) { ifelse(T.h(t) < param$Tmin, 0, ifelse(T.h(t) <= param$Toptr, param$rMax*exp(-((T.h(t)-param$Toptr)^2)/(2*param$sr^2)),
                                                               param$rMax*(1 - ((T.h(t)-param$Toptr)/(param$Toptr-param$Tmaxr))^2))) }
    R0.h <- function(t) { ifelse(T.h(t) < param$Tmin, 0, param$R0Topt*exp(-((T.h(t)-param$ToptR0)^2)/(2*param$sR0^2))) }
    # Life history traits
    b.h <- function(t) { ifelse(T.h(t) < param$Tmin, 0, param$bTopt*exp(-((T.h(t)-param$Toptb)^2)/(2*param$sb^2))) }
    #m.h <- function(t) { param$gTR*(T.h(t)/param$TR)*exp(param$Ag*(1/param$TR-1/T.h(t)))/(1+exp(param$AL*(1/param$TL-1/T.h(t)))+exp(param$AH*(1/param$TH-1/T.h(t)))) }
    m.h <- function(t) { ifelse(T.h(t) < param$Tmin, 0, ifelse(T.h(t) <= param$Toptg, param$gTR*(T.h(t)/param$TR)*exp(param$Ag*(1/param$TR-1/T.h(t)))/(1+exp(param$AL*(1/param$TL-1/T.h(t)))), # if T(t) < Toptg, use monotonic mJ(T)
                                                               ifelse(T.h(t) <= param$Tmaxg, param$gMax, 0))) } # If T(t) < Tmaxg, mJ(T) = gMax; otherwise, mJ(T) = 0
    #m.h <- function(t) { ifelse(T.h(t) < param$Tmin, 0, param$gTR*(T.h(t)/param$TR)*exp(param$Ag*(1/param$TR-1/T.h(t)))) }
    dJ.h <- function(t) { ifelse(T.h(t) < param$Tmin, 0, param$dJTR*exp(param$AdJ*(1/param$TR-1/T.h(t)))) }
    dA.h <- function(t) { ifelse(T.h(t) < param$Tmin, 0, param$dATR*exp(param$AdA*(1/param$TR-1/T.h(t)))) }
    # Lifetime fecundity, survivorship, and recruitment
    f.h <- function(t) { ifelse(T.h(t) < param$Tmin, 0, b.h(t)/dA.h(t)) }
    s.h <- function(t) { ifelse(T.h(t) < param$Tmin, 0, exp(-dJ.h(t)/m.h(t))) }
    #R.h <- function(t) { ifelse(T.h(t) < param$Tmin, 0, m.h(t) * lambertW0(b.h(t)/m.h(t) * exp((dA.h(t)-dJ.h(t))/m.h(t)))) }
    R.h <- function(t) { ifelse(T.h(t) < param$Tmin, 0, b.h(t) * s.h(t)) }
    
    # Integrate across active season
    season.h <- end.h - start.h # season length
    ifelse(daily == TRUE, length <- 0.5, length <- 1)
    for(t in seq(start.h,end.h,length)) { if(T.h(t) < param$Tmin) { season.h <- season.h - length }} # number of days when T(t) > Tmin
    
    # Integration (Note: hcubature must be used with overwintering)
    # fitness, R0, lifetime fecundity, survivorship, and recruitment
    r.TPC.h <- cubintegrate(r.h, lower = start.h, upper = end.h, method = "hcubature")$integral/season.h
    R0.TPC.h <- cubintegrate(R0.h, lower = start.h, upper = end.h, method = "hcubature")$integral/season.h
    f.TPC.h <- cubintegrate(f.h, lower = start.h, upper = end.h, method = "hcubature")$integral/season.h
    s.TPC.h <- cubintegrate(s.h, lower = start.h, upper = end.h, method = "hcubature")$integral/season.h
    R.TPC.h <- cubintegrate(R.h, lower = start.h, upper = end.h, method = "hcubature")$integral/season.h
    # life history traits
    b.TPC.h <- cubintegrate(b.h, lower = start.h, upper = end.h, method = "hcubature")$integral/season.h
    m.TPC.h <- cubintegrate(m.h, lower = start.h, upper = end.h, method = "hcubature")$integral/(end.h-start.h)
    dA.TPC.h <- cubintegrate(dA.h, lower = start.h, upper = end.h, method = "hcubature")$integral/season.h
  }
  
  
  ##################################### TPC: FUTURE CLIMATE ####################################
  # Integrate across life history traits
  if(overw == FALSE) {
    # Fitness and R0
    r.f <- function(t) { ifelse(T.f(t) <= param$Toptr, param$rMax*exp(-((T.f(t)-param$Toptr)^2)/(2*param$sr^2)),
                                                               param$rMax*(1 - ((T.f(t)-param$Toptr)/(param$Toptr-param$Tmaxr))^2)) }
    R0.f <- function(t) { param$R0Topt*exp(-((T.f(t)-param$ToptR0)^2)/(2*param$sR0^2)) }
    # Life history traits
    b.f <- function(t) { param$bTopt*exp(-((T.f(t)-param$Toptb)^2)/(2*param$sb^2)) }
    #m.f <- function(t) { param$gTR*(T.f(t)/param$TR)*exp(param$Ag*(1/param$TR-1/T.f(t)))/(1+exp(param$AL*(1/param$TL-1/T.f(t)))+exp(param$AH*(1/param$TH-1/T.f(t)))) }
    m.f <- function(t) { ifelse(T.f(t) <= param$Toptg, param$gTR*(T.f(t)/param$TR)*exp(param$Ag*(1/param$TR-1/T.f(t)))/(1+exp(param$AL*(1/param$TL-1/T.f(t)))), # if T(t) < Toptg, use monotonic mJ(T)
                                                               ifelse(T.f(t) <= param$Tmaxg, param$gMax, 0)) } # If T(t) < Tmaxg, mJ(T) = gMax; otherwise, mJ(T) = 0
    #m.f <- function(t) { param$gTR*(T.f(t)/param$TR)*exp(param$Ag*(1/param$TR-1/T.f(t))) }
    dJ.f <- function(t) { param$dJTR*exp(param$AdJ*(1/param$TR-1/T.f(t))) }
    dA.f <- function(t) { param$dATR*exp(param$AdA*(1/param$TR-1/T.f(t))) }
    # Lifetime fecundity, survivorship, and recruitment
    f.f <- function(t) { b.f(t)/dA.f(t) }
    s.f <- function(t) { exp(-dJ.f(t)/m.f(t)) }
    #R.f <- function(t) { m.f(t) * lambertW0(b.f(t)/m.f(t) * exp((dA.f(t)-dJ.f(t))/m.f(t))) }
    R.f <- function(t) { b.f(t) * s.f(t) }
    
    # Integration (Note: pcubature is faster but cannot be used with overwintering)
    # fitness, R0, lifetime fecundity, survivorship, and recruitment
    r.TPC.f <- cubintegrate(r.f, lower = start.f, upper = end.f, method = "pcubature")$integral/(end.f-start.f)
    R0.TPC.f <- cubintegrate(R0.f, lower = start.f, upper = end.f, method = "pcubature")$integral/(end.f-start.f)
    f.TPC.f <- cubintegrate(f.f, lower = start.f, upper = end.f, method = "pcubature")$integral/(end.f-start.f)
    s.TPC.f <- cubintegrate(s.f, lower = start.f, upper = end.f, method = "pcubature")$integral/(end.f-start.f)
    R.TPC.f <- cubintegrate(R.f, lower = start.f, upper = end.f, method = "pcubature")$integral/(end.f-start.f)
    # life history traits
    b.TPC.f <- cubintegrate(b.f, lower = start.f, upper = end.f, method = "pcubature")$integral/(end.f-start.f)
    m.TPC.f <- cubintegrate(m.f, lower = start.f, upper = end.f, method = "pcubature")$integral/(end.f-start.f)
    dA.TPC.f <- cubintegrate(dA.f, lower = start.f, upper = end.f, method = "pcubature")$integral/(end.f-start.f)
  }
  if(overw == TRUE) {
    # Fitness and R0
    r.f <- function(t) { ifelse(T.f(t) < param$Tmin, 0, ifelse(T.f(t) <= param$Toptr, param$rMax*exp(-((T.f(t)-param$Toptr)^2)/(2*param$sr^2)),
                                                               param$rMax*(1 - ((T.f(t)-param$Toptr)/(param$Toptr-param$Tmaxr))^2))) }
    R0.f <- function(t) { ifelse(T.f(t) < param$Tmin, 0, param$R0Topt*exp(-((T.f(t)-param$ToptR0)^2)/(2*param$sR0^2))) }
    # Life history traits
    b.f <- function(t) { ifelse(T.f(t) < param$Tmin, 0, param$bTopt*exp(-((T.f(t)-param$Toptb)^2)/(2*param$sb^2))) }
    #m.f <- function(t) { param$gTR*(T.f(t)/param$TR)*exp(param$Ag*(1/param$TR-1/T.f(t)))/(1+exp(param$AL*(1/param$TL-1/T.f(t)))+exp(param$AH*(1/param$TH-1/T.f(t)))) }
    m.f <- function(t) { ifelse(T.f(t) < param$Tmin, 0, ifelse(T.f(t) <= param$Toptg, param$gTR*(T.f(t)/param$TR)*exp(param$Ag*(1/param$TR-1/T.f(t)))/(1+exp(param$AL*(1/param$TL-1/T.f(t)))), # if T(t) < Toptg, use monotonic mJ(T)
                                                               ifelse(T.f(t) <= param$Tmaxg, param$gMax, 0))) } # If T(t) < Tmaxg, mJ(T) = gMax; otherwise, mJ(T) = 0
    #m.f <- function(t) { ifelse(T.f(t) < param$Tmin, 0, param$gTR*(T.f(t)/param$TR)*exp(param$Ag*(1/param$TR-1/T.f(t)))) }
    dJ.f <- function(t) { ifelse(T.f(t) < param$Tmin, 0, param$dJTR*exp(param$AdJ*(1/param$TR-1/T.f(t)))) }
    dA.f <- function(t) { ifelse(T.f(t) < param$Tmin, 0, param$dATR*exp(param$AdA*(1/param$TR-1/T.f(t)))) }
    # Lifetime fecundity, survivorship, and recruitment
    f.f <- function(t) { ifelse(T.f(t) < param$Tmin, 0, b.f(t)/dA.f(t)) }
    s.f <- function(t) { ifelse(T.f(t) < param$Tmin, 0, exp(-dJ.f(t)/m.f(t))) }
    #R.f <- function(t) { ifelse(T.f(t) < param$Tmin, 0, m.f(t) * lambertW0(b.f(t)/m.f(t) * exp((dA.f(t)-dJ.f(t))/m.f(t)))) }
    R.f <- function(t) { ifelse(T.f(t) < param$Tmin, 0, b.f(t) * s.f(t)) }
                                
    # Integrate across active season
    season.f <- end.f - start.f # season length
    ifelse(daily == TRUE, length <- 0.5, length <- 1)
    for(t in seq(start.f,end.f,length)) { if(T.f(t) < param$Tmin) {season.f <- season.f - length }} # number of days when T(t) > Tmin
    
    # Integration (Note: hcubature must be used with overwintering)
    # fitness, R0, lifetime fecundity, survivorship, and recruitment
    r.TPC.f <- cubintegrate(r.f, lower = start.f, upper = end.f, method = "hcubature")$integral/season.f
    R0.TPC.f <- cubintegrate(R0.f, lower = start.f, upper = end.f, method = "hcubature")$integral/season.f
    f.TPC.f <- cubintegrate(f.f, lower = start.f, upper = end.f, method = "hcubature")$integral/season.f
    s.TPC.f <- cubintegrate(s.f, lower = start.f, upper = end.f, method = "hcubature")$integral/season.f
    R.TPC.f <- cubintegrate(R.f, lower = start.f, upper = end.f, method = "hcubature")$integral/season.f
    # life history traits
    b.TPC.f <- cubintegrate(b.f, lower = start.f, upper = end.f, method = "hcubature")$integral/season.f
    m.TPC.f <- cubintegrate(m.f, lower = start.f, upper = end.f, method = "hcubature")$integral/(end.f-start.f)
    dA.TPC.f <- cubintegrate(dA.f, lower = start.f, upper = end.f, method = "hcubature")$integral/season.f
  }


  ################################# MODEL: HISTORICAL CLIMATE ##################################
  # Life history traits
  b.h <- function(t) { param$bTopt*exp(-((T.h(t)-param$Toptb)^2)/(2*param$sb^2)) }
  #m.h <- function(t) { param$gTR*(T.h(t)/param$TR)*exp(param$Ag*(1/param$TR-1/T.h(t)))/(1+exp(param$AL*(1/param$TL-1/T.h(t)))+exp(param$AH*(1/param$TH-1/T.h(t)))) }
  m.h <- function(t) { ifelse(T.h(t) <= param$Toptg, param$gTR*(T.h(t)/param$TR)*exp(param$Ag*(1/param$TR-1/T.h(t)))/(1+exp(param$AL*(1/param$TL-1/T.h(t)))), # if T(t) < Toptg, use monotonic mJ(T)
                              ifelse(T.h(t) <= param$Tmax, param$gTR*(param$Toptg/param$TR)*exp(param$Ag*(1/param$TR-1/param$Toptg))/(1+exp(param$AL*(1/param$TL-1/param$Toptg))), 0)) } # If T(t) < Tmax, mJ(T) = mJ(Toptg); otherwise, mJ(T) = 0
  #m.h <- function(t) { param$gTR*(T.h(t)/param$TR)*exp(param$Ag*(1/param$TR-1/T.h(t))) }
  dJ.h <- function(t) { param$dJTR*exp(param$AdJ*(1/param$TR-1/T.h(t))) }
  dA.h <- function(t) { param$dATR*exp(param$AdA*(1/param$TR-1/T.h(t))) }
  # Add R0, fecundity, and recruitment at each time-step in DDE model
  TS.h$R0 <- b.h(TS.h$Time)/dA.h(TS.h$Time)*TS.h$S
  TS.h$f <- b.h(TS.h$Time)/dA.h(TS.h$Time)
  TS.h$R <- b.h(TS.h$Time - TS.h$tau) * m.h(TS.h$Time)/m.h(TS.h$Time - TS.h$tau) * TS.h$S
  # Add birth rate and adult mortality at each time-step in DDE model
  TS.h$b <- b.h(TS.h$Time)
  TS.h$dA <- dA.h(TS.h$Time)
  
  # Integrate across data from DDE model
  # initial values
  r.model.h <- 0
  r.max.h <- 0
  R0.model.h <- 0
  f.model.h <- 0
  s.model.h <- 0
  R.model.h <- 0
  b.model.h <- 0
  tau.model.h <- 0
  dA.model.h <- 0
  # sum across model time-steps
  count.h <- 0
  for(i in start.h:end.h) {
    if(overw == FALSE || T.h(i) >= param$Tmin || trait == "Development") {
      # Fitness
      if(TS.h$A[i] > 0 && TS.h$A[i-1] > 0) {
        r.max.h <- max(r.max.h, log(TS.h$A[i]/TS.h$A[i-1]))
        r.model.h <- r.model.h + log(TS.h$A[i]/TS.h$A[i-1]) }
      # R0, fecundity, survivorship, and recruitment
      R0.model.h <- R0.model.h + TS.h$R0[i]
      f.model.h <- f.model.h + TS.h$f[i]
      s.model.h <- s.model.h + TS.h$S[i]
      R.model.h <- R.model.h + TS.h$R[i]
      # Life history traits
      b.model.h <- b.model.h + TS.h$b[i]
      tau.model.h <- tau.model.h + TS.h$tau[i]
      dA.model.h <- dA.model.h + TS.h$dA[i]
      # Active season
      count.h <- count.h + 1
    }
  }
  # average across active season
  r.model.h <- r.model.h/count.h
  if((all == FALSE && species == "Brevicoryne brassicae") || (all == TRUE && s == 25)) {r.max.h <- param$rMax} # NOTE: use rMax for Brevicoryne brassicae b/c r.max.h is biased by one extreme value 
  R0.model.h <- R0.model.h/count.h
  f.model.h <- f.model.h/count.h
  s.model.h <- s.model.h/count.h
  R.model.h <- R.model.h/count.h
  b.model.h <- b.model.h/count.h
  tau.model.h <- tau.model.h/(end.h-start.h)
  dA.model.h <- dA.model.h/count.h
  
  
  ################################### MODEL: FUTURE CLIMATE ####################################
  # Life history traits
  b.f <- function(t) { param$bTopt*exp(-((T.f(t)-param$Toptb)^2)/(2*param$sb^2)) }
  #m.f <- function(t) { param$gTR*(T.f(t)/param$TR)*exp(param$Ag*(1/param$TR-1/T.f(t)))/(1+exp(param$AL*(1/param$TL-1/T.f(t)))+exp(param$AH*(1/param$TH-1/T.f(t)))) }
  m.f <- function(t) { ifelse(T.f(t) <= param$Toptg, param$gTR*(T.f(t)/param$TR)*exp(param$Ag*(1/param$TR-1/T.f(t)))/(1+exp(param$AL*(1/param$TL-1/T.f(t)))), # if T(t) < Toptg, use monotonic mJ(T)
                              ifelse(T.f(t) <= param$Tmax, param$gTR*(param$Toptg/param$TR)*exp(param$Ag*(1/param$TR-1/param$Toptg))/(1+exp(param$AL*(1/param$TL-1/param$Toptg))), 0.0000001)) } # If T(t) < Tmax, mJ(T) = mJ(Toptg); otherwise, mJ(T) = 0
  #m.f <- function(t) { param$gTR*(T.f(t)/param$TR)*exp(param$Ag*(1/param$TR-1/T.f(t))) }
  dJ.f <- function(t) { param$dJTR*exp(param$AdJ*(1/param$TR-1/T.f(t))) }
  dA.f <- function(t) { param$dATR*exp(param$AdA*(1/param$TR-1/T.f(t))) }
  # Add R0, fecundity, and recruitment at each time-step in DDE model
  TS.f$R0 <- (b.f(TS.f$Time)/dA.f(TS.f$Time))*TS.f$S
  TS.f$f <- b.f(TS.f$Time)/dA.f(TS.f$Time)
  TS.f$R <- b.f(TS.f$Time - TS.f$tau) * m.f(TS.f$Time)/m.f(TS.f$Time - TS.f$tau) * TS.f$S
  # Add birth rate and adult mortality at each time-step in DDE model
  TS.f$b <- b.f(TS.f$Time)
  TS.f$dA <- dA.f(TS.f$Time)
  
  # Integrate across data from DDE model
  # initial values
  r.model.f <- 0
  r.max.f <- 0
  R0.model.f <- 0
  f.model.f <- 0
  s.model.f <- 0
  R.model.f <- 0
  b.model.f <- 0
  tau.model.f <- 0
  dA.model.f <- 0
  # sum across model time-steps
  count.f <- 0
  for(i in start.f:end.f) {
    if(overw == FALSE || T.f(i) >= param$Tmin || trait == "Development") {
      # Fitness
      if(TS.f$A[i] > 0 && TS.f$A[i-1] > 0) {
        r.max.f <- max(r.max.f, log(TS.f$A[i]/TS.f$A[i-1]))
        r.model.f <- r.model.f + log(TS.f$A[i]/TS.f$A[i-1]) }
      # R0, fecundity, survivorship, and recruitment
      R0.model.f <- R0.model.f + TS.f$R0[i]
      f.model.f <- f.model.f + TS.f$f[i]
      s.model.f <- s.model.f + TS.f$S[i]
      R.model.f <- R.model.f + TS.f$R[i]
      # Life history traits
      b.model.f <- b.model.f + TS.f$b[i]
      tau.model.f <- tau.model.f + TS.f$tau[i]
      dA.model.f <- dA.model.f + TS.f$dA[i]
      # Active season
      count.f <- count.f + 1
    }
  }
  # Average across active season
  r.model.f <- r.model.f/count.f
  if((all == FALSE && species == "Brevicoryne brassicae") || (all == TRUE && s == 25)) {r.max.f <- param$rMax} # NOTE: use rMax for Brevicoryne brassicae b/c r.max.f is biased by one extreme value 
  R0.model.f <- R0.model.f/count.f
  f.model.f <- f.model.f/count.f
  s.model.f <- s.model.f/count.f
  R.model.f <- R.model.f/count.f
  b.model.f <- b.model.f/count.f
  tau.model.f <- tau.model.f/(end.f-start.f)
  dA.model.f <- dA.model.f/count.f

  
  ################################# RECORD AND DISPLAY RESULTS ###############################
  # Calculate optimum and maximum values of key traits for scaling results
  mTopt <- param$gTR*(param$Toptg/param$TR)*exp(param$Ag*(1/param$TR-1/param$Toptg))/(1+exp(param$AL*(1/param$TL-1/param$Toptg)))
  lTopt <- 1/min(param$dATR*exp(param$AdA*(1/param$TR-1/T.h(TS.h$Time))))
  
  # Input results into arrays
  if(all == TRUE) {
    if(trait == "Fitness") {
      results[s,6] <- round(r.TPC.h/param$rMax,3)
      results[s,7] <- round(r.TPC.f/param$rMax,3)
      results[s,8] <- round(r.model.h/r.max.h,3)
      results[s,9] <- round(r.model.f/r.max.h,3)
      results[s,10] <- round((r.TPC.f - r.TPC.h)/param$rMax,3)
      results[s,11] <- round((r.model.f - r.model.h)/r.max.h,3)
    }
    if(trait == "R0") {
      results[s,6] <- round(R0.TPC.h/param$R0Topt,3)
      results[s,7] <- round(R0.TPC.f/param$R0Topt,3)
      results[s,8] <- round(R0.model.h/max(TS.h$R0),3)
      results[s,9] <- round(R0.model.f/max(TS.h$R0),3)
      results[s,10] <- round((R0.TPC.f - R0.TPC.h)/param$R0Topt,3)
      results[s,11] <- round((R0.model.f - R0.model.h)/max(TS.h$R0),3)
    }
    # if(trait == "Fecundity") {
    #   results[s,5] <- f.TPC.h/(param$bTopt*lTopt)
    #   results[s,6] <- f.TPC.f/(param$bTopt*lTopt)
    #   results[s,7] <- f.model.h/(param$bTopt*lTopt)
    #   results[s,8] <- f.model.f/(param$bTopt*lTopt)
    #   results[s,9] <- (f.TPC.f - f.TPC.h)/(param$bTopt*lTopt)
    #   results[s,10] <- (f.model.f - f.model.h)/(param$bTopt*lTopt)
    # }
    if(trait == "Survival") {
      results[s,6] <- round(s.TPC.h,3)
      results[s,7] <- round(s.TPC.f,3)
      results[s,8] <- round(s.model.h,3)
      results[s,9] <- round(s.model.f,3)
      results[s,10] <- round((s.TPC.f - s.TPC.h),3)
      results[s,11] <- round((s.model.f - s.model.h),3)
    }
    if(trait == "Birth") {
      results[s,6] <- round(b.TPC.h/param$bTopt,3)
      results[s,7] <- round(b.TPC.f/param$bTopt,3)
      results[s,8] <- round(b.model.h/param$bTopt,3)
      results[s,9] <- round(b.model.f/param$bTopt,3)
      results[s,10] <- round((b.TPC.f - b.TPC.h)/param$bTopt,3)
      results[s,11] <- round((b.model.f - b.model.h)/param$bTopt,3)
    }
    if(trait == "Development") {
      results[s,6] <- round((1/m.TPC.h)/(1/mTopt),3)
      results[s,7] <- round((1/m.TPC.f)/(1/mTopt),3)
      results[s,8] <- round(tau.model.h/min(TS.h$tau),3)
      results[s,9] <- round(tau.model.f/min(TS.h$tau),3)
      results[s,10] <- round((1/m.TPC.f - 1/m.TPC.h)/(1/mTopt),3)
      results[s,11] <- round((tau.model.f - tau.model.h)/min(TS.h$tau),3)
    }
    if(trait == "Longevity") {
      results[s,6] <- round((1/dA.TPC.h)/lTopt,3)
      results[s,7] <- round((1/dA.TPC.f)/lTopt,3)
      results[s,8] <- round((1/dA.model.h)/(1/min(TS.h$dA)),3)
      results[s,9] <- round((1/dA.model.f)/(1/min(TS.h$dA)),3)
      results[s,10] <- round((1/dA.TPC.f - 1/dA.TPC.h)/lTopt,3)
      results[s,11] <- round((1/dA.model.f - 1/dA.model.h)/(1/min(TS.h$dA)),3)
    }
    # if(trait == "Recruitment") {
    #   results[s,5] <- R.TPC.h/param$bTopt
    #   results[s,6] <- R.TPC.f/param$bTopt
    #   results[s,7] <- R.model.h/param$bTopt
    #   results[s,8] <- R.model.f/param$bTopt
    #   results[s,9] <- (R.TPC.f - R.TPC.h)/param$bTopt
    #   results[s,10] <- (R.model.f - R.model.h)/param$bTopt
    # }
  }

  
  # BREAK FOR LOOP IF ANALYSES ARE RUN FOR A SPECIFIED SPECIES
  if(all == FALSE) { break  }
}


# OUTPUT RESULTS IN CSV FILE
if(output == TRUE && all == TRUE) {
  if(trait == "Fitness") { write_csv(results, "Predictions/Predictions Dev fitness 2.csv") }
  if(trait == "R0") { write_csv(results, "Predictions/Predictions Dev R0 2.csv") }
  #if(trait == "Fecundity") { write_csv(results, "Predictions/Predictions Dev lifetime fecundity.csv") }
  if(trait == "Survival") { write_csv(results, "Predictions/Predictions Dev survival 2.csv") }
  if(trait == "Birth") { write_csv(results, "Predictions/Predictions Dev birth 2.csv") }
  if(trait == "Development") { write_csv(results, "Predictions/Predictions Dev development 2.csv") }
  if(trait == "Longevity") { write_csv(results, "Predictions/Predictions Dev longevity 2.csv") }
  #if(trait == "Recruitment") { write_csv(results, "Predictions/Predictions Dev recruitment.csv") }
}


###################### PLOT TEMPERATURE RESPONSE AND SUMMARIZE RESULTS ######################
# # Read in climate data
# if(all == FALSE) {
#   temp.h <- as.data.frame(read_csv(paste0("Climate data/Historical climate data ",location,".csv")))
#   temp.f <- as.data.frame(read_csv(paste0("Climate data/Future climate data ",location,".csv")))
#   
#   # For maximum temperature, remove daily minimum temperatures (if daily == FALSE)
#   #if(daily == FALSE) { temp.h <- temp.h[temp.h$day %% 1 != 0,] }
#   #if(daily == FALSE) { temp.f <- temp.f[temp.f$day %% 1 != 0,] }
#   
#   # average daily maximum and minimum temperatures (if daily == FALSE)
#   if(daily == FALSE) {
#     # historical
#     temp.h$day <- floor(temp.h$day)
#     temp.h.min <- temp.h[duplicated(temp.h$day),]
#     temp.h.max <- temp.h[duplicated(temp.h$day, fromLast=TRUE),]
#     temp.h <- data.frame(temp.h.min$day, (temp.h.min$T + temp.h.max$T)/2)
#     names(temp.h) <- c("day", "T")
#     # future
#     temp.f$day <- floor(temp.f$day)
#     temp.f.min <- temp.f[duplicated(temp.f$day),]
#     temp.f.max <- temp.f[duplicated(temp.f$day, fromLast=TRUE),]
#     temp.f <- data.frame(temp.f.min$day, (temp.f.min$T + temp.f.max$T)/2)
#     names(temp.f) <- c("day", "T")
#   }
# }
# 
# 
# # PLOT
# if(all == FALSE) {
#   # PLOT OPTIONS
#   Tmin <- round(min(temp.h$T,temp.f$T),0) - 3
#   Tmax <- round(max(temp.h$T,temp.f$T),0) + 3
#   deltaT <- 5 # tick marks for x-axis
#   ymin <- 0
#   ymax1 <- 1
#   if(trait == "Development") {ymin <- 1; ymax1 <- 5}
#   if(trait == "Longevity") {ymin <- 0; ymax1 <- 2}
#   ymax2 <- 0.4 # for temperature histogram
#   Tbin  <- 0.5 # bin size for temperature histogram
#   
#   # FUNCTIONS
#   r <- ifelse(seq(Tmin,Tmax,0.1) <= param$Toptr, param$rMax*exp(-1*((seq(Tmin,Tmax,0.1)-param$Toptr)/(2*param$sr))^2),
#               param$rMax*(1 - ((seq(Tmin,Tmax,0.1)-param$Toptr)/(param$Toptr-param$Tmaxr))^2))
#   R0 <- param$R0Topt*exp(-((seq(Tmin,Tmax,0.1)-param$ToptR0)^2)/(2*param$sR0^2))
#   b <- param$bTopt*exp(-((seq(Tmin,Tmax,0.1)-param$Toptb)^2)/(2*param$sb^2))
#   m <- param$gTR*(seq(Tmin,Tmax,0.1)/param$TR)*exp(param$Ag*(1/param$TR-1/seq(Tmin,Tmax,0.1)))/(1+exp(param$AL*(1/param$TL-1/seq(Tmin,Tmax,0.1)))+exp(param$AH*(1/param$TH-1/seq(Tmin,Tmax,0.1))))
#   dJ <- param$dJTR*exp(param$AdJ*(1/param$TR-1/seq(Tmin,Tmax,0.1)))
#   dA <- param$dATR*exp(param$AdA*(1/param$TR-1/seq(Tmin,Tmax,0.1)))
#   f <- b/dA
#   s <- exp(-dJ/m)
#   R <- b*s
#   
#   # TPC PLOTS
#   dev.new(width=3, height=3, unit="in")
#   plot(-100, xlim=c(Tmin,Tmax), ylim=c(ymin,ymax1), xaxt = "n", xlab="", ylab="", cex.axis=2)
#   axis(1, at=seq(Tmin,Tmax,deltaT), labels=seq(Tmin-273,Tmax-273,deltaT), cex.axis=2)
#   if(trait == "Fitness") {  points(seq(Tmin,Tmax,0.1), r/param$rMax, type="l", lwd=4, col="black") }
#   if(trait == "R0") {  points(seq(Tmin,Tmax,0.1), R0/param$R0Topt, type="l", lwd=4, col="black") }
#   if(trait == "Fecundity") {  points(seq(Tmin,Tmax,0.1), f/(param$bTopt*lTopt), type="l", lwd=4, col="black") }
#   if(trait == "Survival") {  points(seq(Tmin,Tmax,0.1), s, type="l", lwd=4, col="black") }
#   if(trait == "Birth") {  points(seq(Tmin,Tmax,0.1), b/param$bTopt, type="l", lwd=4, col="black") }
#   if(trait == "Development") {  points(seq(Tmin,Tmax,0.1), 1/(m/mTopt), type="l", lwd=4, col="black") }
#   if(trait == "Longevity") {  points(seq(Tmin,Tmax,0.1), (1/dA)/lTopt, type="l", lwd=4, col="black") }
#   if(trait == "Recruitment") {  points(seq(Tmin,Tmax,0.1), R/param$bTopt, type="l", lwd=4, col="black") }
#   # TEMPERATURE PARAMETERS
#   #abline(v = t.param$meanT.h, col="#0072B2", lwd=3, lty=1)
#   #abline(v = t.param$meanT.h + abs(t.param$amplT.h) + abs(t.param$amplD.h), col="#0072B2", lwd=3, lty=2)
#   #abline(v = t.param$meanT.h - abs(t.param$amplT.h) - abs(t.param$amplD.h), col="#0072B2", lwd=3, lty=2)
#   #abline(v = t.param$meanT.f + t.param$delta_mean.f*365*75 , col="#D55E00", lwd=3, lty=1)
#   #abline(v = t.param$meanT.f + t.param$delta_mean.f*365*75 + abs(t.param$amplT.f) + t.param$delta_ampl.f*365*80 + t.param$amplD.f, col="#D55E00", lwd=3, lty=2)
#   #abline(v = t.param$meanT.f + t.param$delta_mean.f*365*75 - abs(t.param$amplT.f) - t.param$delta_ampl.f*365*80 - t.param$amplD.f, col="#D55E00", lwd=3, lty=2)
#   #if(overw == TRUE) { abline(v = param$Tmin, col="black", lwd=3, lty=2) }
#   # AVERAGE VALUES
#   # Fitness (careful on which side of curve for loop starts, may have to swap Tmin and Tmax and add - to seq)
#   if(trait == "Fitness") {
#     for(i in seq(Tmax,Tmin,-0.1)) { if(ifelse(i <= param$Toptr, param$rMax*exp(-1*((i-param$Toptr)/(2*param$sr))^2),
#                                               param$rMax*(1 - ((i-param$Toptr)/(param$Toptr-param$Tmaxr))^2)) > r.TPC.h) { break } }
#     xvalue.h <- i
#     for(i in seq(Tmax,Tmin,-0.1)) { if(ifelse(i <= param$Toptr, param$rMax*exp(-1*((i-param$Toptr)/(2*param$sr))^2),
#                                               param$rMax*(1 - ((i-param$Toptr)/(param$Toptr-param$Tmaxr))^2)) > r.TPC.f) { break } }
#     xvalue.f <- i
#     points(xvalue.h, r.TPC.h/param$rMax, pch=19, cex=3, col="#0072B2")
#     points(xvalue.f, r.TPC.f/param$rMax, pch=19, cex=3, col="#D55E00") }
#   # Fecundity (careful on which side of curve for loop starts, may have to swap Tmin and Tmax and add - to seq)
#   if(trait == "Fecundity") {
#     for(i in seq(Tmax,Tmin,-0.1)) { if(exp(-((i-param$Toptb)^2)/(2*param$sb^2))/(lTopt*param$dATR*exp(param$AdA*(1/param$TR-1/i))) > f.TPC.h/(param$bTopt*lTopt)) { break } }
#     xvalue.h <- i
#     for(i in seq(Tmax,Tmin,-0.1)) { if(exp(-((i-param$Toptb)^2)/(2*param$sb^2))/(lTopt*param$dATR*exp(param$AdA*(1/param$TR-1/i))) > f.TPC.f/(param$bTopt*lTopt)) { break } }
#     xvalue.f <- i
#     points(xvalue.h, f.TPC.h/(param$bTopt*lTopt), pch=19, cex=3, col="#0072B2")
#     points(xvalue.f, f.TPC.f/(param$bTopt*lTopt), pch=19, cex=3, col="#D55E00") }
#   # Survival (careful on which side of curve for loop starts, may have to swap Tmin and Tmax and add - to seq)
#   if(trait == "Survival") {
#     for(i in seq(Tmax,Tmin,-0.1)) { if(exp(-param$dJTR*exp(param$AdJ*(1/param$TR-1/i))/(ifelse(i <= param$Toptg, param$gTR*(i/param$TR)*exp(param$Ag*(1/param$TR-1/i))/(1+exp(param$AL*(1/param$TL-1/i))),
#                                                                                                ifelse(i <= param$Tmax, param$gTR*(param$Toptg/param$TR)*exp(param$Ag*(1/param$TR-1/param$Toptg))/(1+exp(param$AL*(1/param$TL-1/param$Toptg))), 0)))) > s.TPC.h) { break } }
#     xvalue.h <- i
#     for(i in seq(Tmax,Tmin,-0.1)) { if(exp(-param$dJTR*exp(param$AdJ*(1/param$TR-1/i))/(ifelse(i <= param$Toptg, param$gTR*(i/param$TR)*exp(param$Ag*(1/param$TR-1/i))/(1+exp(param$AL*(1/param$TL-1/i))),
#                                                                                                ifelse(i <= param$Tmax, param$gTR*(param$Toptg/param$TR)*exp(param$Ag*(1/param$TR-1/param$Toptg))/(1+exp(param$AL*(1/param$TL-1/param$Toptg))), 0)))) > s.TPC.f) { break } }
#     xvalue.f <- i
#     points(xvalue.h, s.TPC.h, pch=19, cex=3, col="#0072B2")
#     points(xvalue.f, s.TPC.f, pch=19, cex=3, col="#D55E00") }
#   # TEMPERATURE HISTOGRAMS
#   par(new = T)
#   hist(temp.h$T, xlim=c(Tmin,Tmax), ylim=c(ymin,ymax2), axes=F, xlab=NA, ylab=NA, breaks=seq(from=Tmin, to=Tmax, by=Tbin), col=rgb(0,114,178, max = 255, alpha = 80), border=rgb(0,114,178, max = 255, alpha = 80), freq=FALSE, main = NULL)
#   hist(temp.f[temp.f$day>365*65,"T"], xlim=c(Tmin,Tmax), ylim=c(ymin,ymax2), breaks=seq(from=Tmin, to=Tmax, by=Tbin), ylab="r", col=rgb(213,94,0, max = 255, alpha = 80), border=rgb(213,94,0, max = 255, alpha = 80), freq=FALSE, main = NULL, add=TRUE)
#   axis(side = 4, cex.axis=2)
# }


# SUMMARIZE RESULTS
if(trait == "Fitness") { 
  print(paste("r.TPC.h:", r.TPC.h/param$rMax))
  print(paste("r.TPC.f:", r.TPC.f/param$rMax))
  print(paste("r.model.h:", r.model.h/r.max.h))
  print(paste("r.model.f:", r.model.f/r.max.h))
}
if(trait == "R0") { 
  print(paste("R0.TPC.h:", R0.TPC.h/param$R0Topt))
  print(paste("R0.TPC.f:", R0.TPC.f/param$R0Topt))
  print(paste("R0.model.h:", R0.model.h/max(TS.h$R0)))
  print(paste("R0.model.f:", R0.model.f/max(TS.h$R0)))
}
# if(trait == "Fecundity") { 
#   print(paste("f.TPC.h:", f.TPC.h/(param$bTopt*lTopt)))
#   print(paste("f.TPC.f:", f.TPC.f/(param$bTopt*lTopt)))
#   print(paste("f.model.h:", f.model.h/(param$bTopt*lTopt)))
#   print(paste("f.model.f:", f.model.f/(param$bTopt*lTopt)))
# }
if(trait == "Survival") { 
  print(paste("s.TPC.h:", s.TPC.h))
  print(paste("s.TPC.f:", s.TPC.f))
  print(paste("s.model.h:", s.model.h))
  print(paste("s.model.f:", s.model.f))
}
if(trait == "Birth") { 
  print(paste("b.TPC.h:", b.TPC.h/param$bTopt))
  print(paste("b.TPC.f:", b.TPC.f/param$bTopt))
  print(paste("b.model.h:", b.model.h/param$bTopt))
  print(paste("b.model.f:", b.model.f/param$bTopt))
}
if(trait == "Development") { 
  print(paste("tau.TPC.h:", (1/m.TPC.h)/(1/mTopt)))
  print(paste("tau.TPC.f:", (1/m.TPC.f)/(1/mTopt)))
  print(paste("tau.model.h:", tau.model.h/min(TS.h$tau)))
  print(paste("tau.model.f:", tau.model.f/min(TS.h$tau)))
}
if(trait == "Longevity") { 
  print(paste("1/dA.TPC.h:", (1/dA.TPC.h)/lTopt))
  print(paste("1/dA.TPC.f:", (1/dA.TPC.f)/lTopt))
  print(paste("1/dA.model.h:", (1/dA.model.h)/(1/min(TS.h$dA))))
  print(paste("1/dA.model.f:", (1/dA.model.f)/(1/min(TS.h$dA))))
}
# if(trait == "Recruitment") { 
#   print(paste("R.TPC.h:", R.TPC.h/param$bTopt))
#   print(paste("R.TPC.f:", R.TPC.f/param$bTopt))
#   print(paste("R.model.h:", R.model.h/param$bTopt))
#   print(paste("R.model.f:", R.model.f/param$bTopt))
# }
if(all == TRUE) { print(results) }

# Plot changes in life history traits
#barplot(c((b.TPC.f-b.TPC.h), (b.model.f-b.model.h)), col=c("Darkgreen","Orange"), ylim=c(-0.4,0.6), main=expression("Change in r"))
