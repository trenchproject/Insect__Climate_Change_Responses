############################################################################
#### This R script calculates species' traits using TPCs and DDE model #####
############################################################################

# Load packages and set working directory
library(tidyr)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(cubature)
library(lamW)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# USER: choose "Fitness", "R0", "Fecundity", "Survival",
#               "Birth", "Development", "Longevity", or "Recruitment"
trait <- "Fitness"

# USER: enter species and location or set "all" to TRUE to run analysis for all species
species <- "Macrosiphum euphorbiae"
location <- "Canada"
all <- FALSE

# USER: include overwintering? (i.e., do not integrate over temperatures below Tmin)
overw <- TRUE

# USER: include diurnal variation?
daily <- FALSE

# USER: output results in csv (only if all == TRUE)?
output <- TRUE


# READ LIFE HISTORY AND TEMPERATURE PARAMETERS
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
  results <- data.frame(param.all[,1], param.all[,2], param.all[,3], param.all[,4], 1:nrow(param.all), 1:nrow(param.all), 1:nrow(param.all), 1:nrow(param.all), 1:nrow(param.all), 1:nrow(param.all), 1:nrow(param.all), 1:nrow(param.all))
  names(results) <- c("Species","Latitude","Habitat","Subfamily","TPC.h","TPC.f","Model.h","Model.f","max.h","max.f","delta.TPC","delta.model")
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
    ifelse(daily == TRUE, TS.h <- as.data.frame(read_csv(paste0("Time series data/Historical time series ",species," ",location,".csv"))),
           TS.h <- as.data.frame(read_csv(paste0("Time series data DI Tave Dev/Historical time series ",species," ",location,".csv"))))
    ifelse(daily == TRUE, TS.f <- as.data.frame(read_csv(paste0("Time series data/Future time series ",species," ",location,".csv"))),
           TS.f <- as.data.frame(read_csv(paste0("Time series data DI Tave Dev/Future time series ",species," ",location,".csv"))))
  }
  if(all == TRUE) {
    ifelse(daily == TRUE, TS.h <- as.data.frame(read_csv(paste0("Time series data/Historical time series ",param[1],".csv"))),
           TS.h <- as.data.frame(read_csv(paste0("Time series data DI Tave Dev/Historical time series ",param[1],".csv"))))
    ifelse(daily == TRUE, TS.f <- as.data.frame(read_csv(paste0("Time series data/Future time series ",param[1],".csv"))),
           TS.f <- as.data.frame(read_csv(paste0("Time series data DI Tave Dev/Future time series ",param[1],".csv"))))
  }
  
  # Define temperature function and the start and end dates for integration
  init_years <- 65 # from Python DDE model
  # habitat temperature
  T.h <- function(t) { (t.param$meanT.h + t.param$delta_mean.h*(t+init_years*365)) - (t.param$amplT.h + t.param$delta_ampl.h*(t+init_years*365))*cos(2*pi*(t + t.param$shiftT.h)/365) - t.param$amplD.h*cos(2*pi*t) }
  T.f <- function(t) { (t.param$meanT.f + t.param$delta_mean.f*(t+init_years*365)) - (t.param$amplT.f + t.param$delta_ampl.f*(t+init_years*365))*cos(2*pi*(t + t.param$shiftT.f)/365) - t.param$amplD.f*cos(2*pi*t) }
  # start and end times for integration
  end.h <- nrow(TS.h)
  end.f <- nrow(TS.f)
  start.h <- max(end.h - 365*5 + 1, end.h %% 365) # integrate over last 5 years of time-series (max number of years before end date if <5 years in data)
  start.f <- max(end.f - 365*5 + 1, end.f %% 365) # integrate over last 5 years of time-series (max number of years before end date if <5 years in data)
  if(start.h == end.h) { start.h <- 31 } # if < 1 year, set start = 31 to avoid initial transients
  if(start.f == end.f) { start.f <- 31 } # if < 1 year, set start = 31 to avoid initial transients
  
  
  ################################## TPC: HISTORICAL CLIMATE ###################################
  # Integrate across life history traits
  if(overw == FALSE) {
    # Fitness
    r.h <- function(t) { ifelse(T.h(t) <= param$rTopt, param$rMax*exp(-1*((T.h(t)-param$rTopt)/(2*param$rs))^2),
                                param$rMax*(1 - ((T.h(t)-param$rTopt)/(param$rTopt-param$rTmax))^2)) }
    # Life history traits
    b.h <- function(t) { param$bTopt*exp(-((T.h(t)-param$Toptb)^2)/(2*param$sb^2)) }
    #m.h <- function(t) { param$mTR*(T.h(t)/param$TR)*exp(param$AmJ*(1/param$TR-1/T.h(t)))/(1+exp(param$AL*(1/param$TL-1/T.h(t)))+exp(param$AH*(1/param$TH-1/T.h(t)))) }
    m.h <- function(t) { ifelse(T.h(t) <= param$Topt, param$mTR*(T.h(t)/param$TR)*exp(param$AmJ*(1/param$TR-1/T.h(t)))/(1+exp(param$AL*(1/param$TL-1/T.h(t)))), # if T(t) < Topt, use monotonic mJ(T)
                                ifelse(T.h(t) <= param$Tmax, param$mTR*(param$Topt/param$TR)*exp(param$AmJ*(1/param$TR-1/param$Topt))/(1+exp(param$AL*(1/param$TL-1/param$Topt))), 0)) } # If T(t) < Tmax, mJ(T) = mJ(Topt); otherwise, mJ(T) = 0
    dJ.h <- function(t) { param$dJTR*exp(param$AdJ*(1/param$TR-1/T.h(t))) }
    dA.h <- function(t) { param$dATR*exp(param$AdA*(1/param$TR-1/T.h(t))) }
    # R0, lifetime fecundity, survivorship, and recruitment
    R0.h <- function(t) { b.h(t)/dA.h(t) * m.h(t)/(m.h(t)+dJ.h(t)) }
    f.h <- function(t) { b.h(t)/dA.h(t) }
    s.h <- function(t) { exp(-dJ.h(t)/m.h(t)) }
    #R.h <- function(t) { m.h(t) * lambertW0(b.h(t)/m.h(t) * exp((dA.h(t)-dJ.h(t))/m.h(t))) }
    R.h <- function(t) { b.h(t) * s.h(t) }
        
    # Integration (Note: pcubature is faster but cannot be used with overwintering)
    # fitness, R0, lifetime fecundity, survivorship, and recruitment
    r.TPC.h <- cubintegrate(r.h, lower = start.h, upper = end.h, method = "pcubature")$integral/season.h
    R0.TPC.h <- cubintegrate(R0.h, lower = start.h, upper = end.h, method = "pcubature")$integral/season.h
    f.TPC.h <- cubintegrate(f.h, lower = start.h, upper = end.h, method = "pcubature")$integral/season.h
    s.TPC.h <- cubintegrate(s.h, lower = start.h, upper = end.h, method = "pcubature")$integral/season.h
    R.TPC.h <- cubintegrate(R.h, lower = start.h, upper = end.h, method = "pcubature")$integral/season.h
    # life history traits
    b.TPC.h <- cubintegrate(b.h, lower = start.h, upper = end.h, method = "pcubature")$integral/season.h
    m.TPC.h <- cubintegrate(m.h, lower = start.h, upper = end.h, method = "pcubature")$integral/season.h
    dA.TPC.h <- cubintegrate(dA.h, lower = start.h, upper = end.h, method = "pcubature")$integral/season.h
  }
  if(overw == TRUE) {
    # Fitness
    r.h <- function(t) { ifelse(T.h(t) < param$Tmin, 0, ifelse(T.h(t) <= param$rTopt, param$rMax*exp(-1*((T.h(t)-param$rTopt)/(2*param$rs))^2),
                    param$rMax*(1 - ((T.h(t)-param$rTopt)/(param$rTopt-param$rTmax))^2))) }
    # Life history traits
    b.h <- function(t) { ifelse(T.h(t) < param$Tmin, 0, param$bTopt*exp(-((T.h(t)-param$Toptb)^2)/(2*param$sb^2))) }
    m.h <- function(t) { ifelse(T.h(t) < param$Tmin, 0, param$mTR*(T.h(t)/param$TR)*exp(param$AmJ*(1/param$TR-1/T.h(t)))/(1+exp(param$AL*(1/param$TL-1/T.h(t)))+exp(param$AH*(1/param$TH-1/T.h(t))))) }
    dJ.h <- function(t) { ifelse(T.h(t) < param$Tmin, 0, param$dJTR*exp(param$AdJ*(1/param$TR-1/T.h(t)))) }
    dA.h <- function(t) { ifelse(T.h(t) < param$Tmin, 0, param$dATR*exp(param$AdA*(1/param$TR-1/T.h(t)))) }
    # R0, lifetime fecundity, survivorship, and recruitment
    R0.h <- function(t) { ifelse(T.h(t) < param$Tmin, 0, b.h(t)/dA.h(t) * m.h(t)/(m.h(t)+dJ.h(t))) }
    f.h <- function(t) { ifelse(T.h(t) < param$Tmin, 0, b.h(t)/dA.h(t)) }
    s.h <- function(t) { ifelse(T.h(t) < param$Tmin, 0, exp(-dJ.h(t)/m.h(t))) }
    #R.h <- function(t) { ifelse(T.h(t) < param$Tmin, 0, m.h(t) * lambertW0(b.h(t)/m.h(t) * exp((dA.h(t)-dJ.h(t))/m.h(t)))) }
    R.h <- function(t) { ifelse(T.h(t) < param$Tmin, 0, b.h(t) * s.h(t)) }
    
    # Integrate across active season
    season.h <- end.h - start.h # season length
    ifelse(daily == TRUE, length <- 0.5, length <- 1)
    for(t in seq(start.h,end.h,length)) { if(T.h(t) < param$Tmin) {season.h <- season.h - length }} # number of days when T(t) > Tmin
    
    # Integration (Note: hcubature must be used with overwintering)
    # fitness, R0, lifetime fecundity, survivorship, and recruitment
    r.TPC.h <- cubintegrate(r.h, lower = start.h, upper = end.h, method = "hcubature")$integral/season.h
    R0.TPC.h <- cubintegrate(R0.h, lower = start.h, upper = end.h, method = "hcubature")$integral/season.h
    f.TPC.h <- cubintegrate(f.h, lower = start.h, upper = end.h, method = "hcubature")$integral/season.h
    s.TPC.h <- cubintegrate(s.h, lower = start.h, upper = end.h, method = "hcubature")$integral/season.h
    R.TPC.h <- cubintegrate(R.h, lower = start.h, upper = end.h, method = "hcubature")$integral/season.h
    # life history traits
    b.TPC.h <- cubintegrate(b.h, lower = start.h, upper = end.h, method = "hcubature")$integral/season.h
    m.TPC.h <- cubintegrate(m.h, lower = start.h, upper = end.h, method = "hcubature")$integral/season.h
    dA.TPC.h <- cubintegrate(dA.h, lower = start.h, upper = end.h, method = "hcubature")$integral/season.h
  }
  
  
  ##################################### TPC: FUTURE CLIMATE ####################################
  # Integrate across life history traits
  if(overw == FALSE) {
    # Fitness
    r.f <- function(t) { ifelse(T.f(t) <= param$rTopt, param$rMax*exp(-1*((T.f(t)-param$rTopt)/(2*param$rs))^2),
                                param$rMax*(1 - ((T.f(t)-param$rTopt)/(param$rTopt-param$rTmax))^2)) }
    # Life history traits
    b.f <- function(t) { param$bTopt*exp(-((T.f(t)-param$Toptb)^2)/(2*param$sb^2)) }
    #m.f <- function(t) { param$mTR*(T.f(t)/param$TR)*exp(param$AmJ*(1/param$TR-1/T.f(t)))/(1+exp(param$AL*(1/param$TL-1/T.f(t)))+exp(param$AH*(1/param$TH-1/T.f(t)))) }
    m.f <- function(t) { ifelse(T.f(t) <= param$Topt, param$mTR*(T.f(t)/param$TR)*exp(param$AmJ*(1/param$TR-1/T.f(t)))/(1+exp(param$AL*(1/param$TL-1/T.f(t)))), # if T(t) < Topt, use monotonic mJ(T)
                                ifelse(T.f(t) <= param$Tmax, param$mTR*(param$Topt/param$TR)*exp(param$AmJ*(1/param$TR-1/param$Topt))/(1+exp(param$AL*(1/param$TL-1/param$Topt))), 0)) } # If T(t) < Tmax, mJ(T) = mJ(Topt); otherwise, mJ(T) = 0
    dJ.f <- function(t) { param$dJTR*exp(param$AdJ*(1/param$TR-1/T.f(t))) }
    dA.f <- function(t) { param$dATR*exp(param$AdA*(1/param$TR-1/T.f(t))) }
    # R0, lifetime fecundity, survivorship, and recruitment
    R0.f <- function(t) { b.f(t)/dA.f(t) * m.f(t)/(m.f(t)+dJ.f(t)) }
    f.f <- function(t) { b.f(t)/dA.f(t) }
    s.f <- function(t) { exp(-dJ.f(t)/m.f(t)) }
    #R.f <- function(t) { m.f(t) * lambertW0(b.f(t)/m.f(t) * exp((dA.f(t)-dJ.f(t))/m.f(t))) }
    R.f <- function(t) { b.f(t) * s.f(t) }
    
    # Integration (Note: pcubature is faster but cannot be used with overwintering)
    # fitness, R0, lifetime fecundity, survivorship, and recruitment
    r.TPC.f <- cubintegrate(r.f, lower = start.f, upper = end.f, method = "pcubature")$integral/season.f
    R0.TPC.f <- cubintegrate(R0.f, lower = start.f, upper = end.f, method = "pcubature")$integral/season.f
    f.TPC.f <- cubintegrate(f.f, lower = start.f, upper = end.f, method = "pcubature")$integral/season.f
    s.TPC.f <- cubintegrate(s.f, lower = start.f, upper = end.f, method = "pcubature")$integral/season.f
    R.TPC.f <- cubintegrate(R.f, lower = start.f, upper = end.f, method = "pcubature")$integral/season.f
    # life history traits
    b.TPC.f <- cubintegrate(b.f, lower = start.f, upper = end.f, method = "pcubature")$integral/season.f
    m.TPC.f <- cubintegrate(m.f, lower = start.f, upper = end.f, method = "pcubature")$integral/season.f
    dA.TPC.f <- cubintegrate(dA.f, lower = start.f, upper = end.f, method = "pcubature")$integral/season.f
  }
  if(overw == TRUE) {
    # Fitness
    r.f <- function(t) { ifelse(T.f(t) < param$Tmin, 0, ifelse(T.f(t) <= param$rTopt, param$rMax*exp(-1*((T.f(t)-param$rTopt)/(2*param$rs))^2),
                                                               param$rMax*(1 - ((T.f(t)-param$rTopt)/(param$rTopt-param$rTmax))^2))) }
    # Life history traits
    b.f <- function(t) { ifelse(T.f(t) < param$Tmin, 0, param$bTopt*exp(-((T.f(t)-param$Toptb)^2)/(2*param$sb^2))) }
    m.f <- function(t) { ifelse(T.f(t) < param$Tmin, 0, param$mTR*(T.f(t)/param$TR)*exp(param$AmJ*(1/param$TR-1/T.f(t)))/(1+exp(param$AL*(1/param$TL-1/T.f(t)))+exp(param$AH*(1/param$TH-1/T.f(t))))) }
    dJ.f <- function(t) { ifelse(T.f(t) < param$Tmin, 0, param$dJTR*exp(param$AdJ*(1/param$TR-1/T.f(t)))) }
    dA.f <- function(t) { ifelse(T.f(t) < param$Tmin, 0, param$dATR*exp(param$AdA*(1/param$TR-1/T.f(t)))) }
    # R0, lifetime fecundity, survivorship, and recruitment
    R0.f <- function(t) { ifelse(T.f(t) < param$Tmin, 0, b.f(t)/dA.f(t) * m.f(t)/(m.f(t)+dJ.f(t))) }
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
    m.TPC.f <- cubintegrate(m.f, lower = start.f, upper = end.f, method = "hcubature")$integral/season.f
    dA.TPC.f <- cubintegrate(dA.f, lower = start.f, upper = end.f, method = "hcubature")$integral/season.f
  }


  ################################# MODEL: HISTORICAL CLIMATE ##################################
  # Life history traits
  b.h <- function(t) { param$bTopt*exp(-((T.h(t)-param$Toptb)^2)/(2*param$sb^2)) }
  #m.h <- function(t) { param$mTR*(T.h(t)/param$TR)*exp(param$AmJ*(1/param$TR-1/T.h(t)))/(1+exp(param$AL*(1/param$TL-1/T.h(t)))+exp(param$AH*(1/param$TH-1/T.h(t)))) }
  m.h <- function(t) { ifelse(T.h(t) <= param$Topt, param$mTR*(T.h(t)/param$TR)*exp(param$AmJ*(1/param$TR-1/T.h(t)))/(1+exp(param$AL*(1/param$TL-1/T.h(t)))), # if T(t) < Topt, use monotonic mJ(T)
                              ifelse(T.h(t) <= param$Tmax, param$mTR*(param$Topt/param$TR)*exp(param$AmJ*(1/param$TR-1/param$Topt))/(1+exp(param$AL*(1/param$TL-1/param$Topt))), 0)) } # If T(t) < Tmax, mJ(T) = mJ(Topt); otherwise, mJ(T) = 0
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
    if(overw == FALSE || (overw == TRUE && T.h(i) >= param$Tmin)) {
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
  R0.model.h <- R0.model.h/count.h
  f.model.h <- f.model.h/count.h
  s.model.h <- s.model.h/count.h
  R.model.h <- R.model.h/count.h
  b.model.h <- b.model.h/count.h
  tau.model.h <- tau.model.h/count.h
  dA.model.h <- dA.model.h/count.h

  
  ################################### MODEL: FUTURE CLIMATE ####################################
  # Life history traits
  b.f <- function(t) { param$bTopt*exp(-((T.f(t)-param$Toptb)^2)/(2*param$sb^2)) }
  #m.f <- function(t) { param$mTR*(T.f(t)/param$TR)*exp(param$AmJ*(1/param$TR-1/T.f(t)))/(1+exp(param$AL*(1/param$TL-1/T.f(t)))+exp(param$AH*(1/param$TH-1/T.f(t)))) }
  m.f <- function(t) { ifelse(T.f(t) <= param$Topt, param$mTR*(T.f(t)/param$TR)*exp(param$AmJ*(1/param$TR-1/T.f(t)))/(1+exp(param$AL*(1/param$TL-1/T.f(t)))), # if T(t) < Topt, use monotonic mJ(T)
                              ifelse(T.f(t) <= param$Tmax, param$mTR*(param$Topt/param$TR)*exp(param$AmJ*(1/param$TR-1/param$Topt))/(1+exp(param$AL*(1/param$TL-1/param$Topt))), 0)) } # If T(t) < Tmax, mJ(T) = mJ(Topt); otherwise, mJ(T) = 0
  dJ.f <- function(t) { param$dJTR*exp(param$AdJ*(1/param$TR-1/T.f(t))) }
  dA.f <- function(t) { param$dATR*exp(param$AdA*(1/param$TR-1/T.f(t))) }
  # Add R0, fecundity, and recruitment at each time-step in DDE model
  TS.f$R0 <- b.f(TS.f$Time)/dA.f(TS.f$Time)*TS.f$S
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
    if(overw == FALSE || (overw == TRUE && T.f(i) >= param$Tmin)) {
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
  # average across active season
  r.model.f <- r.model.f/count.f
  R0.model.f <- R0.model.f/count.f
  f.model.f <- f.model.f/count.f
  s.model.f <- s.model.f/count.f
  R.model.f <- R.model.f/count.f
  b.model.f <- b.model.f/count.f
  tau.model.f <- tau.model.f/count.f
  dA.model.f <- dA.model.f/count.f

  
  ################################# RECORD AND DISPLAY RESULTS ###############################
  # INPUT RESUTS INTO ARRAY
  if(all == TRUE) {
    if(trait == "Fitness") {
      results[s,5] <- r.TPC.h/param$rMax
      results[s,6] <- r.TPC.f/param$rMax
      results[s,7] <- r.model.h/param$rMax
      results[s,8] <- r.model.f/param$rMax
      results[s,9] <- r.max.h
      results[s,10] <- r.max.f
      results[s,11] <- (r.TPC.f - r.TPC.h)/param$rMax
      results[s,12] <- (r.model.f - r.model.h)/param$rMax
    }
    if(trait == "R0") {
      results[s,5] <- R0.TPC.h
      results[s,6] <- R0.TPC.f
      results[s,7] <- R0.model.h
      results[s,8] <- R0.model.f
      results[s,9] <- max(TS.h[-c(1:start.h),"R0"])
      results[s,10] <- max(TS.f[-c(1:start.f),"R0"])
      results[s,11] <- (R0.TPC.f - R0.TPC.h)/results[s,9]
      results[s,12] <- (R0.model.f - R0.model.h)/results[s,9]
    }
    if(trait == "Fecundity") {
      results[s,5] <- f.TPC.h
      results[s,6] <- f.TPC.f
      results[s,7] <- f.model.h
      results[s,8] <- f.model.f
      results[s,9] <- max(TS.h[-c(1:start.h),"f"])
      results[s,10] <- max(TS.f[-c(1:start.f),"f"])
      results[s,11] <- (f.TPC.f - f.TPC.h)/results[s,9]
      results[s,12] <- (f.model.f - f.model.h)/results[s,9]
    }
    if(trait == "Survival") {
      results[s,5] <- s.TPC.h
      results[s,6] <- s.TPC.f
      results[s,7] <- s.model.h
      results[s,8] <- s.model.f
      results[s,9] <- max(TS.h[-c(1:start.h),"S"])
      results[s,10] <- max(TS.f[-c(1:start.f),"S"])
      results[s,11] <- s.TPC.f - s.TPC.h
      results[s,12] <- s.model.f - s.model.h
    }
    if(trait == "Birth") {
      results[s,5] <- b.TPC.h/param$bTopt
      results[s,6] <- b.TPC.f/param$bTopt
      results[s,7] <- b.model.h/param$bTopt
      results[s,8] <- b.model.f/param$bTopt
      results[s,9] <- max(TS.h[-c(1:start.h),"b"])
      results[s,10] <- max(TS.f[-c(1:start.f),"b"])
      results[s,11] <- (b.TPC.f - b.TPC.h)/param$bTopt
      results[s,12] <- (b.model.f - b.model.h)/param$bTopt
    }
    if(trait == "Development") {
      results[s,5] <- 1/m.TPC.h
      results[s,6] <- 1/m.TPC.f
      results[s,7] <- tau.model.h
      results[s,8] <- tau.model.f
      results[s,9] <- max(TS.h[-c(1:start.h),"tau"])
      results[s,10] <- max(TS.f[-c(1:start.f),"tau"])
      results[s,11] <- (1/m.TPC.f - 1/m.TPC.h)/results[s,9]
      results[s,12] <- (tau.model.f - tau.model.h)/results[s,9]
    }
    if(trait == "Longevity") {
      results[s,5] <- 1/dA.TPC.h
      results[s,6] <- 1/dA.TPC.f
      results[s,7] <- 1/dA.model.h
      results[s,8] <- 1/dA.model.f
      results[s,9] <- 1/min(TS.h[-c(1:start.h),"dA"])
      results[s,10] <- 1/min(TS.f[-c(1:start.f),"dA"])
      results[s,11] <- (1/dA.TPC.f - 1/dA.TPC.h)/results[s,9]
      results[s,12] <- (1/dA.model.f - 1/dA.model.h)/results[s,9]
    }
    if(trait == "Recruitment") {
      results[s,5] <- R.TPC.h/param$bTopt
      results[s,6] <- R.TPC.f/param$bTopt
      results[s,7] <- R.model.h/param$bTopt
      results[s,8] <- R.model.f/param$bTopt
      results[s,9] <- max(TS.h[-c(1:start.h),"R"])
      results[s,10] <- max(TS.f[-c(1:start.f),"R"])
      results[s,11] <- (R.TPC.f - R.TPC.h)/param$bTopt
      results[s,12] <- (R.model.f - R.model.h)/param$bTopt
    }
  }

  
  # BREAK FOR LOOP IF ANALYSES ARE RUN FOR A SPECIFIED SPECIES
  if(all == FALSE) { break  }
}


# OUTPUT RESULTS IN CSV FILE
if(output == TRUE && all == TRUE) {
  if(trait == "Fitness") { write_csv(results, "Predictions/Predictions Dev fitness.csv") }
  if(trait == "R0") { write_csv(results, "Predictions/Predictions Dev R0.csv") }
  if(trait == "Fecundity") { write_csv(results, "Predictions/Predictions Dev fecundity.csv") }
  if(trait == "Survival") { write_csv(results, "Predictions/Predictions Dev survival.csv") }
  if(trait == "Birth") { write_csv(results, "Predictions/Predictions Dev birth.csv") }
  if(trait == "Development") { write_csv(results, "Predictions/Predictions Dev development.csv") }
  if(trait == "Longevity") { write_csv(results, "Predictions/Predictions Dev longevity.csv") }
  if(trait == "Recruitment") { write_csv(results, "Predictions/Predictions Dev recruitment.csv") }
}


###################### PLOT TEMPERATURE RESPONSE AND SUMMARIZE RESULTS ######################
# Read in climate data
if(all == FALSE) {
  temp.h <- as.data.frame(read_csv(paste0("Climate data/Historical climate data ",location,".csv")))
  temp.f <- as.data.frame(read_csv(paste0("Climate data/Future climate data ",location,".csv")))
  
  # For maximum temperature, remove daily minimum temperatures (if daily == FALSE)
  #if(daily == FALSE) { temp.h <- temp.h[temp.h$day %% 1 != 0,] }
  #if(daily == FALSE) { temp.f <- temp.f[temp.f$day %% 1 != 0,] }
  
  # average daily maximum and minimum temperatures (if daily == FALSE)
  if(daily == FALSE) {
    # historical
    temp.h$day <- floor(temp.h$day)
    temp.h.min <- temp.h[duplicated(temp.h$day),]
    temp.h.max <- temp.h[duplicated(temp.h$day, fromLast=TRUE),]
    temp.h <- data.frame(temp.h.min$day, (temp.h.min$T + temp.h.max$T)/2)
    names(temp.h) <- c("day", "T")
    # future
    temp.f$day <- floor(temp.f$day)
    temp.f.min <- temp.f[duplicated(temp.f$day),]
    temp.f.max <- temp.f[duplicated(temp.f$day, fromLast=TRUE),]
    temp.f <- data.frame(temp.f.min$day, (temp.f.min$T + temp.f.max$T)/2)
    names(temp.f) <- c("day", "T")
  }
}

# Plot
if(all == FALSE) {
  # Plot options
  Tmin <- round(min(temp.h$T,temp.f$T),0) - 3
  Tmax <- round(max(temp.h$T,temp.f$T),0) + 3
  ymin <- 0
  if(trait == "Fitness") { ymax1 <- round(param$rMax,1) }
  if(trait == "R0") { ymax1 <- round(param$bTopt/param$dATR,0) }
  if(trait == "Fecundity") { ymax1 <- round(param$bTopt/param$dATR,0) }
  if(trait == "Survival") { ymax1 <- 1 }
  if(trait == "Birth") { ymax1 <- round(param$bTopt,2) + 0.1 }
  if(trait == "Development") { ymax1 <- round(param$mTR,2) + 0.1 }
  if(trait == "Longevity") { ymax1 <- 2*round(1/param$dATR,1) }
  if(trait == "Recruitment") { ymax1 <- round(param$bTopt,2) + 0.1 }
  ymax2 <- 0.1 # for temperature histogram
  # Functions
  r <- ifelse(seq(Tmin,Tmax,1) <= param$rTopt, param$rMax*exp(-1*((seq(Tmin,Tmax,1)-param$rTopt)/(2*param$rs))^2),
              param$rMax*(1 - ((seq(Tmin,Tmax,1)-param$rTopt)/(param$rTopt-param$rTmax))^2))
  b <- param$bTopt*exp(-((seq(Tmin,Tmax,1)-param$Toptb)^2)/(2*param$sb^2))
  m <- param$mTR*(seq(Tmin,Tmax,1)/param$TR)*exp(param$AmJ*(1/param$TR-1/seq(Tmin,Tmax,1)))/(1+exp(param$AL*(1/param$TL-1/seq(Tmin,Tmax,1)))+exp(param$AH*(1/param$TH-1/seq(Tmin,Tmax,1))))
  dJ <- param$dJTR*exp(param$AdJ*(1/param$TR-1/seq(Tmin,Tmax,1)))
  dA <- param$dATR*exp(param$AdA*(1/param$TR-1/seq(Tmin,Tmax,1)))
  R0 <- b/dA * m/(m+dJ)
  f <- b/dA
  s <- exp(-dJ/m)
  #R <- m * lambertW0(b/m * exp((dA-dJ)/m))
  R <- b*s
  # TPC plots
  # Fitness
  plot(seq(Tmin,Tmax,1), r, type="l", lwd=4, col="black", xlim=c(Tmin,Tmax), ylim=c(ymin,ymax1), xlab="T", ylab="r(T)")
  # R0
  #plot(seq(Tmin,Tmax,1), R0, type="l", lwd=4, col="black", xlim=c(Tmin,Tmax), ylim=c(ymin,ymax1), xlab="T", ylab="R0(T)")
  # Fecundity
  #plot(seq(Tmin,Tmax,1), f, type="l", lwd=4, col="black", xlim=c(Tmin,Tmax), ylim=c(ymin,ymax1), xlab="T", ylab="f(T)")
  # Survival
  #plot(seq(Tmin,Tmax,1), s, type="l", lwd=4, col="black", xlim=c(Tmin,Tmax), ylim=c(ymin,ymax1), xlab="T", ylab="s(T)")
  # Birth rate
  #plot(seq(Tmin,Tmax,1), b, type="l", lwd=4, col="black", xlim=c(Tmin,Tmax), ylim=c(ymin,ymax1), xlab="T", ylab="b(T)")
  # Development time (maturation rate)
  #plot(seq(Tmin,Tmax,1), m, type="l", lwd=4, col="black", xlim=c(Tmin,Tmax), ylim=c(ymin,ymax1), xlab="T", ylab="m(T)")
  # Adult longevity
  #plot(seq(Tmin,Tmax,1), 1/dA, type="l", lwd=4, col="black", xlim=c(Tmin,Tmax), ylim=c(ymin,ymax1), xlab="T", ylab="1/dA(T)")
  # Recruitment
  #plot(seq(Tmin,Tmax,1), R, type="l", lwd=4, col="black", xlim=c(Tmin,Tmax), ylim=c(ymin,ymax1), xlab="T", ylab="R(T)")
  abline(v = t.param$meanT.h, col="blue", lwd=3, lty=1)
  abline(v = t.param$meanT.h + abs(t.param$amplT.h) + abs(t.param$amplD.h), col="blue", lwd=3, lty=2)
  abline(v = t.param$meanT.h - abs(t.param$amplT.h) - abs(t.param$amplD.h), col="blue", lwd=3, lty=2)
  abline(v = t.param$meanT.f + t.param$delta_mean.f*365*75 , col="red", lwd=3, lty=1)
  abline(v = t.param$meanT.f + t.param$delta_mean.f*365*75 + abs(t.param$amplT.f) + t.param$delta_ampl.f*365*80 + t.param$amplD.f, col="red", lwd=3, lty=2)
  abline(v = t.param$meanT.f + t.param$delta_mean.f*365*75 - abs(t.param$amplT.f) - t.param$delta_ampl.f*365*80 - t.param$amplD.f, col="red", lwd=3, lty=2)
  if(overw == TRUE) { abline(v = param$Tmin, col="black", lwd=3, lty=2) }
  # Temperature histograms
  par(new = T)
  hist(temp.h$T, xlim=c(Tmin,Tmax), ylim=c(ymin,ymax2), axes=F, xlab=NA, ylab=NA, breaks=seq(from=Tmin, to=Tmax, by=1), col=rgb(0,0,255, max = 255, alpha = 80), border=rgb(0,0,255, max = 255, alpha = 80), freq=FALSE, main = NULL)
  hist(temp.f[temp.f$day>365*65,"T"], xlim=c(Tmin,Tmax), ylim=c(ymin,ymax2), breaks=seq(from=Tmin, to=Tmax, by=1), ylab="r", col=rgb(255,0,0, max = 255, alpha = 80), border=rgb(255,0,0, max = 255, alpha = 80), freq=FALSE, main = NULL, add=TRUE)
  axis(side = 4)
}

# Summarize results
if(trait == "Fitness") { 
  print(paste("r.TPC.h:", r.TPC.h/param$rMax))
  print(paste("r.TPC.f:", r.TPC.f/param$rMax))
  print(paste("r.model.h:", r.model.h/param$rMax))
  print(paste("r.model.f:", r.model.f/param$rMax))
  print(paste("r.max.h:", r.max.h))
  print(paste("r.max.f:", r.max.f))
}
if(trait == "R0") { 
  print(paste("R0.TPC.h:", R0.TPC.h))
  print(paste("R0.TPC.f:", R0.TPC.f))
  print(paste("R0.model.h:", R0.model.h))
  print(paste("R0.model.f:", R0.model.f))
  print(paste("R0.max.h:", max(TS.h[-c(1:start.h),"R0"])))
  print(paste("R0.max.f:", max(TS.f[-c(1:start.f),"R0"])))
}
if(trait == "Fecundity") { 
  print(paste("f.TPC.h:", f.TPC.h))
  print(paste("f.TPC.f:", f.TPC.f))
  print(paste("f.model.h:", f.model.h))
  print(paste("f.model.f:", f.model.f))
  print(paste("f.max.h:", max(TS.h[-c(1:start.h),"f"])))
  print(paste("f.max.f:", max(TS.f[-c(1:start.f),"f"])))
}
if(trait == "Survival") { 
  print(paste("s.TPC.h:", s.TPC.h))
  print(paste("s.TPC.f:", s.TPC.f))
  print(paste("s.model.h:", s.model.h))
  print(paste("s.model.f:", s.model.f))
  print(paste("s.max.h:", max(TS.h[-c(1:start.h),"S"])))
  print(paste("s.max.f:", max(TS.f[-c(1:start.f),"S"])))
}
if(trait == "Birth") { 
  print(paste("b.TPC.h:", b.TPC.h))
  print(paste("b.TPC.f:", b.TPC.f))
  print(paste("b.model.h:", b.model.h))
  print(paste("b.model.f:", b.model.f))
  print(paste("b.max.h:", max(TS.h[-c(1:start.h),"b"])))
  print(paste("b.max.f:", max(TS.f[-c(1:start.f),"b"])))
}
if(trait == "Development") { 
  print(paste("G.TPC.h:", 1/m.TPC.h))
  print(paste("G.TPC.f:", 1/m.TPC.f))
  print(paste("G.model.h:", tau.model.h))
  print(paste("G.model.f:", tau.model.f))
  print(paste("G.max.h:", max(TS.h[-c(1:start.h),"tau"])))
  print(paste("G.max.f:", max(TS.f[-c(1:start.f),"tau"])))
}
if(trait == "Longevity") { 
  print(paste("1/dA.TPC.h:", 1/dA.TPC.h))
  print(paste("1/dA.TPC.f:", 1/dA.TPC.f))
  print(paste("1/dA.model.h:", 1/dA.model.h))
  print(paste("1/dA.model.f:", 1/dA.model.f))
  print(paste("1/dA.max.h:", 1/max(TS.h[-c(1:start.h),"dA"])))
  print(paste("1/dA.max.f:", 1/max(TS.f[-c(1:start.f),"dA"])))
}
if(trait == "Recruitment") { 
  print(paste("R.TPC.h:", R.TPC.h))
  print(paste("R.TPC.f:", R.TPC.f))
  print(paste("R.model.h:", R.model.h))
  print(paste("R.model.f:", R.model.f))
  print(paste("R.max.h:", max(TS.h[-c(1:start.h),"R"])))
  print(paste("R.max.f:", max(TS.f[-c(1:start.f),"R"])))
}
if(all == TRUE) { print(results) }

# Plot changes in life history traits
#barplot(c((b.TPC.f-b.TPC.h), (b.model.f-b.model.h)), col=c("Darkgreen","Orange"), ylim=c(-0.4,0.6), main=expression("Change in r"))
