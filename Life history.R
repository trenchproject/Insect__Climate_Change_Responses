#################################################################################
#### This R script calculates a species' fecundity using TPCs and DDE model #####
#################################################################################

# Load packages and set working directory
library(tidyr)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(cubature)
library(lamW)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# USER: choose "Fecundity" or ""Survival"
trait <- "Fecundity"

# USER: enter species and location or set "all" to TRUE to run analysis for all species
species <- "Macrosiphum euphorbiae"
location <- "Canada"
all <- TRUE

# USER: include overwintering? (i.e., do not integrate over temperatures below Tmin)
overw <- TRUE

# USER: include diurnal variation?
daily <- FALSE


# Read in temperature response and temperature parameters, and temperature response data for selected insect
param <- subset(as.data.frame(read_csv("Temperature response parameters.csv")), Species == paste(species,location))
# Read in temperature parameters
ifelse(daily == TRUE, t.param <- subset(as.data.frame(read_csv("Temperature parameters.csv")), Species == paste(species,location)),
       t.param <- subset(as.data.frame(read_csv("Temperature parameters Tave.csv")), Species == paste(species,location)))

# Select datasets if running the analysis for all species
if(all == TRUE) {
  # Read in temperature response and temperature parameters, and temperature response data
  param.all <- as.data.frame(read_csv("Temperature response parameters.csv"))
  # Read in temperature parameters
  ifelse(daily == TRUE, t.param.all <- as.data.frame(read_csv("Temperature parameters.csv")),
         t.param.all <- as.data.frame(read_csv("Temperature parameters Tave.csv")))
  # Create array for results
  results <- data.frame(param.all[,1], 1:nrow(param.all), 1:nrow(param.all), 1:nrow(param.all), 1:nrow(param.all), 1:nrow(param.all), 1:nrow(param.all), 1:nrow(param.all), 1:nrow(param.all))
  names(results) <- c("Species","TPC.h","TPC.f","Model.h","Model.f","max.h","max.f","delta.TPC","delta.model")
}

# Run analysis for each species
for(s in 1:nrow(param.all)) {
  
  # Select species
  if(all == TRUE) {
    param <- param.all[s,]
    t.param <- t.param.all[s,]
  }
  
  
  ################################## TPC: HISTORICAL CLIMATE ###################################
  # Read in climate data (for plot)
  if(all == FALSE) {
    temp.h <- as.data.frame(read_csv(paste0("Climate data/Historical climate data ",location,".csv")))
    
    # For maximum temperature, remove daily minimum temperatures (if daily == FALSE)
    #if(daily == FALSE) { temp.h <- temp.h[temp.h$day %% 1 != 0,] }
    
    # average daily maximum and minimum temperatures (if daily == FALSE)
    if(daily == FALSE) {
      temp.h$day <- floor(temp.h$day)
      temp.h.min <- temp.h[duplicated(temp.h$day),]
      temp.h.max <- temp.h[duplicated(temp.h$day, fromLast=TRUE),]
      temp.h <- data.frame(temp.h.min$day, (temp.h.min$T + temp.h.max$T)/2)
      names(temp.h) <- c("day", "T") }
  }
  
  # Temperature function
  T.h <- function(t) { (t.param$meanT.h+t.param$delta_mean.h*t) - (t.param$amplT.h+t.param$delta_ampl.h*t)*cos(2*pi*(t + t.param$shiftT.h)/365) - t.param$amplD.h*cos(2*pi*t) }
  start <- 0
  end <- 5*365
  
  # Integrate across life history traits
  if(overw == FALSE) {
    # Life history traits
    b.h <- function(t) { param$bTopt*exp(-((T.h(t)-param$Toptb)^2)/(2*param$sb^2)) / (param$dATR*exp(param$AdA*(1/param$TR-1/T.h(t)))) }
    s.h <- function(t) { exp(-param$dJTR*exp(param$AdJ*(1/param$TR-1/T.h(t))) / (param$mTR*(T.h(t)/param$TR)*exp(param$AmJ*(1/param$TR-1/T.h(t)))/(1+exp(param$AL*(1/param$TL-1/T.h(t)))+exp(param$AH*(1/param$TH-1/T.h(t)))))) }
    
    # Note: pcubature is faster but cannot be used with overwintering
    b.TPC.h <- cubintegrate(b.h, lower = start, upper = end, method = "pcubature")$integral/(end-start)
    s.TPC.h <- cubintegrate(s.h, lower = start, upper = end, method = "pcubature")$integral/(end-start)
  }
  if(overw == TRUE) {
    # Life history traits
    b.h <- function(t) {
      ifelse(T.h(t) < param$Tmin, 0, param$bTopt*exp(-((T.h(t)-param$Toptb)^2)/(2*param$sb^2))) / (param$dATR*exp(param$AdA*(1/param$TR-1/T.h(t)))) }
    s.h <- function(t) {
      ifelse(T.h(t) < param$Tmin, 0, exp(-param$dJTR*exp(param$AdJ*(1/param$TR-1/T.h(t))) / (param$mTR*(T.h(t)/param$TR)*exp(param$AmJ*(1/param$TR-1/T.h(t)))/(1+exp(param$AL*(1/param$TL-1/T.h(t)))+exp(param$AH*(1/param$TH-1/T.h(t))))))) }
    
    # Integrate across active season
    season <- end - start # season length
    ifelse(daily == TRUE, length <- 0.5, length <- 1)
    for(t in seq(start,end,length)) { if(T.h(t) < param$Tmin) {season <- season - length }} # number of days when T(t) > Tmin
    # Note: hcubature must be used with overwintering
    b.TPC.h <- cubintegrate(b.h, lower = start, upper = end, method = "hcubature")$integral/season
    s.TPC.h <- cubintegrate(s.h, lower = start, upper = end, method = "hcubature")$integral/season
  }
  
  
  ##################################### TPC: FUTURE CLIMATE ####################################
  # Read in climate data (for plot)
  if(all == FALSE) {
    temp.f <- as.data.frame(read_csv(paste0("Climate data/Future climate data ",location,".csv")))
    
    # For maximum temperature, remove daily minimum temperatures (if daily == FALSE)
    #if(daily == FALSE) { temp.f <- temp.f[temp.f$day %% 1 != 0,] }
    
    # average daily maximum and minimum temperatures (if daily == FALSE)
    if(daily == FALSE) {
      temp.f$day <- floor(temp.f$day)
      temp.f.min <- temp.f[duplicated(temp.f$day),]
      temp.f.max <- temp.f[duplicated(temp.f$day, fromLast=TRUE),]
      temp.f <- data.frame(temp.f.min$day, (temp.f.min$T + temp.f.max$T)/2)
      names(temp.f) <- c("day", "T") }
  }
  
  # Temperature function
  T.f <- function(t) { (t.param$meanT.f+t.param$delta_mean.f*t) - (t.param$amplT.f+t.param$delta_ampl.f*t)*cos(2*pi*(t + t.param$shiftT.f)/365) - t.param$amplD.f*cos(2*pi*t) }
  start <- 365*70 # start 2095
  end <- 365*75 # end 2100
  
  # Integrate across life history traits
  if(overw == FALSE) {
    # Life history traits
    b.f <- function(t) { param$bTopt*exp(-((T.f(t)-param$Toptb)^2)/(2*param$sb^2)) / (param$dATR*exp(param$AdA*(1/param$TR-1/T.f(t)))) }
    s.f <- function(t) { exp(-param$dJTR*exp(param$AdJ*(1/param$TR-1/T.f(t))) / (param$mTR*(T.f(t)/param$TR)*exp(param$AmJ*(1/param$TR-1/T.f(t)))/(1+exp(param$AL*(1/param$TL-1/T.f(t)))+exp(param$AH*(1/param$TH-1/T.f(t)))))) }
    
    # Note: pcubature is faster but cannot be used with overwintering
    b.TPC.f <- cubintegrate(b.f, lower = start, upper = end, method = "pcubature")$integral/(end-start)
    s.TPC.f <- cubintegrate(s.f, lower = start, upper = end, method = "pcubature")$integral/(end-start)
  }
  if(overw == TRUE) {
    # Life history traits
    b.f <- function(t) { ifelse(T.f(t) < param$Tmin, 0, param$bTopt*exp(-((T.f(t)-param$Toptb)^2)/(2*param$sb^2))) / (param$dATR*exp(param$AdA*(1/param$TR-1/T.f(t)))) }
    s.f <- function(t) { ifelse(T.f(t) < param$Tmin, 0, exp(-param$dJTR*exp(param$AdJ*(1/param$TR-1/T.f(t))) / (param$mTR*(T.f(t)/param$TR)*exp(param$AmJ*(1/param$TR-1/T.f(t)))/(1+exp(param$AL*(1/param$TL-1/T.f(t)))+exp(param$AH*(1/param$TH-1/T.f(t))))))) }
    
    # Integrate across active season
    season <- end - start # season length
    ifelse(daily == TRUE, length <- 0.5, length <- 1)
    for(t in seq(start,end,length)) { if(T.f(t) < param$Tmin) {season <- season - length }} # number of days when T(t) > Tmin
    # Note: hcubature must be used with overwintering
    b.TPC.f <- cubintegrate(b.f, lower = start, upper = end, method = "hcubature")$integral/season
    s.TPC.f <- cubintegrate(s.f, lower = start, upper = end, method = "hcubature")$integral/season
  }

  
  ################################ PLOT LIFE HISTORY TRAIT RESPONSE ############################
  if(all == FALSE) { 
    Tmin <- round(min(temp.h$T,temp.f$T),0) - 3
    Tmax <- round(max(temp.h$T,temp.f$T),0) + 3
    ymin <- 0
    ymax1 <- round(param$bTopt,1) + 0.05
    ymax2 <- 0.2
    # TPC plots
    # fecundity
    plot(seq(Tmin,Tmax,1), param$bTopt*exp(-((seq(Tmin,Tmax,1)-param$Toptb)^2)/(2*param$sb^2)), type="l", lwd=4, col="black", xlim=c(Tmin,Tmax), ylim=c(ymin,ymax1), xlab="T", ylab="b(T)")
    # survival
    #plot(seq(Tmin,Tmax,1), exp(-param$dJTR*exp(AdJ*(1/TR-1/seq(Tmin,Tmax,1))) * mTR*(seq(Tmin,Tmax,1)/param$TR)*exp(param$AmJ*(1/param$TR-1/seq(Tmin,Tmax,1)))/(1+exp(param$AL*(1/param$TL-1/seq(Tmin,Tmax,1)))+exp(param$AH*(1/param$TH-1/seq(Tmin,Tmax,1))))), type="l", lwd=4, col="black", xlim=c(Tmin,Tmax), ylim=c(ymin,ymax1), xlab="T", ylab="s(T)")
    abline(v = t.param$meanT.h, col="blue", lwd=3, lty=1)
    abline(v = t.param$meanT.h + abs(t.param$amplT.h) + abs(t.param$amplD.h), col="blue", lwd=3, lty=2)
    abline(v = t.param$meanT.h - abs(t.param$amplT.h) - abs(t.param$amplD.h), col="blue", lwd=3, lty=2)
    abline(v = t.param$meanT.f + t.param$delta_mean.f*365*75 , col="red", lwd=3, lty=1)
    abline(v = t.param$meanT.f + t.param$delta_mean.f*365*75 + abs(t.param$amplT.f) + t.param$delta_ampl.f*365*80 + t.param$amplD.f, col="red", lwd=3, lty=2)
    abline(v = t.param$meanT.f + t.param$delta_mean.f*365*75 - abs(t.param$amplT.f) - t.param$delta_ampl.f*365*80 - t.param$amplD.f, col="red", lwd=3, lty=2)
    if(overw == TRUE) { abline(v = param$Tmin, col="black", lwd=3, lty=2) }
    # histograms of temperature
    par(new = T)
    hist(temp.h$T, xlim=c(Tmin,Tmax), ylim=c(ymin,ymax2), axes=F, xlab=NA, ylab=NA, breaks=seq(from=Tmin, to=Tmax, by=1), col=rgb(0,0,255, max = 255, alpha = 80), border=rgb(0,0,255, max = 255, alpha = 80), freq=FALSE, main = NULL)
    hist(temp.f[temp.f$day>365*65,"T"], xlim=c(Tmin,Tmax), ylim=c(ymin,ymax2), breaks=seq(from=Tmin, to=Tmax, by=1), ylab="r", col=rgb(255,0,0, max = 255, alpha = 80), border=rgb(255,0,0, max = 255, alpha = 80), freq=FALSE, main = NULL, add=TRUE)
    axis(side = 4)
  }
  
  
  ################################# MODEL: HISTORICAL CLIMATE ##################################
  # Read in climate data and temperature response parameters for selected insect
  if(all == FALSE) {
    ifelse(daily == TRUE, TS.h <- as.data.frame(read_csv(paste0("Time series data/Historical time series ",species," ",location,".csv"))),
           TS.h <- as.data.frame(read_csv(paste0("Time series data Tave/Historical time series ",species," ",location,".csv"))))
  }
  
  # Select species if running the analysis for all species
  if(all == TRUE) {
    ifelse(daily == TRUE, TS.h <- as.data.frame(read_csv(paste0("Time series data/Historical time series ",param[1],".csv"))),
           TS.h <- as.data.frame(read_csv(paste0("Time series data Tave/Historical time series ",param[1],".csv"))))
  }
  
  # Remove rows with NA
  TS.h <- na.omit(TS.h)
  # Set rows with negative survival to zero
  TS.h$S <- pmax(TS.h$S, 0)
  
  # Calculate life history traits at each time-step in model
  init_years <- 0 # from Python DDE model
  T.h <- function(t) { (t.param$meanT.h + t.param$delta_mean.h*t) - (t.param$amplT.h + t.param$delta_ampl.h*t)*cos(2*pi*(t + t.param$shiftT.h)/365) - t.param$amplD.h*cos(2*pi*t) }
  b <- function(t) { param$bTopt*exp(-((T.h(t)-param$Toptb)^2)/(2*param$sb^2)) }
  mJ <- function(t) { param$mTR*(T.h(t)/param$TR)*exp(param$AmJ*(1/param$TR-1/T.h(t)))/(1+exp(param$AL*(1/param$TL-1/T.h(t)))+exp(param$AH*(1/param$TH-1/T.h(t)))) }
  dJ <- function(t) {  param$dJTR*exp(param$AdJ*(1/param$TR-1/T.h(t))) }
  dA <- function(t) {  param$dATR*exp(param$AdA*(1/param$TR-1/T.h(t))) }
  TS.h$b <- b(TS.h$Time) / dA(TS.h$Time)
  
  # Integrate across daily life history traits from DDE model
  b.model.h <- 0
  s.model.h <- 0
  start <- nrow(TS.h) - 365*5 + 1 # integrate over last 5 years of time-series
  end <- nrow(TS.h)
  count <- 0
  for(i in start:end) {
    if(overw == FALSE) {
      b.model.h <- b.model.h + TS.h$b[i]
      s.model.h <- s.model.h + TS.h$S[i]
    }
    if(overw == TRUE) {
      if(T.h(i) >= param$Tmin) {
        b.model.h <- b.model.h + TS.h$b[i]
        s.model.h <- s.model.h + TS.h$S[i]
        count <- count + 1 } # number of days when T(t) > Tmin
    }
  }
  b.model.h <- b.model.h/count
  s.model.h <- s.model.h/count
  
  # Plot trait value over time
  #plot(TS.h[-c(1:start),"Time"],TS.h[-c(1:start),"b"], col="blue")
  #plot(TS.h$Time,TS.h$b, col="blue")
  
  
  ################################### MODEL: FUTURE CLIMATE ####################################
  # Read in climate data and temperature response parameters for selected insect
  if(all == FALSE) {
    ifelse(daily == TRUE, TS.f <- as.data.frame(read_csv(paste0("Time series data/Future time series ",species," ",location,".csv"))),
           TS.f <- as.data.frame(read_csv(paste0("Time series data Tave/Future time series ",species," ",location,".csv"))))
  }
  
  # Select species if running the analysis for all species
  if(all == TRUE) {
    ifelse(daily == TRUE, TS.f <- as.data.frame(read_csv(paste0("Time series data/Future time series ",param[1],".csv"))),
           TS.f <- as.data.frame(read_csv(paste0("Time series data Tave/Future time series ",param[1],".csv"))))
  }
  
  # Remove rows with NA or negative values
  TS.f <- na.omit(TS.f)
  # Set rows with negative survival to zero
  TS.f$S <- pmax(TS.f$S, 0)
  
  # Calculate r at each time-step in model
  init_years <- 0 # from Python DDE model
  T.f <- function(t) { (t.param$meanT.f + t.param$delta_mean.f*t) - (t.param$amplT.f + t.param$delta_ampl.f*t)*cos(2*pi*(t + t.param$shiftT.f)/365) - t.param$amplD.f*cos(2*pi*t) }
  b <- function(t) { param$bTopt*exp(-((T.f(t)-param$Toptb)^2)/(2*param$sb^2)) }
  mJ <- function(t) { param$mTR*(T.f(t)/param$TR)*exp(param$AmJ*(1/param$TR-1/T.f(t)))/(1+exp(param$AL*(1/param$TL-1/T.f(t)))+exp(param$AH*(1/param$TH-1/T.f(t)))) }
  dJ <- function(t) {  param$dJTR*exp(param$AdJ*(1/param$TR-1/T.f(t))) }
  dA <- function(t) {  param$dATR*exp(param$AdA*(1/param$TR-1/T.f(t))) }
  TS.f$b <- b(TS.f$Time) / dA(TS.f$Time)
  
  # Integrate across daily life history traits from DDE model
  b.model.f <- 0
  s.model.f <- 0
  start <- nrow(TS.f) - 365*5 + 1 # integrate over last 5 years of time-series
  end <- nrow(TS.f)
  count <- 0
  for(i in start:end) {
    if(overw == FALSE) {
      b.model.f <- b.model.f + TS.f$b[i]
      s.model.f <- s.model.f + TS.f$S[i]
    }
    if(overw == TRUE) {
      if(T.f(i) >= param$Tmin) {
        b.model.f <- b.model.f + TS.f$b[i]
        s.model.f <- s.model.f + TS.f$S[i]
        count <- count + 1 } # number of days when T(t) > Tmin
    }
  }
  b.model.f <- b.model.f/count
  s.model.f <- s.model.f/count
  
  # Plot life history traits over time
  #plot(TS.f[-c(1:start),"Time"],TS.f[-c(1:start),"b"], col="blue")
  #plot(TS.f$Time,TS.f$b, col="blue")
  
  
  # INPUT RESUTS INTO ARRAY
  if(all == TRUE) {
    if(trait == "Fecundity") {
      results[s,2] <- b.TPC.h
      results[s,3] <- b.TPC.f
      results[s,4] <- b.model.h
      results[s,5] <- b.model.f
      results[s,6] <- max(TS.h[-c(1:start),"b"])
      results[s,7] <- max(TS.f[-c(1:start),"b"])
      results[s,8] <- (b.TPC.f - b.TPC.h)/results[s,6]
      results[s,9] <- (b.model.f - b.model.h)/results[s,6]
    }
    if(trait == "Survival") {
      results[s,2] <- s.TPC.h
      results[s,3] <- s.TPC.f
      results[s,4] <- s.model.h
      results[s,5] <- s.model.f
      results[s,6] <- max(TS.h[-c(1:start),"S"])
      results[s,7] <- max(TS.f[-c(1:start),"S"])
      results[s,8] <- s.TPC.f - s.TPC.h
      results[s,9] <- s.model.f - s.model.h
    }
  }

  
  # BREAK FOR LOOP IF ANALYSES ARE RUN FOR A SPECIFIED SPECIES
  if(all == FALSE) { break  }
}


# OUTPUT RESULTS IN CSV FILE
if(all == TRUE) {
  if(trait == "Fecundity") { write_csv(results, "Fecundity.csv") }
  if(trait == "Survival") { write_csv(results, "Survival.csv") }
}


# SUMMARIZE RESULTS
if(trait == "Fecundity") { 
  print(paste("b.TPC.h:",b.TPC.h))
  print(paste("b.TPC.f:",b.TPC.f))
  print(paste("b.model.h:",b.model.h))
  print(paste("b.model.f:",b.model.f))
  print(paste("b.max.h:",max(TS.h[-c(1:start),"b"])))
  print(paste("b.max.f:",max(TS.f[-c(1:start),"b"])))
}
if(trait == "Survival") { 
  print(paste("s.TPC.h:",s.TPC.h))
  print(paste("s.TPC.f:",s.TPC.f))
  print(paste("s.model.h:",s.model.h))
  print(paste("s.model.f:",s.model.f))
  print(paste("s.max.h:",max(TS.h[-c(1:start),"S"])))
  print(paste("s.max.f:",max(TS.f[-c(1:start),"S"])))
}
if(all == TRUE) { print(results) }


# PLOT CHANGES IN LIFE HISTORY TRAIT
#barplot(c((b.TPC.f-b.TPC.h), (b.model.f-b.model.h)), col=c("Darkgreen","Orange"), ylim=c(-0.4,0.6), main=expression("Change in r"))


# STATISTICS
#model <- lm(results$delta.model ~ results$delta.TPC, data=results)
model <- lm(results$delta.model ~ 0 + results$delta.TPC, data=results)
summary(model) # significant


# PLOT
Xmin <- -0.5
Xmax <- 0.1
Ymin <- -1
Ymax <- 0.2
plot(results$delta.TPC, results$delta.model, pch=21, col="black", bg="black", xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Prop. change (TPC)", ylab="Prop. change (model)")
#points(seq(Xmin,Xmax,0.1), coef(model)[2]*seq(Xmin,Xmax,0.1)+coef(model)[1], type="l", col="black")
points(seq(Xmin,Xmax,0.1), coef(model)[1]*seq(Xmin,Xmax,0.1), type="l", col="black")
abline(0, 1, col="gray")
abline(0, 0, col="gray", lty="longdash")
abline(v = 0, col="gray", lty="longdash")

