##################################################################################
#### This R script calculates a species' fecundity using TPCs and DDE model ######
##################################################################################

# Load packages and set working directory
library(tidyr)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(cubature)
library(lamW)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# USER: enter species and location
species <- "Clavigralla shadabi"
location <- "Benin"

# USER: include overwintering? (i.e., do not integrate over temperatures below Tmin)
overw <- FALSE

# USER: include diurnal variation?
daily <- FALSE


# Read in temperature response and temperature parameters, and temperature response data for selected insect
param <- subset(as.data.frame(read_csv("Temperature response parameters.csv")), Species == paste(species,location))
# Read in temperature parameters
ifelse(daily == TRUE, t.param <- subset(as.data.frame(read_csv("Temperature parameters.csv")), Species == paste(species,location)),
       t.param <- subset(as.data.frame(read_csv("Temperature parameters Tave.csv")), Species == paste(species,location)))


################################## TPC: HISTORICAL CLIMATE ###################################
# Read in climate data
temp.h <- as.data.frame(read_csv(paste0("Climate data/Historical climate data ",location,".csv")))

# Remove daily minimum temperatures (if daily == FALSE)
#if(daily == FALSE) { temp.h <- temp.h[temp.h$day %% 1 != 0,] }

# average daily maximum and minimum temperatures (if daily == FALSE)
if(daily == FALSE) { temp.h$day <- floor(temp.h$day)
  temp.h.min <- temp.h[duplicated(temp.h$day),]
  temp.h.max <- temp.h[duplicated(temp.h$day, fromLast=TRUE),]
  temp.h <- data.frame(temp.h.min$day, (temp.h.min$T + temp.h.max$T)/2)
  names(temp.h) <- c("day", "T") }

# Integrate across b(T(t))
T.h <- function(t) { (t.param$meanT.h+t.param$delta_mean.h*t) - (t.param$amplT.h+t.param$delta_ampl.h*t)*cos(2*pi*(t + t.param$shiftT.h)/365) - t.param$amplD.h*cos(2*pi*t) }
start <- 0
end <- 5*365

if(overw == FALSE) {
  b.h <- function(t) { param$bTopt*exp(-((T.h(t)-param$Toptb)^2)/(2*param$sb^2)) }
  (b.TPC.h <- cubintegrate(b.h, lower = start, upper = end, method = "pcubature")$integral/(end-start)) # pcubature is faster but cannot be used with overwintering
}
if(overw == TRUE) {
  # b during active season
  b.h <- function(t) {
    ifelse(T.h(t) <= param$Tmin + 0*abs(t.param$amplD.h), 0, param$bTopt*exp(-((T.h(t)-param$Toptb)^2)/(2*param$sb^2)))
  }
  # integrate across active season
  season <- end # season length
  ifelse(daily == TRUE, length <- 0.5, length <- 1)
  for(t in seq(0,end,length)) { if(T.h(t) <= param$Tmin + 0*abs(t.param$amplD.h)) {season <- season - length }} # number of days when T(t) > Tmin
  (b.TPC.h <- cubintegrate(b.h, lower = start, upper = end, method = "hcubature")$integral/season)
}


##################################### TPC: FUTURE CLIMATE ####################################
# Read in climate data
temp.f <- as.data.frame(read_csv(paste0("Climate data/Future climate data ",location,".csv")))

# Remove daily minimum temperatures (if daily == FALSE)
#if(daily == FALSE) { temp.f <- temp.f[temp.f$day %% 1 != 0,] }

# average daily maximum and minimum temperatures (if daily == FALSE)
if(daily == FALSE) { temp.f$day <- floor(temp.f$day)
temp.f.min <- temp.f[duplicated(temp.f$day),]
temp.f.max <- temp.f[duplicated(temp.f$day, fromLast=TRUE),]
temp.f <- data.frame(temp.f.min$day, (temp.f.min$T + temp.f.max$T)/2)
names(temp.f) <- c("day", "T") }

# Integrate across b(T(t))
T.f <- function(t) { (t.param$meanT.f+t.param$delta_mean.f*t) - (t.param$amplT.f+t.param$delta_ampl.f*t)*cos(2*pi*(t + t.param$shiftT.f)/365) - t.param$amplD.f*cos(2*pi*t) }
start <- 365*70 # start 2090
end <- 365*75 # end 2100

if(overw == FALSE) {
  b.f <- function(t) { param$bTopt*exp(-((T.f(t)-param$Toptb)^2)/(2*param$sb^2)) }
  (b.TPC.f <- cubintegrate(b.f, lower = start, upper = end, method = "pcubature")$integral/(end-start)) # pcubature is faster but cannot be used with overwintering
}
if(overw == TRUE) {
  # b during active season
  b.f <- function(t) {
    ifelse(T.f(t) <= param$Tmin  + 0*abs(t.param$amplD.f), 0, param$bTopt*exp(-((T.f(t)-param$Toptb)^2)/(2*param$sb^2)))
  }
  # integrate across active season
  season <- end - start # season length
  ifelse(daily == TRUE, length <- 0.5, length <- 1)
  for(t in seq(start,end,length)) { if(T.f(t) <= param$Tmin  + 0*abs(t.param$amplD.f)) {season <- season - length }} # number of days when T(t) > Tmin
  (b.TPC.f <- cubintegrate(b.f, lower = start, upper = end, method = "hcubature")$integral/season)
}



# PLOT
Tmin <- round(min(temp.h$T,temp.f$T),0) - 3
Tmax <- round(max(temp.h$T,temp.f$T),0) + 3
ymin <- 0
ymax1 <- round(param$bTopt,0) + 1
ymax2 <- 0.5
# TPC plots
plot(seq(Tmin,Tmax,1), param$bTopt*exp(-((seq(Tmin,Tmax,1)-param$Toptb)^2)/(2*param$sb^2)), type="l", lwd=4, col="black", xlim=c(Tmin,Tmax), ylim=c(ymin,ymax1), xlab="T", ylab="b(T)")
abline(v = t.param$meanT.h, col="blue", lwd=3, lty=1)
abline(v = t.param$meanT.h + abs(t.param$amplT.h) + abs(t.param$amplD.h), col="blue", lwd=3, lty=2)
abline(v = t.param$meanT.h - abs(t.param$amplT.h) - abs(t.param$amplD.h), col="blue", lwd=3, lty=2)
abline(v = t.param$meanT.f + t.param$delta_mean.f*365*75 , col="red", lwd=3, lty=1)
abline(v = t.param$meanT.f + t.param$delta_mean.f*365*75 + abs(t.param$amplT.f) + t.param$delta_ampl.f*365*80 + t.param$amplD.f, col="red", lwd=3, lty=2)
abline(v = t.param$meanT.f + t.param$delta_mean.f*365*75 - abs(t.param$amplT.f) - t.param$delta_ampl.f*365*80 - t.param$amplD.f, col="red", lwd=3, lty=2)
if(overw == TRUE) { abline(v = param$Tmin  + 0*abs(t.param$amplD.h), col="black", lwd=3, lty=2) }
# histograms of temperature
par(new = T)
hist(temp.h$T, xlim=c(Tmin,Tmax), ylim=c(ymin,ymax2), axes=F, xlab=NA, ylab=NA, breaks=seq(from=Tmin, to=Tmax, by=1), col=rgb(0,0,255, max = 255, alpha = 80), border=rgb(0,0,255, max = 255, alpha = 80), freq=FALSE, main = NULL)
hist(temp.f[temp.f$day>365*65,"T"], xlim=c(Tmin,Tmax), ylim=c(ymin,ymax2), breaks=seq(from=Tmin, to=Tmax, by=1), ylab="r", col=rgb(255,0,0, max = 255, alpha = 80), border=rgb(255,0,0, max = 255, alpha = 80), freq=FALSE, main = NULL, add=TRUE)
axis(side = 4)



################################# MODEL: HISTORICAL CLIMATE ##################################
# Read in climate data and temperature response parameters for selected insect
ifelse(daily == TRUE, TS.h <- as.data.frame(read_csv(paste0("Time series data/Historical time series ",species," ",location,".csv"))),
       TS.h <- as.data.frame(read_csv(paste0("Time series data Tave/Historical time series ",species," ",location,".csv"))))

# Remove rows with NA
TS.h <- na.omit(TS.h)
# Set rows with negative survivorship to zero
TS.h$S <- pmax(TS.h$S, 0)

# Calculate b at each time-step in model
init_years <- 0 # from Python DDE model
T <- function(t) { (t.param$meanT.h + t.param$delta_mean.h*t) - (t.param$amplT.h + t.param$delta_ampl.h*t)*cos(2*pi*(t + t.param$shiftT.h)/365) - t.param$amplD.h*cos(2*pi*t) }
b <- function(t) { param$bTopt*exp(-((T(t)-param$Toptb)^2)/(2*param$sb^2)) }
TS.h$b <- b(TS.h$Time)

# Integrate across daily per capita birth rate from DDE model
b.model.h <- 0
start <- nrow(TS.h) - 365*5 + 1 # integrate over last 5 years of time-series
end <- nrow(TS.h)
count <- end - start
for(i in start:end) { b.model.h <- b.model.h + TS.h$b[i]
if(T(i) < param$Tmin + 0*abs(t.param$amplD.h)) { count <- count - 1 } # number of days when T(t) > Tmin
}
(b.model.h <- b.model.h/count)

# Plot b over time
plot(TS.h[-c(1:start),"Time"],TS.h[-c(1:start),"b"], col="blue")
#plot(TS.h$Time,TS.h$b, col="blue")


################################### MODEL: FUTURE CLIMATE ####################################
# Read in climate data and temperature response parameters for selected insect
ifelse(daily == TRUE, TS.f <- as.data.frame(read_csv(paste0("Time series data/Future time series ",species," ",location,".csv"))),
       TS.f <- as.data.frame(read_csv(paste0("Time series data Tave/Future time series ",species," ",location,".csv"))))

# Remove rows with NA or negative values
TS.f <- na.omit(TS.f)
# Set rows with negative survivorship to zero
TS.f$S <- pmax(TS.f$S, 0)

# Calculate b at each time-step in model
init_years <- 0 # from Python DDE model
T <- function(t) { (t.param$meanT.f + t.param$delta_mean.f*t) - (t.param$amplT.f + t.param$delta_ampl.f*t)*cos(2*pi*(t + t.param$shiftT.f)/365) - t.param$amplD.f*cos(2*pi*t) }
b <- function(t) { param$bTopt*exp(-((T(t)-param$Toptb)^2)/(2*param$sb^2)) }
TS.f$b <- b(TS.f$Time)

# Integrate across daily per capita birth rate from DDE model
b.model.f <- 0
start <- nrow(TS.f) - 365*5 + 1 # integrate over last 5 years of time-series
end <- nrow(TS.f)
count <- end - start
for(i in start:end) { b.model.f <- b.model.f + TS.f$b[i]
if(T(i) < param$Tmin + 0*abs(t.param$amplD.f)) { count <- count - 1 } # number of days when T(t) > Tmin
}
(b.model.f <- b.model.f/count)

# Plot b over time
plot(TS.f[-c(1:start),"Time"],TS.f[-c(1:start),"b"], col="blue")
#plot(TS.f$Time,TS.f$b, col="blue")



# SUMMARIZE RESULTS
b.TPC.h
b.TPC.f
b.model.h
b.model.f

# PLOT CHANGES IN b
#barplot(c((b.TPC.f-b.TPC.h), (b.model.f-b.model.h)), col=c("Darkgreen","Orange"), ylim=c(-0.4,0.6), main=expression("Change in b"))

