##################################################################################
#### This R script calculates a species' intrinsic per capita growth rate (r) ####
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
species <- "Brevicoryne brassicae"
location <- "US Columbia"

# USER: include overwintering? (i.e., do not integrate over temperatures below Tmin)
overw <- TRUE

# USER: include diurnal variation?
daily <- FALSE

# USER: include resource variation due to precipitation?
res <- FALSE

# Read in temperature response and temperature parameters, and temperature response data for selected insect
param <- subset(as.data.frame(read_csv("Temperature response parameters.csv")), Species == paste(species,location))
# Read in temperature parameters
ifelse(daily == TRUE, t.param <- subset(as.data.frame(read_csv("Temperature parameters.csv")), Species == paste(species,location)),
       t.param <- subset(as.data.frame(read_csv("Temperature parameters Tave.csv")), Species == paste(species,location)))


################################## TPC: HISTORICAL CLIMATE ###################################
# Read in climate data
temp.h <- as.data.frame(read_csv(paste0("Climate data/Historical climate data ",location,".csv")))

# Remove daily minimum temperatures (if daily == FALSE)
if(daily == FALSE) { temp.h <- temp.h[temp.h$day %% 1 != 0,] }

# Integrate across r(T(t))
T.h <- function(t) { (t.param$meanT.h+t.param$delta_mean.h*t) - (t.param$amplT.h+t.param$delta_ampl.h*t)*cos(2*pi*(t + t.param$shiftT.h)/365) - t.param$amplD.h*cos(2*pi*t) }
start <- 0
end <- 5*365

if(overw == FALSE) {
  r.h <- function(t) {
           ifelse(T.h(t) <= param$rTopt, param$rMax*exp(-1*((T.h(t)-param$rTopt)/(2*param$rs))^2),
                  param$rMax*(1 - ((T.h(t)-param$rTopt)/(param$rTopt-param$rTmax))^2)) # from Deutsch et al. 2008
  }
  (r.TPC.h <- cubintegrate(r.h, lower = start, upper = end, method = "pcubature")$integral/(end-start)) # pcubature is faster but cannot be used with overwintering
}
if(overw == TRUE) {
  # r during active season
  r.h <- function(t) {
    ifelse(T.h(t) <= param$Tmin + 0*abs(t.param$amplD.h), 0,
           ifelse(T.h(t) <= param$rTopt, param$rMax*exp(-1*((T.h(t)-param$rTopt)/(2*param$rs))^2),
                  param$rMax*(1 - ((T.h(t)-param$rTopt)/(param$rTopt-param$rTmax))^2))) # from Deutsch et al. 2008
  }
  # integrate across active season
  season <- end # season length
  ifelse(daily == TRUE, length <- 0.5, length <- 1)
  for(t in seq(0,end,length)) { if(T.h(t) <= param$Tmin + 0*abs(t.param$amplD.h)) {season <- season - length }} # number of days when T(t) > Tmin
  (r.TPC.h <- cubintegrate(r.h, lower = start, upper = end, method = "hcubature")$integral/season)
}


##################################### TPC: FUTURE CLIMATE ####################################
# Read in climate data
temp.f <- as.data.frame(read_csv(paste0("Climate data/Future climate data ",location,".csv")))

# Remove daily minimum temperatures (if daily == FALSE)
if(daily == FALSE) { temp.f <- temp.f[temp.f$day %% 1 != 0,] }

# Integrate across r(T(t))
T.f <- function(t) { (t.param$meanT.f+t.param$delta_mean.f*t) - (t.param$amplT.f+t.param$delta_ampl.f*t)*cos(2*pi*(t + t.param$shiftT.f)/365) - t.param$amplD.f*cos(2*pi*t) }
start <- 365*70 # start 2090
end <- 365*75 # end 2100

if(overw == FALSE) {
  r.f <- function(t) {
           ifelse(T.f(t) <= param$rTopt, param$rMax*exp(-1*((T.f(t)-param$rTopt)/(2*param$rs))^2),
                  param$rMax*(1 - ((T.f(t)-param$rTopt)/(param$rTopt-param$rTmax))^2)) # from Deutsch et al. 2008
  }
  (r.TPC.f <- cubintegrate(r.f, lower = start, upper = end, method = "pcubature")$integral/(end-start)) # pcubature is faster but cannot be used with overwintering
}
if(overw == TRUE) {
  # r during active season
  r.f <- function(t) {
    ifelse(T.f(t) <= param$Tmin  + 0*abs(t.param$amplD.f), 0,
           ifelse(T.f(t) <= param$rTopt, param$rMax*exp(-1*((T.f(t)-param$rTopt)/(2*param$rs))^2),
                  param$rMax*(1 - ((T.f(t)-param$rTopt)/(param$rTopt-param$rTmax))^2))) # from Deutsch et al. 2008
  }
  # integrate across active season
  season <- end - start # season length
  ifelse(daily == TRUE, length <- 0.5, length <- 1)
  for(t in seq(start,end,length)) { if(T.f(t) <= param$Tmin  + 0*abs(t.param$amplD.f)) {season <- season - length }} # number of days when T(t) > Tmin
  (r.TPC.f <- cubintegrate(r.f, lower = start, upper = end, method = "hcubature")$integral/season)
}



# PLOT
Tmin <- round(min(temp.h$T,temp.f$T),0) - 3
Tmax <- round(max(temp.h$T,temp.f$T),0) + 1
ymin <- 0
ymax <- round(param$rMax,1) + 0.1
hist(temp.h$T, xlim=c(Tmin,Tmax), ylim=c(ymin,ymax), breaks=seq(from=Tmin, to=Tmax, by=1), ylab="r", col=rgb(0,0,255, max = 255, alpha = 80), border=rgb(0,0,255, max = 255, alpha = 80), freq=FALSE, main = NULL)
hist(temp.f[temp.f$day>365*65,"T"], xlim=c(Tmin,Tmax), ylim=c(ymin,ymax), breaks=seq(from=Tmin, to=Tmax, by=1), ylab="r", col=rgb(255,0,0, max = 255, alpha = 80), border=rgb(255,0,0, max = 255, alpha = 80), freq=FALSE, main = NULL, add=TRUE)
#abline(v = mean(param$rTopt), col="gray", lwd=3, lty=1)
#abline(v = mean(param$rTmax), col="gray", lwd=3, lty=2)
points(seq(Tmin,Tmax,1), ifelse(seq(Tmin,Tmax,1) <= param$rTopt, param$rMax*exp(-1*((seq(Tmin,Tmax,1)-param$rTopt)/(2*param$rs))^2),
                                param$rMax*(1 - ((seq(Tmin,Tmax,1)-param$rTopt)/(param$rTopt-param$rTmax))^2)), type="l", lwd=4, col="black")
abline(v = t.param$meanT.h, col="blue", lwd=3, lty=1)
abline(v = t.param$meanT.h + abs(t.param$amplT.h) + abs(t.param$amplD.h), col="blue", lwd=3, lty=2)
abline(v = t.param$meanT.h - abs(t.param$amplT.h) - abs(t.param$amplD.h), col="blue", lwd=3, lty=2)
abline(v = t.param$meanT.f + t.param$delta_mean.f*365*75 , col="red", lwd=3, lty=1)
abline(v = t.param$meanT.f + t.param$delta_mean.f*365*75 + abs(t.param$amplT.f) + t.param$delta_ampl.f*365*80 + t.param$amplD.f, col="red", lwd=3, lty=2)
abline(v = t.param$meanT.f + t.param$delta_mean.f*365*75 - abs(t.param$amplT.f) - t.param$delta_ampl.f*365*80 - t.param$amplD.f, col="red", lwd=3, lty=2)
if(overw == TRUE) { abline(v = param$Tmin  + 0*abs(t.param$amplD.h), col="black", lwd=3, lty=2) }
r.TPC.h
r.TPC.f



################################# MODEL: HISTORICAL CLIMATE ##################################
# Read in climate data and temperature response parameters for selected insect
ifelse(daily == TRUE, TS.h <- as.data.frame(read_csv(paste0("Time series data/Historical time series ",species," ",location,".csv"))),
       TS.h <- as.data.frame(read_csv(paste0("Time series data Tave/Historical time series ",species," ",location,".csv"))))

# Remove rows with NA
TS.h <- na.omit(TS.h)
# Set rows with negative survivorship to zero
TS.h$S <- pmax(TS.h$S, 0)

# Calculate r at each time-step in model
init_years <- 0 # from Python DDE model
T <- function(t) { (t.param$meanT.h + t.param$delta_mean.h*t) - (t.param$amplT.h + t.param$delta_ampl.h*t)*cos(2*pi*(t + t.param$shiftT.h)/365) - t.param$amplD.h*cos(2*pi*t) }
ifelse(res == FALSE, R <- function(t) {1}, R <- function(t) { ifelse(t.param$meanP.h - t.param$amplP.h * cos(2*pi*((t-init_years*365) + shiftP.h)/365) < 0, 0, t.param$meanP.h - t.param$amplP.h * cos(2*pi*((t-init_years*365) + shiftP.h)/365) )})
ifelse(overw == FALSE, M <- function(t) {1}, M <- function(t) { ifelse(T(t) < param$Tmin + 0*abs(t.param$amplD.h), 0, 1) })
b <- function(t) { param$bTopt*exp(-((T(t)-param$Toptb)^2)/(2*param$sb^2)) }
mJ <- function(t) { param$mTR*(T(t)/param$TR)*exp(param$AmJ*(1/param$TR-1/T(t)))/(1+exp(param$AL*(1/param$TL-1/T(t)))+exp(param$AH*(1/param$TH-1/T(t)))) }
dJ <- function(t) {  param$dJTR*exp(param$JdA*(1/param$TR-1/T(t))) }
dA <- function(t) {  param$dATR*exp(param$AdA*(1/param$TR-1/T(t))) }
ifelse(overw == FALSE, TS.h$r <- (lambertW0(R(S.h$Time)*M(TS.h$Time-TS.h$tau)*b(TS.h$Time-TS.h$tau) * mJ(TS.h$Time)/mJ(TS.h$Time-TS.h$tau)*TS.h$S * TS.h$tau * exp(dA(TS.h$Time)*TS.h$tau)) - dA(TS.h$Time) * TS.h$tau) / TS.h$tau,
       TS.h$r <- ifelse(T(TS.h$Time) < param$Tmin + 0*abs(t.param$amplD.h), 0, (lambertW0(R(S.h$Time)*M(TS.h$Time-TS.h$tau)*b(TS.h$Time-TS.h$tau) * mJ(TS.h$Time)/mJ(TS.h$Time-TS.h$tau)*TS.h$S * TS.h$tau * exp(dA(TS.h$Time)*TS.h$tau)) - dA(TS.h$Time) * TS.h$tau) / TS.h$tau))

# Integrate across daily per capita population growth rate from DDE model
r.model.h <- 0
start <- nrow(TS.h) - 365*5 + 1 # integrate over last 5 years of time-series
end <- nrow(TS.h)
count <- end - start
for(i in start:end) { r.model.h <- r.model.h + TS.h$r[i]
  if(T(i) < param$Tmin + 0*abs(t.param$amplD.h)) { count <- count - 1 } # number of days when T(t) > Tmin
}
(r.model.h <- r.model.h/count)

# Plot r over time
plot(TS.h[-c(1:start),"Time"],TS.h[-c(1:start),"r"], col="blue")
#plot(TS.h$Time,TS.h$r, col="blue")


# Integrate model time-series across ln(t/(t-1))
# r.model.h <- 0
# count.h <- 0
# for(i in 3:nrow(TS.h)) {
#   if(TS.h$A[i] >= 0 && TS.h$A[i-1] >= 0 && is.na(TS.h$A[i]) == FALSE && is.na(TS.h$A[i-1]) == FALSE) {
#     r.model.h <- r.model.h + log(TS.h$A[i]/TS.h$A[i-1])
#     count.h <- count.h + 1
# }}
# (r.model.h <- r.model.h/count.h)


################################### MODEL: FUTURE CLIMATE ####################################
# Read in climate data and temperature response parameters for selected insect
ifelse(daily == TRUE, TS.f <- as.data.frame(read_csv(paste0("Time series data/Future time series ",species," ",location,".csv"))),
       TS.f <- as.data.frame(read_csv(paste0("Time series data Tave/Future time series ",species," ",location,".csv"))))

# Remove rows with NA or negative values
TS.f <- na.omit(TS.f)
# Set rows with negative survivorship to zero
TS.f$S <- pmax(TS.f$S, 0)

# Calculate r at each time-step in model
init_years <- 0 # from Python DDE model
T <- function(t) { (t.param$meanT.f + t.param$delta_mean.f*t) - (t.param$amplT.f + t.param$delta_ampl.f*t)*cos(2*pi*(t + t.param$shiftT.f)/365) - t.param$amplD.f*cos(2*pi*t) }
ifelse(res == FALSE, R <- function(t) {1}, R <- function(t) { ifelse(t.param$meanP.f - t.param$amplP.f * cos(2*pi*((t-init_years*365) + shiftP.f)/365) < 0, 0, t.param$meanP.f - t.param$amplP.f * cos(2*pi*((t-init_years*365) + shiftP.f)/365) )})
ifelse(overw == FALSE, M <- function(t) {1}, M <- function(t) { ifelse(T(t) < param$Tmin + 0*abs(t.param$amplD.f), 0, 1) })
b <- function(t) { param$bTopt*exp(-((T(t)-param$Toptb)^2)/(2*param$sb^2)) }
mJ <- function(t) { param$mTR*(T(t)/param$TR)*exp(param$AmJ*(1/param$TR-1/T(t)))/(1+exp(param$AL*(1/param$TL-1/T(t)))+exp(param$AH*(1/param$TH-1/T(t)))) }
dJ <- function(t) {  param$dJTR*exp(param$AdJ*(1/param$TR-1/T(t))) }
dA <- function(t) {  param$dATR*exp(param$AdA*(1/param$TR-1/T(t))) }
ifelse(overw == FALSE, TS.f$r <- (lambertW0(R(S.f$Time)*M(TS.f$Time-TS.f$tau)*b(TS.f$Time-TS.f$tau) * mJ(TS.f$Time)/mJ(TS.f$Time-TS.f$tau)*TS.f$S * TS.f$tau * exp(dA(TS.f$Time)*TS.f$tau)) - dA(TS.f$Time) * TS.f$tau) / TS.f$tau,
       TS.f$r <- ifelse(T(TS.f$Time) < param$Tmin + 0*abs(t.param$amplD.f), 0, (lambertW0(R(S.f$Time)*M(TS.f$Time-TS.f$tau)*b(TS.f$Time-TS.f$tau) * mJ(TS.f$Time)/mJ(TS.f$Time-TS.f$tau)*TS.f$S * TS.f$tau * exp(dA(TS.f$Time)*TS.f$tau)) - dA(TS.f$Time) * TS.f$tau) / TS.f$tau))

# Integrate across daily per capita population growth rate from DDE model
r.model.f <- 0
start <- nrow(TS.f) - 365*5 + 1 # integrate over last 5 years of time-series
end <- nrow(TS.f)
count <- end - start
for(i in start:end) { r.model.f <- r.model.f + TS.f$r[i] #- dA(i)
  if(T(i) < param$Tmin + 0*abs(t.param$amplD.f)) { count <- count - 1 } # number of days when T(t) > Tmin
}
(r.model.f <- r.model.f/count)

# Plot r over time
plot(TS.f[-c(1:start),"Time"],TS.f[-c(1:start),"r"], col="blue")
#plot(TS.f$Time,TS.f$r, col="blue")


# Integrate across ln(t/(t-1))
# r.model.f <- 0
# count.f <- 0
# for(i in 3:nrow(TS.f)) {
#   if(TS.f$A[i] >= 0 && TS.f$A[i-1] >= 0 && is.na(TS.f$A[i]) == FALSE && is.na(TS.f$A[i-1]) == FALSE) {
#     r.model.f <- r.model.f + log(TS.f$A[i]/TS.f$A[i-1])
#     count.f <- count.f + 1
#   }}
# (r.model.f <- r.model.f/count.f)


# SUMMARIZE RESULTS
r.TPC.h
r.TPC.f
r.model.h
r.model.f
#r.TPC.f/r.TPC.h
#r.model.f/r.model.h

# PLOT CHANGES IN r
#barplot(c((r.TPC.f-r.TPC.h)/r.TPC.h, (r.model.f-r.model.h)/r.model.h), col=c("Darkgreen","Orange"), ylim=c(-0.4,0.6), main=expression("Proportional change in r"))
