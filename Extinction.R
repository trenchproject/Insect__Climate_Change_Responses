###################################################################################
#### This R script calculates a species' temperature sensitivity to extinction ####
###################################################################################

# Load packages and set working directory
library(tidyr)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(cubature)
library(lamW)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# USER: run analyses for temperature "mean", "ampl", or "both"?
case <- "ampl"

# USER: include overwintering? (i.e., do not integrate over temperatures below Tmin)
overw <- TRUE

# USER: include diurnal variation?
daily <- FALSE


# READ PARAMETERS AND TEMPERATURE DATA
param.all <- as.data.frame(read_csv("Temperature response parameters.csv"))
ifelse(daily == TRUE, t.param.all <- as.data.frame(read_csv("Temperature parameters.csv")),
       t.param.all <- as.data.frame(read_csv("Temperature parameters Tave.csv")))


# CREATE ARRAY FOR RESULTS
results <- data.frame(param.all[,1], param.all[,2], param.all[,3], param.all[,4], 1:nrow(param.all), 1:nrow(param.all))
names(results) <- c("Species","Latitude","Habitat","Subfamily","TPC","Model")


# RUN ANALYSES FOR EACH SPECIES
for(s in 1:nrow(param.all)) {
  
# Select species
param <- param.all[s,]
t.param <- t.param.all[s,]


############################################ TPC ############################################
# Increase in mean temperature mean and amplitude of temperature fluctuations
delta.mean <- 0
delta.ampl <- 0

# Repeat integration of TPC as mean temperature increase until r < 0
repeat{

# Integrate across r(T(t))
T <- function(t) { (t.param$meanT.h+delta.mean) - (t.param$amplT.h+delta.ampl)*cos(2*pi*(t + t.param$shiftT.h)/365) - t.param$amplD.h*cos(2*pi*t) }
start <- 0
end <- 10*365

if(overw == FALSE) {
  r <- function(t) {
    ifelse(T(t) <= param$rTopt, param$rMax*exp(-1*((T(t)-param$rTopt)/(2*param$rs))^2),
           param$rMax*(1 - ((T(t)-param$rTopt)/(param$rTopt-param$rTmax))^2)) # from Deutsch et al. 2008
  }
  r.TPC <- cubintegrate(r, lower = start, upper = end, method = "pcubature")$integral/(end-start) # pcubature is faster but cannot be used with overwintering
}
if(overw == TRUE) {
  # r during active season
  r <- function(t) {
    ifelse(T(t) <= param$Tmin, 0,
           ifelse(T(t) <= param$rTopt, param$rMax*exp(-1*((T(t)-param$rTopt)/(2*param$rs))^2),
                  param$rMax*(1 - ((T(t)-param$rTopt)/(param$rTopt-param$rTmax))^2))) # from Deutsch et al. 2008
  }
  # integrate across active season
  season <- end # season length
  ifelse(daily == TRUE, length <- 0.5, length <- 1)
  for(t in seq(0,end,length)) { if(T(t) <= param$Tmin) {season <- season - length }} # number of days when T(t) > Tmin
  r.TPC <- cubintegrate(r, lower = start, upper = end, method = "hcubature")$integral/season
}
#print(r.TPC)

# Evaluate whether r < 0 and if so, break repeat
if(r.TPC <= 0){
  if(case == "mean") { results[s,5] <- delta.mean }
  if(case == "ampl") { results[s,5] <- delta.ampl }
  break
}

# Increase in temperature
if(case == "mean") { delta.mean <- delta.mean + 0.1 }
if(case == "ampl") { delta.ampl <- delta.ampl + 0.1 }
}


########################################### MODEL ############################################
if(case == "mean") { TS <- as.data.frame(read_csv(paste0("Time series data Ext meanT/Time series ",param[1],".csv"))) }
if(case == "ampl") { TS <- as.data.frame(read_csv(paste0("Time series data Ext amplT/Time series ",param[1],".csv"))) }

# Time of extinction
ext.time <- TS[(TS$J == 0 & TS$A == 0), "Time"][1]
#print(ext.time)

# Calculate change in temperature at time of extinction
if(case == "mean") { results[s,6] <- 0.1/365*ext.time }
if(case == "ampl") { results[s,6] <- 0.1/365*ext.time }
}


# OUTPUT RESULTS IN CSV FILE
if(case == "mean") { write_csv(results, "Extinction meanT.csv") }
if(case == "ampl") { write_csv(results, "Extinction amplT.csv") }


