###################################################################
#### This R script fits functions to temperature response data ####
###################################################################


# Load packages and set working directory
library(tidyr)
library(ggplot2)
library(dplyr)
library(tidyverse)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Read data
data <- as.data.frame(read_csv("Temperature response data.csv"))

# Select an insect by removing # in front of name and placing # in front of other species
sp.data <- subset(data, Species == "Clavigralla shadabi")
#sp.data <- subset(data, Species == "Clavigralla tomentosicollis")


# Remove columns that do not contain temperature data
sp.data <- sp.data[-c(1:8,12,14,16,18,20,21,23,24,26,27,29,31,32,34,35)]


# Set some option for nls and plots
Tmin <- 288
Tmax <- 313
TR <- 298

################################ FECUNDITY #####################################
# NLS
fec <- nls(Birth_Rate ~ bTopt*exp(-((T_K-Toptb)^2)/(2*sb^2)), data=sp.data,
           start=list(bTopt=2.5, Toptb=TR, sb=8))
summary(fec)

# Plot model fits
plot(sp.data$T_K, sp.data$Birth_Rate)
points(seq(Tmin,Tmax,1),coef(fec)[1]*exp(-((seq(Tmin,Tmax,1)-coef(fec)[2])^2)/(2*coef(fec)[3]^2)), type="l", col="blue")


############################### DEVELOPMENT ####################################
mTR <- sp.data[sp.data$T_K==TR,"Development"]

# NLS (monotonic response)
# NOTE: removed data beyond max development
dev.mon <- nls(Development ~ mTR*T_K/TR*exp(A*(1/TR-1/T_K)), data=sp.data[-c(nrow(sp.data)),],
                      start=list(A=1000))
summary(dev.mon)

# Plot model fits
plot(sp.data$T_K, sp.data$Development)
points(seq(Tmin,Tmax,1),mTR*(seq(Tmin,Tmax,1)/TR)*exp(coef(dev.mon)[1]*(1/TR-1/seq(Tmin,Tmax,1))), type="l", col="blue")


# NLS for development (Sharpe-Schoolfield response)
# estimate AL, TL, and TH
# NOTE: if needed, set AH based on NLS below and run code iteratively until AL and AH do not appreciably change)
dev.SS <- nls(Development ~ mTR*(T_K/TR)*exp(coef(dev.mon)[1]*(1/TR-1/T_K))/(1+exp(AL*(1/TL-1/T_K))+exp(AH*(1/TH-1/T_K))),
              data=sp.data, start=list(AL=-5000, AH=200000, TL=290, TH=310))
summary(dev.SS)

# Plot model fits
plot(sp.data$T_K, sp.data$Development)
points(seq(Tmin,Tmax,1), mTR*(seq(Tmin,Tmax,1)/TR)*exp(coef(dev.mon)[1]*(1/TR-1/seq(Tmin,Tmax,1)))/
         (1+(exp(coef(dev.SS)[1]*(1/coef(dev.SS)[3]-1/seq(Tmin,Tmax,1)))+exp(coef(dev.SS)[2]*(1/coef(dev.SS)[4]-1/seq(Tmin,Tmax,1))))), type="l", col="blue")

# estimate AH separately if needed
#dev.AH <- nls(Development ~ mTR*(T_K/TR)*exp(coef(dev.mon)[1]*(1/TR-1/T_K))/(1+exp(coef(dev.SS)[1]*(1/coef(dev.SS)[2]-1/T_K))+exp(AH*(1/coef(dev.SS)[3]-1/T_K))),
#               data=sp.data, start=list(AH=10000))
#summary(dev.AH)

# Plot model fits
#plot(sp.data$T_K, sp.data$Development)
#points(seq(Tmin,Tmax,1), mTR*(seq(Tmin,Tmax,1)/TR)*exp(coef(dev.mon)[1]*(1/TR-1/seq(Tmin,Tmax,1)))/
#         (1+(exp(coef(dev.SS)[1]*(1/coef(dev.SS)[2]-1/seq(Tmin,Tmax,1)))+exp(coef(dev.AH)[1]*(1/coef(dev.SS)[3]-1/seq(Tmin,Tmax,1))))), type="l", col="blue")


################################ MORTALITY #####################################
# Mortality estimated using data at reference temperature
dJTR <- sp.data[sp.data$T_K==TR,"Juv_Mortality"]
dATR <- sp.data[sp.data$T_K==TR,"Adult_Mortality"]

# NLS for juvenile mortality
# NOTE: removed data with 100% mortality
mort.J <- nls(Juv_Mortality ~ dJTR*exp(A*(1/TR-1/T_K)), data=sp.data[-c(nrow(sp.data)),],
              start=list(A=1000))
summary(mort.J)

# Plot model fits
plot(sp.data$T_K, sp.data$Juv_Mortality)
points(seq(Tmin,Tmax,1), dJTR*exp(coef(mort.J)[1]*(1/TR-1/seq(Tmin,Tmax,1))), type="l", col="blue")


# NLS for adult mortality 
# NOTE: removed data with 100% mortality
mort.A <- nls(Adult_Mortality ~ dATR*exp(A*(1/TR-1/T_K)), data=sp.data,
              start=list(A=1000))
summary(mort.A)

# Plot model fits
plot(sp.data$T_K, sp.data$Adult_Mortality)
points(seq(Tmin,Tmax,1), dATR*exp(coef(mort.A)[1]*(1/TR-1/seq(Tmin,Tmax,1))), type="l", col="blue")



# Mortality estimated using fit at reference temperature
# NLS for juvenile mortality
# NOTE: removed data with 100% mortality
mort.J <- nls(Juv_Mortality ~ xTR*exp(A*(1/TR-1/T_K)), data=sp.data[-c(nrow(sp.data)),],
              start=list(xTR=0.01, A=1000))
summary(mort.J)

# Plot model fits
plot(sp.data$T_K, sp.data$Juv_Mortality)
points(seq(Tmin,Tmax,1), coef(mort.J)[1]*exp(coef(mort.J)[2]*(1/TR-1/seq(Tmin,Tmax,1))), type="l", col="blue")


# NLS for adult mortality
# NOTE: removed data with 100% mortality
mort.A <- nls(Adult_Mortality ~ xTR*exp(A*(1/TR-1/T_K)), data=sp.data,
              start=list(xTR=0.01, A=1000))
summary(mort.A)

# Plot model fits
plot(sp.data$T_K, sp.data$Adult_Mortality)
points(seq(Tmin,Tmax,1), coef(mort.A)[1]*exp(coef(mort.A)[2]*(1/TR-1/seq(Tmin,Tmax,1))), type="l", col="blue")
