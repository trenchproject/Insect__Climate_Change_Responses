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
#sp.data <- subset(data, Species == "Clavigralla tomentosicollis Benin")
#sp.data <- subset(data, Species == "Clavigralla tomentosicollis Nigeria")
#sp.data <- subset(data, Species == "Clavigralla tomentosicollis Burkina Faso")
#sp.data <- subset(data, Species == "Apolygus lucorum")
#sp.data <- subset(data, Species == "Adelphocoris suturalis")
#sp.data <- subset(data, Species == "Macrosiphum euphorbiae Brazil")
#sp.data <- subset(data, Species == "Aulacorthum solani Brazil")
#sp.data <- subset(data, Species == "Uroleucon ambrosiae")
#sp.data <- subset(data, Species == "Lygus lineolaris")
#sp.data <- subset(data, Species == "Pilophorus typicus")
#sp.data <- subset(data, Species == "Macrolophus pygmaeus on Myzus persicae")
#sp.data <- subset(data, Species == "Macrolophus pygmaeus on Trialeurodes vaporariorum")

# Remove columns that do not contain temperature data
sp.data <- sp.data[-c(1:8,12,14,16,18,20,22,24,26,27,29,31,32,34,35,37,39,40,42,44,46,48,49)]


# Set some option for nls and plots
Tmin <- 285
Tmax <- 315
TR <- 298



################################ FECUNDITY #####################################
# NLS
fec <- nls(Birth_Rate ~ bTopt*exp(-((T_K-Toptb)^2)/(2*sb^2)), data=sp.data,
           start=list(bTopt=5, Toptb=TR+2, sb=3))
summary(fec)
# Plot model fits
plot(sp.data$T_K, sp.data$Birth_Rate)
points(seq(Tmin,Tmax,1),coef(fec)[1]*exp(-((seq(Tmin,Tmax,1)-coef(fec)[2])^2)/(2*coef(fec)[3]^2)), type="l", col="blue")



############################### DEVELOPMENT ####################################
# NLS (monotonic response)
dev.mon <- nls(Development ~ xTR*T_K/TR*exp(A*(1/TR-1/T_K)), data=sp.data,
               start=list(xTR=0.01, A=5000))
summary(dev.mon)
# Plot model fits
plot(sp.data$T_K, sp.data$Development)
points(seq(Tmin,Tmax,1), coef(dev.mon)[1]*(seq(Tmin,Tmax,1)/TR)*exp(coef(dev.mon)[2]*(1/TR-1/seq(Tmin,Tmax,1))), type="l", col="blue")


# NLS for development (Sharpe-Schoolfield response)
# estimate all parameters
dev <- nls(Development ~ xTR*(T_K/TR)*exp(A*(1/TR-1/T_K))/(1+exp(AL*(1/TL-1/T_K))+exp(AH*(1/TH-1/T_K))),
              data=sp.data, start=list(xTR=0.03, A=7000, AL=-30000, AH=270000, TL=286, TH=306))
summary(dev)
# Plot model fits
plot(sp.data$T_K, sp.data$Development)
points(seq(Tmin,Tmax,1), coef(dev)[1]*(seq(Tmin,Tmax,1)/TR)*exp(coef(dev)[2]*(1/TR-1/seq(Tmin,Tmax,1)))/
         (1+(exp(coef(dev)[3]*(1/coef(dev)[5]-1/seq(Tmin,Tmax,1)))+exp(coef(dev)[4]*(1/coef(dev)[6]-1/seq(Tmin,Tmax,1))))), type="l", col="blue")

# estimate xTR and A
# NOTE: removed data beyond max development
dev.mon <- nls(Development ~ xTR*T_K/TR*exp(A*(1/TR-1/T_K)), data=sp.data[-c((nrow(sp.data)-1):nrow(sp.data)),],
               start=list(xTR=0.1, A=1000))
summary(dev.mon)
# Plot model fits
plot(sp.data$T_K, sp.data$Development)
points(seq(Tmin,Tmax,1), coef(dev.mon)[1]*(seq(Tmin,Tmax,1)/TR)*exp(coef(dev.mon)[2]*(1/TR-1/seq(Tmin,Tmax,1))), type="l", col="blue")

# estimate AL, AH, TL, and TH
# NOTE: if needed, set AH based on NLS below and run code iteratively until AL and AH do not appreciably change)
dev.SS <- nls(Development ~ coef(dev.mon)[1]*(T_K/TR)*exp(coef(dev.mon)[2]*(1/TR-1/T_K))/(1+exp(AL*(1/TL-1/T_K))+exp(AH*(1/TH-1/T_K))),
              data=sp.data, start=list(AL=-30000, AH=270000, TL=286, TH=306))
summary(dev.SS)
# Plot model fits
plot(sp.data$T_K, sp.data$Development)
points(seq(Tmin,Tmax,1), coef(dev.mon)[1]*(seq(Tmin,Tmax,1)/TR)*exp(coef(dev.mon)[2]*(1/TR-1/seq(Tmin,Tmax,1)))/
         (1+(exp(coef(dev.SS)[1]*(1/coef(dev.SS)[3]-1/seq(Tmin,Tmax,1)))+exp(coef(dev.SS)[2]*(1/coef(dev.SS)[4]-1/seq(Tmin,Tmax,1))))), type="l", col="blue")

# estimate AL and AH separately from TL and TH if needed
kTL <- 285
kTH <- 305
dev.A <- nls(Development ~ coef(dev.mon)[1]*(T_K/TR)*exp(coef(dev.mon)[2]*(1/TR-1/T_K))/(1+exp(AL*(1/kTL-1/T_K))+exp(AH*(1/kTH-1/T_K))),
              data=sp.data, start=list(AL=-50000, AH=50000))
summary(dev.A)
dev.T <- nls(Development ~ coef(dev.mon)[1]*(T_K/TR)*exp(coef(dev.mon)[2]*(1/TR-1/T_K))/(1+exp(coef(dev.A)[1]*(1/TL-1/T_K))+exp(coef(dev.A)[2]*(1/TH-1/T_K))),
             data=sp.data, start=list(TL=kTL, TH=kTH))
summary(dev.T)
# Plot model fits
plot(sp.data$T_K, sp.data$Development)#, ylim=c(0,0.07))
points(seq(Tmin,Tmax,1), coef(dev.mon)[1]*(seq(Tmin,Tmax,1)/TR)*exp(coef(dev.mon)[2]*(1/TR-1/seq(Tmin,Tmax,1)))/
         (1+(exp(coef(dev.A)[1]*(1/coef(dev.T)[1]-1/seq(Tmin,Tmax,1)))+exp(coef(dev.A)[2]*(1/coef(dev.T)[2]-1/seq(Tmin,Tmax,1))))), type="l", col="blue")

# estimate TH and AH separately from TL and AL if needed
kAL <- -100000
dev.H <- nls(Development ~ coef(dev.mon)[1]*(T_K/TR)*exp(coef(dev.mon)[2]*(1/TR-1/T_K))/(1+exp(kAL*(1/kTL-1/T_K))+exp(AH*(1/TH-1/T_K))),
             data=sp.data, start=list(TH=kTH, AH=50000))
summary(dev.H)
dev.AL <- nls(Development ~ coef(dev.mon)[1]*(T_K/TR)*exp(coef(dev.mon)[2]*(1/TR-1/T_K))/(1+exp(AL*(1/kTL-1/T_K))+exp(coef(dev.H)[2]*(1/coef(dev.H)[1]-1/T_K))),
             data=sp.data, start=list(AL=-100000))
summary(dev.AL)
dev.TL <- nls(Development ~ coef(dev.mon)[1]*(T_K/TR)*exp(coef(dev.mon)[2]*(1/TR-1/T_K))/(1+exp(coef(dev.AL)[1]*(1/TL-1/T_K))+exp(coef(dev.H)[2]*(1/coef(dev.H)[1]-1/T_K))),
             data=sp.data, start=list(TL=kTL))
summary(dev.TL)
# Plot model fits
plot(sp.data$T_K, sp.data$Development)
points(seq(Tmin,Tmax,1), coef(dev.mon)[1]*(seq(Tmin,Tmax,1)/TR)*exp(coef(dev.mon)[2]*(1/TR-1/seq(Tmin,Tmax,1)))/
         (1+(exp(coef(dev.AL)[1]*(1/coef(dev.TL)[1]-1/seq(Tmin,Tmax,1)))+exp(coef(dev.H)[2]*(1/coef(dev.H)[1]-1/seq(Tmin,Tmax,1))))), type="l", col="blue")



################################ MORTALITY #####################################
# Mortality estimated using fit at reference temperature
# NLS for juvenile mortality
mort.J <- nls(Juv_Mortality ~ xTR*exp(A*(1/TR-1/T_K)), data=sp.data,
              start=list(xTR=0.01, A=10000))
summary(mort.J)
# Plot model fits
plot(sp.data$T_K, sp.data$Juv_Mortality, ylim=c(0,1))
points(seq(Tmin,Tmax,1), coef(mort.J)[1]*exp(coef(mort.J)[2]*(1/TR-1/seq(Tmin,Tmax,1))), type="l", col="blue")


# NLS for adult mortality
mort.A <- nls(Adult_Mortality ~ xTR*exp(A*(1/TR-1/T_K)), data=sp.data,
              start=list(xTR=0.01, A=10000))
summary(mort.A)
# Plot model fits
plot(sp.data$T_K, sp.data$Adult_Mortality, ylim=c(0,1))
points(seq(Tmin,Tmax,1), coef(mort.A)[1]*exp(coef(mort.A)[2]*(1/TR-1/seq(Tmin,Tmax,1))), type="l", col="blue")


# Mortality estimated using data at reference temperature
dJTR <- sp.data[sp.data$T_K==TR,"Juv_Mortality"]
dATR <- sp.data[sp.data$T_K==TR,"Adult_Mortality"]

# NLS for juvenile mortality
mort.J <- nls(Juv_Mortality ~ dJTR*exp(A*(1/TR-1/T_K)), data=sp.data,
              start=list(A=1000))
summary(mort.J)
# Plot model fits
plot(sp.data$T_K, sp.data$Juv_Mortality, ylim=c(0,0.2))
points(seq(Tmin,Tmax,1), dJTR*exp(coef(mort.J)[1]*(1/TR-1/seq(Tmin,Tmax,1))), type="l", col="blue")


# NLS for adult mortality
mort.A <- nls(Adult_Mortality ~ dATR*exp(A*(1/TR-1/T_K)), data=sp.data,
              start=list(A=1000))
summary(mort.A)
# Plot model fits
plot(sp.data$T_K, sp.data$Adult_Mortality, ylim=c(0,0.1))
points(seq(Tmin,Tmax,1), dATR*exp(coef(mort.A)[1]*(1/TR-1/seq(Tmin,Tmax,1))), type="l", col="blue")



########################### R0 (NET REPRODUCTIVE RATE) ##############################
# NLS
fit.R0 <- nls(R0 ~ R0Topt*exp(-((T_K-ToptR0)^2)/(2*sR0^2)), data=sp.data,
           start=list(R0Topt=100, ToptR0=TR, sR0=10))
summary(fit.R0)
# Plot model fits
plot(sp.data$T_K, sp.data$R0)
points(seq(Tmin,Tmax,1),coef(fit.R0)[1]*exp(-((seq(Tmin,Tmax,1)-coef(fit.R0)[2])^2)/(2*coef(fit.R0)[3]^2)), type="l", col="blue")



########################### r (INTRINSIC GROWTH RATE) ###############################
# estimate all parameters
r <- nls(r ~ xTR*(T_K/TR)*exp(A*(1/TR-1/T_K))/(1+exp(AL*(1/TL-1/T_K))+exp(AH*(1/TH-1/T_K))),
           data=sp.data, start=list(xTR=0.01, A=5000, AL=-150000, AH=200000, TL=290, TH=310))
summary(r)
# Plot model fits
plot(sp.data$T_K, sp.data$r)
points(seq(Tmin,Tmax,1), coef(r)[1]*(seq(Tmin,Tmax,1)/TR)*exp(coef(r)[2]*(1/TR-1/seq(Tmin,Tmax,1)))/
         (1+(exp(coef(r)[3]*(1/coef(r)[5]-1/seq(Tmin,Tmax,1)))+exp(coef(r)[4]*(1/coef(r)[6]-1/seq(Tmin,Tmax,1))))), type="l", col="blue")

# estimate xTR and A
# NOTE: removed data beyond max r
r.mon <- nls(r ~ xTR*T_K/TR*exp(A*(1/TR-1/T_K)), data=sp.data[-c((nrow(sp.data)-2):nrow(sp.data)),],
               start=list(xTR=0.1, A=1000))
summary(r.mon)
# Plot model fits
plot(sp.data$T_K, sp.data$r)
points(seq(Tmin,Tmax,1), coef(r.mon)[1]*(seq(Tmin,Tmax,1)/TR)*exp(coef(r.mon)[2]*(1/TR-1/seq(Tmin,Tmax,1))), type="l", col="blue")

# estimate AL, AH, TL, and TH
# NOTE: if needed, set AH based on NLS below and run code iteratively until AL and AH do not appreciably change)
r.SS <- nls(r ~ coef(r.mon)[1]*(T_K/TR)*exp(coef(r.mon)[2]*(1/TR-1/T_K))/(1+exp(AL*(1/TL-1/T_K))+exp(AH*(1/TH-1/T_K))),
              data=sp.data, start=list(AL=-62000, AH=30000, TL=285.2, TH=308.1))
summary(r.SS)
# Plot model fits
plot(sp.data$T_K, sp.data$r)
points(seq(Tmin,Tmax,1), coef(r.mon)[1]*(seq(Tmin,Tmax,1)/TR)*exp(coef(r.mon)[2]*(1/TR-1/seq(Tmin,Tmax,1)))/
         (1+(exp(coef(r.SS)[1]*(1/coef(r.SS)[3]-1/seq(Tmin,Tmax,1)))+exp(coef(r.SS)[2]*(1/coef(r.SS)[4]-1/seq(Tmin,Tmax,1))))), type="l", col="blue")

# estimate AL and AH separately from TL and TH if needed
kTL <- 292.2
kTH <- 304.7
r.A <- nls(r ~ coef(r.mon)[1]*(T_K/TR)*exp(coef(r.mon)[2]*(1/TR-1/T_K))/(1+exp(AL*(1/kTL-1/T_K))+exp(AH*(1/kTH-1/T_K))),
             data=sp.data, start=list(AL=-50000, AH=230000))
summary(r.A)
r.T <- nls(r ~ coef(r.mon)[1]*(T_K/TR)*exp(coef(r.mon)[2]*(1/TR-1/T_K))/(1+exp(coef(r.A)[1]*(1/TL-1/T_K))+exp(coef(r.A)[2]*(1/TH-1/T_K))),
             data=sp.data, start=list(TL=kTL, TH=kTH))
summary(r.T)
# Plot model fits
plot(sp.data$T_K, sp.data$r)#, ylim=c(0,0.07))
points(seq(Tmin,Tmax,1), coef(r.mon)[1]*(seq(Tmin,Tmax,1)/TR)*exp(coef(r.mon)[2]*(1/TR-1/seq(Tmin,Tmax,1)))/
         (1+(exp(coef(r.A)[1]*(1/coef(r.T)[1]-1/seq(Tmin,Tmax,1)))+exp(coef(r.A)[2]*(1/coef(r.T)[2]-1/seq(Tmin,Tmax,1))))), type="l", col="blue")

# estimate TH and AH separately from TL and AL if needed
kAL <- -100000
r.H <- nls(r ~ coef(r.mon)[1]*(T_K/TR)*exp(coef(r.mon)[2]*(1/TR-1/T_K))/(1+exp(kAL*(1/kTL-1/T_K))+exp(AH*(1/TH-1/T_K))),
             data=sp.data, start=list(TH=kTH, AH=50000))
summary(r.H)
r.AL <- nls(r ~ coef(r.mon)[1]*(T_K/TR)*exp(coef(r.mon)[2]*(1/TR-1/T_K))/(1+exp(AL*(1/kTL-1/T_K))+exp(coef(r.H)[2]*(1/coef(r.H)[1]-1/T_K))),
              data=sp.data, start=list(AL=-100000))
summary(r.AL)
r.TL <- nls(r ~ coef(r.mon)[1]*(T_K/TR)*exp(coef(r.mon)[2]*(1/TR-1/T_K))/(1+exp(coef(r.AL)[1]*(1/TL-1/T_K))+exp(coef(r.H)[2]*(1/coef(r.H)[1]-1/T_K))),
              data=sp.data, start=list(TL=kTL))
summary(r.TL)
# Plot model fits
plot(sp.data$T_K, sp.data$r)
points(seq(Tmin,Tmax,1), coef(r.mon)[1]*(seq(Tmin,Tmax,1)/TR)*exp(coef(r.mon)[2]*(1/TR-1/seq(Tmin,Tmax,1)))/
         (1+(exp(coef(r.AL)[1]*(1/coef(r.TL)[1]-1/seq(Tmin,Tmax,1)))+exp(coef(r.H)[2]*(1/coef(r.H)[1]-1/seq(Tmin,Tmax,1))))), type="l", col="blue")

# Calculate Xmax and Topt for r
Xmax <- 0
Topt <- 0
for(i in 0:100) {
  T <- coef(r.T)[1] + i*(coef(r.T)[2]-coef(r.T)[1])/100
  X.T <- coef(r.mon)[1]*(T/TR)*exp(coef(r.mon)[2]*(1/TR-1/T))/(1+exp(coef(r.A)[1]*(1/coef(r.T)[1]-1/T))+exp(coef(r.A)[2]*(1/coef(r.T)[2]-1/T)))
  if(X.T > Xmax) {Xmax <- X.T
     Topt <- T}}
Xmax
Topt



