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

# USER: enter species name (used in temperature response data.csv)
name <- "Hyperomyzus lactucae"

# Assign species
sp.data <- data[data$Species == name,]
  
  
# Remove columns that do not contain temperature data
sp.data <- sp.data[-c(1:8,12,14,16,18,20,21,23,24,26,27,29,31,32,34,36,37,38)]

# Set some option for nls and plots
Tmin <- 280
Tmax <- 310
TR <- 293



####################################### FECUNDITY ###########################################
# NLS
fec <- nls(Birth_Rate ~ bTopt*exp(-((T_K-Toptb)^2)/(2*sb^2)), data=sp.data,
           start=list(bTopt=5, Toptb=TR, sb=3))
summary(fec)
# Plot model fits
plot(sp.data$T_K, sp.data$Birth_Rate)
points(seq(Tmin,Tmax,1),coef(fec)[1]*exp(-((seq(Tmin,Tmax,1)-coef(fec)[2])^2)/(2*coef(fec)[3]^2)), type="l", col="blue")


###################################### DEVELOPMENT ##########################################
# estimate xTR and A
# NOTE: removed data beyond max development
dev.mon <- nls(Development ~ xTR*T_K/TR*exp(A*(1/TR-1/T_K)), data=sp.data[-c((nrow(sp.data)-0):nrow(sp.data)),],
               start=list(xTR=0.1, A=1000))
summary(dev.mon)
# Plot model fits
plot(sp.data$T_K, sp.data$Development)
points(seq(Tmin,Tmax,1), coef(dev.mon)[1]*(seq(Tmin,Tmax,1)/TR)*exp(coef(dev.mon)[2]*(1/TR-1/seq(Tmin,Tmax,1))), type="l", col="blue")

# estimate AL and AH separately from TL and TH if needed
xTR <- coef(dev.mon)[1]
A <- coef(dev.mon)[2]
kTL <- 280
kTH <- 300
dev.A <- nls(Development ~ xTR*(T_K/TR)*exp(A*(1/TR-1/T_K))/(1+exp(AL*(1/kTL-1/T_K))+exp(AH*(1/kTH-1/T_K))),
              data=sp.data, start=list(AL=-50000, AH=50000))
summary(dev.A)
AL <- coef(dev.A)[1]
AH <- coef(dev.A)[2]
dev.T <- nls(Development ~ xTR*(T_K/TR)*exp(A*(1/TR-1/T_K))/(1+exp(AL*(1/TL-1/T_K))+exp(AH*(1/TH-1/T_K))),
             data=sp.data, start=list(TL=kTL, TH=kTH))
summary(dev.T)
# Plot model fits
plot(sp.data$T_K, sp.data$Development, xlim=c(Tmin,Tmax), ylim=c(0,max(sp.data$Development)+0.05))
points(seq(Tmin,Tmax,1), xTR*(seq(Tmin,Tmax,1)/TR)*exp(A*(1/TR-1/seq(Tmin,Tmax,1)))/
         (1+(exp(AL*(1/coef(dev.T)[1]-1/seq(Tmin,Tmax,1)))+exp(AH*(1/coef(dev.T)[2]-1/seq(Tmin,Tmax,1))))), type="l", col="blue")


# estimate TH and AH separately from TL and AL if needed
kAL <- -50000
dev.H <- nls(Development ~ coef(dev.mon)[1]*(T_K/TR)*exp(coef(dev.mon)[2]*(1/TR-1/T_K))/(1+exp(kAL*(1/kTL-1/T_K))+exp(AH*(1/TH-1/T_K))),
             data=sp.data, start=list(TH=kTH, AH=10000))
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


# Minimum developmental temperature
# NOTE: removed data beyond max development
dev.min <- nls(Development ~ m*T_K+b, data=sp.data[-c((nrow(sp.data)-0):nrow(sp.data)),],
               start=list(m=0.01, b=0))
summary(dev.min)
# Plot model fits
plot(sp.data$T_K, sp.data$Development, xlim=c(280,Tmax))
points(seq(280,Tmax,1), coef(dev.min)[1]*seq(280,Tmax,1)+coef(dev.min)[2], type="l", col="blue")

# Calculate Tmin
Tmin <- (-coef(dev.min)[2]/coef(dev.min)[1])[[1]]
Tmin



#################################### MORTALITY ######################################
# NLS for juvenile mortality
mort.J <- nls(Juv_Mortality ~ xTR*exp(A*(1/TR-1/T_K)), data=sp.data,
              start=list(xTR=0.01, A=10000))
summary(mort.J)
# Plot model fits
plot(sp.data$T_K, sp.data$Juv_Mortality) #, ylim=c(0,0.2))
points(seq(Tmin,Tmax,1), coef(mort.J)[1]*exp(coef(mort.J)[2]*(1/TR-1/seq(Tmin,Tmax,1))), type="l", col="blue")


# NLS for adult mortality
mort.A <- nls(Adult_Mortality ~ xTR*exp(A*(1/TR-1/T_K)), data=sp.data,
              start=list(xTR=0.01, A=10000))
summary(mort.A)
# Plot model fits
plot(sp.data$T_K, sp.data$Adult_Mortality) #, ylim=c(0,0.2))
points(seq(Tmin,Tmax,1), coef(mort.A)[1]*exp(coef(mort.A)[2]*(1/TR-1/seq(Tmin,Tmax,1))), type="l", col="blue")



########################### r (INTRINSIC GROWTH RATE) ###############################
# estimate all parameters
r <- nls(r ~ ifelse(T_K <= Topt, rMax*exp(-1*((T_K-Topt)/(2*sr))^2),
                    rMax*(1 - ((T_K-Topt)/(Topt-Tmax))^2)),
                    data=sp.data, start=list(sr=5, Topt=295, Tmax=300, rMax=1.5))
summary(r)
# Plot model fits
plot(sp.data$T_K, sp.data$r, xlim = c(Tmin,Tmax+5), ylim = c(0,max(sp.data$r)+0.05))
points(seq(Tmin,Tmax+5,1), ifelse(seq(Tmin,Tmax+5,1) <= coef(r)[2], coef(r)[4]*exp(-1*((seq(Tmin,Tmax+5,1)-coef(r)[2])/(2*coef(r)[1]))^2),
                                coef(r)[4]*(1 - ((seq(Tmin,Tmax+5,1)-coef(r)[2])/(coef(r)[2]-coef(r)[3]))^2)), type="l", col="blue")


# set Topt and rMax (NOTE: Topt cannot equal Tmax in nls)
Topt <- sp.data[nrow(sp.data)-1,"T_K"]
rMax <- sp.data[sp.data$T_K==Topt,"r"]
# estimate all other parameters
r <- nls(r ~ ifelse(T_K <= Topt, rMax*exp(-1*((T_K-Topt)/(2*sr))^2),
                    rMax*(1 - ((T_K-Topt)/(Topt-Tmax))^2)),
         data=sp.data, start=list(sr=3, Tmax=303))
summary(r)
# Plot model fits
plot(sp.data$T_K, sp.data$r, xlim=c(Tmin,Tmax), ylim=c(0,rMax))
points(seq(Tmin,Tmax+5,1), ifelse(seq(Tmin,Tmax+5,1) <= Topt, rMax*exp(-1*((seq(Tmin,Tmax+5,1)-Topt)/(2*coef(r)[1]))^2),
                                rMax*(1 - ((seq(Tmin,Tmax+5,1)-Topt)/(Topt-coef(r)[2]))^2)), type="l", col="blue")
rMax
Topt
coef(r)[2]






##################################### EGG DEVELOPMENT ########################################
# estimate xTR and A
# NOTE: removed data beyond max development
dev.mon <- nls(Egg_Dev ~ xTR*T_K/TR*exp(A*(1/TR-1/T_K)), data=sp.data[-c((nrow(sp.data)-0):nrow(sp.data)),],
               start=list(xTR=0.1, A=1000))
summary(dev.mon)
# Plot model fits
plot(sp.data$T_K, sp.data$Egg_Dev)
points(seq(Tmin,Tmax,1), coef(dev.mon)[1]*(seq(Tmin,Tmax,1)/TR)*exp(coef(dev.mon)[2]*(1/TR-1/seq(Tmin,Tmax,1))), type="l", col="blue")

# estimate AL and AH separately from TL and TH if needed
kTL <- 297
kTH <- 311
dev.A <- nls(Egg_Dev ~ coef(dev.mon)[1]*(T_K/TR)*exp(coef(dev.mon)[2]*(1/TR-1/T_K))/(1+exp(AL*(1/kTL-1/T_K))+exp(AH*(1/kTH-1/T_K))),
             data=sp.data, start=list(AL=-50000, AH=50000))
summary(dev.A)
dev.T <- nls(Egg_Dev ~ coef(dev.mon)[1]*(T_K/TR)*exp(coef(dev.mon)[2]*(1/TR-1/T_K))/(1+exp(coef(dev.A)[1]*(1/TL-1/T_K))+exp(coef(dev.A)[2]*(1/TH-1/T_K))),
             data=sp.data, start=list(TL=kTL, TH=kTH))
summary(dev.T)
# Plot model fits
plot(sp.data$T_K, sp.data$Egg_Dev, xlim=c(Tmin,Tmax), ylim=c(0,0.4))
points(seq(Tmin,Tmax,1), coef(dev.mon)[1]*(seq(Tmin,Tmax,1)/TR)*exp(coef(dev.mon)[2]*(1/TR-1/seq(Tmin,Tmax,1)))/
         (1+(exp(coef(dev.A)[1]*(1/coef(dev.T)[1]-1/seq(Tmin,Tmax,1)))+exp(coef(dev.A)[2]*(1/coef(dev.T)[2]-1/seq(Tmin,Tmax,1))))), type="l", col="blue")

# estimate TH and AH separately from TL and AL if needed
kAL <- -100000
dev.H <- nls(Egg_Dev ~ coef(dev.mon)[1]*(T_K/TR)*exp(coef(dev.mon)[2]*(1/TR-1/T_K))/(1+exp(kAL*(1/kTL-1/T_K))+exp(AH*(1/TH-1/T_K))),
             data=sp.data, start=list(TH=kTH, AH=50000))
summary(dev.H)
dev.AL <- nls(Egg_Dev ~ coef(dev.mon)[1]*(T_K/TR)*exp(coef(dev.mon)[2]*(1/TR-1/T_K))/(1+exp(AL*(1/kTL-1/T_K))+exp(coef(dev.H)[2]*(1/coef(dev.H)[1]-1/T_K))),
              data=sp.data, start=list(AL=-100000))
summary(dev.AL)
dev.TL <- nls(Egg_Dev ~ coef(dev.mon)[1]*(T_K/TR)*exp(coef(dev.mon)[2]*(1/TR-1/T_K))/(1+exp(coef(dev.AL)[1]*(1/TL-1/T_K))+exp(coef(dev.H)[2]*(1/coef(dev.H)[1]-1/T_K))),
              data=sp.data, start=list(TL=kTL))
summary(dev.TL)
# Plot model fits
plot(sp.data$T_K, sp.data$Egg_Dev)
points(seq(Tmin,Tmax,1), coef(dev.mon)[1]*(seq(Tmin,Tmax,1)/TR)*exp(coef(dev.mon)[2]*(1/TR-1/seq(Tmin,Tmax,1)))/
         (1+(exp(coef(dev.AL)[1]*(1/coef(dev.TL)[1]-1/seq(Tmin,Tmax,1)))+exp(coef(dev.H)[2]*(1/coef(dev.H)[1]-1/seq(Tmin,Tmax,1))))), type="l", col="blue")



###################################### DEVELOPMENT ##########################################
# Monotonic response
# dev.mon <- nls(Development ~ xTR*T_K/TR*exp(A*(1/TR-1/T_K)), data=sp.data,
#                start=list(xTR=0.01, A=5000))
# summary(dev.mon)
# Plot model fits
# plot(sp.data$T_K, sp.data$Development)
# points(seq(Tmin,Tmax,1), coef(dev.mon)[1]*(seq(Tmin,Tmax,1)/TR)*exp(coef(dev.mon)[2]*(1/TR-1/seq(Tmin,Tmax,1))), type="l", col="blue")


# Sharpe-Schoolfield response
# estimate all parameters
# dev <- nls(Development ~ xTR*(T_K/TR)*exp(A*(1/TR-1/T_K))/(1+exp(AL*(1/TL-1/T_K))+exp(AH*(1/TH-1/T_K))),
#               data=sp.data, start=list(xTR=0.03, A=7000, AL=-30000, AH=270000, TL=286, TH=306))
# summary(dev)
# Plot model fits
# plot(sp.data$T_K, sp.data$Development)
# points(seq(Tmin,Tmax,1), coef(dev)[1]*(seq(Tmin,Tmax,1)/TR)*exp(coef(dev)[2]*(1/TR-1/seq(Tmin,Tmax,1)))/
#          (1+(exp(coef(dev)[3]*(1/coef(dev)[5]-1/seq(Tmin,Tmax,1)))+exp(coef(dev)[4]*(1/coef(dev)[6]-1/seq(Tmin,Tmax,1))))), type="l", col="blue")


# estimate AL, AH, TL, and TH
# NOTE: if needed, set AH based on NLS below and run code iteratively until AL and AH do not appreciably change)
# dev.SS <- nls(Development ~ coef(dev.mon)[1]*(T_K/TR)*exp(coef(dev.mon)[2]*(1/TR-1/T_K))/(1+exp(AL*(1/TL-1/T_K))+exp(AH*(1/TH-1/T_K))),
#               data=sp.data, start=list(AL=-66000, AH=235000, TL=288, TH=309.6))
# summary(dev.SS)
# Plot model fits
# plot(sp.data$T_K, sp.data$Development)
# points(seq(Tmin,Tmax,1), coef(dev.mon)[1]*(seq(Tmin,Tmax,1)/TR)*exp(coef(dev.mon)[2]*(1/TR-1/seq(Tmin,Tmax,1)))/
#          (1+(exp(coef(dev.SS)[1]*(1/coef(dev.SS)[3]-1/seq(Tmin,Tmax,1)))+exp(coef(dev.SS)[2]*(1/coef(dev.SS)[4]-1/seq(Tmin,Tmax,1))))), type="l", col="blue")



########################### LACTIN MODEL OF DEVELOPMENT #############################
# estimate all parameters
#lactin <- nls(Development ~ exp(mTopt*T_K)-exp(mTopt*Tmax-(Tmax-T_K)/mDelta) + mTmin,
#              data=sp.data, start=list(mTopt=0.05, Tmax=305, mDelta=10, mTmin=0.01))
#summary(lactin)
# Plot model fits
#plot(sp.data$T_K, sp.data$Development)
#points(seq(Tmin,Tmax,1), coef(dev.mon)[1]*(seq(Tmin,Tmax,1)/TR)*exp(coef(dev.mon)[2]*(1/TR-1/seq(Tmin,Tmax,1)))/
#         (1+(exp(coef(dev.SS)[1]*(1/coef(dev.SS)[3]-1/seq(Tmin,Tmax,1)))+exp(coef(dev.SS)[2]*(1/coef(dev.SS)[4]-1/seq(Tmin,Tmax,1))))), type="l", col="blue")

# estimate mDelta and mTmin
#lactin.test <- nls(Development ~ exp(0.06*T_K)-exp(0.06*kTH-(kTH-T_K)/mDelta) + mTmin,
#              data=sp.data, start=list(mDelta=1, mTmin=0.01))
#summary(lactin.test)
# Plot model fits
#plot(sp.data$T_K, sp.data$Development)
#points(seq(Tmin,Tmax,1), exp(0.06*seq(Tmin,Tmax,1))-exp(0.06*305-(305-seq(Tmin,Tmax,1))/1) + 0.01, type="l", col="blue")



#################################### EGG MORTALITY ######################################
# Mortality estimated using fit at reference temperature
# NLS for egg mortality
mort.J <- nls(Egg_Mortality ~ xTR*exp(A*(1/TR-1/T_K)), data=sp.data,
              start=list(xTR=0.01, A=10000))
summary(mort.J)
# Plot model fits
plot(sp.data$T_K, sp.data$Egg_Mortality, ylim=c(0,0.2))
points(seq(Tmin,Tmax,1), coef(mort.J)[1]*exp(coef(mort.J)[2]*(1/TR-1/seq(Tmin,Tmax,1))), type="l", col="blue")



# Mortality estimated using data at reference temperature
dJTR <- sp.data[sp.data$T_K==TR,"Juv_Mortality"]
dATR <- sp.data[sp.data$T_K==TR,"Adult_Mortality"]

# NLS for juvenile mortality
mort.J <- nls(Juv_Mortality ~ dJTR*exp(A*(1/TR-1/T_K)), data=sp.data,
              start=list(A=10000))
summary(mort.J)
# Plot model fits
plot(sp.data$T_K, sp.data$Juv_Mortality) #, ylim=c(0,0.2))
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
              start=list(R0Topt=200, ToptR0=TR, sR0=10))
summary(fit.R0)
# Plot model fits
plot(sp.data$T_K, sp.data$R0)
points(seq(Tmin,Tmax,1),coef(fit.R0)[1]*exp(-((seq(Tmin,Tmax,1)-coef(fit.R0)[2])^2)/(2*coef(fit.R0)[3]^2)), type="l", col="blue")

# calculate R0Tmax: temperature at which R0 = 1
R0Tmax <- 0
deltaT <- (5*coef(fit.R0)[3]/1000)[[1]] # temperature range = 5 sigma away from Topt
for(i in 0:1000) {
  if(coef(fit.R0)[1]*exp(-(((coef(fit.R0)[2] + i*deltaT) - coef(fit.R0)[2])^2)/(2*coef(fit.R0)[3])^2) <= 1 & R0Tmax == 0)
  { R0Tmax <- (coef(fit.R0)[2] + i*deltaT)[[1]] }}
R0Tmax


