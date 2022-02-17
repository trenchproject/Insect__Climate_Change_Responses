#############################################################################
#### This R script analyzes thermal performance curve and DDE model data ####
#############################################################################


# Load packages and set working directory
library(tidyr)
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyverse)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# READ IN DATA
# Life history data
r.data <- as.data.frame(read_csv("Predictions/Predictions Dev fitness.csv"))
R0.data <- as.data.frame(read_csv("Predictions/Predictions Dev R0.csv"))
b.data <- as.data.frame(read_csv("Predictions/Predictions Dev birth.csv"))
tau.data <- as.data.frame(read_csv("Predictions/Predictions Dev development.csv"))
s.data <- as.data.frame(read_csv("Predictions/Predictions Dev survival.csv"))
L.data <- as.data.frame(read_csv("Predictions/Predictions Dev longevity.csv"))
f.data <- as.data.frame(read_csv("Predictions/Predictions Dev fecundity.csv"))
R.data <- as.data.frame(read_csv("Predictions/Predictions Dev recruitment.csv"))
# Population dynamics data
pop.data <- as.data.frame(read_csv("Predictions/Predictions population dynamics.csv"))
# Life history trait data
LH.data <- as.data.frame(read_csv("Temperature response parameters.csv"))
# Temperature data
temp.data <- as.data.frame(read_csv("Temperature parameters Tave.csv"))
# Extinction results
# results.m <- as.data.frame(read_csv("Extinction meanT.csv"))
# results.a <- as.data.frame(read_csv("Extinction amplT.csv"))
# results.b <- as.data.frame(read_csv("Extinction both.csv"))


# EXCLUDE DATA
# Macrolophus pygmaeus (only predator and separate thermal responses for each prey)
# Bemisia argentifollii and Diaphorina citri (different suborder)
r.data <- r.data[-c(12,13,15,18),]
R0.data <- R0.data[-c(12,13,15,18),]
b.data <- b.data[-c(12,13,15,18),]
tau.data$active.h <- pop.data$active.h # add activity period data to developmental data
tau.data$active.f <- pop.data$active.f # add activity period data to developmental data
tau.data <- tau.data[-c(12,13,15,18),]
s.data <- s.data[-c(12,13,15,18),]
L.data <- L.data[-c(12,13,15,18),]
f.data <- f.data[-c(12,13,15,18),]
R.data <- R.data[-c(12,13,15,18),]
LH.data <- LH.data[-c(12,13,15,18),]
# Macrosiphum euphorbiae Canada and Brevicoryne brassicae (went extinct)
pop.data <- pop.data[-c(12,13,15,18,20,25),]
# results.m <- results.m[-c(12,13,15,18),]
# results.a <- results.a[-c(12,13,15,18),]
# results.b <- results.b[-c(12,13,15,18),]


# SCALE LOWEST FITNESS CHANGE TO -1
r.data$delta.TPC <- pmax(-1, r.data$delta.TPC)
r.data$delta.model <- pmax(-1, r.data$delta.model)



######################################## STATISTICS #########################################
# FITNESS
# Change in relative fitness
r.delta <- lm(delta.model ~ delta.TPC, data=r.data)
summary(r.delta) # significant
# Model vs Latitude
# linear
r.lat <- lm(delta.model ~ Latitude, data=r.data)
summary(r.lat) # non-significant
# non-linear
#r.lat2 <- nls(delta.model ~ a + b*Latitude + c*Latitude^2, data=r.data, start=list(a=1, b=-0.1, c=1))
#summary(r.lat2) # non-significant

# R0
# Change in relative R0
R0.delta <- lm(delta.model ~ delta.TPC, data=R0.data)
summary(R0.delta) # significant!
# Model vs Latitude
R0.lat <- lm(delta.model ~ Latitude, data=R0.data)
summary(R0.lat) # non-significant

# BIRTH RATE
# Model vs Latitude
b.lat <- lm(delta.model ~ Latitude, data=b.data)
summary(b.lat) # non-significant

# DEVELOPMENT TIME
# Change in development time
tau.delta <- lm(delta.model ~ delta.TPC, data=tau.data)
summary(tau.delta)  # significant!
# Model vs Latitude
tau.lat <- lm(delta.model ~ Latitude, data=tau.data)
summary(tau.lat) # significant!

# SURVIVAL
# Change in survival
s.delta <- lm(delta.model ~ delta.TPC, data=s.data)
summary(s.delta)  # significant!
# Model vs Latitude
s.lat <- lm(delta.model ~ Latitude, data=s.data)
summary(s.lat) # non-significant

# ADULT LONGEVITY
# Model vs Latitude
L.lat <- lm(delta.model ~ Latitude, data=L.data)
summary(L.lat) # marginally-significant

# LIFETIME FECUNDITY
# Model vs Latitude
#f.lat <- lm(delta.model ~ Latitude, data=f.data)
#summary(f.lat) # significant!

# ADULT RECRUITMENT
# Change in recruitment
#R.delta <- lm(delta.model ~ delta.TPC, data=R.data)
#summary(R.delta) # significant!
# Model vs Latitude
#R.lat <- lm(delta.model ~ Latitude, data=R.data)
#summary(R.lat) # non-significant


# POPULATION DYNAMICS
# Mean density vs Latitude
mean.lat <- lm(delta.mean ~ Latitude, data=pop.data)
summary(mean.lat) # non-significant
# CV of density vs Latitude
CV.lat <- lm(delta.CV ~ Latitude, data=pop.data)
summary(CV.lat) # significant!
# Active period vs Latitude (NOTE: non-significant for temperate species only)
active.lat <- lm(delta.active ~ Latitude, data=pop.data) #[pop.data$Habitat == "Temperate",])
summary(active.lat) # significant!


# LIFE HISTORY TRAITS
# BIRTH RATE
# bTopt vs Latitude
bTopt.lat <- lm(bTopt ~ Latitude, data=LH.data)
summary(bTopt.lat) # significant!
plot(LH.data$Latitude,LH.data$bTopt, ylim=c(0,10))
# Toptb vs Latitude
Toptb.lat <- lm(Toptb ~ Latitude, data=LH.data)
summary(Toptb.lat) # significant!
plot(LH.data$Latitude,LH.data$Toptb)
# sb vs Latitude
sb.lat <- lm(sb ~ Latitude, data=LH.data)
summary(sb.lat) # significant!
plot(LH.data$Latitude,LH.data$sb)
# Toptb - meanT vs Latitude
LH.data$delta_b <- LH.data$Toptb - (temp.data$meanT.f + temp.data$delta_mean.f*75*365)
delta.b.lat <- lm(delta_b ~ Latitude, data=LH.data)
summary(delta.b.lat) # significant!
plot(LH.data$Latitude,LH.data$delta_b)

# DEVELOPMENT RATE
# mTR vs Latitude
mTR.lat <- lm(mTR ~ Latitude, data=LH.data)
summary(mTR.lat) # non-significant
plot(LH.data$Latitude,LH.data$mTR)
# AL vs Latitude
AL.lat <- lm(AL ~ Latitude, data=LH.data)
summary(AL.lat) # non-significant
plot(LH.data$Latitude,LH.data$AL)
# AH vs Latitude
AH.lat <- lm(AH ~ Latitude, data=LH.data)
summary(AH.lat) # non-significant
plot(LH.data$Latitude,LH.data$AH)
# TL vs Latitude
TL.lat <- lm(TL ~ Latitude, data=LH.data)
summary(TL.lat) # significant!
plot(LH.data$Latitude,LH.data$TL)
# TH vs Latitude
TH.lat <- lm(TH ~ Latitude, data=LH.data)
summary(TH.lat) # non-significant
plot(LH.data$Latitude,LH.data$TH)
# Tmin vs Latitude
Tmin.lat <- lm(Tmin ~ Latitude, data=LH.data)
summary(Tmin.lat) # significant!
plot(LH.data$Latitude,LH.data$Tmin)
# Topt vs Latitude
Topt.lat <- lm(Topt ~ Latitude, data=LH.data)
summary(Topt.lat) # significant!
plot(LH.data$Latitude,LH.data$Topt)
# Tmax vs Latitude
Tmax.lat <- lm(Tmax ~ Latitude, data=LH.data)
summary(Tmax.lat) # non-significant
plot(LH.data$Latitude,LH.data$Tmax)
# Topt - meanT vs Latitude
LH.data$delta_Topt <- LH.data$Topt - (temp.data$meanT.f  + temp.data$delta_mean.f*75*365)
delta.Topt.lat <- lm(delta_Topt ~ Latitude, data=LH.data)
summary(delta.Topt.lat) # significant!
plot(LH.data$Latitude,LH.data$delta_Topt)
# Topt - meanT vs Latitude
LH.data$delta_Topt <- LH.data$Topt - (temp.data$meanT.f  + temp.data$delta_mean.f*75*365 +
                                        abs(temp.data$amplT.f)  + abs(temp.data$delta_ampl.f*75*365))
delta.Topt.lat <- lm(delta_Topt ~ Latitude, data=LH.data)
summary(delta.Topt.lat) # non-significant
plot(LH.data$Latitude,LH.data$delta_Topt)

# JUVENILE MORTALITY RATE
# dJTR vs Latitude
dJTR.lat <- lm(dJTR ~ Latitude, data=LH.data)
summary(dJTR.lat) # non-significant
plot(LH.data$Latitude,LH.data$dJTR)
# AdJ vs Latitude
AdJ.lat <- lm(AdJ ~ Latitude, data=LH.data)
summary(AdJ.lat) # non-significant
plot(LH.data$Latitude,LH.data$AdJ)
# dJ(Tmean) vs Latitude
LH.data$dJ_Tmean <- LH.data$dJTR*exp(LH.data$AdJ*(1/LH.data$TR-1/(temp.data$meanT.f + temp.data$delta_mean.f*75*365)))
dJ.Tmean.lat <- lm(dJ_Tmean ~ Latitude, data=LH.data)
summary(dJ.Tmean.lat) # significant
plot(LH.data$Latitude,LH.data$dJ_Tmean)

# ADULT MORTALITY RATE
# dATR vs Latitude
dATR.lat <- lm(dATR ~ Latitude, data=LH.data)
summary(dATR.lat) # marginally-significant
plot(LH.data$Latitude,LH.data$dATR)
# AdA vs Latitude
AdA.lat <- lm(AdA ~ Latitude, data=LH.data)
summary(AdA.lat) # non-significant
plot(LH.data$Latitude,LH.data$AdA)
# dA(Tmean) vs Latitude
LH.data$dA_Tmean <- LH.data$dATR*exp(LH.data$AdA*(1/LH.data$TR-1/(temp.data$meanT.f + temp.data$delta_mean.f*75*365)))
dA.Tmean.lat <- lm(dA_Tmean ~ Latitude, data=LH.data)
summary(dA.Tmean.lat) # marginally-significant
plot(LH.data$Latitude,LH.data$dA_Tmean)

# NUMBER OF GENERATIONS
# historical period
LH.data$Gen.h <- (tau.data$active.h / tau.data$Model.h)
Gen.h.lat <- lm(Gen.h ~ Latitude, data=LH.data)
summary(Gen.h.lat) # non-significant
#plot(LH.data$Latitude,LH.data$Gen.h)
# future period
LH.data$Gen.f <- (tau.data$active.f / tau.data$Model.f)
Gen.f.lat <- lm(Gen.f ~ Latitude, data=LH.data)
summary(Gen.f.lat) # non-significant
#plot(LH.data$Latitude,LH.data$Gen.f)
# change in number of generations
LH.data$delta_Gen <- (LH.data$Gen.f - LH.data$Gen.h) / LH.data$Gen.h
delta.Gen.lat <- lm(delta_Gen ~ Latitude, data=LH.data[-c(16,21),])
summary(delta.Gen.lat) # non-significant
plot(LH.data$Latitude[-c(16,21)],LH.data$delta_Gen[-c(16,21)])


# TEMPERATURE
# delta_meanT vs Latitude
Tmean.lat <- lm(delta_mean.f ~ Latitude, data=temp.data)
summary(Tmean.lat) # non-significant
# delta_amplT temperature vs Latitude
Tampl.lat <- lm(delta_ampl.f ~ Latitude, data=temp.data)
summary(Tampl.lat) # significant!
# Change in mean temperature vs Latitude
temp.data$mean_ch <- (temp.data$meanT.f + temp.data$delta_mean.f*75*365)- temp.data$meanT.h
mean.ch.lat <- lm(mean_ch ~ Latitude, data=temp.data)
summary(mean.ch.lat) # non-significant
# Change in ampl temperature vs Latitude
temp.data$ampl_ch <- (temp.data$amplT.f + temp.data$delta_ampl.f*75*365)- temp.data$amplT.h
ampl.ch.lat <- lm(ampl_ch ~ Latitude, data=temp.data)
summary(ampl.ch.lat) # marginally-significant


# EXTINCTION
# Increase in mean temperature
# TPC vs model
#mean.delta <- lm(Model ~ TPC, data=results.m)
#summary(mean.delta) # significant!
# Model vs Latitude
#mean.lat <- lm(Model ~ Latitude, data=results.m)
#summary(mean.lat) # significant!

# Increase in temperature amplitude
# TPC vs model
#ampl.delta <- lm(Model ~ TPC, data=results.a)
#summary(ampl.delta) # significant
# Model vs Latitude
#ampl.lat <- lm(Model ~ Latitude, data=results.a)
#summary(ampl.lat) # non-significant

# Increase in temperature mean and amplitude
# TPC vs model
#both.delta <- lm(Model ~ TPC, data=results.b)
#summary(both.delta) # non-significant
# Model vs Latitude
#both.lat <- lm(Model ~ Latitude, data=results.b)
#summary(both.lat) # non-significant



########################################### PLOTS ###########################################
# RELATIVE FITNESS
# Model vs TPCs
Xmin <- -1
Xmax <- 0.5
Ymin <- -1
Ymax <- 0.5
#dev.new(width=3, height=3, unit="in")
plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="TPC", ylab="Model")
polygon(c(2*Xmin,0,2*Xmax,2*Xmax),c(2*Xmin,0,-2*Xmax,2*Ymin), col = "#E2E2E2", border = NA)
polygon(c(2*Xmin,0,2*Xmax),c(-2*Xmin,0,2*Xmax), col = "#E2E2E2", border = NA)
abline(0, 1, col="gray", lwd=3)
abline(0, -1, col="gray", lwd=1)
abline(0, 0, col="gray", lwd=3, lty="longdash")
abline(v = 0, col="gray", lwd=3, lty="longdash")
points(r.data[r.data$Habitat=="Tropical","delta.TPC"], r.data[r.data$Habitat=="Tropical","delta.model"], pch=19, cex=1.5, col="#FFB000") # orange
points(r.data[r.data$Habitat=="Subtropical","delta.TPC"], r.data[r.data$Habitat=="Subtropical","delta.model"], pch=19, cex=1.5, col="#40B0A6") # teal
points(r.data[r.data$Habitat=="Mediterranean","delta.TPC"], r.data[r.data$Habitat=="Mediterranean","delta.model"], pch=19, cex=1.5, col="#40B0A6") # teal
points(r.data[r.data$Habitat=="Temperate","delta.TPC"], r.data[r.data$Habitat=="Temperate","delta.model"], pch=19, cex=1.5, col="#785EF0") # purple
points(seq(2*Xmin,2*Xmax,0.1), coef(r.delta)[2]*seq(2*Xmin,2*Xmax,0.1)+coef(r.delta)[1], type="l", lwd=3, col="black")

# Model vs latitude
Xmin <- 0
Xmax <- 60
Ymin <- -1
Ymax <- 0.5
#dev.new(width=3, height=3, unit="in")
plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Latitude", ylab="Model")
abline(0, 0, col="gray", lwd=3, lty="longdash")
points(r.data[r.data$Habitat=="Tropical","Latitude"], r.data[r.data$Habitat=="Tropical","delta.model"], pch=19, cex=1.5, col="#FFB000") # orange
points(r.data[r.data$Habitat=="Subtropical","Latitude"], r.data[r.data$Habitat=="Subtropical","delta.model"], pch=19, cex=1.5, col="#6FD012") # green
points(r.data[r.data$Habitat=="Mediterranean","Latitude"], r.data[r.data$Habitat=="Mediterranean","delta.model"], pch=19, cex=1.5, col="#6FD012") # green
points(r.data[r.data$Habitat=="Temperate","Latitude"], r.data[r.data$Habitat=="Temperate","delta.model"], pch=19, cex=1.5, col="#785EF0") # purple
points(seq(2*Xmin,2*Xmax,1), coef(r.lat)[2]*seq(2*Xmin,2*Xmax,1) + coef(r.lat)[1], type="l", lwd=3, col="black", lty="longdash")
#points(seq(Xmin,Xmax,1), coef(r.lat2)[1] + coef(r.lat2)[2]*seq(Xmin,Xmax,1) + coef(r.lat2)[3]*seq(Xmin,Xmax,1)^2, type="l", lwd=3, col="black")


# R0
# Model vs TPCs
Xmin <- -1
Xmax <- 0.2
Ymin <- -1
Ymax <- 0.2
#dev.new(width=3, height=3, unit="in")
plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="TPC", ylab="Model")
polygon(c(2*Xmin,0,2*Xmax,2*Xmax),c(2*Xmin,0,-2*Xmax,2*Ymin), col = "lightgray", border = NA)
polygon(c(2*Xmin,0,2*Xmax),c(-2*Xmin,0,2*Xmax), col = "lightgray", border = NA)
abline(0, 1, col="gray", lwd=3)
abline(0, -1, col="gray", lwd=1)
abline(0, 0, col="gray", lwd=3, lty="longdash")
abline(v = 0, col="gray", lwd=3, lty="longdash")
points(R0.data[R0.data$Habitat=="Tropical","delta.TPC"], R0.data[R0.data$Habitat=="Tropical","delta.model"], pch=19, cex=1.5, col="#FFB000") # orange
points(R0.data[R0.data$Habitat=="Subtropical","delta.TPC"], R0.data[R0.data$Habitat=="Subtropical","delta.model"], pch=19, cex=1.5, col="#6FD012") # green
points(R0.data[R0.data$Habitat=="Mediterranean","delta.TPC"], R0.data[R0.data$Habitat=="Mediterranean","delta.model"], pch=19, cex=1.5, col="#6FD012") # green
points(R0.data[R0.data$Habitat=="Temperate","delta.TPC"], R0.data[R0.data$Habitat=="Temperate","delta.model"], pch=19, cex=1.5, col="#785EF0") # purple
points(seq(2*Xmin,2*Xmax,0.1), coef(R0.delta)[2]*seq(2*Xmin,2*Xmax,0.1)+coef(R0.delta)[1], type="l", lwd=3, col="black")

# Model vs latitude
Xmin <- 0
Xmax <- 60
Ymin <- -0.8
Ymax <- 0.2
#dev.new(width=3, height=3, unit="in")
plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Latitude", ylab="Model")
abline(0, 0, col="gray", lwd=3, lty="longdash")
points(R0.data[R0.data$Habitat=="Tropical","Latitude"], R0.data[R0.data$Habitat=="Tropical","delta.model"], pch=19, cex=1.5, col="#FFB000") # orange
points(R0.data[R0.data$Habitat=="Subtropical","Latitude"], R0.data[R0.data$Habitat=="Subtropical","delta.model"], pch=19, cex=1.5, col="#6FD012") # green
points(R0.data[R0.data$Habitat=="Mediterranean","Latitude"], R0.data[R0.data$Habitat=="Mediterranean","delta.model"], pch=19, cex=1.5, col="#6FD012") # green
points(R0.data[R0.data$Habitat=="Temperate","Latitude"], R0.data[R0.data$Habitat=="Temperate","delta.model"], pch=19, cex=1.5, col="#785EF0") # purple
points(seq(2*Xmin,2*Xmax,1), coef(R0.lat)[2]*seq(2*Xmin,2*Xmax,1) + coef(R0.lat)[1], type="l", lwd=3, col="black", lty="longdash")


# BIRTH RATE
# Model vs TPCs
Xmin <- -0.6
Xmax <- 0.2
Ymin <- -0.6
Ymax <- 0.2
#dev.new(width=3, height=3, unit="in")
plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="TPC", ylab="Model")
polygon(c(2*Xmin,0,2*Xmax,2*Xmax),c(2*Xmin,0,-2*Xmax,2*Ymin), col = "lightgray", border = NA)
polygon(c(2*Xmin,0,2*Xmax),c(-2*Xmin,0,2*Xmax), col = "lightgray", border = NA)
abline(0, 1, col="gray", lwd=3)
abline(0, -1, col="gray", lwd=1)
abline(0, 0, col="gray", lwd=3, lty="longdash")
abline(v = 0, col="gray", lwd=3, lty="longdash")
points(b.data[b.data$Habitat=="Tropical","delta.TPC"], b.data[b.data$Habitat=="Tropical","delta.model"], pch=19, cex=1.5, col="#FFB000") # orange
points(b.data[b.data$Habitat=="Subtropical","delta.TPC"], b.data[b.data$Habitat=="Subtropical","delta.model"], pch=19, cex=1.5, col="#6FD012") # green
points(b.data[b.data$Habitat=="Mediterranean","delta.TPC"], b.data[b.data$Habitat=="Mediterranean","delta.model"], pch=19, cex=1.5, col="#6FD012") # green
points(b.data[b.data$Habitat=="Temperate","delta.TPC"], b.data[b.data$Habitat=="Temperate","delta.model"], pch=19, cex=1.5, col="#785EF0") # purple

# Model vs latitude
Xmin <- 0
Xmax <- 60
Ymin <- -0.6
Ymax <- 0.2
#dev.new(width=3, height=3, unit="in")
plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Latitude", ylab="Model")
abline(0, 0, col="gray", lwd=3, lty="longdash")
points(b.data[b.data$Habitat=="Tropical","Latitude"], b.data[b.data$Habitat=="Tropical","delta.model"], pch=19, cex=1.5, col="#FFB000") # orange
points(b.data[b.data$Habitat=="Subtropical","Latitude"], b.data[b.data$Habitat=="Subtropical","delta.model"], pch=19, cex=1.5, col="#6FD012") # green
points(b.data[b.data$Habitat=="Mediterranean","Latitude"], b.data[b.data$Habitat=="Mediterranean","delta.model"], pch=19, cex=1.5, col="#6FD012") # green
points(b.data[b.data$Habitat=="Temperate","Latitude"], b.data[b.data$Habitat=="Temperate","delta.model"], pch=19, cex=1.5, col="#785EF0") # purple
points(seq(2*Xmin,2*Xmax,1), coef(b.lat)[2]*seq(2*Xmin,2*Xmax,1) + coef(b.lat)[1], type="l", lwd=3, col="black", lty="longdash")


# DEVELOPMENT TIME
# Model vs TPCs
Xmin <- -0.3
Xmax <- 0.2
Ymin <- -0.3
Ymax <- 0.1
#dev.new(width=3, height=3, unit="in")
plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="TPC", ylab="Model")
polygon(c(2*Xmin,0,2*Xmax,2*Xmax),c(2*Xmin,0,-2*Xmax,2*Ymin), col = "lightgray", border = NA)
polygon(c(2*Xmin,0,2*Xmax),c(-2*Xmin,0,2*Xmax), col = "lightgray", border = NA)
abline(0, 1, col="gray", lwd=3)
abline(0, -1, col="gray", lwd=1)
abline(0, 0, col="gray", lwd=3, lty="longdash")
abline(v = 0, col="gray", lwd=3, lty="longdash")
points(tau.data[tau.data$Habitat=="Tropical","delta.TPC"], tau.data[tau.data$Habitat=="Tropical","delta.model"], pch=19, cex=1.5, col="#FFB000") # orange
points(tau.data[tau.data$Habitat=="Subtropical","delta.TPC"], tau.data[tau.data$Habitat=="Subtropical","delta.model"], pch=19, cex=1.5, col="#6FD012") # green
points(tau.data[tau.data$Habitat=="Mediterranean","delta.TPC"], tau.data[tau.data$Habitat=="Mediterranean","delta.model"], pch=19, cex=1.5, col="#6FD012") # green
points(tau.data[tau.data$Habitat=="Temperate","delta.TPC"], tau.data[tau.data$Habitat=="Temperate","delta.model"], pch=19, cex=1.5, col="#785EF0") # purple
points(seq(2*Xmin,2*Xmax,0.1), coef(tau.delta)[2]*seq(2*Xmin,2*Xmax,0.1)+coef(tau.delta)[1], type="l", lwd=3, col="black")

# Model vs latitude
Xmin <- 0
Xmax <- 60
Ymin <- -0.6
Ymax <- 0
#dev.new(width=3, height=3, unit="in")
plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Latitude", ylab="Model")
abline(0, 0, col="gray", lwd=3, lty="longdash")
points(tau.data[tau.data$Habitat=="Tropical","Latitude"], tau.data[tau.data$Habitat=="Tropical","delta.model"], pch=19, cex=1.5, col="#FFB000") # orange
points(tau.data[tau.data$Habitat=="Subtropical","Latitude"], tau.data[tau.data$Habitat=="Subtropical","delta.model"], pch=19, cex=1.5, col="#6FD012") # green
points(tau.data[tau.data$Habitat=="Mediterranean","Latitude"], tau.data[tau.data$Habitat=="Mediterranean","delta.model"], pch=19, cex=1.5, col="#6FD012") # green
points(tau.data[tau.data$Habitat=="Temperate","Latitude"], tau.data[tau.data$Habitat=="Temperate","delta.model"], pch=19, cex=1.5, col="#785EF0") # purple
points(seq(2*Xmin,2*Xmax,1), coef(tau.lat)[2]*seq(2*Xmin,2*Xmax,1) + coef(tau.lat)[1], type="l", lwd=3, col="black")


# SURVIVAL
# Model vs TPCs
Xmin <- -0.6
Xmax <- 0.2
Ymin <- -1
Ymax <- 0.2
#dev.new(width=3, height=3, unit="in")
plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="TPC", ylab="Model")
polygon(c(2*Xmin,0,2*Xmax,2*Xmax),c(2*Xmin,0,-2*Xmax,2*Ymin), col = "lightgray", border = NA)
polygon(c(2*Xmin,0,2*Xmax),c(-2*Xmin,0,2*Xmax), col = "lightgray", border = NA)
abline(0, 1, col="gray", lwd=3)
abline(0, -1, col="gray", lwd=1)
abline(0, 0, col="gray", lwd=3, lty="longdash")
abline(v = 0, col="gray", lwd=3, lty="longdash")
points(s.data[s.data$Habitat=="Tropical","delta.TPC"], s.data[s.data$Habitat=="Tropical","delta.model"], pch=19, cex=1.5, col="#FFB000") # orange
points(s.data[s.data$Habitat=="Subtropical","delta.TPC"], s.data[s.data$Habitat=="Subtropical","delta.model"], pch=19, cex=1.5, col="#6FD012") # green
points(s.data[s.data$Habitat=="Mediterranean","delta.TPC"], s.data[s.data$Habitat=="Mediterranean","delta.model"], pch=19, cex=1.5, col="#6FD012") # green
points(s.data[s.data$Habitat=="Temperate","delta.TPC"], s.data[s.data$Habitat=="Temperate","delta.model"], pch=19, cex=1.5, col="#785EF0") # purple
points(seq(2*Xmin,2*Xmax,0.1), coef(s.delta)[2]*seq(2*Xmin,2*Xmax,0.1)+coef(s.delta)[1], type="l", lwd=3, col="black")

# Model vs latitude
Xmin <- 0
Xmax <- 60
Ymin <- -1
Ymax <- 0.2
#dev.new(width=3, height=3, unit="in")
plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Latitude", ylab="Model")
abline(0, 0, col="gray", lwd=3, lty="longdash")
points(s.data[s.data$Habitat=="Tropical","Latitude"], s.data[s.data$Habitat=="Tropical","delta.model"], pch=19, cex=1.5, col="#FFB000") # orange
points(s.data[s.data$Habitat=="Subtropical","Latitude"], s.data[s.data$Habitat=="Subtropical","delta.model"], pch=19, cex=1.5, col="#6FD012") # green
points(s.data[s.data$Habitat=="Mediterranean","Latitude"], s.data[s.data$Habitat=="Mediterranean","delta.model"], pch=19, cex=1.5, col="#6FD012") # green
points(s.data[s.data$Habitat=="Temperate","Latitude"], s.data[s.data$Habitat=="Temperate","delta.model"], pch=19, cex=1.5, col="#785EF0") # purple
points(seq(2*Xmin,2*Xmax,1), coef(s.lat)[2]*seq(2*Xmin,2*Xmax,1) + coef(s.lat)[1], type="l", lwd=3, col="black", lty="longdash")


# ADULT LONGEVITY
# Model vs TPCs
Xmin <- -0.4
Xmax <- 0.1
Ymin <- -0.4
Ymax <- 0.1
#dev.new(width=3, height=3, unit="in")
plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="TPC", ylab="Model")
polygon(c(2*Xmin,0,2*Xmax,2*Xmax),c(2*Xmin,0,-2*Xmax,2*Ymin), col = "lightgray", border = NA)
polygon(c(2*Xmin,0,2*Xmax),c(-2*Xmin,0,2*Xmax), col = "lightgray", border = NA)
abline(0, 1, col="gray", lwd=3)
abline(0, -1, col="gray", lwd=1)
abline(0, 0, col="gray", lwd=3, lty="longdash")
abline(v = 0, col="gray", lwd=3, lty="longdash")
points(L.data[L.data$Habitat=="Tropical","delta.TPC"], L.data[L.data$Habitat=="Tropical","delta.model"], pch=19, cex=1.5, col="#FFB000") # orange
points(L.data[L.data$Habitat=="Subtropical","delta.TPC"], L.data[L.data$Habitat=="Subtropical","delta.model"], pch=19, cex=1.5, col="#6FD012") # green
points(L.data[L.data$Habitat=="Mediterranean","delta.TPC"], L.data[L.data$Habitat=="Mediterranean","delta.model"], pch=19, cex=1.5, col="#6FD012") # green
points(L.data[L.data$Habitat=="Temperate","delta.TPC"], L.data[L.data$Habitat=="Temperate","delta.model"], pch=19, cex=1.5, col="#785EF0") # purple

# Model vs latitude
Xmin <- 0
Xmax <- 60
Ymin <- -1.5
Ymax <- 0.5
#dev.new(width=3, height=3, unit="in")
plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Latitude", ylab="Model")
abline(0, 0, col="gray", lwd=3, lty="longdash")
points(L.data[L.data$Habitat=="Tropical","Latitude"], L.data[L.data$Habitat=="Tropical","delta.model"], pch=19, cex=1.5, col="#FFB000") # orange
points(L.data[L.data$Habitat=="Subtropical","Latitude"], L.data[L.data$Habitat=="Subtropical","delta.model"], pch=19, cex=1.5, col="#6FD012") # green
points(L.data[L.data$Habitat=="Mediterranean","Latitude"], L.data[L.data$Habitat=="Mediterranean","delta.model"], pch=19, cex=1.5, col="#6FD012") # green
points(L.data[L.data$Habitat=="Temperate","Latitude"], L.data[L.data$Habitat=="Temperate","delta.model"], pch=19, cex=1.5, col="#785EF0") # purple
points(seq(2*Xmin,2*Xmax,1), coef(L.lat)[2]*seq(2*Xmin,2*Xmax,1) + coef(L.lat)[1], type="l", lwd=3, col="black")



# POPULATION DYNAMICS
# Mean density vs latitude
Xmin <- 0
Xmax <- 60
Ymin <- -0.6
Ymax <- 0.6
#dev.new(width=3, height=3, unit="in")
plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Latitude", ylab="Model")
abline(0, 0, col="gray", lwd=3, lty="longdash")
points(seq(2*Xmin,2*Xmax,1), coef(mean.lat)[2]*seq(2*Xmin,2*Xmax,1) + coef(mean.lat)[1], type="l", lwd=3, col="black", lty="longdash")
points(pop.data[pop.data$Habitat=="Tropical","Latitude"], pop.data[pop.data$Habitat=="Tropical","delta.mean"], pch=19, cex=1.5, col="#FFB000") # orange
points(pop.data[pop.data$Habitat=="Subtropical","Latitude"], pop.data[pop.data$Habitat=="Subtropical","delta.mean"], pch=19, cex=1.5, col="#6FD012") # green
points(pop.data[pop.data$Habitat=="Mediterranean","Latitude"], pop.data[pop.data$Habitat=="Mediterranean","delta.mean"], pch=19, cex=1.5, col="#6FD012") # green
points(pop.data[pop.data$Habitat=="Temperate","Latitude"], pop.data[pop.data$Habitat=="Temperate","delta.mean"], pch=19, cex=1.5, col="#785EF0") # purple

# CV of density vs latitude
Xmin <- 0
Xmax <- 60
Ymin <- -0.6
Ymax <- 0.8
#dev.new(width=3, height=3, unit="in")
plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Latitude", ylab="Model")
abline(0, 0, col="gray", lwd=3, lty="longdash")
points(seq(2*Xmin,2*Xmax,1), coef(CV.lat)[2]*seq(2*Xmin,2*Xmax,1) + coef(CV.lat)[1], type="l", lwd=3, col="black")
points(pop.data[pop.data$Habitat=="Tropical","Latitude"], pop.data[pop.data$Habitat=="Tropical","delta.CV"], pch=19, cex=1.5, col="#FFB000") # orange
points(pop.data[pop.data$Habitat=="Subtropical","Latitude"], pop.data[pop.data$Habitat=="Subtropical","delta.CV"], pch=19, cex=1.5, col="#6FD012") # green
points(pop.data[pop.data$Habitat=="Mediterranean","Latitude"], pop.data[pop.data$Habitat=="Mediterranean","delta.CV"], pch=19, cex=1.5, col="#6FD012") # green
points(pop.data[pop.data$Habitat=="Temperate","Latitude"], pop.data[pop.data$Habitat=="Temperate","delta.CV"], pch=19, cex=1.5, col="#785EF0") # purple

# Active period vs latitude
Xmin <- 0
Xmax <- 60
Ymin <- 0
Ymax <- 0.6
#dev.new(width=3, height=3, unit="in")
plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Latitude", ylab="Model")
abline(0, 0, col="gray", lwd=3, lty="longdash")
points(seq(2*Xmin,2*Xmax,1), coef(active.lat)[2]*seq(2*Xmin,2*Xmax,1) + coef(active.lat)[1], type="l", lwd=3, col="black")
points(pop.data[pop.data$Habitat=="Tropical","Latitude"], pop.data[pop.data$Habitat=="Tropical","delta.active"], pch=19, cex=1.5, col="#FFB000") # orange
points(pop.data[pop.data$Habitat=="Subtropical","Latitude"], pop.data[pop.data$Habitat=="Subtropical","delta.active"], pch=19, cex=1.5, col="#6FD012") # green
points(pop.data[pop.data$Habitat=="Mediterranean","Latitude"], pop.data[pop.data$Habitat=="Mediterranean","delta.active"], pch=19, cex=1.5, col="#6FD012") # green
points(pop.data[pop.data$Habitat=="Temperate","Latitude"], pop.data[pop.data$Habitat=="Temperate","delta.active"], pch=19, cex=1.5, col="#785EF0") # purple



# # EXTINCTION
# # Change in mean temperature
# # Model vs TPCs
# Xmin <- 0
# Xmax <- 25
# Ymin <- 0
# Ymax <- 25
# #dev.new(width=3, height=3, unit="in")
# plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="TPC", ylab="Model")
# polygon(c(-1,2*Xmax,2*Xmax),c(-1,-1,2*Xmax), col = "lightgray", border = NA)
# abline(0, 1, col="gray", lwd=3)
# points(results.m[results.m$Habitat=="Tropical","TPC"], results.m[results.m$Habitat=="Tropical","Model"], pch=19, cex=1.5, col="#FFB000") # orange
# points(results.m[results.m$Habitat=="Subtropical","TPC"], results.m[results.m$Habitat=="Subtropical","Model"], pch=19, cex=1.5, col="#6FD012") # green
# points(results.m[results.m$Habitat=="Mediterranean","TPC"], results.m[results.m$Habitat=="Mediterranean","Model"], pch=19, cex=1.5, col="#6FD012") # green
# points(results.m[results.m$Habitat=="Temperate","TPC"], results.m[results.m$Habitat=="Temperate","Model"], pch=19, cex=1.5, col="#785EF0") # purple
# points(seq(2*Xmin,2*Xmax,0.1), coef(mean.delta)[2]*seq(2*Xmin,2*Xmax,0.1)+coef(mean.delta)[1], type="l", lwd=3, col="black")
# 
# # Model vs latitude
# Xmin <- 0
# Xmax <- 60
# Ymin <- 0
# Ymax <- 25
# #dev.new(width=3, height=3, unit="in")
# plot(results.m[results.m$Habitat=="Tropical","Latitude"], results.m[results.m$Habitat=="Tropical","Model"], pch=19, cex=1.5, col="red",
#      xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Latitude", ylab="change in mean temperature")
# points(results.m[results.m$Habitat=="Subtropical","Latitude"], results.m[results.m$Habitat=="Subtropical","Model"], pch=19, cex=1.5, col="#6FD012") # green
# points(results.m[results.m$Habitat=="Mediterranean","Latitude"], results.m[results.m$Habitat=="Mediterranean","Model"], pch=19, cex=1.5, col="#6FD012") # green
# points(results.m[results.m$Habitat=="Temperate","Latitude"], results.m[results.m$Habitat=="Temperate","Model"], pch=19, cex=1.5, col="#785EF0") # purple
# points(seq(2*Xmin,2*Xmax,1), coef(mean.lat)[2]*seq(2*Xmin,2*Xmax,1) + coef(mean.lat)[1], type="l", lwd=3, col="black")
# 
# 
# # Change in temperature amplitude
# # Model vs TPCs
# Xmin <- 0
# Xmax <- 50
# Ymin <- 0
# Ymax <- 50
# #dev.new(width=3, height=3, unit="in")
# plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="TPC", ylab="Model")
# polygon(c(-10,2*Xmax,2*Xmax),c(-10,-10,2*Xmax), col = "lightgray", border = NA)
# abline(0, 1, col="gray", lwd=3)
# points(results.a[results.a$Habitat=="Tropical","TPC"], results.a[results.a$Habitat=="Tropical","Model"], pch=19, cex=1.5, col="#FFB000") # orange
# points(results.a[results.a$Habitat=="Subtropical","TPC"], results.a[results.a$Habitat=="Subtropical","Model"], pch=19, cex=1.5, col="#6FD012") # green
# points(results.a[results.a$Habitat=="Mediterranean","TPC"], results.a[results.a$Habitat=="Mediterranean","Model"], pch=19, cex=1.5, col="#6FD012") # green
# points(results.a[results.a$Habitat=="Temperate","TPC"], results.a[results.a$Habitat=="Temperate","Model"], pch=19, cex=1.5, col="#785EF0") # purple
# points(seq(2*Xmin,2*Xmax,0.1), coef(ampl.delta)[2]*seq(2*Xmin,2*Xmax,0.1)+coef(ampl.delta)[1], type="l", lwd=3, col="black")
# 
# # Model vs latitude
# Xmin <- 0
# Xmax <- 60
# Ymin <- 0
# Ymax <- 50
# #dev.new(width=3, height=3, unit="in")
# plot(results.a[results.a$Habitat=="Tropical","Latitude"], results.a[results.a$Habitat=="Tropical","Model"], pch=19, cex=1.5, col="red",
#      xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Latitude", ylab="increase in amplitude")
# points(results.a[results.a$Habitat=="Subtropical","Latitude"], results.a[results.a$Habitat=="Subtropical","Model"], pch=19, cex=1.5, col="#6FD012") # green
# points(results.a[results.a$Habitat=="Mediterranean","Latitude"], results.a[results.a$Habitat=="Mediterranean","Model"], pch=19, cex=1.5, col="#6FD012") # green
# points(results.a[results.a$Habitat=="Temperate","Latitude"], results.a[results.a$Habitat=="Temperate","Model"], pch=19, cex=1.5, col="#785EF0") # purple
# points(seq(2*Xmin,2*Xmax,1), coef(ampl.lat)[2]*seq(2*Xmin,2*Xmax,1) + coef(ampl.lat)[1], type="l", lwd=3, col="black", lty="longdash")
# 
# 
# # Change in temperature mean and amplitude
# # Model vs TPCs
# Xmin <- 0
# Xmax <- 30
# Ymin <- 0
# Ymax <- 60
# #dev.new(width=3, height=3, unit="in")
# plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="TPC", ylab="Model")
# polygon(c(-10,2*Xmax,2*Xmax),c(-10,-10,2*Xmax), col = "lightgray", border = NA)
# abline(0, 1, col="gray", lwd=3)
# points(results.b[results.b$Habitat=="Tropical","TPC"], results.b[results.b$Habitat=="Tropical","Model"], pch=19, cex=1.5, col="#FFB000") # orange
# points(results.b[results.b$Habitat=="Subtropical","TPC"], results.b[results.b$Habitat=="Subtropical","Model"], pch=19, cex=1.5, col="#6FD012") # green
# points(results.b[results.b$Habitat=="Mediterranean","TPC"], results.b[results.b$Habitat=="Mediterranean","Model"], pch=19, cex=1.5, col="#6FD012") # green
# points(results.b[results.b$Habitat=="Temperate","TPC"], results.b[results.b$Habitat=="Temperate","Model"], pch=19, cex=1.5, col="#785EF0") # purple
# points(seq(2*Xmin,2*Xmax,0.1), coef(both.delta)[2]*seq(2*Xmin,2*Xmax,0.1)+coef(both.delta)[1], type="l", lwd=3, col="black", lty="longdash")
# 
# # Model vs latitude
# Xmin <- 0
# Xmax <- 60
# Ymin <- 0
# Ymax <- 60
# #dev.new(width=3, height=3, unit="in")
# plot(results.b[results.b$Habitat=="Tropical","Latitude"], results.b[results.b$Habitat=="Tropical","Model"], pch=19, cex=1.5, col="red",
#      xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Latitude", ylab="temperature increase")
# points(results.b[results.b$Habitat=="Subtropical","Latitude"], results.b[results.b$Habitat=="Subtropical","Model"], pch=19, cex=1.5, col="#6FD012") # green
# points(results.b[results.b$Habitat=="Mediterranean","Latitude"], results.b[results.b$Habitat=="Mediterranean","Model"], pch=19, cex=1.5, col="#6FD012") # green
# points(results.b[results.b$Habitat=="Temperate","Latitude"], results.b[results.b$Habitat=="Temperate","Model"], pch=19, cex=1.5, col="#785EF0") # purple
# points(seq(2*Xmin,2*Xmax,1), coef(both.lat)[2]*seq(2*Xmin,2*Xmax,1) + coef(both.lat)[1], type="l", lwd=3, col="black", lty="longdash")





################################# TIME SERIES STATISTICS ##################################
# # CALCULATE CHANGES BETWEEN HISTORICAL AND CLIMATE CHANGE PERIODS
# # mean density
# mean.J <- mean(data.model[["J"]])
# mean.J.CC <- mean(data.model.CC[["J"]])
# mean.A <- mean(data.model[["A"]])
# mean.A.CC <- mean(data.model.CC[["A"]])
# d.mean.J <- (mean.J.CC-mean.J)/mean.J
# d.mean.A <- (mean.A.CC-mean.A)/mean.A
# 
# # peak density
# max.J <- max(data.model[["J"]])
# max.J.CC <- max(data.model.CC[["J"]])
# max.A <- max(data.model[["A"]])
# max.A.CC <- max(data.model.CC[["A"]])
# d.max.J <- (max.J.CC-max.J)/max.J
# d.max.A <- (max.A.CC-max.A)/max.A
# 
# # peak phenology
# time.J <- subset(data.model, J==max.J)$Time%%yr
# time.J.CC <- subset(data.model.CC, J==max.J.CC)$Time%%yr
# time.A <- subset(data.model, A==max.A)$Time%%yr
# time.A.CC <- subset(data.model.CC, A==max.A.CC)$Time%%yr
# d.time.J <- (time.J.CC-time.J)
# d.time.A <- (time.A.CC-time.A)
# 
# # minimum density
# min.J <- 0
# min.J.CC <- 0
# min.A <- 0
# min.A.CC <- 0
# if(min(data.model[["J"]]) > 0) { min.J <- round(min(data.model[["J"]]), digits=2) }
# if(min(data.model.CC[["J"]]) > 0) { min.J.CC <- round(min(data.model.CC[["J"]]), digits=2) }
# if(min(data.model[["A"]]) > 0) { min.A <- round(min(data.model[["A"]]), digits=2) }
# if(min(data.model.CC[["A"]]) > 0) { min.A.CC <- round(min(data.model.CC[["A"]]), digits=2) }
# if(min.J != 0) { d.min.J <- (min.J.CC-min.J)/min.J } else {d.min.J <- 0}
# if(min.A != 0) { d.min.A <- (min.A.CC-min.A)/min.A } else {d.min.A <- 0}
# 
# 
# # temperature functions
# temp <- function(t) (temp.data$meanT + temp.data$delta_mean*(t+time.shift))  - (temp.data$amplT + temp.data$delta_ampl*(t+time.shift)) * cos(2*pi*((t+time.shift) + temp.data$shiftT)/yr)
# temp.CC <- function(t) (temp.data$meanT + temp.data$delta_mean*(t+time.shift.CC))  - (temp.data$amplT + temp.data$delta_ampl*(t+time.shift.CC)) * cos(2*pi*((t+time.shift.CC) + temp.data$shiftT)/yr)
# length <- 360 # length of time over which to compare models
# 
# # calculate activity period (T(t) > Tmin)
# for(i in 0:length) {
#   if(temp(i) > sp.data["Tmin"]) {print(i)
#     i <- length} }
# 
# # reproductive activity period (birth rate > 0.1 bTopt)
# b.period <- 0
# b.period.CC <- 0
# for(i in 0:length) {
#   if(sp.data["bTopt"]*exp(-((temp(i)-sp.data["Toptb"])^2)/(2*sp.data["sb"]^2)) > 0.1*sp.data["bTopt"]) {b.period <- b.period + 1} }
# for(i in 0:length) {
#   if(sp.data["bTopt"]*exp(-((temp.CC(i)-sp.data["Toptb"])^2)/(2*sp.data["sb"]^2)) > 0.1*sp.data["bTopt"]) {b.period.CC <- b.period.CC + 1} }
# d.b <- (b.period.CC-b.period)/b.period
# 
# # developmental activity period (development rate > 0.25 Mmax or > dMin from Jahansson et al. 2020)
# # calculate Mmax and Topt
# Mmax <- 0
# Topt <- 0
# for(i in 0:100) {
#   T <- coef(r.T)[1] + i*(sp.data["TH"]-sp.data["TL"])/100
#   M.T <- sp.data["mTR"]*(T/sp.data["TR"])*exp(sp.data["AmJ"]*(1/sp.data["TR"]-1/T))/(1+sp.data["skew"]*exp(sp.data["AL"]*(1/sp.data["TL"]-1/T))+exp(sp.data["AH"]*(1/sp.data["TH"]-1/T)))
#   if(M.T > Mmax) {Mmax <- M.T[[1]]
#   Topt <- T[[1]]}}
# # calculate time above 0.25 Mmax
# m.period <- 0
# m.period.CC <- 0
# for(i in 0:length) {
#   if(sp.data["mTR"]*(temp(i)/sp.data["TR"])*exp(sp.data["AmJ"]*(1/sp.data["TR"]-1/temp(i)))/(1+sp.data["skew"]*exp(sp.data["AL"]*(1/sp.data["TL"]-1/temp(i)))+exp(sp.data["AH"]*(1/sp.data["TH"]-1/temp(i)))) > 0.25*Mmax) {m.period = m.period + 1} }
# for(i in 0:length) {
#   if(sp.data["mTR"]*(temp.CC(i)/sp.data["TR"])*exp(sp.data["AmJ"]*(1/sp.data["TR"]-1/temp.CC(i)))/(1+sp.data["skew"]*exp(sp.data["AL"]*(1/sp.data["TL"]-1/temp.CC(i)))+exp(sp.data["AH"]*(1/sp.data["TH"]-1/temp.CC(i)))) > 0.25*Mmax) {m.period.CC = m.period.CC + 1} }
# d.m <- (m.period.CC-m.period)/m.period
# 
# # average reproductive rate (mean b(T))
# b.sum <- 0
# b.sum.CC <- 0
# for(i in 0:length) {
#   b.sum <- b.sum + (sp.data["bTopt"]*exp(-((temp(i)-sp.data["Toptb"])^2)/(2*sp.data["sb"]^2)))[[1]] }
# for(i in 0:length) {
#   b.sum.CC <- b.sum.CC + (sp.data["bTopt"]*exp(-((temp.CC(i)-sp.data["Toptb"])^2)/(2*sp.data["sb"]^2)))[[1]] }
# d.b.ave <- (b.sum.CC-b.sum)/b.sum
# 
# # average development rate (mean m(T))
# m.sum <- 0
# m.sum.CC <- 0
# for(i in 0:length) {
#   m.sum <- m.sum + (sp.data["mTR"]*(temp(i)/sp.data["TR"])*exp(sp.data["AmJ"]*(1/sp.data["TR"]-1/temp(i)))/(1+sp.data["skew"]*exp(sp.data["AL"]*(1/sp.data["TL"]-1/temp(i)))+exp(sp.data["AH"]*(1/sp.data["TH"]-1/temp(i)))))[[1]] }
# for(i in 0:length) {
#   m.sum.CC <- m.sum.CC + (sp.data["mTR"]*(temp.CC(i)/sp.data["TR"])*exp(sp.data["AmJ"]*(1/sp.data["TR"]-1/temp.CC(i)))/(1+sp.data["skew"]*exp(sp.data["AL"]*(1/sp.data["TL"]-1/temp.CC(i)))+exp(sp.data["AH"]*(1/sp.data["TH"]-1/temp.CC(i)))))[[1]] }
# d.m.ave <- (m.sum.CC-m.sum)/m.sum
# 
# # average adult mortality rate (mean dA(T))
# dA.sum <- 0
# dA.sum.CC <- 0
# for(i in 0:length) {
#   dA.sum <- dA.sum + (sp.data["dATR"]*exp(sp.data["AdA"]*(1/sp.data["TR"]-1/temp(i))))[[1]] }
# for(i in 0:length) {
#   dA.sum.CC <- dA.sum.CC + (sp.data["dATR"]*exp(sp.data["AdA"]*(1/sp.data["TR"]-1/temp.CC(i))))[[1]] }
# d.dA.ave <- (dA.sum.CC-dA.sum)/dA.sum
# 
# # time above optimum range for reproductive rate (T > ToptR0 + sR0)
# R0.period <- 0
# R0.period.CC <- 0
# for(i in 0:length) {
#   if(temp(i) > sp.data["ToptR0"] + sp.data["sR0"]) {R0.period <- R0.period + 1} }
# for(i in 0:length) {
#   if(temp.CC(i) > sp.data["ToptR0"] + sp.data["sR0"]) {R0.period.CC <- R0.period.CC + 1} }
# if(R0.period !=0) { d.R0 <- (R0.period.CC-R0.period)/R0.period } else {d.R0 <- 0 }
# 
# # time above rTmax (T > rTmax)
# r.period <- 0
# r.period.CC <- 0
# for(i in 0:length) {
#   if(temp(i) > sp.data["rTmax"]) {r.period <- r.period + 1} }
# for(i in 0:length) {
#   if(temp.CC(i) > sp.data["rTmax"]) {r.period.CC <- r.period.CC + 1} }
# if(r.period !=0) { d.r <- (r.period.CC-r.period)/r.period } else {d.r <- 0 }
# 
# # mean thermal safety margin (Toptr - T)
# TSM <- 0
# TSM.CC <- 0
# for(i in 0:length) {
#   TSM <- TSM + (sp.data["Toptr"] - temp(i)) }
# for(i in 0:length) {
#   TSM.CC <- TSM.CC + (sp.data["Toptr"] - temp.CC(i)) }
# TSM <- (TSM/length)[[1]]
# TSM.CC <- (TSM.CC/length)[[1]]
# if(TSM !=0) { d.TSM <- (TSM.CC-TSM)/TSM } else {d.TSM <- 0 }
# 
# 
# # PLOTS
# par(mfrow=c(3,5))
# # density metrics
# barplot(c(d.mean.J,d.mean.A), col=c("#d1495b","#30638e"), ylim=c(-0.6,0.6), main=expression("Mean density"))
# barplot(c(d.max.J,d.max.A), col=c("#d1495b","#30638e"), ylim=c(-0.6,0.6), main=expression("Peak density"))
# barplot(c(d.time.J,d.time.A), col=c("#d1495b","#30638e"), ylim=c(-380,20), main=expression("Timing of peak"))
# barplot(c(d.min.J,d.min.A), col=c("#d1495b","#30638e"), ylim=c(-1,0.2), main=expression("Minimum density"))
# # activity periods
# barplot(d.b, col="#30638e", xlim=c(0.2,2), ylim=c(-0.5,0.5), main=expression("Rep. period"))
# barplot(d.m, col="#d1495b", xlim=c(0.2,2), ylim=c(-0.5,0.5), main=expression("Dev. period"))
# # average life history traits
# barplot(d.b.ave, col="#30638e", xlim=c(0.2,2), ylim=c(-0.5,0.5), main=expression("Mean b(T)"))
# barplot(d.m.ave, col="#d1495b", xlim=c(0.2,2), ylim=c(-0.5,0.5), main=expression("Mean m(T)"))
# barplot(d.dA.ave, col="#30638e", xlim=c(0.2,2), ylim=c(-0.5,0.5), main=expression("Mean dA(T)"))
# # thermal performance curves
# barplot(c(R0.period,R0.period.CC), col=c("purple","#FFB000") # orange, ylim=c(0,100), main=expression("above R0"))
# barplot(c(r.period,r.period.CC), col=c("purple","#FFB000") # orange, ylim=c(0,100), main=expression("above r"))
# barplot(c(TSM,TSM.CC), col=c("purple","#FFB000") # orange, ylim=c(-20,20), main=expression("Thermal margin"))
# # habitat temperatures
# barplot((sp.data["ext_meanT"])[[1]], col="#30638e", xlim=c(0.2,2), ylim=c(0,30), main=expression("meanT"))
# barplot((sp.data["ext_amplT"])[[1]], col="#30638e", xlim=c(0.2,2), ylim=c(0,30), main=expression("amplT"))
# barplot((sp.data["ext_mT_aT"])[[1]], col="#30638e", xlim=c(0.2,2), ylim=c(0,30), main=expression("meanT + amplT"))
# par(mfrow=c(1,1))
