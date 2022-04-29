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
f.data <- as.data.frame(read_csv("Predictions/Predictions Dev lifetime fecundity.csv"))
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
# Clavigralla tomentosicollis Nigeria (used data from Burkina Faso to parameterize birth rate and adult mortality)
# Macrolophus pygmaeus (only predator and separate thermal responses for each prey)
# Bemisia argentifollii and Diaphorina citri (different suborder)
r.data <- r.data[-c(3,12,13,15,18),]
R0.data <- R0.data[-c(3,12,13,15,18),]
b.data <- b.data[-c(3,12,13,15,18),]
tau.data$active.h <- pop.data$active.h # add activity period data to developmental data
tau.data$active.f <- pop.data$active.f # add activity period data to developmental data
tau.data <- tau.data[-c(3,12,13,15,18),]
s.data <- s.data[-c(3,12,13,15,18),]
L.data <- L.data[-c(3,12,13,15,18),]
f.data <- f.data[-c(3,12,13,15,18),]
R.data <- R.data[-c(3,12,13,15,18),]
LH.data <- LH.data[-c(3,12,13,15,18),]
pop.data <- pop.data[-c(3,12,13,15,18),]
# results.m <- results.m[-c(3,12,13,15,18),]
# results.a <- results.a[-c(3,12,13,15,18),]
# results.b <- results.b[-c(3,12,13,15,18),]


# SCALE LOWEST FITNESS CHANGE TO -1
r.data$delta.TPC <- pmax(-1, r.data$delta.TPC)
r.data$delta.model <- pmax(-1, r.data$delta.model)



######################################## STATISTICS #########################################
# FITNESS
# Model vs TPC
r.delta <- lm(delta.model ~ delta.TPC, data=r.data)
summary(r.delta) # significant!
# Model vs Latitude
# linear
r.lat <- lm(delta.model ~ Latitude, data=r.data)
summary(r.lat) # significant!
# non-linear
#r.lat2 <- nls(delta.model ~ a + b*Latitude + c*Latitude^2, data=r.data, start=list(a=1, b=-0.1, c=1))
#summary(r.lat2) # non-significant
# exact binomial test
r.data$sign <- sign(r.data$delta.TPC) == sign(r.data$delta.model)
r.data$under <- abs(r.data$delta.TPC) < 0.99*abs(r.data$delta.model)
r.data$over <- abs(r.data$delta.TPC) > 1.01*abs(r.data$delta.model)
binom.test(12, 21, p=0.5, alternative = "two.sided") # underestimated

# R0
# Model vs TPC
R0.delta <- lm(delta.model ~ delta.TPC, data=R0.data)
summary(R0.delta) # significant!
# Model vs Latitude
R0.lat <- lm(delta.model ~ Latitude, data=R0.data)
summary(R0.lat) # marginally-significant
# exact binomial test
R0.data$sign <- sign(R0.data$delta.TPC) == sign(R0.data$delta.model)
R0.data$under <- abs(R0.data$delta.TPC) < 0.99*abs(R0.data$delta.model)
R0.data$over <- abs(R0.data$delta.TPC) > 1.01*abs(R0.data$delta.model)
binom.test(15, 20, p=0.5, alternative = "two.sided") # underestimated

# BIRTH RATE
# Model vs Latitude
b.lat <- lm(delta.model ~ Latitude, data=b.data)
summary(b.lat) # non-significant
# sign
b.data$sign <- b.data$delta.TPC < 0

# DEVELOPMENT TIME
# Model vs TPC
tau.delta <- lm(delta.model ~ delta.TPC, data=tau.data)
summary(tau.delta)  # significant!
# Model vs Latitude
tau.lat <- lm(delta.model ~ Latitude, data=tau.data)
summary(tau.lat) # significant!
# exact binomial test
tau.data$sign <- sign(tau.data$delta.TPC) == sign(tau.data$delta.model)
tau.data$under <- abs(tau.data$delta.TPC) < 0.99*abs(tau.data$delta.model)
tau.data$over <- abs(tau.data$delta.TPC) > 1.01*abs(tau.data$delta.model)
binom.test(18, 20, p=0.5, alternative = "two.sided") # underestimated

# SURVIVAL
# Model vs TPC
s.delta <- lm(delta.model ~ delta.TPC, data=s.data)
summary(s.delta)  # significant!
# Model vs Latitude
s.lat <- lm(delta.model ~ Latitude, data=s.data)
summary(s.lat) # non-significant
# exact binomial test
s.data$sign <- sign(s.data$delta.TPC) == sign(s.data$delta.model)
s.data$under <- abs(s.data$delta.TPC) < 0.99*abs(s.data$delta.model)
s.data$over <- abs(s.data$delta.TPC) > 1.01*abs(s.data$delta.model)
binom.test(4, 12, p=0.5, alternative = "two.sided") # underestimated

# ADULT LONGEVITY
# Model vs Latitude
L.lat <- lm(delta.model ~ Latitude, data=L.data)
summary(L.lat) # significant!
# sign
L.data$sign <- L.data$delta.TPC < 0

# LIFETIME FECUNDITY
# Model vs Latitude
#f.lat <- lm(delta.model ~ Latitude, data=f.data)
#summary(f.lat) # significant!

# ADULT RECRUITMENT
# Model vs TPC
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
CV.lat <- lm(delta.CV ~ Latitude, data=pop.data[-c(15,20),]) # excluding Macrosiphum euphorbiae Canada and Brevicoryne brassicae (went extinct)
summary(CV.lat) # significant!
# Active period vs Latitude (NOTE: non-significant for temperate species only)
active.lat <- lm(delta.active ~ Latitude, data=pop.data) #[pop.data$Habitat == "Temperate",])
summary(active.lat) # non-significant



########################################### PLOTS ###########################################
# RELATIVE FITNESS
# Model vs TPCs
Xmin <- -1
Xmax <- 0.5
Ymin <- -1
Ymax <- 0.5
#dev.new(width=3, height=3, unit="in")
plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="TPC", ylab="Model", cex.axis=2)
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

# Fitness vs latitude
Xmin <- 0
Xmax <- 60
Ymin <- -1
Ymax <- 0.5
#dev.new(width=3, height=3, unit="in")
plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Latitude", ylab="Model", cex.axis=2)
abline(0, 0, col="gray", lwd=3, lty="longdash")
points(r.data[r.data$Habitat=="Tropical","Latitude"], r.data[r.data$Habitat=="Tropical","delta.model"], pch=19, cex=1.5, col="#FFB000") # orange
points(r.data[r.data$Habitat=="Subtropical","Latitude"], r.data[r.data$Habitat=="Subtropical","delta.model"], pch=19, cex=1.5, col="#6FD012") # green
points(r.data[r.data$Habitat=="Mediterranean","Latitude"], r.data[r.data$Habitat=="Mediterranean","delta.model"], pch=19, cex=1.5, col="#6FD012") # green
points(r.data[r.data$Habitat=="Temperate","Latitude"], r.data[r.data$Habitat=="Temperate","delta.model"], pch=19, cex=1.5, col="#785EF0") # purple
points(seq(Xmin,Xmax,1), coef(r.lat)[2]*seq(Xmin,Xmax,1) + coef(r.lat)[1], type="l", lwd=3, col="black")


# R0
# Model vs TPCs
Xmin <- -0.8
Xmax <- 0.2
Ymin <- -0.8
Ymax <- 0.2
#dev.new(width=3, height=3, unit="in")
plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="TPC", ylab="Model", cex.axis=2)
polygon(c(2*Xmin,0,2*Xmax,2*Xmax),c(2*Xmin,0,-2*Xmax,2*Ymin), col = "#E2E2E2", border = NA)
polygon(c(2*Xmin,0,2*Xmax),c(-2*Xmin,0,2*Xmax), col = "#E2E2E2", border = NA)
abline(0, 1, col="gray", lwd=3)
abline(0, -1, col="gray", lwd=1)
abline(0, 0, col="gray", lwd=3, lty="longdash")
abline(v = 0, col="gray", lwd=3, lty="longdash")
points(R0.data[R0.data$Habitat=="Tropical","delta.TPC"], R0.data[R0.data$Habitat=="Tropical","delta.model"], pch=19, cex=1.5, col="#FFB000") # orange
points(R0.data[R0.data$Habitat=="Subtropical","delta.TPC"], R0.data[R0.data$Habitat=="Subtropical","delta.model"], pch=19, cex=1.5, col="#6FD012") # green
points(R0.data[R0.data$Habitat=="Mediterranean","delta.TPC"], R0.data[R0.data$Habitat=="Mediterranean","delta.model"], pch=19, cex=1.5, col="#6FD012") # green
points(R0.data[R0.data$Habitat=="Temperate","delta.TPC"], R0.data[R0.data$Habitat=="Temperate","delta.model"], pch=19, cex=1.5, col="#785EF0") # purple
points(seq(2*Xmin,2*Xmax,0.1), coef(R0.delta)[2]*seq(2*Xmin,2*Xmax,0.1)+coef(R0.delta)[1], type="l", lwd=3, col="black")

# R0 vs latitude
Xmin <- 0
Xmax <- 60
Ymin <- -0.8
Ymax <- 0.2
#dev.new(width=3, height=3, unit="in")
plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Latitude", ylab="Model", cex.axis=2)
abline(0, 0, col="gray", lwd=3, lty="longdash")
points(R0.data[R0.data$Habitat=="Tropical","Latitude"], R0.data[R0.data$Habitat=="Tropical","delta.model"], pch=19, cex=1.5, col="#FFB000") # orange
points(R0.data[R0.data$Habitat=="Subtropical","Latitude"], R0.data[R0.data$Habitat=="Subtropical","delta.model"], pch=19, cex=1.5, col="#6FD012") # green
points(R0.data[R0.data$Habitat=="Mediterranean","Latitude"], R0.data[R0.data$Habitat=="Mediterranean","delta.model"], pch=19, cex=1.5, col="#6FD012") # green
points(R0.data[R0.data$Habitat=="Temperate","Latitude"], R0.data[R0.data$Habitat=="Temperate","delta.model"], pch=19, cex=1.5, col="#785EF0") # purple
points(seq(Xmin,Xmax,1), coef(R0.lat)[2]*seq(Xmin,Xmax,1) + coef(R0.lat)[1], type="l", lwd=3, col="black", lty="longdash")


# BIRTH RATE
# Model vs TPCs
Xmin <- -0.6
Xmax <- 0.2
Ymin <- -0.6
Ymax <- 0.2
#dev.new(width=3, height=3, unit="in")
plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="TPC", ylab="Model", cex.axis=2)
polygon(c(2*Xmin,0,2*Xmax,2*Xmax),c(2*Xmin,0,-2*Xmax,2*Ymin), col = "#E2E2E2", border = NA)
polygon(c(2*Xmin,0,2*Xmax),c(-2*Xmin,0,2*Xmax), col = "#E2E2E2", border = NA)
abline(0, 1, col="gray", lwd=3)
abline(0, -1, col="gray", lwd=1)
abline(0, 0, col="gray", lwd=3, lty="longdash")
abline(v = 0, col="gray", lwd=3, lty="longdash")
points(b.data[b.data$Habitat=="Tropical","delta.TPC"], b.data[b.data$Habitat=="Tropical","delta.model"], pch=19, cex=1.5, col="#FFB000") # orange
points(b.data[b.data$Habitat=="Subtropical","delta.TPC"], b.data[b.data$Habitat=="Subtropical","delta.model"], pch=19, cex=1.5, col="#6FD012") # green
points(b.data[b.data$Habitat=="Mediterranean","delta.TPC"], b.data[b.data$Habitat=="Mediterranean","delta.model"], pch=19, cex=1.5, col="#6FD012") # green
points(b.data[b.data$Habitat=="Temperate","delta.TPC"], b.data[b.data$Habitat=="Temperate","delta.model"], pch=19, cex=1.5, col="#785EF0") # purple

# Birth rate vs latitude
Xmin <- 0
Xmax <- 60
Ymin <- -0.6
Ymax <- 0.2
#dev.new(width=3, height=3, unit="in")
plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Latitude", ylab="Model", cex.axis=2)
abline(0, 0, col="gray", lwd=3, lty="longdash")
points(b.data[b.data$Habitat=="Tropical","Latitude"], b.data[b.data$Habitat=="Tropical","delta.model"], pch=19, cex=1.5, col="#FFB000") # orange
points(b.data[b.data$Habitat=="Subtropical","Latitude"], b.data[b.data$Habitat=="Subtropical","delta.model"], pch=19, cex=1.5, col="#6FD012") # green
points(b.data[b.data$Habitat=="Mediterranean","Latitude"], b.data[b.data$Habitat=="Mediterranean","delta.model"], pch=19, cex=1.5, col="#6FD012") # green
points(b.data[b.data$Habitat=="Temperate","Latitude"], b.data[b.data$Habitat=="Temperate","delta.model"], pch=19, cex=1.5, col="#785EF0") # purple
points(seq(Xmin,Xmax,1), coef(b.lat)[2]*seq(Xmin,Xmax,1) + coef(b.lat)[1], type="l", lwd=3, col="black", lty="longdash")


# DEVELOPMENT TIME
# Model vs TPCs
Xmin <- -2.5
Xmax <- 0
Ymin <- -8
Ymax <- 0
#dev.new(width=3, height=3, unit="in")
plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="TPC", ylab="Model", cex.axis=2)
polygon(c(2*Xmin,0,2*Xmax,2*Xmax),c(2*Xmin,0,-2*Xmax,2*Ymin), col = "#E2E2E2", border = NA)
polygon(c(2*Xmin,0,2*Xmax),c(-2*Xmin,0,2*Xmax), col = "#E2E2E2", border = NA)
abline(0, 1, col="gray", lwd=3)
abline(0, -1, col="gray", lwd=1)
abline(0, 0, col="gray", lwd=3, lty="longdash")
abline(v = 0, col="gray", lwd=3, lty="longdash")
points(tau.data[tau.data$Habitat=="Tropical","delta.TPC"], tau.data[tau.data$Habitat=="Tropical","delta.model"], pch=19, cex=1.5, col="#FFB000") # orange
points(tau.data[tau.data$Habitat=="Subtropical","delta.TPC"], tau.data[tau.data$Habitat=="Subtropical","delta.model"], pch=19, cex=1.5, col="#6FD012") # green
points(tau.data[tau.data$Habitat=="Mediterranean","delta.TPC"], tau.data[tau.data$Habitat=="Mediterranean","delta.model"], pch=19, cex=1.5, col="#6FD012") # green
points(tau.data[tau.data$Habitat=="Temperate","delta.TPC"], tau.data[tau.data$Habitat=="Temperate","delta.model"], pch=19, cex=1.5, col="#785EF0") # purple
points(seq(2*Xmin,2*Xmax,0.1), coef(tau.delta)[2]*seq(2*Xmin,2*Xmax,0.1)+coef(tau.delta)[1], type="l", lwd=3, col="black")

# Development time vs latitude
Xmin <- 0
Xmax <- 60
Ymin <- -8
Ymax <- 0
#dev.new(width=3, height=3, unit="in")
plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Latitude", ylab="Model", cex.axis=2)
abline(0, 0, col="gray", lwd=3, lty="longdash")
points(tau.data[tau.data$Habitat=="Tropical","Latitude"], tau.data[tau.data$Habitat=="Tropical","delta.model"], pch=19, cex=1.5, col="#FFB000") # orange
points(tau.data[tau.data$Habitat=="Subtropical","Latitude"], tau.data[tau.data$Habitat=="Subtropical","delta.model"], pch=19, cex=1.5, col="#6FD012") # green
points(tau.data[tau.data$Habitat=="Mediterranean","Latitude"], tau.data[tau.data$Habitat=="Mediterranean","delta.model"], pch=19, cex=1.5, col="#6FD012") # green
points(tau.data[tau.data$Habitat=="Temperate","Latitude"], tau.data[tau.data$Habitat=="Temperate","delta.model"], pch=19, cex=1.5, col="#785EF0") # purple
points(seq(Xmin,Xmax,1), coef(tau.lat)[2]*seq(Xmin,Xmax,1) + coef(tau.lat)[1], type="l", lwd=3, col="black")


# SURVIVAL
# Model vs TPCs
Xmin <- -0.4
Xmax <- 0.1
Ymin <- -1
Ymax <- 0.2
#dev.new(width=3, height=3, unit="in")
plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="TPC", ylab="Model", cex.axis=2)
polygon(c(2*Xmin,0,2*Xmax,2*Xmax),c(2*Xmin,0,-2*Xmax,2*Ymin), col = "#E2E2E2", border = NA)
polygon(c(2*Xmin,0,2*Xmax),c(-2*Xmin,0,2*Xmax), col = "#E2E2E2", border = NA)
abline(0, 1, col="gray", lwd=3)
abline(0, -1, col="gray", lwd=1)
abline(0, 0, col="gray", lwd=3, lty="longdash")
abline(v = 0, col="gray", lwd=3, lty="longdash")
points(s.data[s.data$Habitat=="Tropical","delta.TPC"], s.data[s.data$Habitat=="Tropical","delta.model"], pch=19, cex=1.5, col="#FFB000") # orange
points(s.data[s.data$Habitat=="Subtropical","delta.TPC"], s.data[s.data$Habitat=="Subtropical","delta.model"], pch=19, cex=1.5, col="#6FD012") # green
points(s.data[s.data$Habitat=="Mediterranean","delta.TPC"], s.data[s.data$Habitat=="Mediterranean","delta.model"], pch=19, cex=1.5, col="#6FD012") # green
points(s.data[s.data$Habitat=="Temperate","delta.TPC"], s.data[s.data$Habitat=="Temperate","delta.model"], pch=19, cex=1.5, col="#785EF0") # purple
points(seq(2*Xmin,2*Xmax,0.1), coef(s.delta)[2]*seq(2*Xmin,2*Xmax,0.1)+coef(s.delta)[1], type="l", lwd=3, col="black")

# Survival vs latitude
Xmin <- 0
Xmax <- 60
Ymin <- -1
Ymax <- 0.2
#dev.new(width=3, height=3, unit="in")
plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Latitude", ylab="Model", cex.axis=2)
abline(0, 0, col="gray", lwd=3, lty="longdash")
points(s.data[s.data$Habitat=="Tropical","Latitude"], s.data[s.data$Habitat=="Tropical","delta.model"], pch=19, cex=1.5, col="#FFB000") # orange
points(s.data[s.data$Habitat=="Subtropical","Latitude"], s.data[s.data$Habitat=="Subtropical","delta.model"], pch=19, cex=1.5, col="#6FD012") # green
points(s.data[s.data$Habitat=="Mediterranean","Latitude"], s.data[s.data$Habitat=="Mediterranean","delta.model"], pch=19, cex=1.5, col="#6FD012") # green
points(s.data[s.data$Habitat=="Temperate","Latitude"], s.data[s.data$Habitat=="Temperate","delta.model"], pch=19, cex=1.5, col="#785EF0") # purple
points(seq(Xmin,Xmax,1), coef(s.lat)[2]*seq(Xmin,Xmax,1) + coef(s.lat)[1], type="l", lwd=3, col="black", lty="longdash")


# ADULT LONGEVITY
# Model vs TPCs
Xmin <- -0.4
Xmax <- 0.1
Ymin <- -0.4
Ymax <- 0.1
#dev.new(width=3, height=3, unit="in")
plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="TPC", ylab="Model", cex.axis=2)
polygon(c(2*Xmin,0,2*Xmax,2*Xmax),c(2*Xmin,0,-2*Xmax,2*Ymin), col = "#E2E2E2", border = NA)
polygon(c(2*Xmin,0,2*Xmax),c(-2*Xmin,0,2*Xmax), col = "#E2E2E2", border = NA)
abline(0, 1, col="gray", lwd=3)
abline(0, -1, col="gray", lwd=1)
abline(0, 0, col="gray", lwd=3, lty="longdash")
abline(v = 0, col="gray", lwd=3, lty="longdash")
points(L.data[L.data$Habitat=="Tropical","delta.TPC"], L.data[L.data$Habitat=="Tropical","delta.model"], pch=19, cex=1.5, col="#FFB000") # orange
points(L.data[L.data$Habitat=="Subtropical","delta.TPC"], L.data[L.data$Habitat=="Subtropical","delta.model"], pch=19, cex=1.5, col="#6FD012") # green
points(L.data[L.data$Habitat=="Mediterranean","delta.TPC"], L.data[L.data$Habitat=="Mediterranean","delta.model"], pch=19, cex=1.5, col="#6FD012") # green
points(L.data[L.data$Habitat=="Temperate","delta.TPC"], L.data[L.data$Habitat=="Temperate","delta.model"], pch=19, cex=1.5, col="#785EF0") # purple

# Adult longevity vs latitude
Xmin <- 0
Xmax <- 60
Ymin <- -0.4
Ymax <- 0.1
#dev.new(width=3, height=3, unit="in")
plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Latitude", ylab="Model", cex.axis=2)
abline(0, 0, col="gray", lwd=3, lty="longdash")
points(L.data[L.data$Habitat=="Tropical","Latitude"], L.data[L.data$Habitat=="Tropical","delta.model"], pch=19, cex=1.5, col="#FFB000") # orange
points(L.data[L.data$Habitat=="Subtropical","Latitude"], L.data[L.data$Habitat=="Subtropical","delta.model"], pch=19, cex=1.5, col="#6FD012") # green
points(L.data[L.data$Habitat=="Mediterranean","Latitude"], L.data[L.data$Habitat=="Mediterranean","delta.model"], pch=19, cex=1.5, col="#6FD012") # green
points(L.data[L.data$Habitat=="Temperate","Latitude"], L.data[L.data$Habitat=="Temperate","delta.model"], pch=19, cex=1.5, col="#785EF0") # purple
points(seq(Xmin,Xmax,1), coef(L.lat)[2]*seq(Xmin,Xmax,1) + coef(L.lat)[1], type="l", lwd=3, col="black")



# POPULATION DYNAMICS
# Mean density vs latitude
Xmin <- 0
Xmax <- 60
Ymin <- -1
Ymax <- 1
#dev.new(width=3, height=3, unit="in")
plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Latitude", ylab="Model", cex.axis=2)
abline(0, 0, col="gray", lwd=3, lty="longdash")
points(seq(2*Xmin,2*Xmax,1), coef(mean.lat)[2]*seq(2*Xmin,2*Xmax,1) + coef(mean.lat)[1], type="l", lwd=3, col="black", lty="longdash")
points(pop.data[pop.data$Habitat=="Tropical","Latitude"], pop.data[pop.data$Habitat=="Tropical","delta.mean"], pch=19, cex=1.5, col="#FFB000") # orange
points(pop.data[pop.data$Habitat=="Mediterranean","Latitude"], pop.data[pop.data$Habitat=="Mediterranean","delta.mean"], pch=19, cex=1.5, col="#6FD012") # green
points(pop.data[pop.data$Habitat=="Temperate","Latitude"], pop.data[pop.data$Habitat=="Temperate","delta.mean"], pch=19, cex=1.5, col="#785EF0") # purple
points(pop.data[pop.data$Habitat=="Subtropical","Latitude"], pop.data[pop.data$Habitat=="Subtropical","delta.mean"], pch=19, cex=1.5, col="#6FD012") # green

# CV of density vs latitude
Xmin <- 0
Xmax <- 60
Ymin <- -1
Ymax <- 1
#dev.new(width=3, height=3, unit="in")
plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Latitude", ylab="Model", cex.axis=2)
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
Ymax <- 0.4
#dev.new(width=3, height=3, unit="in")
plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Latitude", ylab="Model", cex.axis=2)
abline(0, 0, col="gray", lwd=3, lty="longdash")
points(seq(32,58,1), coef(active.lat)[2]*seq(32,58,1) + coef(active.lat)[1], type="l", lwd=3, col="black", lty="longdash")
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
# polygon(c(-1,2*Xmax,2*Xmax),c(-1,-1,2*Xmax), col = "#E2E2E2", border = NA)
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
# polygon(c(-10,2*Xmax,2*Xmax),c(-10,-10,2*Xmax), col = "#E2E2E2", border = NA)
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
# polygon(c(-10,2*Xmax,2*Xmax),c(-10,-10,2*Xmax), col = "#E2E2E2", border = NA)
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




####################################### OTHER ANALYSES #######################################
# LIFE HISTORY TRAITS
# BIRTH RATE
# # bTopt vs Latitude
# bTopt.lat <- lm(bTopt ~ Latitude, data=LH.data)
# summary(bTopt.lat) # significant!
# plot(LH.data$Latitude,LH.data$bTopt, ylim=c(0,10))
# # Toptb vs Latitude
# Toptb.lat <- lm(Toptb ~ Latitude, data=LH.data)
# summary(Toptb.lat) # significant!
# plot(LH.data$Latitude,LH.data$Toptb)
# # sb vs Latitude
# sb.lat <- lm(sb ~ Latitude, data=LH.data)
# summary(sb.lat) # significant!
# plot(LH.data$Latitude,LH.data$sb)
# # Toptb - meanT vs Latitude
# LH.data$delta_b <- LH.data$Toptb - (temp.data$meanT.f + temp.data$delta_mean.f*75*365)
# delta.b.lat <- lm(delta_b ~ Latitude, data=LH.data)
# summary(delta.b.lat) # significant!
# plot(LH.data$Latitude,LH.data$delta_b)

# DEVELOPMENT RATE
# # mTR vs Latitude
# mTR.lat <- lm(mTR ~ Latitude, data=LH.data)
# summary(mTR.lat) # non-significant
# plot(LH.data$Latitude,LH.data$mTR)
# # AL vs Latitude
# AL.lat <- lm(AL ~ Latitude, data=LH.data)
# summary(AL.lat) # non-significant
# plot(LH.data$Latitude,LH.data$AL)
# # AH vs Latitude
# AH.lat <- lm(AH ~ Latitude, data=LH.data)
# summary(AH.lat) # non-significant
# plot(LH.data$Latitude,LH.data$AH)
# # TL vs Latitude
# TL.lat <- lm(TL ~ Latitude, data=LH.data)
# summary(TL.lat) # significant!
# plot(LH.data$Latitude,LH.data$TL)
# # TH vs Latitude
# TH.lat <- lm(TH ~ Latitude, data=LH.data)
# summary(TH.lat) # non-significant
# plot(LH.data$Latitude,LH.data$TH)
# # Tmin vs Latitude
# Tmin.lat <- lm(Tmin ~ Latitude, data=LH.data)
# summary(Tmin.lat) # significant!
# plot(LH.data$Latitude,LH.data$Tmin)
# # Topt vs Latitude
# Topt.lat <- lm(Topt ~ Latitude, data=LH.data)
# summary(Topt.lat) # significant!
# plot(LH.data$Latitude,LH.data$Topt)
# # Tmax vs Latitude
# Tmax.lat <- lm(Tmax ~ Latitude, data=LH.data)
# summary(Tmax.lat) # non-significant
# plot(LH.data$Latitude,LH.data$Tmax)
# # Topt - meanT vs Latitude
# LH.data$delta_Topt <- LH.data$Topt - (temp.data$meanT.f  + temp.data$delta_mean.f*75*365)
# delta.Topt.lat <- lm(delta_Topt ~ Latitude, data=LH.data)
# summary(delta.Topt.lat) # significant!
# plot(LH.data$Latitude,LH.data$delta_Topt)
# # Topt - meanT vs Latitude
# LH.data$delta_Topt <- LH.data$Topt - (temp.data$meanT.f  + temp.data$delta_mean.f*75*365 +
#                                         abs(temp.data$amplT.f)  + abs(temp.data$delta_ampl.f*75*365))
# delta.Topt.lat <- lm(delta_Topt ~ Latitude, data=LH.data)
# summary(delta.Topt.lat) # non-significant
# plot(LH.data$Latitude,LH.data$delta_Topt)

# JUVENILE MORTALITY RATE
# # dJTR vs Latitude
# dJTR.lat <- lm(dJTR ~ Latitude, data=LH.data)
# summary(dJTR.lat) # non-significant
# plot(LH.data$Latitude,LH.data$dJTR)
# # AdJ vs Latitude
# AdJ.lat <- lm(AdJ ~ Latitude, data=LH.data)
# summary(AdJ.lat) # non-significant
# plot(LH.data$Latitude,LH.data$AdJ)
# # dJ(Tmean) vs Latitude
# LH.data$dJ_Tmean <- LH.data$dJTR*exp(LH.data$AdJ*(1/LH.data$TR-1/(temp.data$meanT.f + temp.data$delta_mean.f*75*365)))
# dJ.Tmean.lat <- lm(dJ_Tmean ~ Latitude, data=LH.data)
# summary(dJ.Tmean.lat) # significant
# plot(LH.data$Latitude,LH.data$dJ_Tmean)

# ADULT MORTALITY RATE
# # dATR vs Latitude
# dATR.lat <- lm(dATR ~ Latitude, data=LH.data)
# summary(dATR.lat) # marginally-significant
# plot(LH.data$Latitude,LH.data$dATR)
# # AdA vs Latitude
# AdA.lat <- lm(AdA ~ Latitude, data=LH.data)
# summary(AdA.lat) # non-significant
# plot(LH.data$Latitude,LH.data$AdA)
# # dA(Tmean) vs Latitude
# LH.data$dA_Tmean <- LH.data$dATR*exp(LH.data$AdA*(1/LH.data$TR-1/(temp.data$meanT.f + temp.data$delta_mean.f*75*365)))
# dA.Tmean.lat <- lm(dA_Tmean ~ Latitude, data=LH.data)
# summary(dA.Tmean.lat) # marginally-significant
# plot(LH.data$Latitude,LH.data$dA_Tmean)

# NUMBER OF GENERATIONS
# # historical period
# LH.data$Gen.h <- (tau.data$active.h / tau.data$Model.h)
# Gen.h.lat <- lm(Gen.h ~ Latitude, data=LH.data)
# summary(Gen.h.lat) # non-significant
# #plot(LH.data$Latitude,LH.data$Gen.h)
# # future period
# LH.data$Gen.f <- (tau.data$active.f / tau.data$Model.f)
# Gen.f.lat <- lm(Gen.f ~ Latitude, data=LH.data)
# summary(Gen.f.lat) # non-significant
# #plot(LH.data$Latitude,LH.data$Gen.f)
# # change in number of generations
# LH.data$delta_Gen <- (LH.data$Gen.f - LH.data$Gen.h) / LH.data$Gen.h
# delta.Gen.lat <- lm(delta_Gen ~ Latitude, data=LH.data[-c(16,21),])
# summary(delta.Gen.lat) # non-significant
# plot(LH.data$Latitude[-c(16,21)],LH.data$delta_Gen[-c(16,21)])


# TEMPERATURE
# delta_meanT vs Latitude
# Tmean.lat <- lm(delta_mean.f ~ Latitude, data=temp.data)
# summary(Tmean.lat) # non-significant
# delta_amplT temperature vs Latitude
# Tampl.lat <- lm(delta_ampl.f ~ Latitude, data=temp.data)
# summary(Tampl.lat) # significant!
# Change in mean temperature vs Latitude
# temp.data$mean_ch <- (temp.data$meanT.f + temp.data$delta_mean.f*75*365)- temp.data$meanT.h
# mean.ch.lat <- lm(mean_ch ~ Latitude, data=temp.data)
# summary(mean.ch.lat) # non-significant
# Change in ampl temperature vs Latitude
# temp.data$ampl_ch <- (temp.data$amplT.f + temp.data$delta_ampl.f*75*365)- temp.data$amplT.h
# ampl.ch.lat <- lm(ampl_ch ~ Latitude, data=temp.data)
# summary(ampl.ch.lat) # marginally-significant


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