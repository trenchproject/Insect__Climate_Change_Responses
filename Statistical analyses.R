#############################################################################
#### This R script analyzes the thermal performance curve and model data ####
#############################################################################


# Load packages and set working directory
library(tidyr)
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyverse)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# READ IN DATA
r.data <- as.data.frame(read_csv("Model results DI Tave.csv"))
f.data <- as.data.frame(read_csv("Fecundity.csv"))
s.data <- as.data.frame(read_csv("Survival.csv"))
results.m <- as.data.frame(read_csv("Extinction meanT.csv"))
results.a <- as.data.frame(read_csv("Extinction amplT.csv"))


# EXCLUDE DELOACH 1974 (QUESTIONABLE DATA on r)
r.data <- r.data[-c(25:27),]
f.data <- f.data[-c(25:27),]
s.data <- s.data[-c(25:27),]
results.m <- results.m[-c(25:27),]
results.a <- results.a[-c(25:27),]


######################################## STATISTICS #########################################
# RELATIVE FITNESS
# Change in relative fitness
delta.r <- lm(delta.model ~ delta.TPC, data=r.data) # significant!
summary(delta.r)
# Relative fitness in historical period
r.h <- lm(r.model.h ~ r.TPC.h, data=r.data) # significant!
summary(r.h)
# Relative fitness in future period
r.f <- lm(r.model.f ~ r.TPC.f, data=r.data) # non-significant
summary(r.f)

# # TPC vs Latitude
# # Linear
# TPC <- lm(delta.TPC ~ Latitude, data=r.data)
# summary(TPC) # non-significant
# # Non-linear
# TPC2 <- nls(delta.TPC ~ a + b*Latitude + c*Latitude^2, data=r.data, start=list(a=1, b=-0.1, c=1))
# summary(TPC2) # non-significant

# Model vs Latitude
# linear
r.model <- lm(delta.model ~ Latitude, data=r.data)
summary(r.model) # non-significant
# non-linear
r.model2 <- nls(delta.model ~ a + b*Latitude + c*Latitude^2, data=r.data, start=list(a=1, b=-0.1, c=1))
summary(r.model2) # non-significant


# RELATIVE FECUNDITY
# Model vs Latitude
# Linear
f.model <- lm(delta.model ~ Latitude, data=f.data)
summary(f.model) # significant!
# Non-linear
f.model2 <- nls(delta.model ~ a + b*Latitude + c*Latitude^2, data=f.data, start=list(a=1, b=-0.1, c=1))
summary(f.model2) # non-significant


# SURVIVAL
# Change in survival
delta.s <- lm(delta.model ~ delta.TPC, data=s.data) # significant!
summary(delta.s)
# Survival in historical period
s.h <- lm(Model.h ~ TPC.h, data=s.data) # significant!
summary(s.h)
# Survival in future period
s.f <- lm(Model.f ~ TPC.f, data=s.data) # significant!
summary(s.f)

# Model vs Latitude
# Linear
s.model <- lm(delta.model ~ Latitude, data=s.data)
summary(s.model) # non-significant
# Non-linear
s.model2 <- nls(delta.model ~ a + b*Latitude + c*Latitude^2, data=s.data, start=list(a=1, b=-0.1, c=1))
summary(s.model2) # non-significant


# EXTINCTION (INCREASE IN MEAN TEMPERATURE)
model.m <- lm(Model ~ TPC, data=results.m)
summary(model.m) # marginally significant


# EXTINCTION (INCREASE IN TEMPERATURE AMPLITUDE)
model.a <- lm(Model ~ TPC, data=results.a)
summary(model.a) # marginally significant



########################################### PLOTS ###########################################
# RELATIVE FITNESS
# Model vs TPCs
Xmin <- -1
Xmax <- 0.5
Ymin <- -1
Ymax <- 0.5
#dev.new(width=3, height=3, unit="in")
plot(r.data[r.data$Habitat=="Tropical","delta.TPC"], r.data[r.data$Habitat=="Tropical","delta.model"], pch=19, cex=1.5, col="red",
     xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="TPC", ylab="Model")
points(r.data[r.data$Habitat=="Subtropical","delta.TPC"], r.data[r.data$Habitat=="Subtropical","delta.model"], pch=19, cex=1.5, col="orange")
points(r.data[r.data$Habitat=="Mediterranean","delta.TPC"], r.data[r.data$Habitat=="Mediterranean","delta.model"], pch=19, cex=1.5, col="orange")
points(r.data[r.data$Habitat=="Temperate","delta.TPC"], r.data[r.data$Habitat=="Temperate","delta.model"], pch=19, cex=1.5, col="blue")
points(seq(Xmin,Xmax,0.1), coef(delta.r)[2]*seq(Xmin,Xmax,0.1)+coef(delta.r)[1], type="l", col="black")
abline(0, 1, col="gray")
abline(0, 0, col="gray", lty="longdash")
abline(v = 0, col="gray", lty="longdash")

# Model vs TPCs in historical period
Xmin <- 0
Xmax <- 1
Ymin <- 0
Ymax <- 1
#dev.new(width=3, height=3, unit="in")
plot(r.data[r.data$Habitat=="Tropical","r.TPC.h"], r.data[r.data$Habitat=="Tropical","r.model.h"], pch=19, cex=1.5, col="red",
     xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="TPC", ylab="Model")
points(r.data[r.data$Habitat=="Subtropical","r.TPC.h"], r.data[r.data$Habitat=="Subtropical","r.model.h"], pch=19, cex=1.5, col="orange")
points(r.data[r.data$Habitat=="Mediterranean","r.TPC.h"], r.data[r.data$Habitat=="Mediterranean","r.model.h"], pch=19, cex=1.5, col="orange")
points(r.data[r.data$Habitat=="Temperate","r.TPC.h"], r.data[r.data$Habitat=="Temperate","r.model.h"], pch=19, cex=1.5, col="blue")
points(seq(Xmin,Xmax,0.1), coef(r.h)[2]*seq(Xmin,Xmax,0.1)+coef(r.h)[1], type="l", col="black")
abline(0, 1, col="gray")
#abline(0, 0, col="gray", lty="longdash")
#abline(v = 0, col="gray", lty="longdash")

# Model vs TPCs in future period
Xmin <- -0.5
Xmax <- 1.5
Ymin <- -0.5
Ymax <- 1.5
#dev.new(width=3, height=3, unit="in")
plot(r.data[r.data$Habitat=="Tropical","r.TPC.f"], r.data[r.data$Habitat=="Tropical","r.model.f"], pch=19, cex=1.5, col="red",
     xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="TPC", ylab="Model")
points(r.data[r.data$Habitat=="Subtropical","r.TPC.f"], r.data[r.data$Habitat=="Subtropical","r.model.f"], pch=19, cex=1.5, col="orange")
points(r.data[r.data$Habitat=="Mediterranean","r.TPC.f"], r.data[r.data$Habitat=="Mediterranean","r.model.f"], pch=19, cex=1.5, col="orange")
points(r.data[r.data$Habitat=="Temperate","r.TPC.f"], r.data[r.data$Habitat=="Temperate","r.model.f"], pch=19, cex=1.5, col="blue")
points(seq(Xmin,Xmax,0.1), coef(r.f)[2]*seq(Xmin,Xmax,0.1)+coef(r.f)[1], type="l", col="black")
abline(0, 1, col="gray")
#abline(0, 0, col="gray", lty="longdash")
#abline(v = 0, col="gray", lty="longdash")

# TPC vs latitude
# Xmin <- 0
# Xmax <- 60
# Ymin <- -1
# Ymax <- 0.5
# plot(r.data[r.data$Habitat=="Tropical","Latitude"], r.data[r.data$Habitat=="Tropical","delta.TPC"], pch=19, cex=1.5, col="red",
#      xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Latitude", ylab="change in r")
# points(r.data[r.data$Habitat=="Subtropical","Latitude"], r.data[r.data$Habitat=="Subtropical","delta.TPC"], pch=19, cex=1.5, col="orange")
# points(r.data[r.data$Habitat=="Mediterranean","Latitude"], r.data[r.data$Habitat=="Mediterranean","delta.TPC"], pch=19, cex=1.5, col="orange")
# points(r.data[r.data$Habitat=="Temperate","Latitude"], r.data[r.data$Habitat=="Temperate","delta.TPC"], pch=19, cex=1.5, col="blue")
# points(seq(Xmin,Xmax,1), coef(TPC)[2]*seq(Xmin,Xmax,1) + coef(TPC)[1], type="l", col="black", lty="longdash")
# points(seq(Xmin,Xmax,1), coef(TPC2)[1] + coef(TPC2)[2]*seq(Xmin,Xmax,1) + coef(TPC2)[3]*seq(Xmin,Xmax,1)^2, type="l", col="black")
# abline(0, 0, col="black") #, lty="longdash")

# Model vs latitude
Xmin <- 0
Xmax <- 60
Ymin <- -1
Ymax <- 0.5
#dev.new(width=3, height=3, unit="in")
plot(r.data[r.data$Habitat=="Tropical","Latitude"], r.data[r.data$Habitat=="Tropical","delta.model"], pch=19, cex=1.5, col="red",
     xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Latitude", ylab="change in r")
points(r.data[r.data$Habitat=="Subtropical","Latitude"], r.data[r.data$Habitat=="Subtropical","delta.model"], pch=19, cex=1.5, col="orange")
points(r.data[r.data$Habitat=="Mediterranean","Latitude"], r.data[r.data$Habitat=="Mediterranean","delta.model"], pch=19, cex=1.5, col="orange")
abline(0, 0, col="gray") #, lty="longdash")
points(r.data[r.data$Habitat=="Temperate","Latitude"], r.data[r.data$Habitat=="Temperate","delta.model"], pch=19, cex=1.5, col="blue")
#points(seq(Xmin,Xmax,1), coef(r.model)[2]*seq(Xmin,Xmax,1) + coef(r.model)[1], type="l", col="black", lty="longdash")
#points(seq(Xmin,Xmax,1), coef(r.model2)[1] + coef(r.model2)[2]*seq(Xmin,Xmax,1) + coef(r.model2)[3]*seq(Xmin,Xmax,1)^2, type="l", col="black")


# RELATIVE LIFETIME FECUNDITY
# Model vs TPCs
Xmin <- -0.6
Xmax <- 0.2
Ymin <- -0.6
Ymax <- 0.2
#dev.new(width=3, height=3, unit="in")
plot(f.data[f.data$Habitat=="Tropical","delta.TPC"], f.data[f.data$Habitat=="Tropical","delta.model"], pch=19, cex=1.5, col="red",
     xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="TPC", ylab="Model")
points(f.data[f.data$Habitat=="Subtropical","delta.TPC"], f.data[f.data$Habitat=="Subtropical","delta.model"], pch=19, cex=1.5, col="orange")
points(f.data[f.data$Habitat=="Mediterranean","delta.TPC"], f.data[f.data$Habitat=="Mediterranean","delta.model"], pch=19, cex=1.5, col="orange")
points(f.data[f.data$Habitat=="Temperate","delta.TPC"], f.data[f.data$Habitat=="Temperate","delta.model"], pch=19, cex=1.5, col="blue")
abline(0, 1, col="gray")
abline(0, 0, col="gray", lty="longdash")
abline(v = 0, col="gray", lty="longdash")

# Model vs TPCs in historical period
Xmin <- 0
Xmax <- 1
Ymin <- 0
Ymax <- 1
#dev.new(width=3, height=3, unit="in")
plot(f.data[f.data$Habitat=="Tropical","TPC.h"], f.data[f.data$Habitat=="Tropical","Model.h"], pch=19, cex=1.5, col="red",
     xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="TPC", ylab="Model")
points(f.data[f.data$Habitat=="Subtropical","TPC.h"], f.data[f.data$Habitat=="Subtropical","Model.h"], pch=19, cex=1.5, col="orange")
points(f.data[f.data$Habitat=="Mediterranean","TPC.h"], f.data[f.data$Habitat=="Mediterranean","Model.h"], pch=19, cex=1.5, col="orange")
points(f.data[f.data$Habitat=="Temperate","TPC.h"], f.data[f.data$Habitat=="Temperate","Model.h"], pch=19, cex=1.5, col="blue")
abline(0, 1, col="gray")
#abline(0, 0, col="gray", lty="longdash")
#abline(v = 0, col="gray", lty="longdash")

# Model vs TPCs in future period
Xmin <- 0
Xmax <- 1
Ymin <- 0
Ymax <- 1
#dev.new(width=3, height=3, unit="in")
plot(f.data[f.data$Habitat=="Tropical","TPC.f"], f.data[f.data$Habitat=="Tropical","Model.f"], pch=19, cex=1.5, col="red",
     xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="TPC", ylab="Model")
points(f.data[f.data$Habitat=="Subtropical","TPC.f"], f.data[f.data$Habitat=="Subtropical","Model.f"], pch=19, cex=1.5, col="orange")
points(f.data[f.data$Habitat=="Mediterranean","TPC.f"], f.data[f.data$Habitat=="Mediterranean","Model.f"], pch=19, cex=1.5, col="orange")
points(f.data[f.data$Habitat=="Temperate","TPC.f"], f.data[f.data$Habitat=="Temperate","Model.f"], pch=19, cex=1.5, col="blue")
abline(0, 1, col="gray")
#abline(0, 0, col="gray", lty="longdash")
#abline(v = 0, col="gray", lty="longdash")

# Model vs latitude
Xmin <- 0
Xmax <- 60
Ymin <- -1
Ymax <- 0.5
#dev.new(width=3, height=3, unit="in")
plot(f.data[f.data$Habitat=="Tropical","Latitude"], f.data[f.data$Habitat=="Tropical","delta.model"], pch=19, cex=1.5, col="red",
     xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Latitude", ylab="change in fecundity")
points(f.data[f.data$Habitat=="Subtropical","Latitude"], f.data[f.data$Habitat=="Subtropical","delta.model"], pch=19, cex=1.5, col="orange")
points(f.data[f.data$Habitat=="Mediterranean","Latitude"], f.data[f.data$Habitat=="Mediterranean","delta.model"], pch=19, cex=1.5, col="orange")
points(f.data[f.data$Habitat=="Temperate","Latitude"], f.data[f.data$Habitat=="Temperate","delta.model"], pch=19, cex=1.5, col="blue")
abline(0, 0, col="gray") #, lty="longdash")
points(seq(Xmin,Xmax,1), coef(f.model)[2]*seq(Xmin,Xmax,1) + coef(f.model)[1], type="l", col="black")
#points(seq(Xmin,Xmax,1), coef(f.model2)[1] + coef(f.model2)[2]*seq(Xmin,Xmax,1) + coef(f.model2)[3]*seq(Xmin,Xmax,1)^2, type="l", col="black")


# SURVIVAL
# Model vs TPCs
Xmin <- -0.6
Xmax <- 0.2
Ymin <- -0.6
Ymax <- 0.2
#dev.new(width=3, height=3, unit="in")
plot(s.data[s.data$Habitat=="Tropical","delta.TPC"], s.data[s.data$Habitat=="Tropical","delta.model"], pch=19, cex=1.5, col="red",
     xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="TPC", ylab="Model")
points(s.data[s.data$Habitat=="Subtropical","delta.TPC"], s.data[s.data$Habitat=="Subtropical","delta.model"], pch=19, cex=1.5, col="orange")
points(s.data[s.data$Habitat=="Mediterranean","delta.TPC"], s.data[s.data$Habitat=="Mediterranean","delta.model"], pch=19, cex=1.5, col="orange")
points(s.data[s.data$Habitat=="Temperate","delta.TPC"], s.data[s.data$Habitat=="Temperate","delta.model"], pch=19, cex=1.5, col="blue")
points(seq(Xmin,Xmax,0.1), coef(delta.s)[2]*seq(Xmin,Xmax,0.1)+coef(delta.s)[1], type="l", col="black")
abline(0, 1, col="gray")
abline(0, 0, col="gray", lty="longdash")
abline(v = 0, col="gray", lty="longdash")

# Model vs TPCs in historical period
Xmin <- 0
Xmax <- 1
Ymin <- 0
Ymax <- 1
#dev.new(width=3, height=3, unit="in")
plot(s.data[s.data$Habitat=="Tropical","TPC.h"], s.data[s.data$Habitat=="Tropical","Model.h"], pch=19, cex=1.5, col="red",
     xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="TPC", ylab="Model")
points(s.data[s.data$Habitat=="Subtropical","TPC.h"], s.data[s.data$Habitat=="Subtropical","Model.h"], pch=19, cex=1.5, col="orange")
points(s.data[s.data$Habitat=="Mediterranean","TPC.h"], s.data[s.data$Habitat=="Mediterranean","Model.h"], pch=19, cex=1.5, col="orange")
points(s.data[s.data$Habitat=="Temperate","TPC.h"], s.data[s.data$Habitat=="Temperate","Model.h"], pch=19, cex=1.5, col="blue")
points(seq(Xmin,Xmax,0.1), coef(s.h)[2]*seq(Xmin,Xmax,0.1)+coef(s.h)[1], type="l", col="black")
abline(0, 1, col="gray")
#abline(0, 0, col="gray", lty="longdash")
#abline(v = 0, col="gray", lty="longdash")

# Model vs TPCs in future period
Xmin <- 0
Xmax <- 1
Ymin <- 0
Ymax <- 1
#dev.new(width=3, height=3, unit="in")
plot(s.data[s.data$Habitat=="Tropical","TPC.f"], s.data[s.data$Habitat=="Tropical","Model.f"], pch=19, cex=1.5, col="red",
     xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="TPC", ylab="Model")
points(s.data[s.data$Habitat=="Subtropical","TPC.f"], s.data[s.data$Habitat=="Subtropical","Model.f"], pch=19, cex=1.5, col="orange")
points(s.data[s.data$Habitat=="Mediterranean","TPC.f"], s.data[s.data$Habitat=="Mediterranean","Model.f"], pch=19, cex=1.5, col="orange")
points(s.data[s.data$Habitat=="Temperate","TPC.f"], s.data[s.data$Habitat=="Temperate","Model.f"], pch=19, cex=1.5, col="blue")
points(seq(Xmin,Xmax,0.1), coef(s.f)[2]*seq(Xmin,Xmax,0.1)+coef(s.f)[1], type="l", col="black")
abline(0, 1, col="gray")
#abline(0, 0, col="gray", lty="longdash")
#abline(v = 0, col="gray", lty="longdash")

# Model vs latitude
Xmin <- 0
Xmax <- 60
Ymin <- -1
Ymax <- 0.5
#dev.new(width=3, height=3, unit="in")
plot(s.data[s.data$Habitat=="Tropical","Latitude"], s.data[s.data$Habitat=="Tropical","delta.model"], pch=19, cex=1.5, col="red",
     xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Latitude", ylab="change in survival")
points(s.data[s.data$Habitat=="Subtropical","Latitude"], s.data[s.data$Habitat=="Subtropical","delta.model"], pch=19, cex=1.5, col="orange")
points(s.data[s.data$Habitat=="Mediterranean","Latitude"], s.data[s.data$Habitat=="Mediterranean","delta.model"], pch=19, cex=1.5, col="orange")
points(s.data[s.data$Habitat=="Temperate","Latitude"], s.data[s.data$Habitat=="Temperate","delta.model"], pch=19, cex=1.5, col="blue")
#points(seq(Xmin,Xmax,1), coef(s.model)[2]*seq(Xmin,Xmax,1) + coef(s.model)[1], type="l", col="black")
#points(seq(Xmin,Xmax,1), coef(s.model2)[1] + coef(s.model2)[2]*seq(Xmin,Xmax,1) + coef(s.model2)[3]*seq(Xmin,Xmax,1)^2, type="l", col="black")
abline(0, 0, col="gray") #, lty="longdash")


# EXTINCTION
# Change in mean temperature
# Model vs TPCs
Xmin <- 0
Xmax <- 25
Ymin <- 0
Ymax <- 25
#dev.new(width=3, height=3, unit="in")
plot(results.m[results.m$Habitat=="Tropical","TPC"], results.m[results.m$Habitat=="Tropical","Model"], pch=19, cex=1.5, col="red",
     xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="TPC", ylab="Model")
points(results.m[results.m$Habitat=="Subtropical","TPC"], results.m[results.m$Habitat=="Subtropical","Model"], pch=19, cex=1.5, col="orange")
points(results.m[results.m$Habitat=="Mediterranean","TPC"], results.m[results.m$Habitat=="Mediterranean","Model"], pch=19, cex=1.5, col="orange")
points(results.m[results.m$Habitat=="Temperate","TPC"], results.m[results.m$Habitat=="Temperate","Model"], pch=19, cex=1.5, col="blue")
points(seq(Xmin,Xmax,0.1), coef(model.m)[2]*seq(Xmin,Xmax,0.1)+coef(model.m)[1], type="l", col="black", lty="longdash")
abline(0, 1, col="gray")
#abline(0, 0, col="gray", lty="longdash")
#abline(v = 0, col="gray", lty="longdash")

# Change in temperature amplitude
# Model vs TPCs
Xmin <- 0
Xmax <- 50
Ymin <- 0
Ymax <- 50
#dev.new(width=3, height=3, unit="in")
plot(results.a[results.a$Habitat=="Tropical","TPC"], results.a[results.a$Habitat=="Tropical","Model"], pch=19, cex=1.5, col="red",
     xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="TPC", ylab="Model")
points(results.a[results.a$Habitat=="Subtropical","TPC"], results.a[results.a$Habitat=="Subtropical","Model"], pch=19, cex=1.5, col="orange")
points(results.a[results.a$Habitat=="Mediterranean","TPC"], results.a[results.a$Habitat=="Mediterranean","Model"], pch=19, cex=1.5, col="orange")
points(results.a[results.a$Habitat=="Temperate","TPC"], results.a[results.a$Habitat=="Temperate","Model"], pch=19, cex=1.5, col="blue")
points(seq(Xmin,Xmax,0.1), coef(model.a)[2]*seq(Xmin,Xmax,0.1)+coef(model.a)[1], type="l", col="black", lty="longdash")
abline(0, 1, col="gray")
#abline(0, 0, col="gray", lty="longdash")
#abline(v = 0, col="gray", lty="longdash")





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
# barplot(c(R0.period,R0.period.CC), col=c("purple","red"), ylim=c(0,100), main=expression("above R0"))
# barplot(c(r.period,r.period.CC), col=c("purple","red"), ylim=c(0,100), main=expression("above r"))
# barplot(c(TSM,TSM.CC), col=c("purple","red"), ylim=c(-20,20), main=expression("Thermal margin"))
# # habitat temperatures
# barplot((sp.data["ext_meanT"])[[1]], col="#30638e", xlim=c(0.2,2), ylim=c(0,30), main=expression("meanT"))
# barplot((sp.data["ext_amplT"])[[1]], col="#30638e", xlim=c(0.2,2), ylim=c(0,30), main=expression("amplT"))
# barplot((sp.data["ext_mT_aT"])[[1]], col="#30638e", xlim=c(0.2,2), ylim=c(0,30), main=expression("meanT + amplT"))
# par(mfrow=c(1,1))
