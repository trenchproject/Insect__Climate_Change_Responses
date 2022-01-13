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
# r
#r.data <- as.data.frame(read_csv("Model results Diurnal.csv"))
r.data <- as.data.frame(read_csv("Model results Tave.csv"))

######################################## STATISTICS #########################################
# COMPARING CHANGE IN PER CAPITA GROWTH RATE BETWEEN MODEL AND TPC
# Compare model vs TPC
# change in r (model vs TPCs)
delta <- lm(delta.model ~ delta.TPC, data=r.data) # significant!
summary(delta)
# r in historical period (model vs TPCs)
r.h <- lm(r.model.h ~ r.TPC.h, data=r.data) # significant!
summary(r.h)
# r in future period (model vs TPCs)
r.f <- lm(r.model.f ~ r.TPC.f, data=r.data) # significant!
summary(r.f)


# TPC vs Latitude
# linear
TPC <- lm(delta.TPC ~ Latitude, data=r.data)
summary(TPC) # non-significant
# non-linear
TPC2 <- nls(delta.TPC ~ a + b*Latitude + c*Latitude^2, data=r.data, start=list(a=1, b=-0.1, c=1))
summary(TPC2) # non-significant

# Model vs Latitude
# linear
model <- lm(delta.model ~ Latitude, data=r.data)
summary(model) # non-significant
# non-linear
model2 <- nls(delta.model ~ a + b*Latitude + c*Latitude^2, data=r.data, start=list(a=1, b=-0.1, c=1))
summary(model2) # non-significant


# # COMPARING PROPORTIONAL CHANGE IN r
# # TPC
# # linear
# TPC.prop <- lm(delta.prop.TPC ~ Latitude, data=r.data)
# summary(TPC.prop) # non-significant
# # non-linear
# TPC2.prop <- nls(delta.prop.TPC ~ a + b*Latitude + c*Latitude^2, data=r.data, start=list(a=1, b=-0.1, c=1))
# summary(TPC2.prop) # non-significant
# 
# # Model
# # linear
# model.prop <- lm(delta.prop.model ~ Latitude, data=r.data)
# summary(model.prop) # non-significant
# # non-linear
# model2.prop <- nls(delta.prop.model ~ a + b*Latitude + c*Latitude^2, data=r.data, start=list(a=1, b=-0.1, c=1))
# summary(model2.prop) # non-significant
# 
# # Compare model vs TPC
# # proportional change in r (model vs TPCs)
# d <- lm(delta.prop.model ~ delta.prop.TPC, data=r.data) # non-significant
# summary(d)


########################################### PLOTS ###########################################
# CHANGE IN r
# Model vs TPCs
Xmin <- -2
Xmax <- 0.5
Ymin <- -1.5
Ymax <- 0.5
plot(r.data[r.data$Habitat=="Tropical","delta.TPC"], r.data[r.data$Habitat=="Tropical","delta.model"], pch=21, col="red", bg="red",
     xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="TPC", ylab="Model")
points(r.data[r.data$Habitat=="Subtropical","delta.TPC"], r.data[r.data$Habitat=="Subtropical","delta.model"], pch=21, col="orange", bg="orange")
points(r.data[r.data$Habitat=="Mediterranean","delta.TPC"], r.data[r.data$Habitat=="Mediterranean","delta.model"], pch=21, col="orange", bg="orange")
points(r.data[r.data$Habitat=="Temperate","delta.TPC"], r.data[r.data$Habitat=="Temperate","delta.model"], pch=21, col="blue", bg="blue")
points(seq(Xmin,Xmax,0.1), coef(delta)[2]*seq(Xmin,Xmax,0.1)+coef(delta)[1], type="l", col="black", lty="longdash")
abline(0, 1, col="gray")
#abline(v = 0, col="gray")

# Model vs TPCs in historical period
Xmin <- 0 #-1
Xmax <- 1 #2.5
Ymin <- 0 #-1
Ymax <- 1 #2.5
plot(r.data[r.data$Habitat=="Tropical","r.TPC.h"], r.data[r.data$Habitat=="Tropical","r.model.h"], pch=21, col="red", bg="red",
     xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="TPC", ylab="Model")
points(r.data[r.data$Habitat=="Subtropical","r.TPC.h"], r.data[r.data$Habitat=="Subtropical","r.model.h"], pch=21, col="orange", bg="orange")
points(r.data[r.data$Habitat=="Mediterranean","r.TPC.h"], r.data[r.data$Habitat=="Mediterranean","r.model.h"], pch=21, col="orange", bg="orange")
points(r.data[r.data$Habitat=="Temperate","r.TPC.h"], r.data[r.data$Habitat=="Temperate","r.model.h"], pch=21, col="blue", bg="blue")
points(seq(Xmin,Xmax,0.1), coef(r.h)[2]*seq(Xmin,Xmax,0.1)+coef(r.h)[1], type="l", col="black", lty="longdash")
abline(0, 1, col="gray")
#abline(v = 0, col="gray")

# Model vs TPCs in future period
Xmin <- -1.5 #-4
Xmax <- 1.5 #2
Ymin <- -1.5 #-4
Ymax <- 1.5 #2
plot(r.data[r.data$Habitat=="Tropical","r.TPC.f"], r.data[r.data$Habitat=="Tropical","r.model.f"], pch=21, col="red", bg="red",
     xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="TPC", ylab="Model")
points(r.data[r.data$Habitat=="Subtropical","r.TPC.f"], r.data[r.data$Habitat=="Subtropical","r.model.f"], pch=21, col="orange", bg="orange")
points(r.data[r.data$Habitat=="Mediterranean","r.TPC.f"], r.data[r.data$Habitat=="Mediterranean","r.model.f"], pch=21, col="orange", bg="orange")
points(r.data[r.data$Habitat=="Temperate","r.TPC.f"], r.data[r.data$Habitat=="Temperate","r.model.f"], pch=21, col="blue", bg="blue")
points(seq(Xmin,Xmax,0.1), coef(r.f)[2]*seq(Xmin,Xmax,0.1)+coef(r.f)[1], type="l", col="black", lty="longdash")
abline(0, 1, col="gray")
#abline(v = 0, col="gray")


# TPC vs latitude
Xmin <- 0
Xmax <- 60
Ymin <- -2 #-4
Ymax <- 0.5 #1
plot(r.data[r.data$Habitat=="Tropical","Latitude"], r.data[r.data$Habitat=="Tropical","delta.TPC"], pch=21, col="red", bg="red",
     xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Latitude", ylab="change in r")
points(r.data[r.data$Habitat=="Subtropical","Latitude"], r.data[r.data$Habitat=="Subtropical","delta.TPC"], pch=21, col="orange", bg="orange")
points(r.data[r.data$Habitat=="Mediterranean","Latitude"], r.data[r.data$Habitat=="Mediterranean","delta.TPC"], pch=21, col="orange", bg="orange")
points(r.data[r.data$Habitat=="Temperate","Latitude"], r.data[r.data$Habitat=="Temperate","delta.TPC"], pch=21, col="blue", bg="blue")
#points(seq(Xmin,Xmax,1), coef(TPC)[2]*seq(Xmin,Xmax,1) + coef(TPC)[1], type="l", col="black", lty="longdash")
#points(seq(Xmin,Xmax,1), coef(TPC2)[1] + coef(TPC2)[2]*seq(Xmin,Xmax,1) + coef(TPC2)[3]*seq(Xmin,Xmax,1)^2, type="l", col="black")
abline(0, 0, col="black") #, lty="longdash")

# Model vs latitude
Xmin <- 0
Xmax <- 60
Ymin <- -1 #-2.5
Ymax <- 0.5
plot(r.data[r.data$Habitat=="Tropical","Latitude"], r.data[r.data$Habitat=="Tropical","delta.model"], pch=21, col="red", bg="red",
     xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Latitude", ylab="change in r")
points(r.data[r.data$Habitat=="Subtropical","Latitude"], r.data[r.data$Habitat=="Subtropical","delta.model"], pch=21, col="orange", bg="orange")
points(r.data[r.data$Habitat=="Mediterranean","Latitude"], r.data[r.data$Habitat=="Mediterranean","delta.model"], pch=21, col="orange", bg="orange")
points(r.data[r.data$Habitat=="Temperate","Latitude"], r.data[r.data$Habitat=="Temperate","delta.model"], pch=21, col="blue", bg="blue")
#points(seq(Xmin,Xmax,1), coef(model)[2]*seq(Xmin,Xmax,1) + coef(model)[1], type="l", col="black", lty="longdash")
#points(seq(Xmin,Xmax,1), coef(model2)[1] + coef(model2)[2]*seq(Xmin,Xmax,1) + coef(model2)[3]*seq(Xmin,Xmax,1)^2, type="l", col="black")
abline(0, 0, col="black") #, lty="longdash")



##################################### TPC METRICS ###########################################
# # Temperature response parameters
# data <- as.data.frame(read_csv("Temperature response parameters.csv"))
# 
# 
# 
# # EXTINCTION THRESHOLD (increase in meanT that leads to extinction in the model)
# # Versus intrinsic growth rate (r)
# # all parameters: rTopt significant only in full model
# r <- lm(ext_meanT ~ rTopt + Toptr + rTmax, data=data)
# summary (r)
# r.Topt <- lm(ext_meanT ~ rTopt, data=data)
# summary (r.Topt)
# 
# # thermal safety margin: non significant
# tsm <- lm(ext_meanT ~ TSM, data=data)
# summary (tsm)
# 
# # warming tolerance: significant
# wt <- lm(ext_meanT ~ WT, data=data)
# #wt <- lm(ext_meanT ~ 0 + WT, data=data) # remove intercept
# summary (wt)
# 
# # warming tolerance (active period): significant
# wt.active <- lm(ext_meanT ~ WT_active, data=data)
# #wt.active <- lm(ext_meanT ~ 0 + WT_active, data=data) # remove intercept
# summary (wt.active)
# 
# 
# # Versus net reproductive rate (R0)
# # all parameters: sR0 significant
# R0 <- lm(ext_meanT ~ R0Topt + ToptR0 + sR0, data=data)
# summary (R0)
# 
# # R0: Tmax (temperature at which R0 = 1): marginally significant
# R0.Tmax <- lm(ext_meanT ~ R0Tmax, data=data)
# summary (R0.Tmax)
# 
# # R0 thermal safety margin: significant
# R0.tsm <- lm(ext_meanT ~ R0TSM, data=data)
# R0.tsm <- lm(ext_meanT ~ 0 + R0TSM, data=data) # remove intercept
# summary (R0.tsm)
# 
# # R0 warming tolerance: significant
# R0.wt <- lm(ext_meanT ~ R0WT, data=data)
# R0.wt <- lm(ext_meanT ~ 0 + R0WT, data=data) # remove intercept
# summary (R0.wt)
# 
# # warming tolerance (active period): significant
# R0.wt.active <- lm(ext_meanT ~ R0WT_active, data=data)
# R0.wt.active <- lm(ext_meanT ~ 0 + R0WT_active, data=data) # remove intercept
# summary (R0.wt.active)
# 
# 
# # Extinction risk across latitude
# lat <- lm(ext_meanT ~ Latitude, data=data)
# #lat <- lm(ext_meanT ~ 0 + Latitude, data=data) # remove intercept
# summary (lat)
# 
# 
# # PLOTS
# # Metrics based on intrinsic growth rate (r)
# # thermal safety margin: non significant
# Xmin <- 0
# Xmax <- 12
# Ymin <- 0
# Ymax <- 15
# plot(data[data$Habitat=="Tropical","TSM"], data[data$Habitat=="Tropical","ext_meanT"], pch=21, col="red", bg="red",
#      xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Thermal safety margin", ylab="Extinction threshold (mean T)")
# points(data[data$Habitat=="Mediterranean","TSM"], data[data$Habitat=="Mediterranean","ext_meanT"], pch=21, col="purple", bg="purple")
# points(data[data$Habitat=="Temperate","TSM"], data[data$Habitat=="Temperate","ext_meanT"], pch=21, col="blue", bg="blue")
# #points(seq(Xmin,Xmax,1), coef(tsm)[1]*seq(Xmin,Xmax,1), type="l", col="blue")
# points(seq(Xmin,Xmax,1), coef(tsm)[2]*seq(Xmin,Xmax,1)+coef(tsm)[1], type="l", col="black", lty="longdash")
# #abline(0, 1, col="black", lty="longdash")
# 
# # warming tolerance: significant (slope = 0.81)
# Xmin <- 0
# Xmax <- 20
# Ymin <- 0
# Ymax <- 15
# plot(data[data$Habitat=="Tropical","WT"], data[data$Habitat=="Tropical","ext_meanT"], pch=21, col="red", bg="red",
#      xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Warming tolerance", ylab="Extinction threshold (mean T)")
# points(data[data$Habitat=="Mediterranean","WT"], data[data$Habitat=="Mediterranean","ext_meanT"], pch=21, col="purple", bg="purple")
# points(data[data$Habitat=="Temperate","WT"], data[data$Habitat=="Temperate","ext_meanT"], pch=21, col="blue", bg="blue")
# #points(seq(Xmin,Xmax,1), coef(wt)[1]*seq(Xmin,Xmax,1), type="l", col="blue")
# points(seq(Xmin,Xmax,1), coef(wt)[2]*seq(Xmin,Xmax,1)+coef(wt)[1], type="l", col="black")
# abline(0, 1, col="black", lty="longdash")
# 
# # warming tolerance (active period): significant (slope = 1.01)
# Xmin <- 0
# Xmax <- 20
# Ymin <- 0
# Ymax <- 15
# plot(data[data$Habitat=="Tropical","WT_active"], data[data$Habitat=="Tropical","ext_meanT"], pch=21, col="red", bg="red",
#      xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Warming tolerance (active period)", ylab="Extinction threshold (mean T)")
# points(data[data$Habitat=="Mediterranean","WT_active"], data[data$Habitat=="Mediterranean","ext_meanT"], pch=21, col="purple", bg="purple")
# points(data[data$Habitat=="Temperate","WT_active"], data[data$Habitat=="Temperate","ext_meanT"], pch=21, col="blue", bg="blue")
# #points(seq(Xmin,Xmax,1), coef(wt.active)[1]*seq(Xmin,Xmax,1), type="l", col="blue")
# points(seq(Xmin,Xmax,1), coef(wt.active)[2]*seq(Xmin,Xmax,1)+coef(wt.active)[1], type="l", col="black")
# abline(0, 1, col="black", lty="longdash")
# 
# 
# 
# # Net reproductive rate (R0)
# # R0 thermal safety margin: significant (slope = 0.46)
# Xmin <- 10
# Xmax <- 30
# Ymin <- 0
# Ymax <- 15
# plot(data$R0TSM, data$ext_meanT, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), main=expression("extinction threshold (meanT) vs. R0 TSM"))
# points(seq(Xmin,Xmax,1), coef(R0.tsm)[1]*seq(Xmin,Xmax,1), type="l", col="blue")
# abline(0, 1)
# 
# # warming tolerance: significant (slope = 0.38)
# Xmin <- 0
# Xmax <- 40
# Ymin <- 0
# Ymax <- 15
# plot(data$R0WT, data$ext_meanT, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), main=expression("extinction threshold (meanT) vs. R0 WT"))
# points(seq(Xmin,Xmax,1), coef(R0.wt)[1]*seq(Xmin,Xmax,1), type="l", col="blue")
# abline(0, 1)
# 
# # warming tolerance (active period): significant (slope = 0.5)
# Xmin <- 0
# Xmax <- 30
# Ymin <- 0
# Ymax <- 15
# plot(data$R0WT_active, data$ext_meanT, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), main=expression("extinction threshold (meanT) vs. R0 WT (active)"))
# points(seq(Xmin,Xmax,1), coef(R0.wt.active)[1]*seq(Xmin,Xmax,1), type="l", col="blue")
# abline(0, 1)
# 
# 
# 
# # Extinction risk across latitude
# Xmin <- 0
# Xmax <- 45
# Ymin <- 0
# Ymax <- 15
# lat.plot <- ggplot(data, aes(Latitude, ext_meanT)) +
#   geom_point(size=3, aes(color=Habitat)) +
#   scale_colour_manual(values=c("purple", "blue", "red")) +
#   geom_function(fun = function(t) (coef(lat)[2]*t+coef(lat)[1]),
#                 size=0.8, linetype="longdash", color="black") +
#   labs(x="", y="") + #labs(x="Absolute latitude", y="Extinction threshold (mean T)") +
#   scale_x_continuous(limits=c(Xmin, Xmax)) +
#   scale_y_continuous(limits=c(Ymin, Ymax)) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
#         axis.line = element_line(colour = "black"), legend.position = "none", 
#         axis.text = element_text(size=13), axis.title = element_text(size=20))
# ggdraw()  + draw_plot(lat.plot, x = 0, y = 0, width = 1, height = 0.4)



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
