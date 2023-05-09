###############################################################################
##### This R script analyzes thermal performance curve and DDE model data #####
###############################################################################

# Load packages and set working directory
library(tidyverse)
library(ggplot2)
library(cowplot)

# Set working directory (if necessary)
#setwd() # enter working directory of main downloaded file (containing R project file)

# Read in predictions from "Model predictions" folder
r.data <- as.data.frame(read_csv("Model predictions/Predictions rm.csv"))
R0.data <- as.data.frame(read_csv("Model predictions/Predictions R0.csv"))
b.data <- as.data.frame(read_csv("Model predictions/Predictions birth.csv"))
g.data <- as.data.frame(read_csv("Model predictions/Predictions development.csv"))
s.data <- as.data.frame(read_csv("Model predictions/Predictions survival.csv"))
l.data <- as.data.frame(read_csv("Model predictions/Predictions longevity.csv"))
pop.data <- as.data.frame(read_csv("Model predictions/Predictions population dynamics.csv"))


################################# STATISTICS ###################################
# INTRINSTIC GROWTH RATE, r_m
# DDE model predictions vs latitude (Fig. 3a)
r.lat <- lm(delta.model ~ Latitude, data=r.data)
summary(r.lat)
# DDE model predictions vs directly integrating TPCs (Fig. 4a)
r.delta <- lm(delta.TPC ~ delta.model, data=r.data)
summary(r.delta)

# NET REPRODUCTIVE RATE, R0
# DDE model predictions vs latitude (Fig. 3b)
R0.lat <- lm(delta.model ~ Latitude, data=R0.data)
summary(R0.lat) # significant
# DDE model predictions vs directly integrating TPCs (Fig. 4b)
R0.delta <- lm(delta.TPC ~ delta.model, data=R0.data)
summary(R0.delta)

# SURVIVAL TO REPRODUCTION
# DDE model predictions vs latitude (Fig. 3c)
s.lat <- lm(delta.model ~ Latitude, data=s.data)
summary(s.lat)
# DDE model predictions vs directly integrating TPCs (Fig. 4c)
s.delta <- lm(delta.TPC ~ delta.model, data=s.data)
summary(s.delta)

# BIRTH RATE
# DDE model predictions vs latitude (Fig. 3d)
b.lat <- lm(delta.model ~ Latitude, data=b.data)
summary(b.lat)

# DEVELOPMENT TIME
# DDE model predictions vs latitude (Fig. 3e)
g.lat <- lm(delta.model ~ Latitude, data=g.data)
summary(g.lat)
# DDE model predictions vs directly integrating TPCs (Fig. 4d)
g.delta <- lm(delta.TPC ~ delta.model, data=g.data)
summary(g.delta)

# ADULT LONGEVITY
# DDE model predictions vs latitude (Fig. 3f)
l.lat <- lm(delta.model ~ Latitude, data=l.data)
summary(l.lat)

# MEAN ADULT DENSITY
# DDE model predictions vs latitude (Fig. 5a)
mean.lat <- lm(delta.mean ~ Latitude, data=pop.data)
summary(mean.lat)

# POPULATION VARIATION
# DDE model predictions vs latitude (Fig. 5b)
CV.lat <- lm(delta.CV ~ Latitude, data=pop.data[-c(14,15,18:20),]) # excluding populations that went extinct
summary(CV.lat)



################################### PLOTS ######################################
# INTRINSIC GROWTH RATE, r_m
# DDE model prediction vs latitude (Fig. 3a)
Xmin <- 0
Xmax <- 60
Ymin <- -1.05
Ymax <- 0.18
#dev.new(width=3, height=3, unit="in")
plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Latitude", ylab="Model", cex.axis=2)
abline(0, 0, col="gray", lwd=3, lty="longdash")
points(r.data[r.data$Habitat=="Tropical","Latitude"], r.data[r.data$Habitat=="Tropical","delta.model"], pch=19, cex=1.5, col="#FFB000")
points(r.data[r.data$Habitat=="Subtropical","Latitude"], r.data[r.data$Habitat=="Subtropical","delta.model"], pch=19, cex=1.5, col="#6FD012")
points(r.data[r.data$Habitat=="Temperate","Latitude"], r.data[r.data$Habitat=="Temperate","delta.model"], pch=19, cex=1.5, col="#785EF0")
points(seq(Xmin,Xmax,1), coef(r.lat)[2]*seq(Xmin,Xmax,1) + coef(r.lat)[1], type="l", lwd=3, col="black")

# DDE model predictions vs directly integrating TPCs (Fig. 4a)
Xmin <- -1.2
Xmax <- 0.2
Ymin <- -1.2
Ymax <- 0.2
#dev.new(width=3, height=3, unit="in")
plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Model", ylab="TPC", cex.axis=2)
polygon(c(2*Xmin,2*Xmin,2*Xmax),c(2*Ymin,2*Ymax,2*Ymax), col = "#E2E2E2", border = NA)
abline(0, 1, col="gray", lwd=1)
abline(0, 0, col="gray", lwd=3)
abline(v = 0, col="gray", lwd=3)
points(r.data[r.data$Habitat=="Tropical","delta.model"], r.data[r.data$Habitat=="Tropical","delta.TPC"], pch=19, cex=1.5, col="#FFB000")
points(r.data[r.data$Habitat=="Subtropical","delta.model"], r.data[r.data$Habitat=="Subtropical","delta.TPC"], pch=19, cex=1.5, col="#40B0A6")
points(r.data[r.data$Habitat=="Temperate","delta.model"], r.data[r.data$Habitat=="Temperate","delta.TPC"], pch=19, cex=1.5, col="#785EF0")
points(seq(2*Xmin,2*Xmax,0.1), coef(r.delta)[2]*seq(2*Xmin,2*Xmax,0.1)+coef(r.delta)[1], type="l", lwd=3, col="black")


# NET REPRODUCTIVE RATE, R0
# DDE model prediction vs latitude (Fig. 3b)
Xmin <- 0
Xmax <- 60
Ymin <- -0.8
Ymax <- 0.2
#dev.new(width=3, height=3, unit="in")
plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Latitude", ylab="Model", cex.axis=2)
abline(0, 0, col="gray", lwd=3, lty="longdash")
points(R0.data[R0.data$Habitat=="Tropical","Latitude"], R0.data[R0.data$Habitat=="Tropical","delta.model"], pch=19, cex=1.5, col="#FFB000")
points(R0.data[R0.data$Habitat=="Subtropical","Latitude"], R0.data[R0.data$Habitat=="Subtropical","delta.model"], pch=19, cex=1.5, col="#6FD012")
points(R0.data[R0.data$Habitat=="Temperate","Latitude"], R0.data[R0.data$Habitat=="Temperate","delta.model"], pch=19, cex=1.5, col="#785EF0")
points(seq(Xmin,Xmax,1), coef(R0.lat)[2]*seq(Xmin,Xmax,1) + coef(R0.lat)[1], type="l", lwd=3, col="black")

# DDE model predictions vs directly integrating TPCs (Fig. 4b)
Xmin <- -0.8
Xmax <- 0.2
Ymin <- -0.8
Ymax <- 0.2
#dev.new(width=3, height=3, unit="in")
plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Model", ylab="TPC", cex.axis=2)
polygon(c(2*Xmin,2*Xmin,2*Xmax),c(2*Ymin,2*Ymax,2*Ymax), col = "#E2E2E2", border = NA)
abline(0, 1, col="gray", lwd=1)
abline(0, 0, col="gray", lwd=3)
abline(v = 0, col="gray", lwd=3)
points(R0.data[R0.data$Habitat=="Tropical","delta.model"], R0.data[R0.data$Habitat=="Tropical","delta.TPC"], pch=19, cex=1.5, col="#FFB000")
points(R0.data[R0.data$Habitat=="Subtropical","delta.model"], R0.data[R0.data$Habitat=="Subtropical","delta.TPC"], pch=19, cex=1.5, col="#40B0A6")
points(R0.data[R0.data$Habitat=="Temperate","delta.model"], R0.data[R0.data$Habitat=="Temperate","delta.TPC"], pch=19, cex=1.5, col="#785EF0")
points(seq(2*Xmin,2*Xmax,0.1), coef(R0.delta)[2]*seq(2*Xmin,2*Xmax,0.1)+coef(R0.delta)[1], type="l", lwd=3, col="black")


# SURVIVAL TO REPRODUCTION
# DDE model prediction vs latitude (Fig. 3c)
Xmin <- 0
Xmax <- 60
Ymin <- -1
Ymax <- 0.2
#dev.new(width=3, height=3, unit="in")
plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Latitude", ylab="Model", cex.axis=2)
abline(0, 0, col="gray", lwd=3, lty="longdash")
points(s.data[s.data$Habitat=="Tropical","Latitude"], s.data[s.data$Habitat=="Tropical","delta.model"], pch=19, cex=1.5, col="#FFB000")
points(s.data[s.data$Habitat=="Subtropical","Latitude"], s.data[s.data$Habitat=="Subtropical","delta.model"], pch=19, cex=1.5, col="#6FD012")
points(s.data[s.data$Habitat=="Temperate","Latitude"], s.data[s.data$Habitat=="Temperate","delta.model"], pch=19, cex=1.5, col="#785EF0")
points(seq(Xmin,Xmax,1), coef(s.lat)[2]*seq(Xmin,Xmax,1) + coef(s.lat)[1], type="l", lwd=3, col="black")

# DDE model predictions vs directly integrating TPCs (Fig. 4c)
Xmin <- -1
Xmax <- 0.2
Ymin <- -1
Ymax <- 0.2
#dev.new(width=3, height=3, unit="in")
plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Model", ylab="TPC", cex.axis=2)
polygon(c(2*Xmin,2*Xmin,2*Xmax),c(2*Ymin,2*Ymax,2*Ymax), col = "#E2E2E2", border = NA)
abline(0, 1, col="gray", lwd=1)
abline(0, 0, col="gray", lwd=3)
abline(v = 0, col="gray", lwd=3)
points(s.data[s.data$Habitat=="Tropical","delta.model"], s.data[s.data$Habitat=="Tropical","delta.TPC"], pch=19, cex=1.5, col="#FFB000")
points(s.data[s.data$Habitat=="Subtropical","delta.model"], s.data[s.data$Habitat=="Subtropical","delta.TPC"], pch=19, cex=1.5, col="#40B0A6")
points(s.data[s.data$Habitat=="Temperate","delta.model"], s.data[s.data$Habitat=="Temperate","delta.TPC"], pch=19, cex=1.5, col="#785EF0")
points(seq(2*Xmin,2*Xmax,0.1), coef(s.delta)[2]*seq(2*Xmin,2*Xmax,0.1)+coef(s.delta)[1], type="l", lwd=3, col="black")


# PER CAPITA BIRTH RATE
# DDE model prediction vs latitude (Fig. 3d)
Xmin <- 0
Xmax <- 60
Ymin <- -0.6
Ymax <- 0.2
#dev.new(width=3, height=3, unit="in")
plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Latitude", ylab="Model", cex.axis=2)
abline(0, 0, col="gray", lwd=3, lty="longdash")
points(b.data[b.data$Habitat=="Tropical","Latitude"], b.data[b.data$Habitat=="Tropical","delta.model"], pch=19, cex=1.5, col="#FFB000")
points(b.data[b.data$Habitat=="Subtropical","Latitude"], b.data[b.data$Habitat=="Subtropical","delta.model"], pch=19, cex=1.5, col="#6FD012")
points(b.data[b.data$Habitat=="Temperate","Latitude"], b.data[b.data$Habitat=="Temperate","delta.model"], pch=19, cex=1.5, col="#785EF0")
points(seq(Xmin,Xmax,1), coef(b.lat)[2]*seq(Xmin,Xmax,1) + coef(b.lat)[1], type="l", lwd=3, col="black")


# DEVELOPMENT TIME
# DDE model prediction vs latitude (Fig. 3e)
Xmin <- 0
Xmax <- 60
Ymin <- -5
Ymax <- 0
#dev.new(width=3, height=3, unit="in")
plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Latitude", ylab="Model", cex.axis=2)
abline(0, 0, col="gray", lwd=3, lty="longdash")
points(g.data[g.data$Habitat=="Tropical","Latitude"], g.data[g.data$Habitat=="Tropical","delta.model"], pch=19, cex=1.5, col="#FFB000")
points(g.data[g.data$Habitat=="Subtropical","Latitude"], g.data[g.data$Habitat=="Subtropical","delta.model"], pch=19, cex=1.5, col="#6FD012")
points(g.data[g.data$Habitat=="Temperate","Latitude"], g.data[g.data$Habitat=="Temperate","delta.model"], pch=19, cex=1.5, col="#785EF0")
points(seq(Xmin,Xmax,1), coef(g.lat)[2]*seq(Xmin,Xmax,1) + coef(g.lat)[1], type="l", lwd=3, col="black")

# DDE model predictions vs directly integrating TPCs (Fig. 4d)
Xmin <- -5
Xmax <- 0
Ymin <- -5
Ymax <- 0
#dev.new(width=3, height=3, unit="in")
plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Model", ylab="TPC", cex.axis=2)
polygon(c(2*Xmin,2*Xmin,1),c(2*Xmin,1,1), col = "#E2E2E2", border = NA)
abline(0, 1, col="gray", lwd=1)
abline(0, 0, col="gray", lwd=3)
abline(v = 0, col="gray", lwd=3)
points(g.data[g.data$Habitat=="Tropical","delta.model"], g.data[g.data$Habitat=="Tropical","delta.TPC"], pch=19, cex=1.5, col="#FFB000")
points(g.data[g.data$Habitat=="Subtropical","delta.model"], g.data[g.data$Habitat=="Subtropical","delta.TPC"], pch=19, cex=1.5, col="#40B0A6")
points(g.data[g.data$Habitat=="Temperate","delta.model"], g.data[g.data$Habitat=="Temperate","delta.TPC"], pch=19, cex=1.5, col="#785EF0")
points(seq(2*Xmin,1,0.1), coef(g.delta)[2]*seq(2*Xmin,1,0.1)+coef(g.delta)[1], type="l", lwd=3, col="black")


# ADULT LONGEVITY
# DDE model prediction vs latitude (Fig. 3f)
Xmin <- 0
Xmax <- 60
Ymin <- -0.4
Ymax <- 0.1
#dev.new(width=3, height=3, unit="in")
plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Latitude", ylab="Model", cex.axis=2)
abline(0, 0, col="gray", lwd=3, lty="longdash")
points(l.data[l.data$Habitat=="Tropical","Latitude"], l.data[l.data$Habitat=="Tropical","delta.model"], pch=19, cex=1.5, col="#FFB000")
points(l.data[l.data$Habitat=="Subtropical","Latitude"], l.data[l.data$Habitat=="Subtropical","delta.model"], pch=19, cex=1.5, col="#6FD012")
points(l.data[l.data$Habitat=="Temperate","Latitude"], l.data[l.data$Habitat=="Temperate","delta.model"], pch=19, cex=1.5, col="#785EF0")
points(seq(Xmin,Xmax,1), coef(l.lat)[2]*seq(Xmin,Xmax,1) + coef(l.lat)[1], type="l", lwd=3, col="black")


# MEAN ADULT DENSITY
# DDE model prediction vs latitude (Fig. 5a)
Xmin <- 0
Xmax <- 60
Ymin <- -1
Ymax <- 0.5
#dev.new(width=3, height=3, unit="in")
plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Latitude", ylab="Model", cex.axis=2)
abline(0, 0, col="gray", lwd=3, lty="longdash")
points(pop.data[pop.data$Habitat=="Tropical","Latitude"], pop.data[pop.data$Habitat=="Tropical","delta.mean"], pch=19, cex=1.5, col="#FFB000")
points(pop.data[pop.data$Habitat=="Temperate","Latitude"], pop.data[pop.data$Habitat=="Temperate","delta.mean"], pch=19, cex=1.5, col="#785EF0")
points(pop.data[pop.data$Habitat=="Subtropical","Latitude"], pop.data[pop.data$Habitat=="Subtropical","delta.mean"], pch=19, cex=1.5, col="#6FD012")


# POPULATION VARIATION
# DDE model prediction vs latitude (Fig. 5b)
Xmin <- 0
Xmax <- 60
Ymin <- -1.1
Ymax <- 1.1
#dev.new(width=3, height=3, unit="in")
plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Latitude", ylab="Model", cex.axis=2)
abline(0, 0, col="gray", lwd=3, lty="longdash")
points(seq(2*Xmin,2*Xmax,1), coef(CV.lat)[2]*seq(2*Xmin,2*Xmax,1) + coef(CV.lat)[1], type="l", lwd=3, col="black")
points(pop.data[pop.data$Habitat=="Tropical","Latitude"], pop.data[pop.data$Habitat=="Tropical","delta.CV"], pch=19, cex=1.5, col="#FFB000")
points(pop.data[pop.data$Habitat=="Subtropical","Latitude"], pop.data[pop.data$Habitat=="Subtropical","delta.CV"], pch=19, cex=1.5, col="#6FD012")
points(pop.data[pop.data$Habitat=="Temperate" & pop.data$CV.f != 0,"Latitude"], pop.data[pop.data$Habitat=="Temperate" & pop.data$CV.f != 0,"delta.CV"], pch=19, cex=1.5, col="#785EF0") # excluding population that went extinct


# ACTIVE PERIOD
# DDE model prediction vs latitude (Fig. 5c)
Xmin <- 0
Xmax <- 60
Ymin <- 0
Ymax <- 0.4
#dev.new(width=3, height=3, unit="in")
plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Latitude", ylab="Model", cex.axis=2)
abline(0, 0, col="gray", lwd=3, lty="longdash")
points(pop.data[pop.data$Habitat=="Tropical","Latitude"], pop.data[pop.data$Habitat=="Tropical","delta.active"], pch=19, cex=1.5, col="#FFB000")
points(pop.data[pop.data$Habitat=="Subtropical","Latitude"], pop.data[pop.data$Habitat=="Subtropical","delta.active"], pch=19, cex=1.5, col="#6FD012")
points(pop.data[pop.data$Habitat=="Temperate","Latitude"], pop.data[pop.data$Habitat=="Temperate","delta.active"], pch=19, cex=1.5, col="#785EF0")