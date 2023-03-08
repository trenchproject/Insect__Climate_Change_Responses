#############################################################################
#### This R script analyzes thermal performance curve and DDE model data ####
#############################################################################


# Load packages and set working directory
library(tidyverse)
library(ggplot2)
library(cowplot)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# READ IN DATA
# Fitness metrics and components
# r.data <- as.data.frame(read_csv("Predictions/Predictions Dev fitness 2.csv"))
# R0.data <- as.data.frame(read_csv("Predictions/Predictions Dev R0 2.csv"))
# b.data <- as.data.frame(read_csv("Predictions/Predictions Dev birth 2.csv"))
# tau.data <- as.data.frame(read_csv("Predictions/Predictions Dev development 2.csv"))
# s.data <- as.data.frame(read_csv("Predictions/Predictions Dev survival 2.csv"))
# L.data <- as.data.frame(read_csv("Predictions/Predictions Dev longevity 2.csv"))

r.data <- as.data.frame(read_csv("Predictions new/Predictions rm.csv"))
R0.data <- as.data.frame(read_csv("Predictions new/Predictions R0.csv"))
b.data <- as.data.frame(read_csv("Predictions new/Predictions birth.csv"))
tau.data <- as.data.frame(read_csv("Predictions new/Predictions development.csv"))
s.data <- as.data.frame(read_csv("Predictions new/Predictions survival.csv"))
L.data <- as.data.frame(read_csv("Predictions new/Predictions longevity.csv"))
# Population dynamics
pop.data <- as.data.frame(read_csv("Predictions/Predictions population dynamics.csv"))


# COMPARE OUTPUTS
old <- as.data.frame(read_csv("Predictions/Predictions Dev longevity.csv"))
old <- old[-c(3,12,13,15,18),]
old <- old[,-(2:4)]
rnd <- 2
old$TPC.h <- round(old$TPC.h,rnd)
old$TPC.f <- round(old$TPC.f,rnd)
old$Model.h <- round(old$Model.h,rnd)
old$Model.f <- round(old$Model.f,rnd)
old$delta.TPC <- round(old$delta.TPC,rnd)
old$delta.model <- round(old$delta.model,rnd)
new <- as.data.frame(read_csv("Predictions new/Predictions longevity.csv"))
new <- new[,-(2:5)]
new$TPC.h <- round(new$TPC.h,rnd)
new$TPC.f <- round(new$TPC.f,rnd)
new$Model.h <- round(new$Model.h,rnd)
new$Model.f <- round(new$Model.f,rnd)
new$delta.TPC <- round(new$delta.TPC,rnd)
new$delta.model <- round(new$delta.model,rnd)
summary(arsenal::comparedf(old,new))


# EXCLUDE DATA
# Clavigralla tomentosicollis Nigeria (used data from Burkina Faso to parameterize birth rate and adult mortality)
# Macrolophus pygmaeus (only predator and separate thermal responses for each prey)
# Bemisia argentifollii and Diaphorina citri (different suborder)
# r.data <- r.data[-c(3,12,13,15,18),]
# R0.data <- R0.data[-c(3,12,13,15,18),]
# b.data <- b.data[-c(3,12,13,15,18),]
# tau.data$active.h <- pop.data$active.h # add activity period data to developmental data
# tau.data$active.f <- pop.data$active.f # add activity period data to developmental data
# tau.data <- tau.data[-c(3,12,13,15,18),]
# s.data <- s.data[-c(3,12,13,15,18),]
# L.data <- L.data[-c(3,12,13,15,18),]
# pop.data <- pop.data[-c(3,12,13,15,18),]


# SCALE LOWEST FITNESS CHANGE TO -1
#r.data$delta.TPC <- pmax(-1, r.data$delta.TPC)
#r.data$delta.model <- pmax(-1, r.data$delta.model)


# QUANTIFY SIGNS AND LOG RATIOS OF MODEL PREDICTIONS TO TPC ESTIMATES
# func <- function(TPC,model) { ifelse(model >= 0,
#                                      ifelse(TPC < 0.95*model || TPC > 1.05*model, model - TPC, 0),
#                                      ifelse(TPC > 0.95*model || TPC < 1.05*model, model - TPC, 0)) }
# # Fitness
# for(i in seq(1,nrow(r.data),1)) { r.data$sign[i] <- sign(func(r.data$delta.TPC[i],r.data$delta.model[i])) }
# r.data$ratio <- abs(r.data$delta.TPC/r.data$delta.model)
# r.data$log.ratio <- log(r.data$ratio)
# # R0
# for(i in seq(1,nrow(R0.data),1)) { R0.data$sign[i] <- sign(func(R0.data$delta.TPC[i],R0.data$delta.model[i])) }
# R0.data$ratio <- abs(R0.data$delta.TPC/R0.data$delta.model)
# R0.data$log.ratio <- log(R0.data$ratio)
# # Survival
# for(i in seq(1,nrow(s.data),1)) { s.data$sign[i] <- sign(func(s.data$delta.TPC[i],s.data$delta.model[i])) }
# s.data$ratio <- abs(s.data$delta.TPC/s.data$delta.model)
# s.data$log.ratio <- log(s.data$ratio)
# # Development
# for(i in seq(1,nrow(tau.data),1)) { tau.data$sign[i] <- sign(func(tau.data$delta.TPC[i],tau.data$delta.model[i])) }
# tau.data$ratio <- abs(tau.data$delta.TPC/tau.data$delta.model)
# tau.data$log.ratio <- log(tau.data$ratio)



######################################## STATISTICS #########################################
# FITNESS
# Model vs Latitude (Fig. 3a)
r.lat <- lm(delta.model ~ Latitude, data=r.data)
summary(r.lat) # significant!
# Model vs TPC (Fig. 4a)
# Log ratio of model to TPC
t.test(r.data$log.ratio, mu=0) # non-significant
# Exact binomial test (fraction underestimated)
#binom.test(11, 22, p=0.5, alternative = "two.sided") # non-significant
# Correlation
r.delta <- lm(delta.TPC ~ delta.model, data=r.data)
summary(r.delta) # significant!

# R0
# Model vs Latitude (Fig. 3b)
R0.lat <- lm(delta.model ~ Latitude, data=R0.data)
summary(R0.lat) # significant!
# Model vs TPC (Fig. 4b)
# Log ratio of model to TPC
t.test(R0.data$log.ratio, mu=0) # non-significant
# Exact binomial test (fraction underestimated)
#binom.test(19, 22, p=0.5, alternative = "two.sided") # significant
# Correlation
R0.delta <- lm(delta.TPC ~ delta.model, data=R0.data)
summary(R0.delta) # significant!

# SURVIVAL
# Model vs Latitude (Fig. 3c)
s.lat <- lm(delta.model ~ Latitude, data=s.data)
summary(s.lat) # marginally-significant
# Model vs TPC (Fig. 4c)
# Log ratio of model to TPC
t.test(s.data$log.ratio, mu=0) # non-significant
# Exact binomial test (fraction underestimated)
#binom.test(9, 22, p=0.5, alternative = "two.sided") # non-significant
# Correlation
s.delta <- lm(delta.TPC ~ delta.model, data=s.data)
summary(s.delta)  # significant!

# BIRTH RATE
# Model vs Latitude (Fig. 3d)
b.lat <- lm(delta.model ~ Latitude, data=b.data)
summary(b.lat) # marginally-significant

# DEVELOPMENT TIME
# Model vs Latitude (Fig. 3e)
tau.lat <- lm(delta.model ~ Latitude, data=tau.data)
summary(tau.lat) # significant!
# Model vs TPC (Fig. 4e)
# Log ratio of model to TPC
t.test(tau.data$log.ratio, mu=0) # significant!
# Exact binomial test (fraction underestimated)
#binom.test(18, 22, p=0.5, alternative = "two.sided") # significant!
# Correlation
tau.delta <- lm(delta.TPC ~ delta.model, data=tau.data)
summary(tau.delta)  # significant!

# ADULT LONGEVITY
# Model vs Latitude (Fig. 3f)
L.lat <- lm(delta.model ~ Latitude, data=L.data)
summary(L.lat) # significant!


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
# Latitude (Fig. 3a)
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

# Model vs TPCs (Fig. 4a)
Xmin <- -1
Xmax <- 0.5
Ymin <- -1
Ymax <- 0.5
#dev.new(width=3, height=3, unit="in")
plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Model", ylab="TPC", cex.axis=2)
polygon(c(2*Xmin,2*Xmin,2*Xmax),c(2*Ymin,2*Ymax,2*Ymax), col = "#E2E2E2", border = NA)
abline(0, 1, col="gray", lwd=1)
abline(0, 0, col="gray", lwd=3)
abline(v = 0, col="gray", lwd=3)
points(r.data[r.data$Habitat=="Tropical","delta.model"], r.data[r.data$Habitat=="Tropical","delta.TPC"], pch=19, cex=1.5, col="#FFB000") # orange
points(r.data[r.data$Habitat=="Subtropical","delta.model"], r.data[r.data$Habitat=="Subtropical","delta.TPC"], pch=19, cex=1.5, col="#40B0A6") # teal
points(r.data[r.data$Habitat=="Mediterranean","delta.model"], r.data[r.data$Habitat=="Mediterranean","delta.TPC"], pch=19, cex=1.5, col="#40B0A6") # teal
points(r.data[r.data$Habitat=="Temperate","delta.model"], r.data[r.data$Habitat=="Temperate","delta.TPC"], pch=19, cex=1.5, col="#785EF0") # purple
points(seq(2*Xmin,2*Xmax,0.1), coef(r.delta)[2]*seq(2*Xmin,2*Xmax,0.1)+coef(r.delta)[1], type="l", lwd=3, col="black")


# R0
# Latitude (Fig. 3b)
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
points(seq(Xmin,Xmax,1), coef(R0.lat)[2]*seq(Xmin,Xmax,1) + coef(R0.lat)[1], type="l", lwd=3, col="black")

# Model vs TPCs (Fig. 4b)
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
points(R0.data[R0.data$Habitat=="Tropical","delta.model"], R0.data[R0.data$Habitat=="Tropical","delta.TPC"], pch=19, cex=1.5, col="#FFB000") # orange
points(R0.data[R0.data$Habitat=="Subtropical","delta.model"], R0.data[R0.data$Habitat=="Subtropical","delta.TPC"], pch=19, cex=1.5, col="#40B0A6") # teal
points(R0.data[R0.data$Habitat=="Mediterranean","delta.model"], R0.data[R0.data$Habitat=="Mediterranean","delta.TPC"], pch=19, cex=1.5, col="#40B0A6") # teal
points(R0.data[R0.data$Habitat=="Temperate","delta.model"], R0.data[R0.data$Habitat=="Temperate","delta.TPC"], pch=19, cex=1.5, col="#785EF0") # purple
points(seq(2*Xmin,2*Xmax,0.1), coef(R0.delta)[2]*seq(2*Xmin,2*Xmax,0.1)+coef(R0.delta)[1], type="l", lwd=3, col="black")


# SURVIVAL
# Latitude (Fig. 3c)
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

# Model vs TPCs (Fig. 4c)
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
points(s.data[s.data$Habitat=="Tropical","delta.model"], s.data[s.data$Habitat=="Tropical","delta.TPC"], pch=19, cex=1.5, col="#FFB000") # orange
points(s.data[s.data$Habitat=="Subtropical","delta.model"], s.data[s.data$Habitat=="Subtropical","delta.TPC"], pch=19, cex=1.5, col="#40B0A6") # teal
points(s.data[s.data$Habitat=="Mediterranean","delta.model"], s.data[s.data$Habitat=="Mediterranean","delta.TPC"], pch=19, cex=1.5, col="#40B0A6") # teal
points(s.data[s.data$Habitat=="Temperate","delta.model"], s.data[s.data$Habitat=="Temperate","delta.TPC"], pch=19, cex=1.5, col="#785EF0") # purple
points(seq(2*Xmin,2*Xmax,0.1), coef(s.delta)[2]*seq(2*Xmin,2*Xmax,0.1)+coef(s.delta)[1], type="l", lwd=3, col="black")


# BIRTH RATE
# Latitude (Fig. 3d)
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

# Model vs TPCs
Xmin <- -0.6
Xmax <- 0.2
Ymin <- -0.6
Ymax <- 0.2
#dev.new(width=3, height=3, unit="in")
plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Model", ylab="TPC", cex.axis=2)
polygon(c(2*Xmin,2*Xmin,2*Xmax),c(2*Ymin,2*Ymax,2*Ymax), col = "#E2E2E2", border = NA)
abline(0, 1, col="gray", lwd=1)
abline(0, 0, col="gray", lwd=3)
abline(v = 0, col="gray", lwd=3)
points(b.data[b.data$Habitat=="Tropical","delta.model"], b.data[b.data$Habitat=="Tropical","delta.TPC"], pch=19, cex=1.5, col="#FFB000") # orange
points(b.data[b.data$Habitat=="Subtropical","delta.model"], b.data[b.data$Habitat=="Subtropical","delta.TPC"], pch=19, cex=1.5, col="#40B0A6") # teal
points(b.data[b.data$Habitat=="Mediterranean","delta.model"], b.data[b.data$Habitat=="Mediterranean","delta.TPC"], pch=19, cex=1.5, col="#40B0A6") # teal
points(b.data[b.data$Habitat=="Temperate","delta.model"], b.data[b.data$Habitat=="Temperate","delta.TPC"], pch=19, cex=1.5, col="#785EF0") # purple


# DEVELOPMENT TIME
# Latitude (Fig. 3e)
Xmin <- 0
Xmax <- 60
Ymin <- -5
Ymax <- 0
#dev.new(width=3, height=3, unit="in")
plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Latitude", ylab="Model", cex.axis=2)
abline(0, 0, col="gray", lwd=3, lty="longdash")
points(tau.data[tau.data$Habitat=="Tropical","Latitude"], tau.data[tau.data$Habitat=="Tropical","delta.model"], pch=19, cex=1.5, col="#FFB000") # orange
points(tau.data[tau.data$Habitat=="Subtropical","Latitude"], tau.data[tau.data$Habitat=="Subtropical","delta.model"], pch=19, cex=1.5, col="#6FD012") # green
points(tau.data[tau.data$Habitat=="Mediterranean","Latitude"], tau.data[tau.data$Habitat=="Mediterranean","delta.model"], pch=19, cex=1.5, col="#6FD012") # green
points(tau.data[tau.data$Habitat=="Temperate","Latitude"], tau.data[tau.data$Habitat=="Temperate","delta.model"], pch=19, cex=1.5, col="#785EF0") # purple
points(seq(Xmin,Xmax,1), coef(tau.lat)[2]*seq(Xmin,Xmax,1) + coef(tau.lat)[1], type="l", lwd=3, col="black")

# Model vs TPCs (Fig. 4e)
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
points(tau.data[tau.data$Habitat=="Tropical","delta.model"], tau.data[tau.data$Habitat=="Tropical","delta.TPC"], pch=19, cex=1.5, col="#FFB000") # orange
points(tau.data[tau.data$Habitat=="Subtropical","delta.model"], tau.data[tau.data$Habitat=="Subtropical","delta.TPC"], pch=19, cex=1.5, col="#40B0A6") # teal
points(tau.data[tau.data$Habitat=="Mediterranean","delta.model"], tau.data[tau.data$Habitat=="Mediterranean","delta.TPC"], pch=19, cex=1.5, col="#40B0A6") # teal
points(tau.data[tau.data$Habitat=="Temperate","delta.model"], tau.data[tau.data$Habitat=="Temperate","delta.TPC"], pch=19, cex=1.5, col="#785EF0") # purple
points(seq(2*Xmin,1,0.1), coef(tau.delta)[2]*seq(2*Xmin,1,0.1)+coef(tau.delta)[1], type="l", lwd=3, col="black")


# ADULT LONGEVITY
# Latitude (Fig. 3f)
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

# Model vs TPCs
Xmin <- -0.4
Xmax <- 0.1
Ymin <- -0.4
Ymax <- 0.1
#dev.new(width=3, height=3, unit="in")
plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Model", ylab="TPC", cex.axis=2)
polygon(c(2*Xmin,2*Xmin,2*Xmax),c(2*Ymin,2*Ymax,2*Ymax), col = "#E2E2E2", border = NA)
abline(0, 1, col="gray", lwd=1)
abline(0, 0, col="gray", lwd=3)
abline(v = 0, col="gray", lwd=3)
points(L.data[L.data$Habitat=="Tropical","delta.model"], L.data[L.data$Habitat=="Tropical","delta.TPC"], pch=19, cex=1.5, col="#FFB000") # orange
points(L.data[L.data$Habitat=="Subtropical","delta.model"], L.data[L.data$Habitat=="Subtropical","delta.TPC"], pch=19, cex=1.5, col="#40B0A6") # teal
points(L.data[L.data$Habitat=="Mediterranean","delta.model"], L.data[L.data$Habitat=="Mediterranean","delta.TPC"], pch=19, cex=1.5, col="#40B0A6") # teal
points(L.data[L.data$Habitat=="Temperate","delta.model"], L.data[L.data$Habitat=="Temperate","delta.TPC"], pch=19, cex=1.5, col="#785EF0") # purple


# POPULATION DYNAMICS
# Mean density vs latitude
Xmin <- 0
Xmax <- 60
Ymin <- -1
Ymax <- 1
#dev.new(width=3, height=3, unit="in")
plot(-100, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Latitude", ylab="Model", cex.axis=2)
abline(0, 0, col="gray", lwd=3, lty="longdash")
#points(seq(2*Xmin,2*Xmax,1), coef(mean.lat)[2]*seq(2*Xmin,2*Xmax,1) + coef(mean.lat)[1], type="l", lwd=3, col="black", lty="longdash")
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
points(pop.data[pop.data$Habitat=="Tropical","Latitude"], pop.data[pop.data$Habitat=="Tropical","delta.active"], pch=19, cex=1.5, col="#FFB000") # orange
points(pop.data[pop.data$Habitat=="Subtropical","Latitude"], pop.data[pop.data$Habitat=="Subtropical","delta.active"], pch=19, cex=1.5, col="#6FD012") # green
points(pop.data[pop.data$Habitat=="Mediterranean","Latitude"], pop.data[pop.data$Habitat=="Mediterranean","delta.active"], pch=19, cex=1.5, col="#6FD012") # green
points(pop.data[pop.data$Habitat=="Temperate","Latitude"], pop.data[pop.data$Habitat=="Temperate","delta.active"], pch=19, cex=1.5, col="#785EF0") # purple

