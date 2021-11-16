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
r.data <- as.data.frame(read_csv("Model results.csv"))


# STATISTICS
# proportional change in r vs latitude
TPC <- lm(delta.TPC ~ Latitude, data=r.data) # non-significant
summary(TPC)

# proportional change in r vs latitude
model <- lm(delta.model ~ Latitude, data=r.data) # non-significant
summary(model)

# proportional change in r (model vs TPCs)
d <- lm(delta.model ~ delta.TPC, data=r.data) # non-significant
summary(d)

# PLOTS
# proportional change in r (TPC) vs latitude
Xmin <- 0
Xmax <- 40
Ymin <- -1.5
Ymax <- 1.5
plot(r.data[r.data$Habitat=="Tropical","Latitude"], r.data[r.data$Habitat=="Tropical","delta.TPC"], pch=21, col="red", bg="red",
     xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Latitude", ylab="change in r")
points(r.data[r.data$Habitat=="Mediterranean","Latitude"], r.data[r.data$Habitat=="Mediterranean","delta.TPC"], pch=21, col="purple", bg="purple")
points(r.data[r.data$Habitat=="Temperate","Latitude"], r.data[r.data$Habitat=="Temperate","delta.TPC"], pch=21, col="blue", bg="blue")
abline(0, 0, col="black", lty="longdash")

# proportional change in r (model) vs latitude
Xmin <- 0
Xmax <- 40
Ymin <- -1.5
Ymax <- 1.5
plot(r.data[r.data$Habitat=="Tropical","Latitude"], r.data[r.data$Habitat=="Tropical","delta.model"], pch=21, col="red", bg="red",
     xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Latitude", ylab="change in r")
points(r.data[r.data$Habitat=="Mediterranean","Latitude"], r.data[r.data$Habitat=="Mediterranean","delta.model"], pch=21, col="purple", bg="purple")
points(r.data[r.data$Habitat=="Temperate","Latitude"], r.data[r.data$Habitat=="Temperate","delta.model"], pch=21, col="blue", bg="blue")
abline(0, 0, col="black", lty="longdash")

# proportional change in r (model vs TPCs)
Xmin <- -1.5
Xmax <- 1.5
Ymin <- -1.5
Ymax <- 1.5
plot(r.data[r.data$Habitat=="Tropical","delta.TPC"], r.data[r.data$Habitat=="Tropical","delta.model"], pch=21, col="red", bg="red",
     xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="TPC", ylab="Model")
points(r.data[r.data$Habitat=="Mediterranean","delta.TPC"], r.data[r.data$Habitat=="Mediterranean","delta.model"], pch=21, col="purple", bg="purple")
points(r.data[r.data$Habitat=="Temperate","delta.TPC"], r.data[r.data$Habitat=="Temperate","delta.model"], pch=21, col="blue", bg="blue")
#points(seq(Xmin,Xmax,1), coef(d)[1]*seq(Xmin,Xmax,1), type="l", col="blue")
points(seq(Xmin,Xmax,1), coef(d)[2]*seq(Xmin,Xmax,1)+coef(d)[1], type="l", col="black", lty="longdash")
#abline(0, 1, col="black", lty="longdash")






# Temperature response parameters
data <- as.data.frame(read_csv("Temperature response parameters.csv"))



# EXTINCTION THRESHOLD (increase in meanT that leads to extinction in the model)
# Versus intrinsic growth rate (r)
# all parameters: rTopt significant only in full model
r <- lm(ext_meanT ~ rTopt + Toptr + rTmax, data=data)
summary (r)
r.Topt <- lm(ext_meanT ~ rTopt, data=data)
summary (r.Topt)

# thermal safety margin: non significant
tsm <- lm(ext_meanT ~ TSM, data=data)
summary (tsm)

# warming tolerance: significant
wt <- lm(ext_meanT ~ WT, data=data)
#wt <- lm(ext_meanT ~ 0 + WT, data=data) # remove intercept
summary (wt)

# warming tolerance (active period): significant
wt.active <- lm(ext_meanT ~ WT_active, data=data)
#wt.active <- lm(ext_meanT ~ 0 + WT_active, data=data) # remove intercept
summary (wt.active)


# Versus net reproductive rate (R0)
# all parameters: sR0 significant
R0 <- lm(ext_meanT ~ R0Topt + ToptR0 + sR0, data=data)
summary (R0)

# R0: Tmax (temperature at which R0 = 1): marginally significant
R0.Tmax <- lm(ext_meanT ~ R0Tmax, data=data)
summary (R0.Tmax)

# R0 thermal safety margin: significant
R0.tsm <- lm(ext_meanT ~ R0TSM, data=data)
R0.tsm <- lm(ext_meanT ~ 0 + R0TSM, data=data) # remove intercept
summary (R0.tsm)

# R0 warming tolerance: significant
R0.wt <- lm(ext_meanT ~ R0WT, data=data)
R0.wt <- lm(ext_meanT ~ 0 + R0WT, data=data) # remove intercept
summary (R0.wt)

# warming tolerance (active period): significant
R0.wt.active <- lm(ext_meanT ~ R0WT_active, data=data)
R0.wt.active <- lm(ext_meanT ~ 0 + R0WT_active, data=data) # remove intercept
summary (R0.wt.active)


# Extinction risk across latitude
lat <- lm(ext_meanT ~ Latitude, data=data)
#lat <- lm(ext_meanT ~ 0 + Latitude, data=data) # remove intercept
summary (lat)


# PLOTS
# Metrics based on intrinsic growth rate (r)
# thermal safety margin: non significant
Xmin <- 0
Xmax <- 12
Ymin <- 0
Ymax <- 15
plot(data[data$Habitat=="Tropical","TSM"], data[data$Habitat=="Tropical","ext_meanT"], pch=21, col="red", bg="red",
     xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Thermal safety margin", ylab="Extinction threshold (mean T)")
points(data[data$Habitat=="Mediterranean","TSM"], data[data$Habitat=="Mediterranean","ext_meanT"], pch=21, col="purple", bg="purple")
points(data[data$Habitat=="Temperate","TSM"], data[data$Habitat=="Temperate","ext_meanT"], pch=21, col="blue", bg="blue")
#points(seq(Xmin,Xmax,1), coef(tsm)[1]*seq(Xmin,Xmax,1), type="l", col="blue")
points(seq(Xmin,Xmax,1), coef(tsm)[2]*seq(Xmin,Xmax,1)+coef(tsm)[1], type="l", col="black", lty="longdash")
#abline(0, 1, col="black", lty="longdash")

# warming tolerance: significant (slope = 0.81)
Xmin <- 0
Xmax <- 20
Ymin <- 0
Ymax <- 15
plot(data[data$Habitat=="Tropical","WT"], data[data$Habitat=="Tropical","ext_meanT"], pch=21, col="red", bg="red",
     xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Warming tolerance", ylab="Extinction threshold (mean T)")
points(data[data$Habitat=="Mediterranean","WT"], data[data$Habitat=="Mediterranean","ext_meanT"], pch=21, col="purple", bg="purple")
points(data[data$Habitat=="Temperate","WT"], data[data$Habitat=="Temperate","ext_meanT"], pch=21, col="blue", bg="blue")
#points(seq(Xmin,Xmax,1), coef(wt)[1]*seq(Xmin,Xmax,1), type="l", col="blue")
points(seq(Xmin,Xmax,1), coef(wt)[2]*seq(Xmin,Xmax,1)+coef(wt)[1], type="l", col="black")
abline(0, 1, col="black", lty="longdash")

# warming tolerance (active period): significant (slope = 1.01)
Xmin <- 0
Xmax <- 20
Ymin <- 0
Ymax <- 15
plot(data[data$Habitat=="Tropical","WT_active"], data[data$Habitat=="Tropical","ext_meanT"], pch=21, col="red", bg="red",
     xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), xlab="Warming tolerance (active period)", ylab="Extinction threshold (mean T)")
points(data[data$Habitat=="Mediterranean","WT_active"], data[data$Habitat=="Mediterranean","ext_meanT"], pch=21, col="purple", bg="purple")
points(data[data$Habitat=="Temperate","WT_active"], data[data$Habitat=="Temperate","ext_meanT"], pch=21, col="blue", bg="blue")
#points(seq(Xmin,Xmax,1), coef(wt.active)[1]*seq(Xmin,Xmax,1), type="l", col="blue")
points(seq(Xmin,Xmax,1), coef(wt.active)[2]*seq(Xmin,Xmax,1)+coef(wt.active)[1], type="l", col="black")
abline(0, 1, col="black", lty="longdash")



# Net reproductive rate (R0)
# R0 thermal safety margin: significant (slope = 0.46)
Xmin <- 10
Xmax <- 30
Ymin <- 0
Ymax <- 15
plot(data$R0TSM, data$ext_meanT, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), main=expression("extinction threshold (meanT) vs. R0 TSM"))
points(seq(Xmin,Xmax,1), coef(R0.tsm)[1]*seq(Xmin,Xmax,1), type="l", col="blue")
abline(0, 1)

# warming tolerance: significant (slope = 0.38)
Xmin <- 0
Xmax <- 40
Ymin <- 0
Ymax <- 15
plot(data$R0WT, data$ext_meanT, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), main=expression("extinction threshold (meanT) vs. R0 WT"))
points(seq(Xmin,Xmax,1), coef(R0.wt)[1]*seq(Xmin,Xmax,1), type="l", col="blue")
abline(0, 1)

# warming tolerance (active period): significant (slope = 0.5)
Xmin <- 0
Xmax <- 30
Ymin <- 0
Ymax <- 15
plot(data$R0WT_active, data$ext_meanT, xlim=c(Xmin,Xmax), ylim=c(Ymin,Ymax), main=expression("extinction threshold (meanT) vs. R0 WT (active)"))
points(seq(Xmin,Xmax,1), coef(R0.wt.active)[1]*seq(Xmin,Xmax,1), type="l", col="blue")
abline(0, 1)



# Extinction risk across latitude
Xmin <- 0
Xmax <- 45
Ymin <- 0
Ymax <- 15
lat.plot <- ggplot(data, aes(Latitude, ext_meanT)) +
  geom_point(size=3, aes(color=Habitat)) +
  scale_colour_manual(values=c("purple", "blue", "red")) +
  geom_function(fun = function(t) (coef(lat)[2]*t+coef(lat)[1]),
                size=0.8, linetype="longdash", color="black") +
  labs(x="", y="") + #labs(x="Absolute latitude", y="Extinction threshold (mean T)") +
  scale_x_continuous(limits=c(Xmin, Xmax)) +
  scale_y_continuous(limits=c(Ymin, Ymax)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20))
ggdraw()  + draw_plot(lat.plot, x = 0, y = 0, width = 1, height = 0.4)



