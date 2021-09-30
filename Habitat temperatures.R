#############################################################################
#### This R script fits sinusoidal functions to habitat temperature data ####
#############################################################################


# Load packages and set working directory
library(tidyr)
library(ggplot2)
library(dplyr)
library(tidyverse)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# USER: enter location, time period, and insect species
location <- "Benin"

# INPUT DATA
# Select a location by removing # in front of name and placing # in front of other locations
data.h <- as.data.frame(read_csv(paste0("Historical climate data ",location,".csv")))
data.f <- as.data.frame(read_csv(paste0("Future climate data ",location,".csv")))


#################################### HISTORICAL CLIMATE #####################################
# Fit sinusoidal function with  annual temperature variation to climate data
fit.h <- nls(T ~ meanT + amplT * sin(2*pi*(day + shiftT)/365), data = data.h,
             start = list(meanT = 300, amplT = 1, shiftT = 30))
summary(fit.h)

# Fit sinusoidal function with  annual and diurnal temperature variation to climate data
fit.h <- nls(T ~ meanT + amplT * sin(2*pi*(day + shiftT)/365) + amplD * sin(2*pi*(day + shiftD)),
             data = data.h, start = list(meanT = 300, amplT = 1, shiftT = 30, amplD = 0.1, shiftD = 0))
summary(fit.h)

# Plot
xmin <- 0
xmax <- nrow(data.h)
ggplot(data.h, aes(x=days, y=T_K)) + 
  geom_point(size=3, color="black") +
  #geom_line(size=0.8) +
  geom_function(fun = function(t){coef(fit.h)[1] + coef(fit.h)[2] * sin(2*pi*(t + coef(fit.h)[3])/360)},
                size=0.8, color="#d1495b") +
  labs(x="Time (days)", y="Mean Temperature (K)") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(coef(fit.h)[1] - abs(coef(fit.h)[2]) - 5, coef(fit.h)[1] + abs(coef(fit.h)[2]) + 5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20))


###################################### FUTURE CLIMATE #######################################
# Fit sinusoidal function with  annual temperature variation to climate data
fit.f <- nls(T ~ meanT + amplT * sin(2*pi*(day + shiftT)/365), data = data.f,
             start = list(meanT = 300, amplT = 1, shiftT = 30))
summary(fit.f)

# Fit sinusoidal function with  annual and diurnal temperature variation to climate data
fit.f <- nls(T ~ meanT + amplT * sin(2*pi*(day + shiftT)/365) + amplD * sin(2*pi*(day + shiftD)),
             data = data.f, start = list(meanT = 300, amplT = 1, shiftT = 30, amplD = 0.1, shiftD = 0))
summary(fit.f)

# Plot
xmin <- 0
xmax <- nrow(data.h)
ggplot(data.h, aes(x=days, y=T_K)) + 
  geom_point(size=3, color="black") +
  #geom_line(size=0.8) +
  geom_function(fun = function(t){coef(fit.h)[1] + coef(fit.h)[2] * sin(2*pi*(t + coef(fit.h)[3])/360)},
                size=0.8, color="#d1495b") +
  labs(x="Time (days)", y="Mean Temperature (K)") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(coef(fit.h)[1] - abs(coef(fit.h)[2]) - 5, coef(fit.h)[1] + abs(coef(fit.h)[2]) + 5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20))




# MEAN TEMPERATURE DURING ACTIVE PERIOD (mean T when T > Tmin)
# input temperature response parameters
TR.data <- as.data.frame(read_csv("Temperature response parameters.csv"))

# Select an insect by removing # in front of name and placing # in front of other species
#sp.data <- subset(TR.data, Species == "Clavigralla shadabi Benin")
#sp.data <- subset(TR.data, Species == "Clavigralla tomentosicollis Benin")
#sp.data <- subset(TR.data, Species == "Clavigralla tomentosicollis Burkina Faso")
#sp.data <- subset(TR.data, Species == "Apolygus lucorum China Dafeng")
#sp.data <- subset(TR.data, Species == "Adelphocoris suturalis China Dafeng")
#sp.data <- subset(TR.data, Species == "Apolygus lucorum China Langfang")
#sp.data <- subset(TR.data, Species == "Adelphocoris suturalis China Xinxiang")
#sp.data <- subset(TR.data, Species == "Macrosiphum euphorbiae Brazil")
#sp.data <- subset(TR.data, Species == "Aulacorthum solani Brazil")
#sp.data <- subset(TR.data, Species == "Uroleucon ambrosiae Brazil")
#sp.data <- subset(TR.data, Species == "Macrolophus pygmaeus on Myzus persicae Greece")
sp.data <- subset(TR.data, Species == "Macrolophus pygmaeus on Trialeurodes vaporariorum Greece")


# Calculate active period in data (T > Tmin)
T.sum <- 0 # sum of temperatures during active period
length <- 0 # length of active period
for(i in 1:nrow(data.f)) {if(data.f[i,"T_K"] > sp.data["Tmin"][[1]])
  {T.sum <- T.sum + data.f[i,"T_K"]
  length <- length + 1}
}
meanT.active <- T.sum/length
meanT.active


# Calculate active period in model (T(t) > Tmin)
yr <- 360 # days in a year (using 360 for simplicity)
# temperature function
T <- function(t) (coef(fit.f)[1] + coef(fit.f)[2] * sin(2*pi*(t + coef(fit.f)[3])/yr))[[1]]
tStart <- 0 # start of active period
tEnd <- 0 # end of active period
T.active <- 0 # sum of temperatures during active period
for(i in 0:yr) {
  if(T(i) > sp.data["Tmin"] & tStart == 0) {tStart <- i}
  if(T(i) > sp.data["Tmin"]) {tEnd <- i
    T.active <- T.active + T(i)}
}
meanT.active <- T.active/(tEnd-tStart)
tStart
tEnd
meanT.active
coef(fit.f)[1][[1]] # meanT during future period


