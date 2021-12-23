########################################################################
#### This R script fits model outputs to empirical time-series data ####
########################################################################


# Load packages and set working directory
library(tidyr)
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyverse)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# USER: enter species, location, if egg stage is modeled, and if diurnal variation is included
species <- "Apolygus lucorum"
location <- "China Dafeng"
egg <- FALSE
daily <- FALSE

# READ IN TEMPERATURE RESPONSE PARAMETERS AND TIME-SERIES DATA
# Temperature response parameters
data <- as.data.frame(read_csv("Temperature response parameters.csv"))
# Temperature parameters
temp.data <- as.data.frame(read_csv("Temperature parameters.csv"))
diurnal.h <- temp.data[temp.data$Species == paste(species,location),"amplD.h"]
diurnal.f <- temp.data[temp.data$Species == paste(species,location),"amplD.f"]
# diurnal variation
if(daily == FALSE) { temp.data <- as.data.frame(read_csv("Temperature parameters Tmax.csv")) }

# Time-series field data
# Select an file by removing # in front of name and placing # in front of other files
#data.density <- read_csv("Population data Nigeria (egg).csv")
data.density <- read_csv("Population data China.csv")

# USER: Select time-series data
#data.TS <- data.density[data.density$Plot=="C",] # select plot A, B, or C
data.TS <- data.density[data.density$location=="Dafeng" & data.density$species=="Apolygus_lucorum",]
#data.TS <- data.density[data.density$location=="Dafeng" & data.density$species=="Adelphocoris_suturalis",]
#data.TS <- data.density[data.density$location=="Langfang" & data.density$species=="Apolygus_lucorum",]
#data.TS <- data.density[data.density$location=="Xinxiang" & data.density$species=="Adelphocoris_suturalis",]

# Data transformation (if needed)
# Convert from log to linear scale
if(location == "Nigeria") {
if(egg == TRUE) {
  data.TS$E <- (10^data.TS$E) - 1
  data.TS$E_SE_L <- (10^data.TS$E_SE_L) - 1
  data.TS$E_SE_H <- (10^data.TS$E_SE_H) - 1
  }
  data.TS$J <- (10^data.TS$J) - 1
  data.TS$J_SE_L <- (10^data.TS$J_SE_L) - 1
  data.TS$J_SE_H <- (10^data.TS$J_SE_H) - 1
  data.TS$A <- (10^data.TS$A) - 1
  data.TS$A_SE_L <- (10^data.TS$A_SE_L) - 1
  data.TS$A_SE_H <- (10^data.TS$A_SE_H) - 1
}

# READ IN TEMPERATURE RESPONSE PARAMETERS, DDE MODEL DYNAMICS, AND TEMPERATURE PARAMETERS
sp.data <- data[data$Species == paste(species,location),]
data.model <- as.data.frame(read_csv(paste0("Time series data/Historical Time Series ",species," ",location,".csv")))
data.model.CC <- as.data.frame(read_csv(paste0("Time series data/Future Time Series ",species," ",location,".csv")))
temp.data <- temp.data[temp.data$Species == paste(species,location),]


# USER: SET PLOT OPTIONS
# Default: plot last 2 year of model
# for historical time period
xmin <- 0
xmax <- 730
ymin <- 0
ymax <- 100
# for climate change time period
xmin.CC <- xmin
xmax.CC <- xmax
ymin.CC <- ymin
ymax.CC <- ymax
# for temperature function
temp.min <- 275 #290
temp.max <- 310
yr <- 365 # days in a year
init_yrs <- 1 # number of years to initiate the model (from Python DDE model)
TS.length <- xmax - xmin # length of time-series data
end <- nrow(data.model)
end.CC <- nrow(data.model.CC)

# FORMAT MODEL OUTPUT TO ALIGN WITH TIME-SERIES DATA
# Remove all rows before time-series data starts
data.model <- data.model[c(-1:-(init_yrs*yr + xmin)), ]

# Remove all rows after xmax days
data.model <- data.model[c(-xmax+1:-end), ]

# climate change period (remove all but last 2 years of data)
data.model.CC <- data.model.CC[c(-1:-(end.CC-1*yr + xmin.CC)), ]

# Re-scale time to start at xmin
# historical period
time.shift <- data.model[[1,1]] - xmin
ifelse(egg == FALSE, data.model <- sweep(data.model, 2, c(time.shift,0,0,0,0)), data.model <- sweep(data.model, 2, c(time.shift,0,0,0,0,0,0,0,0)))
# climate change period
time.shift.CC <- data.model.CC[[1,1]] + xmin.CC
ifelse(egg == FALSE, data.model.CC <- sweep(data.model.CC, 2, c(time.shift.CC,0,0,0,0)), data.model.CC <- sweep(data.model.CC, 2, c(time.shift.CC,0,0,0,0,0,0,0,0)))

# Data transformation (if needed)
# Convert from linear to log scale
# if(egg == TRUE) { data.model$E <- log(data.model$E + 1, 10) }
# data.model$J <- log(data.model$J + 1, 10)
# data.model$A <- log(data.model$A + 1, 10)
# if(egg == TRUE) { data.model.CC$E <- log(data.model.CC$E + 1, 10) }
# data.model.CC$J <- log(data.model.CC$J + 1, 10)
# data.model.CC$A <- log(data.model.CC$A + 1, 10)



# PLOT TIME-SERIES DATA
# Egg density
if (egg == TRUE) {
plot.E = ggplot(data.TS, aes(x=time, y=E, ymin=E_SE_L, ymax=E_SE_H)) + 
  geom_pointrange(size=0.5, color="black") +
  labs(x="Time", y="Density") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(ymin, ymax)) +
  #scale_y_log10(limits=c(ymin, ymax)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20))}

# Juvenile density
plot.J = ggplot(data.TS, aes(x=time, y=J, ymin=J_SE_L, ymax=J_SE_H)) + 
  geom_pointrange(size=0.5, color="#d1495b") + # red color
  labs(x="Time", y="Density") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(ymin, ymax)) +
  #scale_y_log10(limits=c(ymin, ymax)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20))

# Adult density
plot.A = ggplot(data.TS, aes(x=time, y=A, ymin=A_SE_L, ymax=A_SE_H)) + 
  geom_pointrange(size=0.5, color="black") + # blue color
  labs(x="Time", y="Density") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(ymin, ymax)) +
  #scale_y_log10(limits=c(ymin, ymax)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20)) 

# Insect density (Juveniles + Adults)
plot.I = ggplot(data.TS, aes(x=time, y=A, ymin=A_SE_L, ymax=A_SE_H)) + # NOTE: data under column "A", but are for all insect stages
  geom_pointrange(size=0.5, color="black") +
  labs(x="Time", y="Density") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(ymin, ymax)) +
  #scale_y_log10(limits=c(ymin, ymax)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20)) 



# PLOT MODEL DATA
# Historical time period
# Egg density
if (egg == TRUE) {
model.E = ggplot(data.model, aes(x=Time, y=E)) + 
  geom_line(size=0.8, color="#30638e") +
  labs(x="Time", y="Density") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(ymin, ymax)) +
  #scale_y_log10(limits=c(ymin, ymax)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20))}

# Juvenile density
model.J = ggplot(data.model, aes(x=Time, y=J)) + 
  geom_line(size=0.8, color="#d1495b") + # red color
  labs(x="Time", y="Density") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(ymin, ymax)) +
  #scale_y_log10(limits=c(ymin, ymax)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20))

# Adult density
model.A = ggplot(data.model, aes(x=Time, y=A)) + 
  geom_line(size=0.8, color="#30638e") + # blue color
  labs(x="Time", y="Density") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(ymin, ymax)) +
  #scale_y_log10(limits=c(ymin, ymax)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20)) 

# Insect density (juveniles + adults)
model.I = ggplot(data.model, aes(x=Time, y=J+A)) + 
  geom_line(size=0.8, color="black") +
  labs(x="Time", y="Density") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(ymin, ymax)) +
  #scale_y_log10(limits=c(ymin, ymax)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20)) 


# Climate change time period
# Egg density
if (egg == TRUE) {
model.E.CC = ggplot(data.model.CC, aes(x=Time, y=E)) + 
  geom_line(size=0.8, color="black", linetype="longdash") +
  labs(x="Time", y="Density") +
  scale_x_continuous(limits=c(xmin.CC, xmax.CC)) +
  scale_y_continuous(limits=c(ymin.CC, ymax.CC)) +
  #scale_y_log10(limits=c(ymin.CC, ymax.CC)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20))}

# Juvenile density
model.J.CC = ggplot(data.model.CC, aes(x=Time, y=J)) + 
  geom_line(size=0.8, color="#d1495b", linetype="longdash") + # red color
  labs(x="Time", y="Density") +
  scale_x_continuous(limits=c(xmin.CC, xmax.CC)) +
  scale_y_continuous(limits=c(ymin.CC, ymax.CC)) +
  #scale_y_log10(limits=c(ymin.CC, ymax.CC)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20))

# Adult density
model.A.CC = ggplot(data.model.CC, aes(x=Time, y=A)) + 
  #geom_line(size=0.8, color="#30638e", linetype="longdash") + # blue color
  geom_line(size=0.8, color="red") +
  labs(x="Time", y="Density") +
  scale_x_continuous(limits=c(xmin.CC, xmax.CC)) +
  scale_y_continuous(limits=c(ymin.CC, ymax.CC)) +
  #scale_y_log10(limits=c(ymin.CC, ymax.CC)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20)) 

# Insect density (juveniles + adults)
model.I.CC = ggplot(data.model.CC, aes(x=Time, y=J+A)) + 
  geom_line(size=0.8, color="black", linetype="longdash") +
  labs(x="Time", y="Density") +
  scale_x_continuous(limits=c(xmin.CC, xmax.CC)) +
  scale_y_continuous(limits=c(ymin.CC, ymax.CC)) +
  #scale_y_log10(limits=c(ymin.CC, ymax.CC)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20))



# PLOT HABITAT TEMPERATURE FUNCTION
# Historical time period
# data table from Tmin and Tmax functions
temp.fun.h <- data.frame(t=c(xmin:xmax))
temp.fun.h <- data.frame(t=c(xmin:xmax),
                       fun.min = sapply(temp.fun.h$t, FUN = function(t) (temp.data$meanT.h + temp.data$delta_mean.h*(t+time.shift))  - (temp.data$amplT.h + temp.data$delta_ampl.h*(t+time.shift)) * cos(2*pi*((t+time.shift) + temp.data$shiftT.h)/yr) - abs(temp.data$amplD.h)),
                       fun.max = sapply(temp.fun.h$t, FUN = function(t) (temp.data$meanT.h + temp.data$delta_mean.h*(t+time.shift))  - (temp.data$amplT.h + temp.data$delta_ampl.h*(t+time.shift)) * cos(2*pi*((t+time.shift) + temp.data$shiftT.h)/yr) + abs(temp.data$amplD.h)))

# day at which habitat temperature exceeds Tmin
day1.h <- temp.fun.h[temp.fun.h$fun.min >= sp.data$Tmin + abs(diurnal.h), "t"][1]

# plot
plot.temp <- ggplot(temp.fun.h, aes(x=t, y=fun.max)) +
  # Daily average temperature
  geom_function(fun = function(t) (temp.data$meanT.h + temp.data$delta_mean.h*(t+time.shift))  - (temp.data$amplT.h + temp.data$delta_ampl.h*(t+time.shift)) * cos(2*pi*((t+time.shift) + temp.data$shiftT.h)/yr),
                size=0.8, color="blue") +
  # Daily minimum temperature
  #geom_function(fun = function(t) (temp.data$meanT.h + temp.data$delta_mean.h*(t+time.shift))  - (temp.data$amplT.h + temp.data$delta_ampl.h*(t+time.shift)) * cos(2*pi*((t+time.shift) + temp.data$shiftT.h)/yr) - temp.data$amplD.h,
  #              size=0.8, linetype="longdash", color="blue") +
  # Daily maximum temperature
  #geom_function(fun = function(t) (temp.data$meanT.h + temp.data$delta_mean.h*(t+time.shift))  - (temp.data$amplT.h + temp.data$delta_ampl.h*(t+time.shift)) * cos(2*pi*((t+time.shift) + temp.data$shiftT.h)/yr) + temp.data$amplD.h,
  #              size=0.8, linetype="longdash", color="blue") +
  geom_ribbon(aes(ymin = fun.min, ymax = fun.max), fill = "blue", alpha = 0.2) +
  # Minimum developmental temperature
  geom_function(fun = function(t) (sp.data$Tmin + abs(diurnal.h)), size=0.8, color="black") +
  labs(x="Time", y="T(K)") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(temp.min, temp.max)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20))


# Climate change time period
# data table from Tmin and Tmax functions
temp.fun.f <- data.frame(t=c(xmin.CC:xmax.CC))
temp.fun.f <- data.frame(t=c(xmin.CC:xmax.CC),
                       fun.min = sapply(temp.fun.f$t, FUN = function(t) (temp.data$meanT.f + temp.data$delta_mean.f*(t+time.shift.CC))  - (temp.data$amplT.f + temp.data$delta_ampl.f*(t+time.shift.CC)) * cos(2*pi*((t+time.shift.CC) + temp.data$shiftT.f)/yr) - abs(temp.data$amplD.f)),
                       fun.max = sapply(temp.fun.f$t, FUN = function(t) (temp.data$meanT.f + temp.data$delta_mean.f*(t+time.shift.CC))  - (temp.data$amplT.f + temp.data$delta_ampl.f*(t+time.shift.CC)) * cos(2*pi*((t+time.shift.CC) + temp.data$shiftT.f)/yr) + abs(temp.data$amplD.f)))

# day at which habitat temperature exceeds Tmin
day1.f <- temp.fun.f[temp.fun.f$fun.min >= sp.data$Tmin + abs(diurnal.f), "t"][1]

# plot
plot.temp.CC <- ggplot(temp.fun.f, aes(x=t, y=fun.max)) +
  # Daily average temperature
  geom_function(fun = function(t) (temp.data$meanT.f + temp.data$delta_mean.f*(t+time.shift.CC))  - (temp.data$amplT.f + temp.data$delta_ampl.f*(t+time.shift.CC)) * cos(2*pi*((t+time.shift.CC) + temp.data$shiftT.f)/yr),
                size=0.8, color="red") +
  # Daily minimum temperature
  #geom_function(fun = function(t) (temp.data$meanT.f + temp.data$delta_mean.f*(t+time.shift.CC))  - (temp.data$amplT.f + temp.data$delta_ampl.f*(t+time.shift.CC)) * cos(2*pi*((t+time.shift.CC) + temp.data$shiftT.f)/yr) - temp.data$amplD.f,
  #              size=0.8, linetype="longdash", color="red") +
  # Daily maximum temperature
  #geom_function(fun = function(t) (temp.data$meanT.f + temp.data$delta_mean.f*(t+time.shift.CC))  - (temp.data$amplT.f + temp.data$delta_ampl.f*(t+time.shift.CC)) * cos(2*pi*((t+time.shift.CC) + temp.data$shiftT.f)/yr) + temp.data$amplD.f,
  #              size=0.8, linetype="longdash", color="red") +
  geom_ribbon(aes(ymin = fun.min, ymax = fun.max), fill = "red", alpha = 0.2) +
  # Minimum developmental temperature
  geom_function(fun = function(t) (sp.data$Tmin + abs(diurnal.f)), size=0.8, color="black") +
  labs(x="Time", y="T(K)") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(temp.min, temp.max)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20))



# DRAW FINAL PLOTS
# Temperature plots
plot.climate <- ggdraw()  +
  draw_plot(plot.temp, x = 0, y = 0, width = 1, height = 0.3) +
  draw_plot(plot.temp.CC, x = 0, y = 0, width = 1, height = 0.3)
#plot.climate

# Historical time period
# by life stages
ifelse(egg == TRUE,
plot <- ggdraw()  +
  draw_plot(plot.temp, x = 0, y = 0, width = 1, height = 0.3) +
  draw_plot(plot.E, x = 0, y = 0.3, width = 1, height = 0.7) +
  draw_plot(plot.J, x = 0, y = 0.3, width = 1, height = 0.7) +
  draw_plot(plot.A, x = 0, y = 0.3, width = 1, height = 0.7) +
  draw_plot(model.E, x = 0, y = 0.3, width = 1, height = 0.7) +
  draw_plot(model.J, x = 0, y = 0.3, width = 1, height = 0.7) +
  draw_plot(model.A, x = 0, y = 0.3, width = 1, height = 0.7),
plot <- ggdraw()  +
   draw_plot(plot.temp, x = 0, y = 0, width = 1, height = 0.3) +
   #draw_plot(plot.J, x = 0, y = 0.3, width = 1, height = 0.7) +
   draw_plot(plot.A, x = 0, y = 0.3, width = 1, height = 0.7) +
   #draw_plot(model.J, x = 0, y = 0.3, width = 1, height = 0.7) +
   draw_plot(model.A, x = 0, y = 0.3, width = 1, height = 0.7))
#plot

# total insects
(plot <- ggdraw()  +
 draw_plot(plot.temp, x = 0, y = 0, width = 1, height = 0.3) +
 draw_plot(plot.I, x = 0, y = 0.3, width = 1, height = 0.7) +
 draw_plot(model.J, x = 0, y = 0.3, width = 1, height = 0.7) +
 draw_plot(model.A, x = 0, y = 0.3, width = 1, height = 0.7) +
 draw_plot(model.I, x = 0, y = 0.3, width = 1, height = 0.7))
plot

day1.h
day1.f

# Climate change time period
# by life stages
#plot.CC <- ggdraw()  +
#  draw_plot(plot.temp.CC, x = 0, y = 0, width = 1, height = 0.3) +
#  draw_plot(model.J.CC, x = 0, y = 0.3, width = 1, height = 0.7) +
#  draw_plot(model.A.CC, x = 0, y = 0.3, width = 1, height = 0.7)
#plot.CC

# total insects
# plot.CC <- ggdraw()  +
#  draw_plot(plot.temp.CC, x = 0, y = 0, width = 1, height = 0.3) +
#  draw_plot(model.J.CC, x = 0, y = 0.3, width = 1, height = 0.7) +
#  draw_plot(model.A.CC, x = 0, y = 0.3, width = 1, height = 0.7) +
#  draw_plot(model.I.CC, x = 0, y = 0.3, width = 1, height = 0.7)
#plot.CC

# Compare historical and climate change periods
# by life stages
plot.compare <- ggdraw()  +
   draw_plot(plot.temp, x = 0, y = 0, width = 1, height = 0.3) +
   draw_plot(plot.temp.CC, x = 0, y = 0, width = 1, height = 0.3) +
#   draw_plot(model.J, x = 0, y = 0.3, width = 1, height = 0.7) +
   draw_plot(model.A, x = 0, y = 0.3, width = 1, height = 0.7) +
#   draw_plot(model.J.CC, x = 0, y = 0.3, width = 1, height = 0.7) +
   draw_plot(model.A.CC, x = 0, y = 0.3, width = 1, height = 0.7)
#plot.compare
#plot.compare

# total insects
# plot.compare <- ggdraw()  +
#  draw_plot(plot.temp, x = 0, y = 0, width = 1, height = 0.3) +
#  draw_plot(plot.temp.CC, x = 0, y = 0, width = 1, height = 0.3) +
#  draw_plot(model.I, x = 0, y = 0.3, width = 1, height = 0.7) +
#  draw_plot(model.I.CC, x = 0, y = 0.3, width = 1, height = 0.7)
# plot.compare
# plot.compare

