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


############# NOTE: RE-RUN ENTIRE SCRIPT IF ELAPSED TIME LIMIT ERROR OCCURS ##################


# USER: enter species and location
species <- "Clavigralla tomentosicollis"
location <- "Nigeria"
field_plot <- "A" # for Nigeria, must specify plot "A", "B", or "C" 
                  # densities excluded during dry pods in plot A, B and harmattan in plot C

# USER: include diurnal variation?
daily <- FALSE

# USER: Use left-skewed function for development?
left_skew <- TRUE # if FALSE, development plateaus between Topt and Tmax before going to zero above Tmax

# USER: SET PLOT OPTIONS
xmin <- 0 # start date
xmax <- 730 # end date
ymin <- 0 # min density
ymax <- 200 # max density
temp.min <- 0 # min temperature
temp.max <- 40 # max temperature


# READ IN TEMPERATURE RESPONSE PARAMETERS AND TIME-SERIES DATA
# Temperature response parameters
data <- as.data.frame(read_csv("Temperature response parameters.csv"))
# Temperature parameters
ifelse(daily == TRUE, temp.data <- as.data.frame(read_csv("Temperature parameters.csv")),
       temp.data <- as.data.frame(read_csv("Temperature parameters Tave.csv")))

# Time-series field data
if(str_split(location, boundary("word"), simplify = T)[,1] == "China") { data.density <- read_csv("Population data China.csv") }
if(location == "Nigeria") { data.density <- read_csv("Population data Nigeria.csv") }

# Select time-series data
if(str_split(location, boundary("word"), simplify = T)[,1] == "China") { data.TS <- data.density[data.density$location==str_split(location, boundary("word"), simplify = T)[,2] & data.density$species==species,] }
if(location == "Nigeria") { data.TS <- data.density[data.density$Plot==field_plot,] }


# DATA TRANSFORMATIONS (IF NECESSARY)
# Convert from log to linear scale
if(location == "Nigeria") {
  data.TS$J <- (10^data.TS$J) - 1
  data.TS$J_SE_L <- (10^data.TS$J_SE_L) - 1
  data.TS$J_SE_H <- (10^data.TS$J_SE_H) - 1
  data.TS$A <- (10^data.TS$A) - 1
  data.TS$A_SE_L <- (10^data.TS$A_SE_L) - 1
  data.TS$A_SE_H <- (10^data.TS$A_SE_H) - 1
}


# READ IN TEMPERATURE RESPONSE PARAMETERS, TEMPERATURE PARAMETERS, AND DDE MODEL DYNAMICS
sp.data <- data[data$Species == paste(species,location),]
temp.data <- temp.data[temp.data$Species == paste(species,location),]
if((species == "Clavigralla shadabi" && location == "Nigeria") || (species == "Apolygus lucorum" && str_split(location, boundary("word"), simplify = T)[,1] == "China")) {
  # For species with census data
  if(daily == TRUE) {
    data.model <- as.data.frame(read_csv(paste0("Time series data Census/Historical Time Series ",species," ",location,".csv")))
    data.model.CC <- as.data.frame(read_csv(paste0("Time series data Census/Future Time Series ",species," ",location,".csv"))) }
  if(daily == FALSE && left_skew == TRUE) {
    data.model <- as.data.frame(read_csv(paste0("Time series data Census/Historical Time Series Tave ",species," ",location,".csv")))
    data.model.CC <- as.data.frame(read_csv(paste0("Time series data Census/Future Time Series Tave ",species," ",location,".csv"))) }
  if(daily == FALSE && left_skew == FALSE) {
    data.model <- as.data.frame(read_csv(paste0("Time series data Census/Historical Time Series Tave Dev ",species," ",location,".csv")))
    data.model.CC <- as.data.frame(read_csv(paste0("Time series data Census/Future Time Series Tave Dev ",species," ",location," .csv"))) }
} else{
  # For species without census data
  if(daily == TRUE) {
    data.model <- as.data.frame(read_csv(paste0("Time series data/Historical Time Series ",species," ",location,".csv")))
    data.model.CC <- as.data.frame(read_csv(paste0("Time series data/Future Time Series ",species," ",location,".csv"))) }
  if(daily == FALSE && left_skew == TRUE) {
    data.model <- as.data.frame(read_csv(paste0("Time series data Tave/Historical Time Series ",species," ",location,".csv")))
    data.model.CC <- as.data.frame(read_csv(paste0("Time series data Tave/Future Time Series ",species," ",location,".csv"))) }
  if(daily == FALSE && left_skew == FALSE) {
    data.model <- as.data.frame(read_csv(paste0("Time series data Tave Dev/Historical Time Series ",species," ",location,".csv")))
    data.model.CC <- as.data.frame(read_csv(paste0("Time series data Tave Dev/Future Time Series ",species," ",location," .csv"))) }
}


# PLOT OPTIONS
# for climate change time period
xmin.CC <- xmin
xmax.CC <- xmax
ymin.CC <- ymin
ymax.CC <- ymax
yr <- 365 # days in a year
init_yrs <- 8 # number of years before taking time-series data
TS.length <- xmax - xmin # length of time-series data
end <- nrow(data.model)
end.CC <- nrow(data.model.CC)


# FORMAT MODEL OUTPUT TO ALIGN WITH TIME-SERIES DATA
# Remove all rows before time-series data starts
data.model <- data.model[-c(1:(init_yrs*yr + xmin)), ]

# Remove all rows after xmax days
if(xmax < end) { data.model <- data.model[-c(xmax+1:end), ] }

# climate change period (remove all but last 2 years of data)
if(xmax.CC < end.CC) { data.model.CC <- data.model.CC[-c(1:(end.CC-2*yr + xmin.CC)), ] }

# Re-scale time to start at xmin
# historical period
time.shift <- data.model[[1,1]] - xmin
data.model <- sweep(data.model, 2, c(time.shift,0,0,0,0))
# climate change period
time.shift.CC <- data.model.CC[[1,1]] + xmin.CC
data.model.CC <- sweep(data.model.CC, 2, c(time.shift.CC,0,0,0,0))

# Convert from linear to log scale (if necessary)
# data.model$J <- log(data.model$J + 1, 10)
# data.model$A <- log(data.model$A + 1, 10)
# data.model.CC$J <- log(data.model.CC$J + 1, 10)
# data.model.CC$A <- log(data.model.CC$A + 1, 10)

# Quantify total insects
data.model$I <- data.model$J + data.model$A
data.model.CC$I <- data.model.CC$J + data.model.CC$A


########################################## PLOTS #############################################
# TIME-SERIES DATA
# Juvenile density
plot.J = ggplot(data.TS, aes(x=time, y=J, ymin=J_SE_L, ymax=J_SE_H)) + 
  geom_pointrange(size=0.5, color="#E69F00") + # orange color
  labs(x="", y="") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(ymin, ymax)) +
  #scale_y_log10(limits=c(ymin, ymax)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(color = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=1.50))

# Adult density
plot.A = ggplot(data.TS, aes(x=time, y=A, ymin=A_SE_L, ymax=A_SE_H)) + 
  geom_pointrange(size=0.5, color="#009E73") + # green color
  labs(x="", y="") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(ymin, ymax)) +
  #scale_y_log10(limits=c(ymin, ymax)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(color = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=1.50)) 

# Insect density (Juveniles + Adults)
plot.I = ggplot(data.TS, aes(x=time, y=A, ymin=A_SE_L, ymax=A_SE_H)) + # NOTE: data labelled "A", but are for all insect stages
  geom_pointrange(size=0.5, color="black") +
  labs(x="", y="") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(ymin, ymax)) +
  #scale_y_log10(limits=c(ymin, ymax)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(color = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=1.50)) 


# DDE MODEL DATA
# Historical time period
# Juvenile density
model.J = ggplot(data.model, aes(x=Time, y=J)) + 
  geom_line(size=1.5, color="#0072B2") + # blue color
  labs(x="", y="") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(ymin, ymax)) +
  #scale_y_log10(limits=c(ymin, ymax)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(color = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=1.50))

# Adult density
model.A = ggplot(data.model, aes(x=Time, y=A)) + 
  geom_line(size=1.5, color="#0072B2") + # blue color
  labs(x="", y="") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(ymin, ymax)) +
  #scale_y_log10(limits=c(ymin, ymax)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(color = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=1.50)) 

# Insect density (juveniles + adults)
model.I = ggplot(data.model, aes(x=Time, y=I)) + 
  geom_line(size=1.5, color="black") +
  labs(x="", y="") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(ymin, ymax)) +
  #scale_y_log10(limits=c(ymin, ymax)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(color = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=1.50))

# Future time period
# Juvenile density
model.J.CC = ggplot(data.model.CC, aes(x=Time, y=J)) + 
  geom_line(size=1.5, color="#D55E00") + # red color
  labs(x="", y="") +
  scale_x_continuous(limits=c(xmin.CC, xmax.CC)) +
  scale_y_continuous(limits=c(ymin.CC, ymax.CC)) +
  #scale_y_log10(limits=c(ymin.CC, ymax.CC)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(color = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=1.50))

# Adult density
model.A.CC = ggplot(data.model.CC, aes(x=Time, y=A)) + 
  geom_line(size=1.5, color="#D55E00") + # red color
  labs(x="", y="") +
  scale_x_continuous(limits=c(xmin.CC, xmax.CC)) +
  scale_y_continuous(limits=c(ymin.CC, ymax.CC)) +
  #scale_y_log10(limits=c(ymin.CC, ymax.CC)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(color = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=1.50)) 

# Insect density (juveniles + adults)
model.I.CC = ggplot(data.model.CC, aes(x=Time, y=I)) + 
  geom_line(size=1.5, color="black") +
  labs(x="", y="") +
  scale_x_continuous(limits=c(xmin.CC, xmax.CC)) +
  scale_y_continuous(limits=c(ymin.CC, ymax.CC)) +
  #scale_y_log10(limits=c(ymin.CC, ymax.CC)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(color = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=1.50))


# HABITAT TEMPERATURE
# Historical time period
# data table from Tmin and Tmax functions
temp.fun.h <- data.frame(t=c(xmin:xmax))
temp.fun.h$fun.min <- sapply(temp.fun.h$t, FUN = function(t) { (temp.data$meanT.h - 273.15 + temp.data$delta_mean.h*(t+time.shift))  - (temp.data$amplT.h + temp.data$delta_ampl.h*(t+time.shift)) * cos(2*pi*((t+time.shift) + temp.data$shiftT.h)/yr) - abs(temp.data$amplD.h) })
temp.fun.h$fun.max <- sapply(temp.fun.h$t, FUN = function(t) { (temp.data$meanT.h - 273.15 + temp.data$delta_mean.h*(t+time.shift))  - (temp.data$amplT.h + temp.data$delta_ampl.h*(t+time.shift)) * cos(2*pi*((t+time.shift) + temp.data$shiftT.h)/yr) + abs(temp.data$amplD.h) })

# day at which habitat temperature exceeds Tmin
day1.h <- temp.fun.h[temp.fun.h$fun.min >= sp.data$Tmin, "t"][1]

# plot
plot.temp <- ggplot(temp.fun.h, aes(x=t, y=fun.max)) +
  # Daily average temperature
  geom_function(fun = function(t) (temp.data$meanT.h - 273.15 + temp.data$delta_mean.h*(t+time.shift))  - (temp.data$amplT.h + temp.data$delta_ampl.h*(t+time.shift)) * cos(2*pi*((t+time.shift) + temp.data$shiftT.h)/yr),
                size=1.5, color="#0072B2") + # blue color
  # Daily minimum temperature
  #geom_function(fun = function(t) (temp.data$meanT.h - 273.15 + temp.data$delta_mean.h*(t+time.shift))  - (temp.data$amplT.h + temp.data$delta_ampl.h*(t+time.shift)) * cos(2*pi*((t+time.shift) + temp.data$shiftT.h)/yr) - temp.data$amplD.h,
  #              size=1.5, linetype="dashed", color="#0072B2") +
  # Daily maximum temperature
  #geom_function(fun = function(t) (temp.data$meanT.h - 273.15 + temp.data$delta_mean.h*(t+time.shift))  - (temp.data$amplT.h + temp.data$delta_ampl.h*(t+time.shift)) * cos(2*pi*((t+time.shift) + temp.data$shiftT.h)/yr) + temp.data$amplD.h,
  #              size=1.5, linetype="dashed", color="#0072B2") +
  geom_ribbon(aes(ymin = fun.min, ymax = fun.max), fill = "#0072B2", alpha = 0.2) +
  # Minimum developmental temperature
  geom_function(fun = function(t) (sp.data$Tmin), size=1.5, color="black") +
  labs(x="", y="") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(temp.min, temp.max)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(color = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=1.50))

# Future time period
# data table from Tmin and Tmax functions
temp.fun.f <- data.frame(t=c(xmin:xmax))
temp.fun.f$fun.min <- sapply(temp.fun.f$t, FUN = function(t) { (temp.data$meanT.f - 273.15 + temp.data$delta_mean.f*(t+time.shift))  - (temp.data$amplT.f + temp.data$delta_ampl.f*(t+time.shift)) * cos(2*pi*((t+time.shift) + temp.data$shiftT.f)/yr) - abs(temp.data$amplD.f) })
temp.fun.f$fun.max <- sapply(temp.fun.f$t, FUN = function(t) { (temp.data$meanT.f - 273.15 + temp.data$delta_mean.f*(t+time.shift))  - (temp.data$amplT.f + temp.data$delta_ampl.f*(t+time.shift)) * cos(2*pi*((t+time.shift) + temp.data$shiftT.f)/yr) + abs(temp.data$amplD.f) })

# day at which habitat temperature exceeds Tmin
day1.f <- temp.fun.f[temp.fun.f$fun.min >= sp.data$Tmin, "t"][1]

# plot
plot.temp.CC <- ggplot(temp.fun.f, aes(x=t, y=fun.max)) +
  # Daily average temperature
  geom_function(fun = function(t) (temp.data$meanT.f - 273.15 + temp.data$delta_mean.f*(t+time.shift.CC))  - (temp.data$amplT.f + temp.data$delta_ampl.f*(t+time.shift.CC)) * cos(2*pi*((t+time.shift.CC) + temp.data$shiftT.f)/yr),
                size=1.5, color="#D55E00") + # red color
  # Daily minimum temperature
  #geom_function(fun = function(t) (temp.data$meanT.f - 273.15 + temp.data$delta_mean.f*(t+time.shift.CC))  - (temp.data$amplT.f + temp.data$delta_ampl.f*(t+time.shift.CC)) * cos(2*pi*((t+time.shift.CC) + temp.data$shiftT.f)/yr) - temp.data$amplD.f,
  #              size=1.5, linetype="dashed", color="#D55E00") +
  # Daily maximum temperature
  #geom_function(fun = function(t) (temp.data$meanT.f - 273.15 + temp.data$delta_mean.f*(t+time.shift.CC))  - (temp.data$amplT.f + temp.data$delta_ampl.f*(t+time.shift.CC)) * cos(2*pi*((t+time.shift.CC) + temp.data$shiftT.f)/yr) + temp.data$amplD.f,
  #              size=1.5, linetype="dashed", color="#D55E00") +
  geom_ribbon(aes(ymin = fun.min, ymax = fun.max), fill = "#D55E00", alpha = 0.2) +
  # Minimum developmental temperature
  geom_function(fun = function(t) (sp.data$Tmin), size=1.5, color="black") +
  labs(x="", y="") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(temp.min, temp.max)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(color = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=1.50))



#################################### DRAW FINAL PLOTS #######################################
# COMPILE PLOTS
# Temperature plots
plot.climate <- ggdraw()  +
  draw_plot(plot.temp, x = 0, y = 0, width = 1, height = 0.3) +
  draw_plot(plot.temp.CC, x = 0, y = 0, width = 1, height = 0.3)
#plot.climate

# Historical time period
if(location == "Nigeria") {
  plot <- ggdraw()  +
     draw_plot(plot.temp, x = 0, y = 0, width = 1, height = 0.3) +
     #draw_plot(plot.J, x = 0, y = 0.3, width = 1, height = 0.7) +
     draw_plot(plot.A, x = 0, y = 0.3, width = 1, height = 0.7) +
     #draw_plot(model.J, x = 0, y = 0.3, width = 1, height = 0.7) +
     draw_plot(model.A, x = 0, y = 0.3, width = 1, height = 0.7)
}
if(str_split(location, boundary("word"), simplify = T)[,1] == "China") { 
  plot <- ggdraw()  +
   draw_plot(plot.temp, x = 0, y = 0, width = 1, height = 0.3) +
   draw_plot(plot.I, x = 0, y = 0.3, width = 1, height = 0.7) +
   #draw_plot(model.J, x = 0, y = 0.3, width = 1, height = 0.7) +
   #draw_plot(model.A, x = 0, y = 0.3, width = 1, height = 0.7) +
   draw_plot(model.I, x = 0, y = 0.3, width = 1, height = 0.7)
}

# Future time period
if(location == "Nigeria") {
  plot.CC <- ggdraw()  +
    draw_plot(plot.temp.CC, x = 0, y = 0, width = 1, height = 0.3) +
    draw_plot(model.J.CC, x = 0, y = 0.3, width = 1, height = 0.7) +
    draw_plot(model.A.CC, x = 0, y = 0.3, width = 1, height = 0.7)
}
if(str_split(location, boundary("word"), simplify = T)[,1] == "China") { 
  plot.CC <- ggdraw()  +
    draw_plot(plot.temp.CC, x = 0, y = 0, width = 1, height = 0.3) +
    draw_plot(model.J.CC, x = 0, y = 0.3, width = 1, height = 0.7) +
    draw_plot(model.A.CC, x = 0, y = 0.3, width = 1, height = 0.7) +
    draw_plot(model.I.CC, x = 0, y = 0.3, width = 1, height = 0.7)
}

# Compare historical and future time periods
if(location == "Nigeria") {
  plot.compare <- ggdraw()  +
     draw_plot(plot.temp, x = 0, y = 0, width = 1, height = 0.3) +
     draw_plot(plot.temp.CC, x = 0, y = 0, width = 1, height = 0.3) +
     draw_plot(model.J, x = 0, y = 0.3, width = 1, height = 0.7) +
     draw_plot(model.A, x = 0, y = 0.3, width = 1, height = 0.7) +
     draw_plot(model.J.CC, x = 0, y = 0.3, width = 1, height = 0.7) +
     draw_plot(model.A.CC, x = 0, y = 0.3, width = 1, height = 0.7)
}
if(str_split(location, boundary("word"), simplify = T)[,1] == "China") { 
  plot.compare <- ggdraw()  +
    draw_plot(plot.temp, x = 0, y = 0, width = 1, height = 0.3) +
    draw_plot(plot.temp.CC, x = 0, y = 0, width = 1, height = 0.3) +
    draw_plot(model.J, x = 0, y = 0.3, width = 1, height = 0.7) +
    draw_plot(model.A, x = 0, y = 0.3, width = 1, height = 0.7) +
    #draw_plot(model.I, x = 0, y = 0.3, width = 1, height = 0.7) +
    draw_plot(model.J.CC, x = 0, y = 0.3, width = 1, height = 0.7) +
    draw_plot(model.A.CC, x = 0, y = 0.3, width = 1, height = 0.7)
    #draw_plot(model.I.CC, x = 0, y = 0.3, width = 1, height = 0.7)
}

# PRINT PLOTS
# View plot in RStudio
plot
#plot.CC
plot.compare


# OUTPUT PLOTS
#dev.new()
# Time-series plots
# Nigeria
# ggdraw() +
#   draw_plot(plot.A, width = 1, height = 0.4) +
#   draw_plot(model.A, width = 1, height = 0.4)
# China
# ggdraw() +
#   draw_plot(plot.I, width = 1, height = 0.4) +
#   draw_plot(model.A, width = 1, height = 0.4) +
#   draw_plot(model.I, width = 1, height = 0.4)

# Climate change plots
# ggdraw()  +
#   draw_plot(model.J, width = 1, height = 0.4) +
#   draw_plot(model.J.CC, width = 1, height = 0.4)
# ggdraw()  +
#    draw_plot(model.A, width = 1, height = 0.45) +
#    draw_plot(model.A.CC, width = 1, height = 0.45)

# Temperature plots
# ggdraw()  +
#    draw_plot(plot.temp, x = 0, y = 0, width = 1, height = 0.3) +
#    draw_plot(plot.temp.CC, x = 0, y = 0, width = 1, height = 0.3)

