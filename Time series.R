########################################################################
#### This R script fits model outputs to empirical time-series data ####
########################################################################


# Load packages and set working directory
library(tidyr)
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyverse)
library(cvequality)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# NOTE: CLEAR GLOBAL ENVIRONMENT IF ELAPSED TIME LIMIT ERROR OCCURS OR LINES ARE MISSING FROM PLOTS


# USER: enter species and location
species <- "Clavigralla tomentosicollis"
location <- "Nigeria"
field_plot <- "Mean" # for Nigeria, must specify plot "A", "B", "C", "Mean", or "All"

# USER: Include diurnal variation?
daily <- FALSE

# USER: Use left-skewed function for development?
left_skew <- FALSE # if FALSE, development plateaus between Topt and Tmax before going to zero above Tmax

# USER: Fit to census data?
census = TRUE

# USER: SET PLOT OPTIONS
xmin <- 0 # start date
xmax <- 2*365 # end date
num_yrs <- (xmax - xmin)/365
ymin <- 0 # min density
ymax <- 100 # max density
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


# READ IN TEMPERATURE RESPONSE PARAMETERS, TEMPERATURE PARAMETERS, AND DDE MODEL DYNAMICS
sp.data <- data[data$Species == paste(species,location),]
temp.data <- temp.data[temp.data$Species == paste(species,location),]
if(census == TRUE) {
  # For species with census data
  if(daily == TRUE) {
    data.model <- as.data.frame(read_csv(paste0("Time series data Census/Historical time series Diurnal ",species," ",location,".csv")))
    data.model.census <- as.data.frame(read_csv(paste0("Time series data Census/Historical time series Diurnal Census ",species," ",location,".csv")))
    data.model.CC <- as.data.frame(read_csv(paste0("Time series data Census/Future time series Diurnal ",species," ",location,".csv"))) }
  if(daily == FALSE && left_skew == TRUE) {
    data.model <- as.data.frame(read_csv(paste0("Time series data Census/Historical time series Tave ",species," ",location,".csv")))
    data.model.census <- as.data.frame(read_csv(paste0("Time series data Census/Historical time series Tave Census ",species," ",location,".csv")))
    data.model.CC <- as.data.frame(read_csv(paste0("Time series data Census/Future time series Tave ",species," ",location,".csv"))) }
  if(daily == FALSE && left_skew == FALSE) {
    data.model <- as.data.frame(read_csv(paste0("Time series data Census/Historical time series Tave Dev ",species," ",location,".csv")))
    data.model.census <- as.data.frame(read_csv(paste0("Time series data Census/Historical time series Tave Dev Census ",species," ",location,".csv")))
    data.model.CC <- as.data.frame(read_csv(paste0("Time series data Census/Future time series Tave Dev ",species," ",location,".csv"))) }
} else{
  # For species without census data
  if(daily == TRUE) {
    data.model <- as.data.frame(read_csv(paste0("Time series data Diurnal/Historical time series ",species," ",location,".csv")))
    data.model.CC <- as.data.frame(read_csv(paste0("Time series data/Future time series ",species," ",location,".csv"))) }
  if(daily == FALSE && left_skew == TRUE) {
    data.model <- as.data.frame(read_csv(paste0("Time series data Tave/Historical time series ",species," ",location,".csv")))
    data.model.CC <- as.data.frame(read_csv(paste0("Time series data Tave/Future time series ",species," ",location,".csv"))) }
  if(daily == FALSE && left_skew == FALSE) {
    data.model <- as.data.frame(read_csv(paste0("Time series data Tave Dev/Historical time series ",species," ",location,".csv")))
    data.model.CC <- as.data.frame(read_csv(paste0("Time series data Tave Dev/Future time series ",species," ",location,".csv"))) }
}


# PLOT OPTIONS
# for climate change time period
xmin.CC <- xmin
xmax.CC <- xmax
ymin.CC <- ymin
ymax.CC <- ymax
yr <- 365 # days in a year
init_yrs <- 65 # years to initialize model (see init_yrs in "DDE population dynamics.py")
start_date <- 10 - num_yrs # number of years before taking time-series data
TS.length <- xmax - xmin # length of time-series data
end <- nrow(data.model)
end.CC <- nrow(data.model.CC)


# FORMAT MODEL OUTPUT TO ALIGN WITH TIME-SERIES DATA
# Remove all rows before time-series data starts
data.model <- data.model[-c(1:(start_date*yr + xmin)), ]
if(census == TRUE) { data.model.census <- data.model.census[-c(1:(start_date*yr + xmin)), ] }

# Remove all rows after xmax days
if(xmax < end) { data.model <- data.model[-c(xmax+1:end), ] }
if(census == TRUE & xmax < end) { data.model.census <- data.model.census[-c(xmax+1:end), ] }

# climate change period (remove all but last num_yrs years of data)
if(xmax.CC < end.CC) { data.model.CC <- data.model.CC[-c(1:(end.CC-num_yrs*yr + xmin.CC)), ] }

# Re-scale time to start at xmin
# historical period
time.shift <- data.model[[1,1]] - xmin
data.model <- sweep(data.model, 2, c(time.shift,0,0,0,0))
if(census == TRUE) {
  time.shift.census <- data.model.census[[1,1]] - xmin
  data.model.census <- sweep(data.model.census, 2, c(time.shift.census,0,0,0,0)) }

# climate change period
time.shift.CC <- data.model.CC[[1,1]] + xmin.CC
data.model.CC <- sweep(data.model.CC, 2, c(time.shift.CC,0,0,0,0))


########################################## PLOTS #############################################
# HABITAT TEMPERATURE
# Historical time period
# data table from Tmin and Tmax functions
temp.fun.h <- data.frame(t=c(xmin:xmax))
temp.fun.h$fun.min <- sapply(temp.fun.h$t, FUN = function(t) { (temp.data$meanT.h - 273.15 + temp.data$delta_mean.h*(t+time.shift+init_yrs*yr))  - (temp.data$amplT.h + temp.data$delta_ampl.h*(t+time.shift+init_yrs*yr)) * cos(2*pi*((t+time.shift+init_yrs*yr) + temp.data$shiftT.h)/yr) - abs(temp.data$amplD.h) })
temp.fun.h$fun.max <- sapply(temp.fun.h$t, FUN = function(t) { (temp.data$meanT.h - 273.15 + temp.data$delta_mean.h*(t+time.shift+init_yrs*yr))  - (temp.data$amplT.h + temp.data$delta_ampl.h*(t+time.shift+init_yrs*yr)) * cos(2*pi*((t+time.shift+init_yrs*yr) + temp.data$shiftT.h)/yr) + abs(temp.data$amplD.h) })

# active season (when habitat temperature = Tmin)
start.h <- temp.fun.h[temp.fun.h$fun.min + 273.15 >= sp.data$Tmin, "t"][1]
end.h <- tail(temp.fun.h[temp.fun.h$fun.min + 273.15 >= sp.data$Tmin & temp.fun.h$t <= 365, "t"], n=1)

# plot
plot.temp <- ggplot(temp.fun.h, aes(x=t, y=fun.max)) +
  # Daily average temperature
  geom_function(fun = function(t) (temp.data$meanT.h - 273.15 + temp.data$delta_mean.h*(t+time.shift+init_yrs*yr))  - (temp.data$amplT.h + temp.data$delta_ampl.h*(t+time.shift+init_yrs*yr)) * cos(2*pi*((t+time.shift+init_yrs*yr) + temp.data$shiftT.h)/yr),
                size=1.5, color="#0072B2") + # blue color
  # Daily minimum temperature
  #geom_function(fun = function(t) (temp.data$meanT.h - 273.15 + temp.data$delta_mean.h*(t+time.shift+init_yrs*yr))  - (temp.data$amplT.h + temp.data$delta_ampl.h*(t+time.shift+init_yrs*yr)) * cos(2*pi*((t+time.shift+init_yrs*yr) + temp.data$shiftT.h)/yr) - temp.data$amplD.h,
  #              size=1.5, linetype="dashed", color="#0072B2") +
  # Daily maximum temperature
  #geom_function(fun = function(t) (temp.data$meanT.h - 273.15 + temp.data$delta_mean.h*(t+time.shift+init_yrs*yr))  - (temp.data$amplT.h + temp.data$delta_ampl.h*(t+time.shift+init_yrs*yr)) * cos(2*pi*((t+time.shift+init_yrs*yr) + temp.data$shiftT.h)/yr) + temp.data$amplD.h,
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
temp.fun.f$fun.min <- sapply(temp.fun.f$t, FUN = function(t) { (temp.data$meanT.f - 273.15 + temp.data$delta_mean.f*(t+time.shift.CC+init_yrs*yr))  - (temp.data$amplT.f + temp.data$delta_ampl.f*(t+time.shift.CC+init_yrs*yr)) * cos(2*pi*((t+time.shift.CC+init_yrs*yr) + temp.data$shiftT.f)/yr) - abs(temp.data$amplD.f) })
temp.fun.f$fun.max <- sapply(temp.fun.f$t, FUN = function(t) { (temp.data$meanT.f - 273.15 + temp.data$delta_mean.f*(t+time.shift.CC+init_yrs*yr))  - (temp.data$amplT.f + temp.data$delta_ampl.f*(t+time.shift.CC+init_yrs*yr)) * cos(2*pi*((t+time.shift.CC+init_yrs*yr) + temp.data$shiftT.f)/yr) + abs(temp.data$amplD.f) })

# active season (when habitat temperature = Tmin)
start.f <- temp.fun.f[temp.fun.f$fun.min + 273.15 >= sp.data$Tmin, "t"][1]
end.f <- tail(temp.fun.f[temp.fun.f$fun.min + 273.15 >= sp.data$Tmin & temp.fun.f$t <= 365, "t"], n=1)

# plot
plot.temp.CC <- ggplot(temp.fun.f, aes(x=t, y=fun.max)) +
  # Daily average temperature
  geom_function(fun = function(t) (temp.data$meanT.f - 273.15 + temp.data$delta_mean.f*(t+time.shift.CC+init_yrs*yr))  - (temp.data$amplT.f + temp.data$delta_ampl.f*(t+time.shift.CC+init_yrs*yr)) * cos(2*pi*((t+time.shift.CC+init_yrs*yr) + temp.data$shiftT.f)/yr),
                size=1.5, color="#D55E00") + # red color
  # Daily minimum temperature
  #geom_function(fun = function(t) (temp.data$meanT.f - 273.15 + temp.data$delta_mean.f*(t+time.shift.CC+init_yrs*yr))  - (temp.data$amplT.f + temp.data$delta_ampl.f*(t+time.shift.CC+init_yrs*yr)) * cos(2*pi*((t+time.shift.CC+init_yrs*yr) + temp.data$shiftT.f)/yr) - temp.data$amplD.f,
  #              size=1.5, linetype="dashed", color="#D55E00") +
  # Daily maximum temperature
  #geom_function(fun = function(t) (temp.data$meanT.f - 273.15 + temp.data$delta_mean.f*(t+time.shift.CC+init_yrs*yr))  - (temp.data$amplT.f + temp.data$delta_ampl.f*(t+time.shift.CC+init_yrs*yr)) * cos(2*pi*((t+time.shift.CC+init_yrs*yr) + temp.data$shiftT.f)/yr) + temp.data$amplD.f,
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


# PLOT BANDS
# Hot season (warmest half of year) in tropical habitats
# for temperature plots
band.hot.temp = ggplot() +
  geom_rect(aes(xmin=xmin, xmax=xmin + 3*yr/4 - temp.data$shiftT.h, ymin=0, ymax=Inf), fill="grey", alpha=0.5) +
  geom_rect(aes(xmin=xmin + 5*yr/4 - temp.data$shiftT.h, xmax=xmin + 7*yr/4 - temp.data$shiftT.h, ymin=0, ymax=Inf), fill="grey", alpha=0.5) +
  geom_rect(aes(xmin=xmin + 9*yr/4 - temp.data$shiftT.h, xmax=xmax, ymin=0, ymax=Inf), fill="grey", alpha=0.5) +
  labs(x="", y="") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(temp.min, temp.max)) +
  #scale_y_log10(limits=c(ymin, ymax)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(color = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=1.50))
# for density plots
band.hot.density = ggplot() +
  geom_rect(aes(xmin=xmin, xmax=xmin + 3*yr/4 - temp.data$shiftT.h, ymin=0, ymax=Inf), fill="grey", alpha=0.5) +
  geom_rect(aes(xmin=xmin + 5*yr/4 - temp.data$shiftT.h, xmax=xmin + 7*yr/4 - temp.data$shiftT.h, ymin=0, ymax=Inf), fill="grey", alpha=0.5) +
  geom_rect(aes(xmin=xmin + 9*yr/4 - temp.data$shiftT.h, xmax=xmax, ymin=0, ymax=Inf), fill="grey", alpha=0.5) +
  labs(x="", y="") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(ymin, ymax)) +
  #scale_y_log10(limits=c(ymin, ymax)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(color = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=1.50))
# Summer (warmest half of year) in temperate habitats
# for temperature plots
band.summer.temp = ggplot() +
  geom_rect(aes(xmin=xmin + yr/4, xmax=xmin + 3*yr/4, ymin=0, ymax=Inf), fill="grey", alpha=0.5) +
  geom_rect(aes(xmin=xmin + 5*yr/4, xmax=xmin + 7*yr/4, ymin=0, ymax=Inf), fill="grey", alpha=0.5) +
  labs(x="", y="") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(temp.min, temp.max)) +
  #scale_y_log10(limits=c(ymin, ymax)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(color = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=1.50))
# for density plots
band.summer.density = ggplot() +
  geom_rect(aes(xmin=xmin + yr/4, xmax=xmin + 3*yr/4, ymin=0, ymax=Inf), fill="grey", alpha=0.5) +
  geom_rect(aes(xmin=xmin + 5*yr/4, xmax=xmin + 7*yr/4, ymin=0, ymax=Inf), fill="grey", alpha=0.5) +
  labs(x="", y="") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(ymin, ymax)) +
  #scale_y_log10(limits=c(ymin, ymax)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(color = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=1.50))


# TIME-SERIES DATA
# Juvenile density
plot.J = ggplot(data.TS, aes(x=time, y=J, ymin=J_SE_L, ymax=J_SE_H)) + 
  geom_pointrange(size=0.5, color="#0072B2") + # blue color
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
  geom_pointrange(size=0.5, color="#0072B2") + # blue color
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

# Adult density: census
if(census == TRUE) {
  model.census = ggplot(data.model.census, aes(x=Time, y=A)) + 
    geom_line(size=1.5, color="#0072B2") + #, linetype = "longdash") + # blue color
    labs(x="", y="") +
    scale_x_continuous(limits=c(xmin, xmax)) +
    scale_y_continuous(limits=c(ymin, ymax)) +
    #scale_y_log10(limits=c(ymin, ymax)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
          axis.line = element_line(color = "black"), legend.position = "none", 
          axis.text = element_text(size=13), axis.title = element_text(size=1.50)) }

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


#################################### DRAW FINAL PLOTS #######################################
# COMPILE PLOTS
# Temperature plots
if(location == "Nigeria" | sp.data$Habitat == "Tropical") {
  plot.climate <- ggdraw()  +
    draw_plot(band.hot.temp, x = 0, y = 0, width = 1, height = 0.3) +
    draw_plot(plot.temp, x = 0, y = 0, width = 1, height = 0.3) +
    draw_plot(plot.temp.CC, x = 0, y = 0, width = 1, height = 0.3) }
if(str_split(location, boundary("word"), simplify = T)[,1] == "China"  | sp.data$Habitat != "Tropical") { 
  plot.climate <- ggdraw()  +
    draw_plot(band.summer.temp, x = 0, y = 0, width = 1, height = 0.3) +
    draw_plot(plot.temp, x = 0, y = 0, width = 1, height = 0.3) +
    draw_plot(plot.temp.CC, x = 0, y = 0, width = 1, height = 0.3) }

# Historical time period
# if(location == "Nigeria" | sp.data$Habitat == "Tropical") {
#   plot <- ggdraw()  +
#     draw_plot(band.hot.temp, x = 0, y = 0, width = 1, height = 0.3) +
#     draw_plot(plot.temp, x = 0, y = 0, width = 1, height = 0.3) +
#     draw_plot(band.hot.density, x = 0, y = 0.3, width = 1, height = 0.7) +
#     #draw_plot(plot.J, x = 0, y = 0.3, width = 1, height = 0.7) +
#     #draw_plot(model.J, x = 0, y = 0.3, width = 1, height = 0.7) +
#     draw_plot(model.A, x = 0, y = 0.3, width = 1, height = 0.7) +
#     draw_plot(model.census, x = 0, y = 0.3, width = 1, height = 0.7) +
#     draw_plot(plot.A, x = 0, y = 0.3, width = 1, height = 0.7) }
# if(str_split(location, boundary("word"), simplify = T)[,1] == "China" | sp.data$Habitat != "Tropical") { 
#   plot <- ggdraw()  +
#     draw_plot(band.summer.temp, x = 0, y = 0, width = 1, height = 0.3) +
#     draw_plot(plot.temp, x = 0, y = 0, width = 1, height = 0.3) +
#     draw_plot(band.summer.density, x = 0, y = 0.3, width = 1, height = 0.7) +
#     #draw_plot(model.J, x = 0, y = 0.3, width = 1, height = 0.7) +
#     draw_plot(model.A, x = 0, y = 0.3, width = 1, height = 0.7) +
#     draw_plot(model.census, x = 0, y = 0.3, width = 1, height = 0.7) +
#     draw_plot(plot.A, x = 0, y = 0.3, width = 1, height = 0.7) }

# Future time period
# if(location == "Nigeria" | sp.data$Habitat == "Tropical") {
#   plot.CC <- ggdraw()  +
#     draw_plot(band.hot.temp, x = 0, y = 0, width = 1, height = 0.3) +
#     draw_plot(plot.temp.CC, x = 0, y = 0, width = 1, height = 0.3) +
#     draw_plot(band.hot.density, x = 0, y = 0.3, width = 1, height = 0.7) +
#     draw_plot(model.J.CC, x = 0, y = 0.3, width = 1, height = 0.7) +
#     draw_plot(model.A.CC, x = 0, y = 0.3, width = 1, height = 0.7)
# }
# if(str_split(location, boundary("word"), simplify = T)[,1] == "China" | sp.data$Habitat != "Tropical") { 
#   plot.CC <- ggdraw()  +
#     draw_plot(band.summer.temp, x = 0, y = 0, width = 1, height = 0.3) +
#     draw_plot(plot.temp.CC, x = 0, y = 0, width = 1, height = 0.3) +
#     draw_plot(band.summer.density, x = 0, y = 0.3, width = 1, height = 0.7) +
#     draw_plot(model.J.CC, x = 0, y = 0.3, width = 1, height = 0.7) +
#     draw_plot(model.A.CC, x = 0, y = 0.3, width = 1, height = 0.7)
# }

# Compare historical and future time periods
if(location == "Nigeria" | sp.data$Habitat == "Tropical") {
  plot.compare <- ggdraw()  +
    draw_plot(band.hot.temp, x = 0, y = 0, width = 1, height = 0.3) +
    draw_plot(plot.temp, x = 0, y = 0, width = 1, height = 0.3) +
    draw_plot(plot.temp.CC, x = 0, y = 0, width = 1, height = 0.3) +
    draw_plot(band.hot.density, x = 0, y = 0.3, width = 1, height = 0.7) +
    #draw_plot(model.J, x = 0, y = 0.3, width = 1, height = 0.7) +
    #draw_plot(model.A, x = 0, y = 0.3, width = 1, height = 0.7) +
    draw_plot(model.census, x = 0, y = 0.3, width = 1, height = 0.7) +
    #draw_plot(model.J.CC, x = 0, y = 0.3, width = 1, height = 0.7) +
    #draw_plot(model.A.CC, x = 0, y = 0.3, width = 1, height = 0.7) +
    draw_plot(plot.A, x = 0, y = 0.3, width = 1, height = 0.7)
}
if(str_split(location, boundary("word"), simplify = T)[,1] == "China" | sp.data$Habitat != "Tropical") { 
  plot.compare <- ggdraw()  +
    draw_plot(band.summer.temp, x = 0, y = 0, width = 1, height = 0.3) +
    draw_plot(plot.temp, x = 0, y = 0, width = 1, height = 0.3) +
    draw_plot(plot.temp.CC, x = 0, y = 0, width = 1, height = 0.3) +
    draw_plot(band.summer.density, x = 0, y = 0.3, width = 1, height = 0.7) +
    #draw_plot(model.J, x = 0, y = 0.3, width = 1, height = 0.7) +
    draw_plot(model.A, x = 0, y = 0.3, width = 1, height = 0.7) +
    #draw_plot(model.census, x = 0, y = 0.3, width = 1, height = 0.7) +
    #draw_plot(model.J.CC, x = 0, y = 0.3, width = 1, height = 0.7) +
    draw_plot(model.A.CC, x = 0, y = 0.3, width = 1, height = 0.7) #+
    #draw_plot(plot.A, x = 0, y = 0.3, width = 1, height = 0.7)
}


# VIEW IN RSTUDIO
#plot
plot.compare


# STATISTICS
# linear regression of field data and model data
stats.data <- data.frame(Time = data.TS$time, Census = data.TS$A, Model = NA)
for(i in seq(0:nrow(stats.data))) {
  stats.data$Model[i] = data.model.census[data.model.census$Time == stats.data$Time[i],"A"]
}
summary(lm(stats.data$Model ~ 0 + stats.data$Census))
plot(stats.data$Census,stats.data$Model)

# mean and CV of historical versus future period
density.data <- rbind(data.frame(period = rep("hist",nrow(data.model)), J = data.model$J, A = data.model$A),
                      data.frame(period = rep("fut",nrow(data.model.CC)), J = data.model.CC$J, A = data.model.CC$A))
# juvenile density
#aggregate(density.data$J, list(density.data$period), function(x) c(mean = mean(x), sd = sd(x), cv = sd(x)/mean(x)))
#summary(aov(J ~ period, data=density.data)) # mean
#with(density.data, asymptotic_test(J, period)) # CV (asymptotic test)
#with(density.data, mslr_test(nr = 1000, J, period)) # CV (modified signed-likelihood ratio test)
# adult density
aggregate(density.data$A, list(density.data$period), function(x) c(mean = mean(x), sd = sd(x), cv = sd(x)/mean(x)))
summary(aov(A ~ period, data=density.data)) # mean
#with(density.data, asymptotic_test(A, period)) # CV (asymptotic test)
with(density.data, mslr_test(nr = 1000, A, period)) # CV (modified signed-likelihood ratio test)


# OUTPUT PLOTS
#dev.new()
# Temperature plots
# if(location == "Nigeria" | sp.data$Habitat == "Tropical") {
#   ggdraw()  +
#     draw_plot(band.hot.temp, x = 0, y = 0, width = 1, height = 0.3) +
#     draw_plot(plot.temp, x = 0, y = 0, width = 1, height = 0.3) +
#     draw_plot(plot.temp.CC, x = 0, y = 0, width = 1, height = 0.3) }
# if(str_split(location, boundary("word"), simplify = T)[,1] == "China" | sp.data$Habitat != "Tropical") {
#   ggdraw()  +
#     draw_plot(band.summer.temp, x = 0, y = 0, width = 1, height = 0.3) +
#     draw_plot(plot.temp, x = 0, y = 0, width = 1, height = 0.3) +
#     draw_plot(plot.temp.CC, x = 0, y = 0, width = 1, height = 0.3) }

# Climate change plots
# if(location == "Nigeria" | sp.data$Habitat == "Tropical") {
#   ggdraw()  +
#     draw_plot(band.hot.density, width = 1, height = 0.4) +
#     draw_plot(model.J, width = 1, height = 0.4) +
#     draw_plot(model.J.CC, width = 1, height = 0.4) }
# if(location == "Nigeria" | sp.data$Habitat == "Tropical") {
#   ggdraw()  +
#     draw_plot(band.hot.density, width = 1, height = 0.4) +
#     draw_plot(model.A, width = 1, height = 0.4) +
#     draw_plot(model.census, width = 1, height = 0.4) +
#     draw_plot(model.A.CC, width = 1, height = 0.4) +
#     draw_plot(plot.A,width = 1, height = 0.4) }
# if(str_split(location, boundary("word"), simplify = T)[,1] == "China" | sp.data$Habitat != "Tropical") {
#   ggdraw()  +
#     draw_plot(band.summer.density, width = 1, height = 0.4) +
#     draw_plot(model.J, width = 1, height = 0.4) +
#     draw_plot(model.J.CC, width = 1, height = 0.4) }
# if(str_split(location, boundary("word"), simplify = T)[,1] == "China" | sp.data$Habitat != "Tropical") {
#   ggdraw()  +
#     draw_plot(band.summer.density, width = 1, height = 0.4) +
#     draw_plot(model.A, width = 1, height = 0.4) +
#     draw_plot(model.census, width = 1, height = 0.4) +
#     draw_plot(model.A.CC, width = 1, height = 0.4) +
#     draw_plot(plot.A,width = 1, height = 0.4) }
