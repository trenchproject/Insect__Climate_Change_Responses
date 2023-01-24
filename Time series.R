########################################################################
#### This R script fits model outputs to empirical time-series data ####
########################################################################

# Load packages and set working directory
library(tidyverse)
library(ggplot2) # for plotting
library(cowplot) # for combining plots
library(stringr) # for manipulating strings

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#setwd('..')

# NOTE: CLEAR GLOBAL ENVIRONMENT IF ELAPSED TIME LIMIT ERROR OCCURS OR LINES ARE MISSING FROM PLOTS

#old <- as.data.frame(read_csv("Time series data Census/Historical time series Clavigralla tomentosicollis Benin.csv"))
#new <- as.data.frame(read_csv("Time series data Census new/Historical time series Clavigralla tomentosicollis Benin.csv"))
#identical(old, new)


# USER: enter species and location
species <- "Clavigralla tomentosicollis"
location <- "Benin"

# USER: Fit to census data?
census <- TRUE


# READ IN HABITAT TEMPERATURE PARAMETERS, FIELD CENSUS DATA, AND MODEL TIME-SERIES DATA
# Habitat temperature parameters
temp.param <- as.data.frame(read_csv("Model parameters/Habitat temperature parameters.csv"))
temp.param <- temp.param[temp.param$Species == paste(species,location),]

# Field census data
if(census == TRUE && location == "Benin") { data.census <- as.data.frame(read_csv("Biological data/Census data Benin.csv")) }
if(location == "China Dafeng") {
  data.census <- as.data.frame(read_csv("Biological data/Census data China.csv"))
  data.census <- data.census %>% filter(Location = "Dafeng", species = "Apolygus lucorum") }

# Model time-series data
if(census == TRUE) {
  model.his <- as.data.frame(read_csv(paste0("Time series data Census new/Census time series ",species," ",location,".csv")))
} else {
  model.his <- as.data.frame(read_csv(paste0("Time series data new/Historical time series ",species," ",location,".csv")))
}
model.fut <- as.data.frame(read_csv(paste0("Time series data new/Future time series ",species," ",location,".csv")))


# SET PARAMETERS AND PLOT OPTIONS
# Parameters for time-series data
yr <- 365 # days per year
# historical time-series
start.his <- 0 # how many years into climate change period to start population dynamics model (see "start_yr" in "DDE population dynamics.py")
num.his <- 10 # how long to run model (see "num_yrs" in "DDE population dynamics.py")
plot.his <- 2 # how many years to plot (i.e., last X years)
start.his <- start.his + num.his - plot.his # recalculate start year for extracting time-series data
end.his <- min(start.his + plot.his, start.his + num.his) # calculate end year for extracting time-series data
length.his <- nrow(model.his)
# future time-series
start.fut <- 0 # how many years into climate change period to start population dynamics model (see "start_yr" in "DDE population dynamics.py")
num.fut <- 75 # how long to run model (see "num_yrs" in "DDE population dynamics.py")
plot.fut <- 2 # how many years to plot (i.e., last X years)
start.fut <- start.fut + num.fut - plot.fut # recalculate start year for extracting time-series data
end.fut <- min(start.fut + plot.fut, start.fut + num.fut) # calculate end year for extracting time-series data
length.fut <- nrow(model.fut)

# Plotting parameters
xmin.his <- 0 # start time for plots (may differ from "start.his" above)
xmax.his <- plot.his*yr # end time for plots (may differ from "end.his" above)
xmin.fut <- 0 # start time for plots (may differ from "start.fut" above)
xmax.fut <- plot.fut*yr # end time for plots (may differ from "end.fut" above)
ymin <- 0 # min density for plots (same for historical and future periods)
ymax <- 120 # max density for plots (same for historical and future periods)
temp.min <- 0 # min temperature for plots (same for historical and future periods)
temp.max <- 40 # max temperature for plots (same for historical and future periods)
TS.length <- xmax.his - xmin.his # length of time-series data for plots


# FORMAT MODEL OUTPUT TO ALIGN WITH TIME-SERIES DATA
# Extract time-series data for plotting
model.his <- model.his[(start.his*yr+1):(end.his*yr),]
model.fut <- model.fut[(start.fut*yr+1):(end.fut*yr),]

# Re-scale time in model time-series data to start at xmin.his
# historical period
time.shift.his <- max(model.his[[1,1]] - xmin.his, 0) # number of days to "shift" time-series
model.his <- sweep(model.his, 2, c(time.shift.his,0,0,0,0)) # shift "Time" column to align with x-axis of plots

# future period
time.shift.fut <- max(model.fut[[1,1]] + xmin.fut, 0) # number of days to "shift" time-series
model.fut <- sweep(model.fut, 2, c(time.shift.fut,0,0,0,0)) # shift "Time" column to align with x-axis of plots


################################ CONSTRUCT PLOTS ################################
# HABITAT TEMPERATURE
# Historical time period
temp.his <- ggplot() +
  geom_function(fun = function(t) (temp.param$meanT.h - 273.15 + temp.param$delta_mean.h*(t+start.his*yr))  - (temp.param$amplT.h + temp.param$delta_ampl.h*(t+start.his*yr)) * cos(2*pi*((t+start.his*yr) + temp.param$shiftT.h)/yr),
                linewidth=1.5, color="#0072B2") +
  labs(x="", y="") +
  scale_x_continuous(limits=c(xmin.his, xmax.his)) +
  scale_y_continuous(limits=c(temp.min, temp.max)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(color = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=1.50))

# Future time period
temp.fut <- ggplot() +
  geom_function(fun = function(t) (temp.param$meanT.f - 273.15 + temp.param$delta_mean.f*(t+start.fut*yr))  - (temp.param$amplT.f + temp.param$delta_ampl.f*(t+start.fut*yr)) * cos(2*pi*((t+start.fut*yr) + temp.param$shiftT.f)/yr),
                linewidth=1.5, color="#D55E00") +
  labs(x="", y="") +
  scale_x_continuous(limits=c(xmin.fut, xmax.fut)) +
  scale_y_continuous(limits=c(temp.min, temp.max)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(color = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=1.50))


# PLOT TEMPERATURE BANDS AROUND WARMEST HALF OF YEAR
# Tropical habitats
# for temperature plots
band.hot.temp <- ggplot() +
  geom_rect(aes(xmin=xmin.his, xmax=xmin.his + 3*yr/4 - temp.param$shiftT.h, ymin=0, ymax=Inf), fill="grey", alpha=0.5) +
  geom_rect(aes(xmin=xmin.his + 5*yr/4 - temp.param$shiftT.h, xmax=xmin.his + 7*yr/4 - temp.param$shiftT.h, ymin=0, ymax=Inf), fill="grey", alpha=0.5) +
  geom_rect(aes(xmin=xmin.his + 9*yr/4 - temp.param$shiftT.h, xmax=xmax.his, ymin=0, ymax=Inf), fill="grey", alpha=0.5) +
  labs(x="", y="") +
  scale_x_continuous(limits=c(xmin.his, xmax.his)) +
  scale_y_continuous(limits=c(temp.min, temp.max)) +
  #scale_y_log10(limits=c(ymin, ymax)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(color = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=1.50))

# for insect density plots
band.hot.density <- ggplot() +
  geom_rect(aes(xmin=xmin.his, xmax=xmin.his + 3*yr/4 - temp.param$shiftT.h, ymin=0, ymax=Inf), fill="grey", alpha=0.5) +
  geom_rect(aes(xmin=xmin.his + 5*yr/4 - temp.param$shiftT.h, xmax=xmin.his + 7*yr/4 - temp.param$shiftT.h, ymin=0, ymax=Inf), fill="grey", alpha=0.5) +
  geom_rect(aes(xmin=xmin.his + 9*yr/4 - temp.param$shiftT.h, xmax=xmax.his, ymin=0, ymax=Inf), fill="grey", alpha=0.5) +
  labs(x="", y="") +
  scale_x_continuous(limits=c(xmin.his, xmax.his)) +
  scale_y_continuous(limits=c(ymin, ymax)) +
  #scale_y_log10(limits=c(ymin, ymax)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(color = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=1.50))

# Temperate habitats
# for temperature plots
band.summer.temp <- ggplot() +
  geom_rect(aes(xmin=xmin.his + yr/4, xmax=xmin.his + 3*yr/4, ymin=0, ymax=Inf), fill="grey", alpha=0.5) +
  geom_rect(aes(xmin=xmin.his + 5*yr/4, xmax=xmin.his + 7*yr/4, ymin=0, ymax=Inf), fill="grey", alpha=0.5) +
  labs(x="", y="") +
  scale_x_continuous(limits=c(xmin.his, xmax.his)) +
  scale_y_continuous(limits=c(temp.min, temp.max)) +
  #scale_y_log10(limits=c(ymin, ymax)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(color = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=1.50))

# for insect density plots
band.summer.density <- ggplot() +
  geom_rect(aes(xmin=xmin.his + yr/4, xmax=xmin.his + 3*yr/4, ymin=0, ymax=Inf), fill="grey", alpha=0.5) +
  geom_rect(aes(xmin=xmin.his + 5*yr/4, xmax=xmin.his + 7*yr/4, ymin=0, ymax=Inf), fill="grey", alpha=0.5) +
  labs(x="", y="") +
  scale_x_continuous(limits=c(xmin.his, xmax.his)) +
  scale_y_continuous(limits=c(ymin, ymax)) +
  #scale_y_log10(limits=c(ymin, ymax)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(color = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=1.50))


# FIELD CENSUS DATA
plot.census <- ggplot(data.census, aes(x=day, y=Adults, ymin=A_SE_L, ymax=A_SE_H)) + 
  geom_pointrange(size=0.5, color="#0072B2") + # blue color
  labs(x="", y="") +
  scale_x_continuous(limits=c(xmin.his, xmax.his)) +
  scale_y_continuous(limits=c(ymin, ymax)) +
  #scale_y_log10(limits=c(ymin, ymax)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(color = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=1.50)) 


# DDE POPULATION MODEL DATA
# Historical time period
# Juvenile density
model.J.his <- ggplot(model.his, aes(x=Time, y=J)) + 
  geom_line(size=1.5, color="#0072B2") + # blue color
  labs(x="", y="") +
  scale_x_continuous(limits=c(xmin.his, xmax.his)) +
  scale_y_continuous(limits=c(ymin, ymax)) +
  #scale_y_log10(limits=c(ymin, ymax)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(color = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=1.50))

# Adult density
model.A.his <- ggplot(model.his, aes(x=Time, y=A)) + 
  geom_line(size=1.5, color="#0072B2") + # blue color
  labs(x="", y="") +
  scale_x_continuous(limits=c(xmin.his, xmax.his)) +
  scale_y_continuous(limits=c(ymin, ymax)) +
  #scale_y_log10(limits=c(ymin, ymax)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(color = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=1.50))

# Future time period
# Juvenile density
model.J.fut <- ggplot(model.fut, aes(x=Time, y=J)) + 
  geom_line(size=1.5, color="#D55E00") + # red color
  labs(x="", y="") +
  scale_x_continuous(limits=c(xmin.fut, xmax.fut)) +
  scale_y_continuous(limits=c(ymin, ymax)) +
  #scale_y_log10(limits=c(ymin, ymax)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(color = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=1.50))

# Adult density
model.A.fut <- ggplot(model.fut, aes(x=Time, y=A)) + 
  geom_line(size=1.5, color="#D55E00") + # red color
  labs(x="", y="") +
  scale_x_continuous(limits=c(xmin.fut, xmax.fut)) +
  scale_y_continuous(limits=c(ymin, ymax)) +
  #scale_y_log10(limits=c(ymin, ymax)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(color = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=1.50)) 


################################# DRAW PLOTS ####################################
# Habitat temperature plot (Fig. 1A,B)
if(temp.param$Habitat == "Tropical") {
  ggdraw()  +
    draw_plot(band.hot.temp, x = 0, y = 0, width = 1, height = 0.3) +
    draw_plot(temp.his, x = 0, y = 0, width = 1, height = 0.3) +
    draw_plot(temp.fut, x = 0, y = 0, width = 1, height = 0.3)
} else { 
  ggdraw()  +
    draw_plot(band.summer.temp, x = 0, y = 0, width = 1, height = 0.3) +
    draw_plot(temp.his, x = 0, y = 0, width = 1, height = 0.3) +
    draw_plot(temp.fut, x = 0, y = 0, width = 1, height = 0.3) }

# Model and field census data plot (Fig. 1C,D)
if(census == TRUE) {
  if(temp.param$Habitat == "Tropical") {
    ggdraw()  +
      draw_plot(band.hot.density, x = 0, y = 0, width = 1, height = 0.4) +
      draw_plot(model.A.his, x = 0, y = 0, width = 1, height = 0.4) +
      draw_plot(plot.census, x = 0, y = 0, width = 1, height = 0.4)
  } else {
    ggdraw()  +
      draw_plot(band.summer.density, x = 0, y = 0, width = 1, height = 0.4) +
      draw_plot(model.A.his, x = 0, y = 0, width = 1, height = 0.4) +
      draw_plot(plot.census, x = 0, y = 0, width = 1, height = 0.4) }
}



# # Compare historical and future time periods
# if(location == "Nigeria" || temp.param$Habitat == "Tropical") {
#   plot.compare <- ggdraw()  +
#     draw_plot(band.hot.temp, x = 0, y = 0, width = 1, height = 0.3) +
#     draw_plot(plot.temp, x = 0, y = 0, width = 1, height = 0.3) +
#     draw_plot(plot.temp.fut, x = 0, y = 0, width = 1, height = 0.3) +
#     draw_plot(band.hot.density, x = 0, y = 0.3, width = 1, height = 0.7) +
#     #draw_plot(model.J, x = 0, y = 0.3, width = 1, height = 0.7) +
#     #draw_plot(model.A, x = 0, y = 0.3, width = 1, height = 0.7) +
#     draw_plot(model.census, x = 0, y = 0.3, width = 1, height = 0.7) +
#     #draw_plot(model.J.fut, x = 0, y = 0.3, width = 1, height = 0.7) +
#     #draw_plot(model.A.fut, x = 0, y = 0.3, width = 1, height = 0.7) +
#     draw_plot(plot.census, x = 0, y = 0.3, width = 1, height = 0.7)
# }
# if(str_split(location, boundary("word"), simplify = T)[,1] == "China" || temp.param$Habitat != "Tropical") { 
#   plot.compare <- ggdraw()  +
#     draw_plot(band.summer.temp, x = 0, y = 0, width = 1, height = 0.3) +
#     draw_plot(plot.temp, x = 0, y = 0, width = 1, height = 0.3) +
#     draw_plot(plot.temp.fut, x = 0, y = 0, width = 1, height = 0.3) +
#     draw_plot(band.summer.density, x = 0, y = 0.3, width = 1, height = 0.7) +
#     #draw_plot(model.J, x = 0, y = 0.3, width = 1, height = 0.7) +
#     draw_plot(model.A, x = 0, y = 0.3, width = 1, height = 0.7) +
#     #draw_plot(model.census, x = 0, y = 0.3, width = 1, height = 0.7) +
#     #draw_plot(model.J.fut, x = 0, y = 0.3, width = 1, height = 0.7) +
#     draw_plot(model.A.fut, x = 0, y = 0.3, width = 1, height = 0.7) #+
#     #draw_plot(plot.census, x = 0, y = 0.3, width = 1, height = 0.7)
# }


# # VIEW IN RSTUDIO
# #plot
# plot.compare
# 
# 
# # STATISTICS
# #linear regression of field data and model data
# stats.data <- data.frame(Time = data.census$time, Census = data.census$A, Model = NA)
# for(i in seq(1:nrow(stats.data))) {
#   stats.data$Model[i] = model.his.census[model.his.census$Time == stats.data$Time[i],"A"]
# }
# summary(lm(stats.data$Model ~ 0 + stats.data$Census))
# plot(stats.data$Census,stats.data$Model, xlim = c(0,100), ylim = c(0,100))


# OUTPUT PLOTS
#dev.new()
# Temperature plots
# if(location == "Nigeria" || temp.param$Habitat == "Tropical") {
#   ggdraw()  +
#     draw_plot(band.hot.temp, x = 0, y = 0, width = 1, height = 0.3) +
#     draw_plot(plot.temp, x = 0, y = 0, width = 1, height = 0.3) +
#     draw_plot(plot.temp.fut, x = 0, y = 0, width = 1, height = 0.3) }
# if(str_split(location, boundary("word"), simplify = T)[,1] == "China" || temp.param$Habitat != "Tropical") {
#   ggdraw()  +
#     draw_plot(band.summer.temp, x = 0, y = 0, width = 1, height = 0.3) +
#     draw_plot(plot.temp, x = 0, y = 0, width = 1, height = 0.3) +
#     draw_plot(plot.temp.fut, x = 0, y = 0, width = 1, height = 0.3) }

# Climate change plots
# Tropical species
# if(location == "Nigeria" || temp.param$Habitat == "Tropical") {
# # juveniles
#   ggdraw()  +
#     draw_plot(band.hot.density, width = 1, height = 0.4) +
#     draw_plot(model.J, width = 1, height = 0.4) +
#     draw_plot(model.J.fut, width = 1, height = 0.4) }
#if(location == "Nigeria" || temp.param$Habitat == "Tropical") {
## adults
  # ggdraw()  +
  #   draw_plot(band.hot.density, width = 1, height = 0.4) +
  #   draw_plot(model.A, width = 1, height = 0.4) +
  #   #draw_plot(model.census, width = 1, height = 0.4) +
  #   draw_plot(model.A.fut, width = 1, height = 0.4) } #+
  #   #draw_plot(plot.census,width = 1, height = 0.4) }

# Temperate species
# if(str_split(location, boundary("word"), simplify = T)[,1] == "China" || temp.param$Habitat != "Tropical") {
## juveniles
#   ggdraw()  +
#     draw_plot(band.summer.density, width = 1, height = 0.4) +
#     draw_plot(model.J, width = 1, height = 0.4) +
#     draw_plot(model.J.fut, width = 1, height = 0.4) }
# if(str_split(location, boundary("word"), simplify = T)[,1] == "China" || temp.param$Habitat != "Tropical") {
## adults
#   ggdraw()  +
#     draw_plot(band.summer.density, width = 1, height = 0.4) +
#     draw_plot(model.A, width = 1, height = 0.4) +
#     draw_plot(model.census, width = 1, height = 0.4) +
#     draw_plot(model.A.fut, width = 1, height = 0.4) +
#     draw_plot(plot.census,width = 1, height = 0.4) }