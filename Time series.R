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

# USER: enter species and location
species <- "Apolygus lucorum"
location <- "China Dafeng"
species <- "Clavigralla shadabi"
location <- "Benin"

# USER: Fit model prediction to field census data?
census <- TRUE


# READ IN HABITAT TEMPERATURE PARAMETERS, FIELD CENSUS DATA, AND MODEL TIME-SERIES DATA
# Habitat temperature parameters
temp.param <- as.data.frame(read_csv("Model parameters/Habitat temperature parameters.csv"))
temp.param <- temp.param[temp.param$Species == paste(species,location),]

# Field census data
if(census == TRUE && location == "Benin") {
  data.census <- as.data.frame(read_csv("Biological data/Census data Benin.csv"))
  start.his <- 1982 }
if(census == TRUE && location == "China Dafeng") {
  data.census <- as.data.frame(read_csv("Biological data/Census data China.csv"))
  data.census <- data.census %>% filter(location == "Dafeng", species == "Apolygus lucorum")
  start.his <- 1996 }

# Model time-series data
if(census == TRUE) {
  model.his <- as.data.frame(read_csv(paste0("Time series data new/Census time series ",species," ",location,".csv")))
} else { model.his <- as.data.frame(read_csv(paste0("Time series data new/Historical time series ",species," ",location,".csv"))) }
model.fut <- as.data.frame(read_csv(paste0("Time series data new/Future time series ",species," ",location,".csv")))


# SET PARAMETERS AND PLOT OPTIONS
# Parameters for time-series data
yr <- 365 # days per year
if(census == FALSE) { start.his <- 0 } # how many years before starting population dynamics model (see "start_yr" in "DDE population dynamics.py")
start.fut <- 65 # how many years into climate change period to start population dynamics model (see "start_yr" in "DDE population dynamics.py")

# Plotting parameters
xmin <- 0 # always plot model starting at time zero
xmax <- 2*yr # plot 2 years of population dynamics
ymin <- 0 # min density for plots
ymax <- 125 # max density for plots
temp.min <- 0 # min temperature in C for plots
temp.max <- 40 # max temperature in C for plots


# FORMAT MODEL OUTPUT TO ALIGN WITH TIME-SERIES DATA
# Extract time-series data for plotting
model.his <- model.his[(nrow(model.his) - 2*yr + 1):nrow(model.his),]
model.fut <- model.fut[(nrow(model.fut) - 2*yr + 1):nrow(model.fut),]

# Re-scale time in model time-series data to start at xmin
model.his <- sweep(model.his, 2, c(model.his[[1,1]],0,0,0,0)) # shift "Time" column to align with x-axis of plots
model.fut <- sweep(model.fut, 2, c(model.fut[[1,1]],0,0,0,0)) # shift "Time" column to align with x-axis of plots


################################ CONSTRUCT PLOTS ################################
# HABITAT TEMPERATURE
# Historical time period
temp.his <- ggplot() +
  geom_function(fun = function(t) (temp.param$meanT.h - 273.15 + temp.param$delta_mean.h*(t+start.his*yr))  - (temp.param$amplT.h + temp.param$delta_ampl.h*(t+start.his*yr)) * cos(2*pi*((t+start.his*yr) + temp.param$shiftT.h)/yr),
                linewidth=1.5, color="#0072B2") +
  labs(x="", y="") +
  scale_x_continuous(limits=c(xmin, xmax)) +
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
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(temp.min, temp.max)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(color = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=1.50))


# PLOT TEMPERATURE BANDS AROUND WARMEST HALF OF YEAR
# Tropical habitats
# for temperature plots
band.hot.temp <- ggplot() +
  geom_rect(aes(xmin=xmin, xmax=xmin + 3*yr/4 - temp.param$shiftT.h, ymin=0, ymax=Inf), fill="grey", alpha=0.5) +
  geom_rect(aes(xmin=xmin + 5*yr/4 - temp.param$shiftT.h, xmax=xmin + 7*yr/4 - temp.param$shiftT.h, ymin=0, ymax=Inf), fill="grey", alpha=0.5) +
  geom_rect(aes(xmin=xmin + 9*yr/4 - temp.param$shiftT.h, xmax=xmax, ymin=0, ymax=Inf), fill="grey", alpha=0.5) +
  labs(x="", y="") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(temp.min, temp.max)) +
  #scale_y_log10(limits=c(ymin, ymax)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(color = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=1.50))

# for insect density plots
band.hot.density <- ggplot() +
  geom_rect(aes(xmin=xmin, xmax=xmin + 3*yr/4 - temp.param$shiftT.h, ymin=0, ymax=Inf), fill="grey", alpha=0.5) +
  geom_rect(aes(xmin=xmin + 5*yr/4 - temp.param$shiftT.h, xmax=xmin + 7*yr/4 - temp.param$shiftT.h, ymin=0, ymax=Inf), fill="grey", alpha=0.5) +
  geom_rect(aes(xmin=xmin + 9*yr/4 - temp.param$shiftT.h, xmax=xmax, ymin=0, ymax=Inf), fill="grey", alpha=0.5) +
  labs(x="", y="") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(ymin, ymax)) +
  #scale_y_log10(limits=c(ymin, ymax)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(color = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=1.50))

# Temperate habitats
# for temperature plots
band.summer.temp <- ggplot() +
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

# for insect density plots
band.summer.density <- ggplot() +
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


# FIELD CENSUS DATA
plot.census <- ggplot(data.census, aes(x=day, y=Adults, ymin=A_SE_L, ymax=A_SE_H)) + 
  geom_pointrange(size=0.5, color="#0072B2") + # blue color
  labs(x="", y="") +
  scale_x_continuous(limits=c(xmin, xmax)) +
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
  scale_x_continuous(limits=c(xmin, xmax)) +
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
  scale_x_continuous(limits=c(xmin, xmax)) +
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
  scale_x_continuous(limits=c(xmin, xmax)) +
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
  scale_x_continuous(limits=c(xmin, xmax)) +
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