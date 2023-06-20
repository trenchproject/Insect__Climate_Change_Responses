#################################################################################
## This R script plots model population dynamics (including field census data) ##
#################################################################################

# Load packages
library(tidyverse)
library(ggplot2) # for plotting
library(cowplot) # for combining plots
library(stringr) # for manipulating strings

# Set working directory
#setwd() # enter working directory of main downloaded file (containing R project file)


# NOTE: CLEAR GLOBAL ENVIRONMENT IF ELAPSED TIME LIMIT ERROR OCCURS OR LINES ARE MISSING FROM PLOTS


# USER: enter species and location
species <- "Clavigralla shadabi"
location <- "Benin"

# USER: Plot model population dynamics with field census data?
census <- TRUE


# READ IN HABITAT TEMPERATURE PARAMETERS, FIELD CENSUS DATA, AND MODEL TIME-SERIES DATA
# Habitat temperature parameters
t.param <- as.data.frame(read_csv("Model parameters/Habitat temperature parameters.csv"))
t.param <- t.param[t.param$Population == paste(species,location),]

# Field census data
if(census == TRUE && location == "Benin") {
  data.census <- as.data.frame(read_csv("Biological data/Census data Benin.csv"))
  start.r <- 1982 }
if(census == TRUE && location == "China Dafeng") {
  data.census <- as.data.frame(read_csv("Biological data/Census data China.csv"))
  start.r <- 1996 }

# Model time-series data
if(census == TRUE) {
  model.r <- as.data.frame(read_csv(paste0("Time series data/Census time series ",species," ",location,".csv")))
} else { model.r <- as.data.frame(read_csv(paste0("Time series data/Recent time series ",species," ",location,".csv"))) }
model.f <- as.data.frame(read_csv(paste0("Time series data/Future time series ",species," ",location,".csv")))


# SET PARAMETERS AND PLOT OPTIONS
# Parameters for time-series data
yr <- 365 # days per year
if(census == FALSE) { start.r <- 0 } # how many years before starting population dynamics model (see "start_yr" in "DDE population dynamics.py")
start.f <- 65 # how many years into climate change period to start population dynamics model (see "start_yr" in "DDE population dynamics.py")

# Plotting parameters
xmin <- 0 # always plot model starting at time zero
xmax <- 2*yr # plot 2 years of population dynamics
ymin <- 0 # min density for plots
ymax <- 125 # max density for plots
temp.min <- 0 # min temperature in C for plots
temp.max <- 40 # max temperature in C for plots


# FORMAT MODEL OUTPUT TO ALIGN WITH CENSUS DATA
# Extract time-series data for plotting
model.r <- model.r[(nrow(model.r) - 2*yr + 1):nrow(model.r),]
model.f <- model.f[(nrow(model.f) - 2*yr + 1):nrow(model.f),]

# Re-scale time in model time-series data to start at xmin
model.r <- sweep(model.r, 2, c(model.r[[1,1]],0,0,0,0)) # shift "Time" column to align with x-axis of plots
model.f <- sweep(model.f, 2, c(model.f[[1,1]],0,0,0,0)) # shift "Time" column to align with x-axis of plots


################################ CONSTRUCT PLOTS ################################
# HABITAT TEMPERATURE
# Recent time period
temp.r <- ggplot() +
  geom_function(fun = function(t) (t.param$meanT.r - 273.15 + t.param$delta_mean.r*(t+start.r*yr))  - (t.param$amplT.r + t.param$delta_ampl.r*(t+start.r*yr)) * cos(2*pi*((t+start.r*yr) + t.param$shiftT.r)/yr),
                linewidth=1.5, color="#0072B2") +
  labs(x="", y="") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(temp.min, temp.max)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(color = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=1.50))

# Future time period
temp.f <- ggplot() +
  geom_function(fun = function(t) (t.param$meanT.f - 273.15 + t.param$delta_mean.f*(t+start.f*yr))  - (t.param$amplT.f + t.param$delta_ampl.f*(t+start.f*yr)) * cos(2*pi*((t+start.f*yr) + t.param$shiftT.f)/yr),
                linewidth=1.5, color="#D55E00") +
  labs(x="", y="") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(temp.min, temp.max)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(color = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=1.50))


# PLOT TEMPERATURE BANDS AROUND THE WARMEST HALF OF YEAR
# Tropical/subtropical habitats
if(t.param$Habitat == "Tropical" || t.param$Habitat == "Subropical") {
  # for habitat temperature plots
  band.temp <- ggplot() +
    geom_rect(aes(xmin=xmin, xmax=xmin + 3*yr/4 - t.param$shiftT.r, ymin=0, ymax=Inf), fill="grey", alpha=0.5) +
    geom_rect(aes(xmin=xmin + 5*yr/4 - t.param$shiftT.r, xmax=xmin + 7*yr/4 - t.param$shiftT.r, ymin=0, ymax=Inf), fill="grey", alpha=0.5) +
    geom_rect(aes(xmin=xmin + 9*yr/4 - t.param$shiftT.r, xmax=xmax, ymin=0, ymax=Inf), fill="grey", alpha=0.5) +
    labs(x="", y="") +
    scale_x_continuous(limits=c(xmin, xmax)) +
    scale_y_continuous(limits=c(temp.min, temp.max)) +
    #scale_y_log10(limits=c(ymin, ymax)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
          axis.line = element_line(color = "black"), legend.position = "none", 
          axis.text = element_text(size=13), axis.title = element_text(size=1.50))
  
  # for insect density plots
  band.density <- ggplot() +
    geom_rect(aes(xmin=xmin, xmax=xmin + 3*yr/4 - t.param$shiftT.r, ymin=0, ymax=Inf), fill="grey", alpha=0.5) +
    geom_rect(aes(xmin=xmin + 5*yr/4 - t.param$shiftT.r, xmax=xmin + 7*yr/4 - t.param$shiftT.r, ymin=0, ymax=Inf), fill="grey", alpha=0.5) +
    geom_rect(aes(xmin=xmin + 9*yr/4 - t.param$shiftT.r, xmax=xmax, ymin=0, ymax=Inf), fill="grey", alpha=0.5) +
    labs(x="", y="") +
    scale_x_continuous(limits=c(xmin, xmax)) +
    scale_y_continuous(limits=c(ymin, ymax)) +
    #scale_y_log10(limits=c(ymin, ymax)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
          axis.line = element_line(color = "black"), legend.position = "none", 
          axis.text = element_text(size=13), axis.title = element_text(size=1.50))

} else {
  # Temperate habitats
  # for habitat temperature plots
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
}


# FIELD CENSUS DATA
if(census == TRUE) {
  plot.census <- ggplot(data.census, aes(x=Time, y=A, ymin=A_SE_L, ymax=A_SE_H)) + 
    geom_pointrange(size=0.5, color="#0072B2") + # blue color
    labs(x="", y="") +
    scale_x_continuous(limits=c(xmin, xmax)) +
    scale_y_continuous(limits=c(ymin, ymax)) +
    #scale_y_log10(limits=c(ymin, ymax)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
          axis.line = element_line(color = "black"), legend.position = "none", 
          axis.text = element_text(size=13), axis.title = element_text(size=1.50))
}


# DDE POPULATION MODEL DATA
# Recent time period
# Juvenile density
model.J.r <- ggplot(model.r, aes(x=Time, y=J)) + 
  geom_line(linewidth=1.5, color="#0072B2") + # blue color
  labs(x="", y="") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(ymin, ymax)) +
  #scale_y_log10(limits=c(ymin, ymax)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(color = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=1.50))

# Adult density
model.A.r <- ggplot(model.r, aes(x=Time, y=A)) + 
  geom_line(linewidth=1.5, color="#0072B2") + # blue color
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
model.J.f <- ggplot(model.f, aes(x=Time, y=J)) + 
  geom_line(linewidth=1.5, color="#D55E00") + # red color
  labs(x="", y="") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(ymin, ymax)) +
  #scale_y_log10(limits=c(ymin, ymax)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(color = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=1.50))

# Adult density
model.A.f <- ggplot(model.f, aes(x=Time, y=A)) + 
  geom_line(linewidth=1.5, color="#D55E00") + # red color
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
#dev.new()
if(t.param$Habitat == "Tropical") {
  ggdraw()  +
    draw_plot(band.temp, x = 0, y = 0, width = 1, height = 0.3) +
    draw_plot(temp.r, x = 0, y = 0, width = 1, height = 0.3) +
    draw_plot(temp.f, x = 0, y = 0, width = 1, height = 0.3)
} else { 
  ggdraw()  +
    draw_plot(band.summer.temp, x = 0, y = 0, width = 1, height = 0.3) +
    draw_plot(temp.r, x = 0, y = 0, width = 1, height = 0.3) +
    draw_plot(temp.f, x = 0, y = 0, width = 1, height = 0.3)
}

# Model and field census plot (Fig. 1C,D)
if(census == TRUE) {
  if(t.param$Habitat == "Tropical") {
    ggdraw()  +
      draw_plot(band.density, x = 0, y = 0, width = 1, height = 0.4) +
      draw_plot(model.A.r, x = 0, y = 0, width = 1, height = 0.4) +
      draw_plot(plot.census, x = 0, y = 0, width = 1, height = 0.4)
  } else {
    ggdraw()  +
      draw_plot(band.summer.density, x = 0, y = 0, width = 1, height = 0.4) +
      draw_plot(model.A.r, x = 0, y = 0, width = 1, height = 0.4) +
      draw_plot(plot.census, x = 0, y = 0, width = 1, height = 0.4) }
}

# Juvenile density plot (Fig. 1E,F)
if(t.param$Habitat == "Tropical") {
  ggdraw()  +
    draw_plot(band.density, x = 0, y = 0, width = 1, height = 0.4) +
    draw_plot(model.J.r, x = 0, y = 0, width = 1, height = 0.4) +
    draw_plot(model.J.f, x = 0, y = 0, width = 1, height = 0.4)
} else { 
  ggdraw()  +
    draw_plot(band.summer.density, x = 0, y = 0, width = 1, height = 0.4) +
    draw_plot(model.J.r, x = 0, y = 0, width = 1, height = 0.4) +
    draw_plot(model.J.f, x = 0, y = 0, width = 1, height = 0.4)
}

# Adult density plot (Fig. 1G,H)
if(t.param$Habitat == "Tropical") {
  ggdraw()  +
    draw_plot(band.density, x = 0, y = 0, width = 1, height = 0.4) +
    draw_plot(model.A.r, x = 0, y = 0, width = 1, height = 0.4) +
    draw_plot(model.A.f, x = 0, y = 0, width = 1, height = 0.4)
} else { 
  ggdraw()  +
    draw_plot(band.summer.density, x = 0, y = 0, width = 1, height = 0.4) +
    draw_plot(model.A.r, x = 0, y = 0, width = 1, height = 0.4) +
    draw_plot(model.A.f, x = 0, y = 0, width = 1, height = 0.4)
}


# STATISTICS
#linear regression of field data and model data
stats.data <- data.frame(Time = data.census$Time, Census = data.census$A, Model = NA)
for(i in seq(1:nrow(stats.data))) {
  stats.data$Model[i] = model.r[model.r$Time == stats.data$Time[i],"A"]
}
summary(lm(stats.data$Model ~ 0 + stats.data$Census))
#plot(stats.data$Census,stats.data$Model, xlim = c(0,100), ylim = c(0,100))