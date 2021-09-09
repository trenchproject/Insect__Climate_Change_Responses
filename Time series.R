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


# READ IN TEMPERATURE RESPONSE PARAMETERS AND TIME-SERIES DATA
# Temperature response parameters
data <- as.data.frame(read_csv("Temperature response parameters.csv"))

# Time-series data
# Select an file by removing # in front of name and placing # in front of other files
#data.density <- read_csv("Population data Nigeria.csv")
#data.density <- read_csv("Population data China.csv")

# Select time-series data
#data.TS <- subset(data.density, Plot=="B") # select plot A, B, or C
#data.TS <- subset(data.density, location=="Dafeng" & species=="Apolygus_lucorum") # select location and species
#data.TS <- subset(data.density, location=="Dafeng" & species=="Adelphocoris_suturalis") # select location and species
#data.TS <- subset(data.density, location=="Langfang" & species=="Apolygus_lucorum") # select location and species
#data.TS <- subset(data.density, location=="Xinxiang" & species=="Adelphocoris_suturalis") # select location and species

# Data transformation (if needed)
# Convert from log to linear scale
#data.TS$J <- (10^data.TS$J) - 1
#data.TS$J_SE <- sqrt(data.TS$J)
#data.TS$A <- (10^data.TS$A) - 1
#data.TS$A_SE <- sqrt(data.TS$A)


# READ IN MODEL OUTPUT
# From DDE population dynamics.py 
# Select an insect by removing # in front of name and placing # in front of other species
#data.model <- as.data.frame(read_csv("Time Series Clavigralla shadabi Benin.csv"))
#sp.data <- subset(data, Species == "Clavigralla shadabi Benin")
#data.model <- as.data.frame(read_csv("Time Series Clavigralla tomentosicollis Benin.csv"))
#sp.data <- subset(data, Species == "Clavigralla tomentosicollis Benin")
#data.model <- as.data.frame(read_csv("Time Series Clavigralla tomentosicollis Nigeria B.csv"))
#sp.data <- subset(data, Species == "Clavigralla tomentosicollis Nigeria B")
#data.model <- as.data.frame(read_csv("Time Series Clavigralla tomentosicollis Burkina Faso.csv"))
#sp.data <- subset(data, Species == "Clavigralla tomentosicollis Burkina Faso")
data.model <- as.data.frame(read_csv("Time Series Apolygus lucorum China Dafeng.csv"))
sp.data <- subset(data, Species == "Apolygus lucorum China Dafeng")
#data.model <- as.data.frame(read_csv("Time Series Adelphocoris suturalis China Dafeng.csv"))
#sp.data <- subset(data, Species == "Adelphocoris suturalis China Dafeng")
#data.model <- as.data.frame(read_csv("Time Series Apolygus lucorum China Langfang.csv"))
#sp.data <- subset(data, Species == "Apolygus lucorum China Langfang")
#data.model <- as.data.frame(read_csv("Time Series Adelphocoris suturalis China Xinxiang.csv"))
#sp.data <- subset(data, Species == "Adelphocoris suturalis China Xinxiang")
#data.model <- as.data.frame(read_csv("Time Series Macrosiphum euphorbiae Brazil.csv"))
#sp.data <- subset(data, Species == "Macrosiphum euphorbiae Brazil")
#data.model <- as.data.frame(read_csv("Time Series Aulacorthum solani Brazil.csv"))
#sp.data <- subset(data, Species == "Aulacorthum solani Brazil")
#data.model <- as.data.frame(read_csv("Time Series Uroleucon ambrosiae Brazil.csv"))
#sp.data <- subset(data, Species == "Uroleucon ambrosiae Brazil")
#data.model <- as.data.frame(read_csv("Time Series Lygus lineolaris Mississippi.csv"))
#sp.data <- subset(data, Species == "Lygus lineolaris Mississippi")
#data.model <- as.data.frame(read_csv("Time Series Pilophorus typicus Japan.csv"))
#sp.data <- subset(data, Species == "Pilophorus typicus Japan")
#data.model <- as.data.frame(read_csv("Time Series Macrolophus pygmaeus on Myzus persicae Greece.csv"))
#sp.data <- subset(data, Species == "Macrolophus pygmaeus on Myzus persicae Greece")
#data.model <- as.data.frame(read_csv("Time Series Macrolophus pygmaeus on Trialeurodes vaporariorum Greece.csv"))
#sp.data <- subset(data, Species == "Macrolophus pygmaeus on Trialeurodes vaporariorum Greece")


# Population dynamics with climate change
data.model.CC <- data.model


# SET PLOT OPTIONS
# Default: plot last 2 year of model
# for historical time period
xmin <- 0
xmax <- 720
ymin <- 0
ymax <- 1200
# for climate change time period
xmin.CC <- 0
xmax.CC <- 720
ymin.CC <- 0
ymax.CC <- 1200
# for temperature function
temp.min <- 275
temp.max <- 310
yr <- 360 # days in a year (using 360 for simplicity)
init_yrs <- 10 # number of years to initiate the model (from Python DDE model)
TS.length <- xmax - xmin # length of time-series data
end <- nrow(data.model)


# FORMAT MODEL OUTPUT TO ALIGN WITH TIME-SERIES DATA
# Remove all rows before time-series data starts
# for Clavigralla tomentosicollis in Nigeria
#data.model <- data.model[c(-1:-((1972-1961+init_yrs)*yr + xmin)), ]
# for Apolygus lucorum in China
#data.model <- data.model[c(-1:-((2005-1961+init_yrs)*yr + xmin)), ]
# for other species
data.model <- data.model[c(-1:-((2020-1961+init_yrs)*yr + xmin)), ]

# Remove all rows after xmax days
data.model <- data.model[c(-(xmax+1):-end), ]

# climate change period (remove all but last 2 years of data)
data.model.CC <- data.model.CC[c(-1:-(end-2*yr + xmin.CC)), ]

# Re-scale time to start at xmin
# historical period
time.shift <- data.model[[1,1]]-xmin
data.model <- sweep(data.model, 2, c(time.shift,0,0,0,0))
# climate change period
time.shift.CC <- data.model.CC[[1,1]] + xmin.CC
data.model.CC <- sweep(data.model.CC, 2, c(time.shift.CC,0,0,0,0))

# Data transformation (if needed)
# Convert from linear to log scale
#data.model$J <- log(data.model$J, 10)
#data.model$A <- log(data.model$A, 10)
#data.model.CC$J <- log(data.model.CC$J, 10)
#data.model.CC$A <- log(data.model.CC$A, 10)



# PLOT TIME-SERIES DATA
# Juvenile density
plot.J = ggplot(data.TS, aes(x=time, y=J, ymin=J-J_SE, ymax=J+J_SE)) + 
  geom_pointrange(size=0.8, color="#d1495b") +
  #geom_line(size=0.8, linetype="longdash", color="#d1495b") +
  labs(x="Time", y="Density") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(ymin, ymax)) +
  #scale_y_log10(limits=c(ymin, ymax)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20))

# Adult density
plot.A = ggplot(data.TS, aes(x=time, y=A, ymin=A-A_SE, ymax=A+A_SE)) + 
  geom_pointrange(size=0.8, color="#30638e") +
  #geom_line(size=0.8, linetype="longdash", color="#30638e") +
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
# Juvenile density
model.J = ggplot(data.model, aes(x=Time, y=J)) + 
  geom_line(size=0.8, color="#d1495b") +
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
  geom_line(size=0.8, color="#30638e") +
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
# Juvenile density
model.J.CC = ggplot(data.model.CC, aes(x=Time, y=J)) + 
  geom_line(size=0.8, color="#d1495b", linetype="longdash") +
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
  geom_line(size=0.8, color="#30638e", linetype="longdash") +
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
plot.temp <- ggplot() +
  geom_function(fun = function(t) (sp.data$meanT + sp.data$delta_mean*(t+time.shift))  + (sp.data$amplT + sp.data$delta_ampl*(t+time.shift)) * sin(2*pi*((t+time.shift) + sp.data$shiftT)/yr),
                size=0.8, color="red") +
  #geom_function(fun = function(t) (299.08 + 0.000127*(t+time.shift) + 1.84*cos(2*pi*30/(30*yr)*((t+time.shift)-15) - 1.86) + 1.51*cos(2*pi*60/(30*yr)*((t+time.shift)-15) - 3.09) + 0.359*cos(2*pi*90/(30*yr)*((t+time.shift)-15) - 2.57)),
  #              size=0.8, color="red") +
  geom_function(fun = function(t) (sp.data$Tmin), size=0.8, color="black") +
  labs(x="Time", y="T (K)") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(temp.min, temp.max)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20))

# climate change time period
plot.temp.CC <- ggplot() +
  geom_function(fun = function(t) (sp.data$meanT + sp.data$delta_mean*(t+time.shift.CC))  + (sp.data$amplT + sp.data$delta_ampl*(t+time.shift.CC)) * sin(2*pi*((t+time.shift.CC) + sp.data$shiftT)/yr),
                size=0.8, linetype="longdash", color="red") +
  #geom_function(fun = function(t) (299.08 + 0.000127*(t+time.shift.CC) + 1.84*cos(2*pi*30/(30*yr)*((t+time.shift.CC)-15) - 1.86) + 1.51*cos(2*pi*60/(30*yr)*((t+time.shift.CC)-15) - 3.09) + 0.359*cos(2*pi*90/(30*yr)*((t+time.shift.CC)-15) - 2.57)),
  #              size=0.8, linetype="longdash", color="red") +
  geom_function(fun = function(t) (sp.data$Tmin), size=0.8, color="black") +
  labs(x="Time", y="T (K)") +
  scale_x_continuous(limits=c(xmin, xmax)) +
  scale_y_continuous(limits=c(temp.min, temp.max)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour = "black"), legend.position = "none", 
        axis.text = element_text(size=13), axis.title = element_text(size=20))


# DRAW FINAL PLOTS
# Historical time period
# Juveniles and adults
plot <- ggdraw()  +
  draw_plot(plot.temp, x = 0, y = 0, width = 1, height = 0.3) +
  #draw_plot(plot.J, x = 0, y = 0.3, width = 1, height = 0.7) +
  draw_plot(plot.A, x = 0, y = 0.3, width = 1, height = 0.7) +
  #draw_plot(model.J, x = 0, y = 0.3, width = 1, height = 0.7) +
  draw_plot(model.A, x = 0, y = 0.3, width = 1, height = 0.7)
plot

# Total insects
#plot <- ggdraw()  +
#  draw_plot(plot.temp, x = 0, y = 0, width = 1, height = 0.3) +
#  draw_plot(plot.A, x = 0, y = 0.3, width = 1, height = 0.7) +
#  draw_plot(model.J, x = 0, y = 0.3, width = 1, height = 0.7) +
#  draw_plot(model.A, x = 0, y = 0.3, width = 1, height = 0.7) +
#  draw_plot(model.I, x = 0, y = 0.3, width = 1, height = 0.7)
#plot

# Climate change time period
# Juveniles and adults
#plot.CC <- ggdraw()  +
#  draw_plot(plot.temp.CC, x = 0, y = 0, width = 1, height = 0.3) +
#  draw_plot(model.J.CC, x = 0, y = 0.3, width = 1, height = 0.7) +
#  draw_plot(model.A.CC, x = 0, y = 0.3, width = 1, height = 0.7)
#plot.CC

# Total insects
#plot.CC <- ggdraw()  +
#  draw_plot(plot.temp.CC, x = 0, y = 0, width = 1, height = 0.3) +
#  draw_plot(model.J.CC, x = 0, y = 0.3, width = 1, height = 0.7) +
#  draw_plot(model.A.CC, x = 0, y = 0.3, width = 1, height = 0.7) +
#  draw_plot(model.I.CC, x = 0, y = 0.3, width = 1, height = 0.7)
#plot.CC

# Compare historical and climate change periods
# Juveniles and adults
plot.compare <- ggdraw()  +
  draw_plot(plot.temp, x = 0, y = 0, width = 1, height = 0.3) +
  draw_plot(plot.temp.CC, x = 0, y = 0, width = 1, height = 0.3) +
  draw_plot(model.J, x = 0, y = 0.3, width = 1, height = 0.7) +
  draw_plot(model.A, x = 0, y = 0.3, width = 1, height = 0.7) +
  draw_plot(model.J.CC, x = 0, y = 0.3, width = 1, height = 0.7) +
  draw_plot(model.A.CC, x = 0, y = 0.3, width = 1, height = 0.7)
plot.compare


# Total insects
#plot.compare <- ggdraw()  +
#  draw_plot(plot.temp, x = 0, y = 0, width = 1, height = 0.3) +
#  draw_plot(plot.temp.CC, x = 0, y = 0, width = 1, height = 0.3) +
#  draw_plot(model.I, x = 0, y = 0.3, width = 1, height = 0.7) +
#  draw_plot(model.I.CC, x = 0, y = 0.3, width = 1, height = 0.7)
#plot.compare





# CALCULATE CHANGES BETWEEN HISTORICAL AND CLIMATE CHANGE PERIODS
# mean density
mean.J <- mean(data.model[["J"]])
mean.J.CC <- mean(data.model.CC[["J"]])
mean.A <- mean(data.model[["A"]])
mean.A.CC <- mean(data.model.CC[["A"]])
d.mean.J <- (mean.J.CC-mean.J)/mean.J
d.mean.A <- (mean.A.CC-mean.A)/mean.A

# peak density
max.J <- max(data.model[["J"]])
max.J.CC <- max(data.model.CC[["J"]])
max.A <- max(data.model[["A"]])
max.A.CC <- max(data.model.CC[["A"]])
d.max.J <- (max.J.CC-max.J)/max.J
d.max.A <- (max.A.CC-max.A)/max.A

# peak phenology
time.J <- subset(data.model, J==max.J)$Time%%yr
time.J.CC <- subset(data.model.CC, J==max.J.CC)$Time%%yr
time.A <- subset(data.model, A==max.A)$Time%%yr
time.A.CC <- subset(data.model.CC, A==max.A.CC)$Time%%yr
d.time.J <- (time.J.CC-time.J)
d.time.A <- (time.A.CC-time.A)

# minimum density
min.J <- 0
min.J.CC <- 0
min.A <- 0
min.A.CC <- 0
if(min(data.model[["J"]]) > 0) { min.J <- round(min(data.model[["J"]]), digits=2) }
if(min(data.model.CC[["J"]]) > 0) { min.J.CC <- round(min(data.model.CC[["J"]]), digits=2) }
if(min(data.model[["A"]]) > 0) { min.A <- round(min(data.model[["A"]]), digits=2) }
if(min(data.model.CC[["A"]]) > 0) { min.A.CC <- round(min(data.model.CC[["A"]]), digits=2) }
if(min.J != 0) { d.min.J <- (min.J.CC-min.J)/min.J } else {d.min.J <- 0}
if(min.A != 0) { d.min.A <- (min.A.CC-min.A)/min.A } else {d.min.A <- 0}


# temperature functions
temp <- function(t) (sp.data$meanT + sp.data$delta_mean*(t+time.shift))  + (sp.data$amplT + sp.data$delta_ampl*(t+time.shift)) * sin(2*pi*((t+time.shift) + sp.data$shiftT)/yr)
temp.CC <- function(t) (sp.data$meanT + sp.data$delta_mean*(t+time.shift.CC))  + (sp.data$amplT + sp.data$delta_ampl*(t+time.shift.CC)) * sin(2*pi*((t+time.shift.CC) + sp.data$shiftT)/yr)
length <- 360 # length of time over which to compare models

# calculate activity period (T(t) > Tmin)
for(i in 0:length) {
  if(temp(i) > sp.data["Tmin"]) {print(i)
    i <- length} }

# reproductive activity period (birth rate > 0.1 bTopt)
b.period <- 0
b.period.CC <- 0
for(i in 0:length) {
  if(sp.data["bTopt"]*exp(-((temp(i)-sp.data["Toptb"])^2)/(2*sp.data["sb"]^2)) > 0.1*sp.data["bTopt"]) {b.period <- b.period + 1} }
for(i in 0:length) {
  if(sp.data["bTopt"]*exp(-((temp.CC(i)-sp.data["Toptb"])^2)/(2*sp.data["sb"]^2)) > 0.1*sp.data["bTopt"]) {b.period.CC <- b.period.CC + 1} }
d.b <- (b.period.CC-b.period)/b.period

# developmental activity period (development rate > 0.25 Mmax or > dMin from Jahansson et al. 2020)
# calculate Mmax and Topt
Mmax <- 0
Topt <- 0
for(i in 0:100) {
  T <- coef(r.T)[1] + i*(sp.data["TH"]-sp.data["TL"])/100
  M.T <- sp.data["mTR"]*(T/sp.data["TR"])*exp(sp.data["AmJ"]*(1/sp.data["TR"]-1/T))/(1+sp.data["skew"]*exp(sp.data["AL"]*(1/sp.data["TL"]-1/T))+exp(sp.data["AH"]*(1/sp.data["TH"]-1/T)))
  if(M.T > Mmax) {Mmax <- M.T[[1]]
  Topt <- T[[1]]}}
# calculate time above 0.25 Mmax
m.period <- 0
m.period.CC <- 0
for(i in 0:length) {
  if(sp.data["mTR"]*(temp(i)/sp.data["TR"])*exp(sp.data["AmJ"]*(1/sp.data["TR"]-1/temp(i)))/(1+sp.data["skew"]*exp(sp.data["AL"]*(1/sp.data["TL"]-1/temp(i)))+exp(sp.data["AH"]*(1/sp.data["TH"]-1/temp(i)))) > 0.25*Mmax) {m.period = m.period + 1} }
for(i in 0:length) {
  if(sp.data["mTR"]*(temp.CC(i)/sp.data["TR"])*exp(sp.data["AmJ"]*(1/sp.data["TR"]-1/temp.CC(i)))/(1+sp.data["skew"]*exp(sp.data["AL"]*(1/sp.data["TL"]-1/temp.CC(i)))+exp(sp.data["AH"]*(1/sp.data["TH"]-1/temp.CC(i)))) > 0.25*Mmax) {m.period.CC = m.period.CC + 1} }
d.m <- (m.period.CC-m.period)/m.period

# average reproductive rate (mean b(T))
b.sum <- 0
b.sum.CC <- 0
for(i in 0:length) {
  b.sum <- b.sum + (sp.data["bTopt"]*exp(-((temp(i)-sp.data["Toptb"])^2)/(2*sp.data["sb"]^2)))[[1]] }
for(i in 0:length) {
  b.sum.CC <- b.sum.CC + (sp.data["bTopt"]*exp(-((temp.CC(i)-sp.data["Toptb"])^2)/(2*sp.data["sb"]^2)))[[1]] }
d.b.ave <- (b.sum.CC-b.sum)/b.sum

# average development rate (mean m(T))
m.sum <- 0
m.sum.CC <- 0
for(i in 0:length) {
  m.sum <- m.sum + (sp.data["mTR"]*(temp(i)/sp.data["TR"])*exp(sp.data["AmJ"]*(1/sp.data["TR"]-1/temp(i)))/(1+sp.data["skew"]*exp(sp.data["AL"]*(1/sp.data["TL"]-1/temp(i)))+exp(sp.data["AH"]*(1/sp.data["TH"]-1/temp(i)))))[[1]] }
for(i in 0:length) {
  m.sum.CC <- m.sum.CC + (sp.data["mTR"]*(temp.CC(i)/sp.data["TR"])*exp(sp.data["AmJ"]*(1/sp.data["TR"]-1/temp.CC(i)))/(1+sp.data["skew"]*exp(sp.data["AL"]*(1/sp.data["TL"]-1/temp.CC(i)))+exp(sp.data["AH"]*(1/sp.data["TH"]-1/temp.CC(i)))))[[1]] }
d.m.ave <- (m.sum.CC-m.sum)/m.sum

# average adult mortality rate (mean dA(T))
dA.sum <- 0
dA.sum.CC <- 0
for(i in 0:length) {
  dA.sum <- dA.sum + (sp.data["dATR"]*exp(sp.data["AdA"]*(1/sp.data["TR"]-1/temp(i))))[[1]] }
for(i in 0:length) {
  dA.sum.CC <- dA.sum.CC + (sp.data["dATR"]*exp(sp.data["AdA"]*(1/sp.data["TR"]-1/temp.CC(i))))[[1]] }
d.dA.ave <- (dA.sum.CC-dA.sum)/dA.sum

# time above optimum range for reproductive rate (T > ToptR0 + sR0)
R0.period <- 0
R0.period.CC <- 0
for(i in 0:length) {
  if(temp(i) > sp.data["ToptR0"] + sp.data["sR0"]) {R0.period <- R0.period + 1} }
for(i in 0:length) {
  if(temp.CC(i) > sp.data["ToptR0"] + sp.data["sR0"]) {R0.period.CC <- R0.period.CC + 1} }
if(R0.period !=0) { d.R0 <- (R0.period.CC-R0.period)/R0.period } else {d.R0 <- 0 }

# time above rTmax (T > rTmax)
r.period <- 0
r.period.CC <- 0
for(i in 0:length) {
  if(temp(i) > sp.data["rTmax"]) {r.period <- r.period + 1} }
for(i in 0:length) {
  if(temp.CC(i) > sp.data["rTmax"]) {r.period.CC <- r.period.CC + 1} }
if(r.period !=0) { d.r <- (r.period.CC-r.period)/r.period } else {d.r <- 0 }

# mean thermal safety margin (Toptr - T)
TSM <- 0
TSM.CC <- 0
for(i in 0:length) {
  TSM <- TSM + (sp.data["Toptr"] - temp(i)) }
for(i in 0:length) {
  TSM.CC <- TSM.CC + (sp.data["Toptr"] - temp.CC(i)) }
TSM <- (TSM/length)[[1]]
TSM.CC <- (TSM.CC/length)[[1]]
if(TSM !=0) { d.TSM <- (TSM.CC-TSM)/TSM } else {d.TSM <- 0 }


# PLOTS
par(mfrow=c(3,5))
# density metrics
barplot(c(d.mean.J,d.mean.A), col=c("#d1495b","#30638e"), ylim=c(-0.6,0.6), main=expression("Mean density"))
barplot(c(d.max.J,d.max.A), col=c("#d1495b","#30638e"), ylim=c(-0.6,0.6), main=expression("Peak density"))
barplot(c(d.time.J,d.time.A), col=c("#d1495b","#30638e"), ylim=c(-380,20), main=expression("Timing of peak"))
barplot(c(d.min.J,d.min.A), col=c("#d1495b","#30638e"), ylim=c(-1,0.2), main=expression("Minimum density"))
# activity periods
barplot(d.b, col="#30638e", xlim=c(0.2,2), ylim=c(-0.5,0.5), main=expression("Rep. period"))
barplot(d.m, col="#d1495b", xlim=c(0.2,2), ylim=c(-0.5,0.5), main=expression("Dev. period"))
# average life history traits
barplot(d.b.ave, col="#30638e", xlim=c(0.2,2), ylim=c(-0.5,0.5), main=expression("Mean b(T)"))
barplot(d.m.ave, col="#d1495b", xlim=c(0.2,2), ylim=c(-0.5,0.5), main=expression("Mean m(T)"))
barplot(d.dA.ave, col="#30638e", xlim=c(0.2,2), ylim=c(-0.5,0.5), main=expression("Mean dA(T)"))
# thermal performance curves
barplot(c(R0.period,R0.period.CC), col=c("purple","red"), ylim=c(0,100), main=expression("above R0"))
barplot(c(r.period,r.period.CC), col=c("purple","red"), ylim=c(0,100), main=expression("above r"))
barplot(c(TSM,TSM.CC), col=c("purple","red"), ylim=c(-20,20), main=expression("Thermal margin"))
# habitat temperatures
barplot((sp.data["ext_meanT"])[[1]], col="#30638e", xlim=c(0.2,2), ylim=c(0,30), main=expression("meanT"))
barplot((sp.data["ext_amplT"])[[1]], col="#30638e", xlim=c(0.2,2), ylim=c(0,30), main=expression("amplT"))
barplot((sp.data["ext_mT_aT"])[[1]], col="#30638e", xlim=c(0.2,2), ylim=c(0,30), main=expression("meanT + amplT"))
par(mfrow=c(1,1))
