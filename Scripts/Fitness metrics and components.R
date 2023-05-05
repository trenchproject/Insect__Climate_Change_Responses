################################################################################
### This R script plots the fitness metrics and fitness components in Fig. 1 ###
################################################################################

# Load packages
library(tidyverse)

# Set working directory (if necessary)
#setwd() # enter working directory of main downloaded file (containing R project file)

# USER: select fitness metric or fitness component by removing "#" in front of desired response
#       (note: make sure there is a "#" in front of all other responses)
response <- "r_m"
#response <- "R0"
#response <- "development"
#response <- "survival"
#response <- "longevity"
#response <- "birth"

# Read in files
param <- as.data.frame(read_csv("Model parameters/Temperature response parameters.csv"))
r.data <- as.data.frame(read_csv("Model predictions/Predictions rm.csv"))
R0.data <- as.data.frame(read_csv("Model predictions/Predictions R0.csv"))

# Select data for Clavigralla shadabi in Benin
param <- param[1,]
r.data <- r.data[1,]
R0.data <- R0.data[1,]

# Quantify daily mean temperatures (for histograms in Fig. 1A,B)
# historical climate
temp.r <- as.data.frame(read_csv(paste0("Climate data/Recent climate data Benin.csv")))
temp.r$day <- floor(temp.r$day)
temp.r.min <- temp.r[duplicated(temp.r$day),]
temp.r.max <- temp.r[duplicated(temp.r$day, fromLast=TRUE),]
temp.r <- data.frame(temp.r.min$day, (temp.r.min$T + temp.r.max$T)/2)
names(temp.r) <- c("day", "T")
# future climate
temp.f <- as.data.frame(read_csv(paste0("Climate data/Future climate data Benin.csv")))
temp.f$day <- floor(temp.f$day)
temp.f.min <- temp.f[duplicated(temp.f$day),]
temp.f.max <- temp.f[duplicated(temp.f$day, fromLast=TRUE),]
temp.f <- data.frame(temp.f.min$day, (temp.f.min$T + temp.f.max$T)/2)
names(temp.f) <- c("day", "T")

# Plot options
Tmin <- round(min(temp.r$T,temp.f$T),0) - 3
Tmax <- round(max(temp.r$T,temp.f$T),0) + 3
temp <- seq(Tmin,Tmax,0.1)
deltaT <- 5 # tick marks for x-axis
ymin <- 0
ymax1 <- 1
if(response == "development") {ymin <- 1; ymax1 <- 4}
if(response == "longevity") {ymin <- 0; ymax1 <- 60}
ymax2 <- 0.4 # for temperature histogram
Tbin  <- 0.5 # bin size for temperature histogram

# Define temperature responses
r <- function(T) { ifelse(T <= param$Toptr, param$rMax*exp(-((T-param$Toptr)^2)/(2*param$sr^2)), param$rMax*(1 - ((T-param$Toptr)/(param$Toptr-param$Tmaxr))^2)) }
R0 <- function(T) { param$R0Topt*exp(-((T-param$ToptR0)^2)/(2*param$sR0^2)) }
dev <- function(T) { ifelse(T <= param$Toptg, 1/(param$gTR*(T/param$TR)*exp(param$Ag*(1/param$TR-1/T))/(1+exp(param$AL*(1/param$TL-1/T)))), ifelse(T <= param$Tmaxg, 1/param$gMax, 1000)) } # using 1000 instead of "Inf" for plotting purposes
long <- function(T) { 1/(param$dATR*exp(param$AdA*(1/param$TR-1/T))) }
surv <- function(T) { exp(-dev(T)*(param$dJTR*exp(param$AdJ*(1/param$TR-1/T))+0.1*exp(-param$AdJ*(1/param$TR-1/T)))) } # note: incorporated increased mortality at low temperatures for illustrative purposes
birth <- function(T) { param$bTopt*exp(-((T-param$Toptb)^2)/(2*param$sb^2)) }

# Plot TPC curves
#dev.new(width=3, height=3, unit="in")
plot(-100, xlim=c(Tmin,Tmax), ylim=c(ymin,ymax1), xaxt = "n", xlab="", ylab="", cex.axis=2) # make blank plot
axis(1, at=seq(Tmin,Tmax,deltaT), labels=seq(Tmin-273,Tmax-273,deltaT), cex.axis=2) # convert x-axis to C
if(response == "r_m") { points(temp, r(temp)/param$rMax, type="l", lwd=4, col="black") }
if(response == "R0") { points(temp, R0(temp)/param$R0Topt, type="l", lwd=4, col="black") }
if(response == "development") { points(temp, dev(temp)*param$gMax, type="l", lwd=4, col="black") }
if(response == "survival") { points(temp, surv(temp), type="l", lwd=4, col="black") }
if(response == "longevity") { points(temp, long(temp), type="l", lwd=4, col="black") }
if(response == "birth") { points(temp, birth(temp)/param$bTopt, type="l", lwd=4, col="black") }

# Quantify mean values of r_m in recent and future climates (Fig. 1A)
if(response == "r_m") {
  # find the temperature at which the TPC is equal to the mean value of r_m in the recent and future climates
  for(i in seq(Tmin,Tmax,0.1)) { if(r(i)/param$rMax >= r.data$TPC.r) { break } }
  for(j in seq(Tmax,Tmin,-0.1)) { if(r(j)/param$rMax >= r.data$TPC.f) { break } }

  # plot mean values (blue and orange point in Fig. 1A)
  points(i, r.data$TPC.r, pch=19, cex=3, col="#0072B2")
  points(j, r.data$TPC.f, pch=19, cex=3, col="#D55E00")
}

# Quantify mean values of R0 in recent and future climates (Fig. 1B)
if(response == "R0") {
  # find the temperature at which the TPC is equal to the mean value of R0 in the recent and future climates
  for(i in seq(Tmax,Tmin,-0.1)) { if(R0(i)/param$R0Topt >= R0.data$TPC.r) { break } }
  for(j in seq(Tmax,Tmin,-0.1)) { if(R0(j)/param$R0Topt >= R0.data$TPC.f) { break } }
  
  # plot mean values (blue and orange point in Fig. 1B)
  points(i, R0.data$TPC.r, pch=19, cex=3, col="#0072B2")
  points(j, R0.data$TPC.f, pch=19, cex=3, col="#D55E00")
}

# Plot temperature histograms for recent (1980-1985) and future (2095-2100) climates (Fig. 1A,B)
if(response == "r_m" || response == "R0") {
  start.r <- (1980-1941)*365-1 # start data for histogram in recent climate (note that climate data starts on Jan 2, 1941; thus, the -1)
  par(new = T)
  hist(temp.r[start.r:(start.r+5*365),"T"], xlim=c(Tmin,Tmax), ylim=c(ymin,ymax2), axes=F, xlab=NA, ylab=NA, breaks=seq(from=Tmin, to=Tmax, by=Tbin), col=rgb(0,114,178, max = 255, alpha = 80), border=rgb(0,114,178, max = 255, alpha = 80), freq=FALSE, main = NULL)
  hist(temp.f[temp.f$day>365*70,"T"], xlim=c(Tmin,Tmax), ylim=c(ymin,ymax2), breaks=seq(from=Tmin, to=Tmax, by=Tbin), ylab="r", col=rgb(213,94,0, max = 255, alpha = 80), border=rgb(213,94,0, max = 255, alpha = 80), freq=FALSE, main = NULL, add=TRUE)
  axis(side = 4, cex.axis=2)
}