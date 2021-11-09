#####################################################################################
#### This R script reads climate data from netCDF files and converts them to CSV ####
#####################################################################################

# Load packages and set working directory
library(tidyr)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(ncdf4)
library(lubridate)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# USER: enter location and start date YYYY-MM-DD (see "Climate station data.xlxs")
loc <- "Benin"
date <- yday("1941-01-02")


################################## HISTORICAL CLIMATE DATA ##################################
# Open netCDF files
nc.max <- nc_open(paste0("Historical Tmax ",loc,".nc"))
nc.min <- nc_open(paste0("Historical Tmin ",loc,".nc"))

# Save metadata
#{
#  sink(paste0("Historical data ",loc,".txt"))
#  print(nc_data)
#  sink()
#}

# Get variables and data from netCDF files
t.max <- ncvar_get(nc.max, "time") # time in Tmax data frame
t.min <- ncvar_get(nc.min, "time") # time in Tmin data frame
Tmax <- ncvar_get(nc.max, "TMAX")
Tmin <- ncvar_get(nc.min, "TMIN")
data.max <- data.frame(date + t.max + 0.5,  273.15 + Tmax) # offset Tmax 0.5 days from Tmin
data.min <- data.frame(date + t.min,  273.15 + Tmin)
names(data.max) <- c("day", "T")
names(data.min) <- c("day", "T")

# Enter data into R data frame
data <- as.data.frame(0) # create empty data frame
data <- data[FALSE,]
data <- rbind(data.min, data.max)
data <- data[order(data$day),]

# Remove rows with NA
data <- na.omit(data)

# Save data in CSV file
write.csv(data, paste0("Historical climate data ",loc,".csv"), row.names = FALSE)

# Close netCDF files
nc_close(nc.max)
nc_close(nc.min)


#################################### FUTURE CLIMATE DATA ####################################
# Get variables and data from CSV files created by climate data.py
data.max <- as.data.frame(read_csv(paste0("Future Tmax ",loc,".csv")))
data.min <- as.data.frame(read_csv(paste0("Future Tmin ",loc,".csv")))
names(data.max) <- c("day", "time", "latitude", "longitude", "T")
names(data.min) <- c("day", "time", "latitude", "longitude", "T")
data.max$day <- data.max$day + 0.5  # offset Tmax 0.5 days from Tmin

# Combine data into R data frame
data <- as.data.frame(0) # create empty data frame
data <- data[FALSE,]
data <- rbind(data.min, data.max)
data <- data[order(data$day),]

# Remove rows with NA
data <- na.omit(data)

# Save data in CSV file
write.csv(data, paste0("Future climate data ",loc,".csv"), row.names = FALSE)

# Remove Tmax and Tmin CSV files
file.remove(paste0("Future Tmax ",loc,".csv"))
file.remove(paste0("Future Tmin ",loc,".csv"))



# # SUPPLEMENTARY CODES
# # Get data
# data <- as.data.frame(read_csv(paste0("Historical climate data ",loc,".csv")))
# 
# # Visualize data
# hist(data$T)
# shapiro.test(data[1:5000,]$T) # Is distribution significantly different from normal?
# 
# 
# # Interpolated function of data (not working)
# T <- approxfun(data[1:730,], rule = 2) # rule = 2 sets any extrapolated point to the closest data extreme
# 
# # Plot data and interpolated function
# xmin <- 0
# xmax <- 730
# ymin <- round(min(data$T),0)
# ymax <- round(max(data$T),0)+1
# ggplot(data, aes(x=day, y=T)) +
#    geom_point(size=0.8, color="red") +
#    geom_function(fun = T, size=0.8, color="black") +
#    labs(x="Time (days)", y="Temperature") +
#    scale_x_continuous(limits=c(xmin, xmax)) +
#    scale_y_continuous(limits=c(ymin, ymax)) +
#    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#          panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
#          axis.line = element_line(colour = "black"), legend.position = "none",
#          axis.text = element_text(size=13), axis.title = element_text(size=20))
# 
# 
# # Fourier transformation (under-estimating amplitude b/c there are so many Fourier terms)
# s <- data$T # temperature series
# n <- length(s)
# 
# # de-trend series
# s.mean <- mean(s)
# s.dt <- s - s.mean
# 
# # Fourier transform
# f <- fft(s.dt)
# amp <- 2*abs(f)/n
# phase <- Arg(f)
# plot(amp, xlim=c(0,400)) # dominant frequencies at 2 and 66
# 
# # build function
# T <- function(t) { s.mean + amp[2]*cos(2*pi*t+phase[2]) + amp[65]*cos(2*pi*t/365+phase[66]) }
# 
# # Plot data and Fourier function
# xmin <- 0
# xmax <- 730
# ymin <- round(min(data$T),0)
# ymax <- round(max(data$T),0)+1
# ggplot(data, aes(x=day, y=T)) +
#   geom_point(size=0.8, color="red") +
#   geom_function(fun = T, size=0.8, color="black") +
#   labs(x="Time (days)", y="Temperature") +
#   scale_x_continuous(limits=c(xmin, xmax)) +
#   scale_y_continuous(limits=c(ymin, ymax)) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_rect(fill="transparent"), plot.background = element_rect(fill="transparent"),
#         axis.line = element_line(colour = "black"), legend.position = "none",
#         axis.text = element_text(size=13), axis.title = element_text(size=20))