#####################################################################################
#### This R script reads climate data from netCDF files and converts them to CSV ####
#####################################################################################

# Load packages
library(readxl)
library(lubridate)
library(ncdf4)
library(tidyverse)

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('..')

# Read climate station data
# NOTE: To analyze populations beyond the manuscript, update line 17 to "Extensions beyond manuscript/Climate station data extension.xlsx"
clim.data <- read_xlsx("Climate data/Climate station data.xlsx")

# USER: Enter location of climate data (see "Climate station data.xlsx")
loc <- "Benin"

# USER: Save climate data?
save <- FALSE

# USER: Save remove climate netCDF files?
remove <- FALSE

# Find information for selected species in "Climate station data.xlsx"
sp.num <- -1
for(i in 1:nrow(clim.data)) {
  if(clim.data[i,]$Location == loc) {
    sp.num <- i
    break }
  }
if(sp.num == -1) {
  print("Location not found in Climate station data.xlsx, please update line 20")
  }

# Assign day number for the first date in which there is climate data
date <- yday(paste0(clim.data[sp.num,]$Start_yr,"-",
                    if(clim.data[sp.num,]$Start_mo < 10) {"0"}, clim.data[sp.num,]$Start_mo,"-", # add 0 before months < 10
                    if(clim.data[sp.num,]$Start_day < 10) {"0"}, clim.data[sp.num,]$Start_day)) # add 0 before days < 10


################################## HISTORICAL CLIMATE DATA ##################################
# Open netCDF files
# NOTE: Must first run "Read download climate nc files, please read "ReadME download historical data.docx"
nc.max <- nc_open(paste0("Climate data/Historical Tmax ",loc,".nc"))
nc.min <- nc_open(paste0("Climate data/Historical Tmin ",loc,".nc"))

# Get variables and data from netCDF files
t.max <- ncvar_get(nc.max, "time") # time in Tmax data frame
t.min <- ncvar_get(nc.min, "time") # time in Tmin data frame
Tmax <- ncvar_get(nc.max, "TMAX")
Tmin <- ncvar_get(nc.min, "TMIN")
data.max <- data.frame(date + t.max + 0.5,  273.15 + Tmax) # offset Tmax 0.5 days from Tmin
data.min <- data.frame(date + t.min,  273.15 + Tmin)
names(data.max) <- c("day", "T")
names(data.min) <- c("day", "T")

# Enter climate data into R data frame
data <- rbind(data.min, data.max)
data <- data[order(data$day),]

# Remove rows with NA
data <- na.omit(data)

# Save data in CSV file
if(save == TRUE) { write.csv(data, paste0("Climate Data/Historical climate data ",loc,".csv"), row.names = FALSE) }

# Close and delete netCDF files
nc_close(nc.max)
nc_close(nc.min)
if(remove == TRUE) {
  file.remove(paste0("Climate Data/Historical Tmax ",loc,".nc"))
  file.remove(paste0("Climate Data/Historical Tmin ",loc,".nc"))
}


#################################### FUTURE CLIMATE DATA ####################################
# Get variables and data from CSV files created by climate data.py
data.max <- as.data.frame(read_csv(paste0("Climate Data/Future Tmax ",loc,".csv")))
data.min <- as.data.frame(read_csv(paste0("Climate Data/Future Tmin ",loc,".csv")))
names(data.max) <- c("day", "time", "latitude", "longitude", "T")
names(data.min) <- c("day", "time", "latitude", "longitude", "T")
data.max$day <- data.max$day + 0.5  # offset Tmax 0.5 days from Tmin

# Combine data into R data frame
data <- rbind(data.min, data.max)
data <- data[order(data$day),]

# Remove rows with NA
data <- na.omit(data)

# Save data in CSV file
if(save == TRUE) { write.csv(data, paste0("Climate Data/Future climate data ",loc,".csv"), row.names = FALSE) }

# Remove Tmax and Tmin CSV files
if(remove == TRUE) {
  file.remove(paste0("Climate Data/Future Tmax ",loc,".csv"))
  file.remove(paste0("Climate Data/Future Tmin ",loc,".csv"))
}