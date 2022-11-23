#####################################################################################
#### This R script reads climate data from netCDF files and converts them to CSV ####
#####################################################################################

# Load packages and set working directory
library(readxl)
library(lubridate)
# need to check if libraries below are needed in this script
library(tidyr)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(ncdf4)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Read climate station datasheet and present user with list of species
clim.data <- read_xlsx("Climate station data extension.xlsx")
clim.data$Species

# USER: enter species and start date YYYY-MM-DD (see "Climate station data.xlxs")
found <- F
while(found == F) {
  sp.num <- readline("Please enter number associated with species and location: ") # must run lines 17-18 to see options
  #sp.num <- 1 # this line is only for checking the code, remove when everything is working
  if(sp.num >= 1 & sp.num <= nrow(clim.data)) {
    found <- T
  } else { print("Species not found, please try again")}
}
sp.num <- as.numeric(sp.num)
print(paste("Species selected is: ", sp.num, clim.data[sp.num,]$Species))

# Set location (loc) and first day of year (date) for climate record for user-selected species above
loc <- clim.data[sp.num,]$Location
date <- yday(paste0(clim.data[sp.num,]$Start_yr,"-",
                    if(clim.data[sp.num,]$Start_mo < 10) {"0"}, clim.data[sp.num,]$Start_mo,"-", # add 0 before months < 10
                    if(clim.data[sp.num,]$Start_day < 10) {"0"}, clim.data[sp.num,]$Start_day)) # add 0 before days < 10


################################## HISTORICAL CLIMATE DATA ##################################
# Open netCDF files
###### NOTE: Must first run "Read download climate nc files, please read "Project overview.docx" under the "Protocols" section
nc.max <- nc_open(paste0("Historical Tmax ",loc,".nc"))
nc.min <- nc_open(paste0("Historical Tmin ",loc,".nc"))

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

# Close and delete netCDF files
nc_close(nc.max)
nc_close(nc.min)
file.remove(paste0("Historical Tmax ",loc,".nc"))
file.remove(paste0("Historical Tmin ",loc,".nc"))


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



