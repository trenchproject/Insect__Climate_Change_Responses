#####################################################################################
#### This R script reads climate data from netCDF files and converts them to CSV ####
#####################################################################################

# Load packages
library(lubridate)
library(ncdf4)
library(tidyverse)

# Set working directory (if necessary)
#setwd() # enter working directory of main downloaded file (containing R project file)

# Read climate station data
clim.data <- read_csv("Climate data/Climate station data.csv")

# USER: Enter location of climate data (see "Climate station data.csv")
location <- "Benin"

# USER: Save climate data?
save <- FALSE

# USER: Save remove climate netCDF files?
remove <- FALSE

# Find information for selected location in "Climate station data.csv"
sp.num <- -1
for(i in 1:nrow(clim.data)) {
  if(clim.data[i,]$Location == location) {
    sp.num <- i
    break }
  }
if(sp.num == -1) {
  print("Location not found in Climate station data.csv, please update location")
  }

# Assign day number for the first date in which there is climate data
start.date <- yday(paste0(clim.data[sp.num,]$Start_yr,"-",
                    if(clim.data[sp.num,]$Start_mo < 10) {"0"}, clim.data[sp.num,]$Start_mo,"-", # add 0 before months < 10
                    if(clim.data[sp.num,]$Start_day < 10) {"0"}, clim.data[sp.num,]$Start_day)) - 1 # add 0 before days < 10


################################## RECENT CLIMATE DATA ##################################
# Open netCDF files
# NOTE: Must have first downloaded climate nc files, please read "ReadME download recent data.docx"
nc.max <- nc_open(paste0("Climate data/Recent Tmax ",location,".nc"))
nc.min <- nc_open(paste0("Climate data/Recent Tmin ",location,".nc"))

# Get variables and data from netCDF files
Tmax.time <- ncvar_get(nc.max, "time") # time in Tmax data frame
Tmin.time <- ncvar_get(nc.min, "time") # time in Tmin data frame
Tmax <- ncvar_get(nc.max, "TMAX")
Tmin <- ncvar_get(nc.min, "TMIN")
if(is.na(Tmax[1]) == FALSE) { data.max <- data.frame(start.date + Tmax.time + 0.5,  273.15 + Tmax) # add start.date to Tmax.time and offset Tmax by 0.5 days from Tmin (e.g., on day 10, the temperature is at Tmin at day 10.0 and at Tmax at day 10.5)
} else { data.max <- data.frame(Tmax.time + 0.5,  273.15 + Tmax) } # offset Tmax.time by 0.5 days from Tmin.time (e.g., on day 10, the temperature is at Tmin at day 10.0 and at Tmax at day 10.5)
if(is.na(Tmin[1]) == FALSE) { data.min <- data.frame(start.date + Tmin.time,  273.15 + Tmin) # add start.date to Tmin.time
} else { data.min <- data.frame(Tmin.time,  273.15 + Tmin) }
names(data.max) <- c("day", "T")
names(data.min) <- c("day", "T")

# Enter climate data into R data frame
data <- rbind(data.min, data.max)
data <- data[order(data$day),]

# Remove rows with NA
data <- na.omit(data)

# Save data in CSV file
if(save == TRUE) { write.csv(data, paste0("Climate Data/Recent climate data ",location,".csv"), row.names = FALSE) }

# Close and delete netCDF files
nc_close(nc.max)
nc_close(nc.min)
if(remove == TRUE) {
  file.remove(paste0("Climate Data/Recent Tmax ",location,".nc"))
  file.remove(paste0("Climate Data/Recent Tmin ",location,".nc"))
}


#################################### FUTURE CLIMATE DATA ####################################
# Get variables and data from CSV files created by climate data.py
data.max <- as.data.frame(read_csv(paste0("Climate Data/Future Tmax ",location,".csv")))
data.min <- as.data.frame(read_csv(paste0("Climate Data/Future Tmin ",location,".csv")))
names(data.max) <- c("day", "time", "latitude", "longitude", "T")
names(data.min) <- c("day", "time", "latitude", "longitude", "T")

# Remove any data before Jan 1, 2025 and combine data.max and data.min
for(i in nrow(data.max):0) {
  if(i == 0) { break } # end for loop if i gets to zero (i.e., no data to remove)
  if(year(data.max[i,]$time) < 2025) { break }  # find row number of Dec 31, 2024
}
for(j in nrow(data.min):0) {
  if(j == 0) { break } # end for loop if i gets to zero (i.e., no data to remove)
  if(year(data.min[j,]$time) < 2025) { break }  # find row number of Dec 31, 2024
}
if(i != 0) { data.max <- data.max[-(1:i),] } # remove data before row Jan 1, 2025
if(j != 0) { data.min <- data.min[-(1:j),] } # remove data before row Jan 1, 2025
data.max$day <- data.max$day - i + 0.5  # set first day to zero and offset Tmax 0.5 days from Tmin
data.min$day <- data.min$day - i  # set first day to zero

# Combine data into R data frame
data <- rbind(data.min, data.max)
data <- data[order(data$day),]

# Remove rows with NA
data <- na.omit(data)

# Save data in CSV file
if(save == TRUE) { write.csv(data, paste0("Climate Data/Future climate data ",location,".csv"), row.names = FALSE) }

# Remove Tmax and Tmin CSV files
if(remove == TRUE) {
  file.remove(paste0("Climate Data/Future Tmax ",location,".csv"))
  file.remove(paste0("Climate Data/Future Tmin ",location,".csv"))
}