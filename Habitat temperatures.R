#############################################################################
#### This R script fits sinusoidal functions to habitat temperature data ####
#############################################################################


# Load packages and set working directory
library(tidyr)
library(ggplot2)
library(dplyr)
library(tidyverse)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Read data
read.data <- read_csv("Habitat Temperatures.csv")
temp.data <- as.data.frame(read.data)

# Select an insect by removing # in front of name and placing # in front of other species
sp.data <- subset(temp.data, Species == "Clavigralla shadabi")
#sp.data <- subset(temp.data, Species == "Clavigralla tomentosicollis Benin")

# Remove columns that do not contain temperature data
sp.data <- sp.data[-1:-7]
rownames(sp.data) <- c("T_K")

# Make data table with time and mean monthly temperature
data <- matrix(nrow = length(sp.data), ncol = 1)
for(i in 1:length(data)) {data[i] <- 30*(i-1)+15} # mid-point of month
colnames(data) <- c("day")
data <- cbind(data, t(sp.data)) # mean monthly temperature data

# Remove blank cells
data <- data[rowSums(is.na(data)) == 0,]
data <- as.data.frame(data)

# Fit sinusoidal function to habitat temperature data
fit <- nls(T_K ~ meanT + amplT*sin(2*pi*(day + shiftT)/365), data = data,
           start = list(meanT = 293, amplT = 2, shiftT = 30))
summary(fit)
