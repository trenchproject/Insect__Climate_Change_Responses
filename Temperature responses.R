###################################################################
#### This R script fits functions to temperature response data ####
###################################################################


# Load packages and set working directory
library(tidyr)
library(ggplot2)
library(dplyr)
library(tidyverse)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Read data
read.data <- read_csv("Habitat temperatures.csv")
temp.data <- as.data.frame(read.data)

# Select an insect by removing # in front of name and placing # in front of other species
#sp.data <- subset(temp.data, Species == "Clavigralla shadabi")
#sp.data <- subset(temp.data, Species == "Clavigralla tomentosicollis Benin")
sp.data <- subset(temp.data, Species == "Clavigralla tomentosicollis Nigeria A")

# Remove columns that do not contain temperature data
sp.data <- sp.data[-1:-7]
rownames(sp.data) <- c("T_K")