################################################################################
########## This R script fit functions to temperature response data ###########
################################################################################

# Load packages
library(tidyverse)

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('..')

# USER: enter species name and location (used in "temperature response data.csv") OR set "all" to TRUE to run analysis for all species
name <- "Clavigralla shadabi Benin"
all <- TRUE

# USER: Save model fit?
save <- FALSE

# Read data and model parameters file
data <- as.data.frame(read_csv("Biological data/Temperature response data.csv"))
params <- as.data.frame(read_csv("Model parameters/Temperature response parameters.csv"))
  
# Run analyses for all species
for(s in 3:3) { #nrow(params)) {

  # Find model parameters for selected species (line 13) or move iteratively through rows of “Temperature response parameters.csv”
  if(all == FALSE) {
    sp.num <- -1
    for(i in 1:nrow(params)) {
      if(params[i,]$Species == name) {
        sp.num <- i
        break }
    }
    if(sp.num == -1) {
      print("Location not found in Habitat temperature parameters.csv, please update line 13")
    }
  } else { sp.num <- s }
  
  # Assign species
  sp.data <- data[data$Species == params[sp.num,1],]
  
  # Remove columns that do not contain temperature data
  sp.data <- sp.data[-c(1:8,12,14,16,18,20,21,23,24,26,27,29,31,32,34,36,37,38)]
  
  # Set minimum, maximum, and reference temperature (TR) for nls and plotting
  if(params[sp.num,]$Habitat == "Tropical") {
    Tmin <- 285
    Tmax <- 315
  } else {
    Tmin <- 275
    Tmax <- 305
  }
  if(params[sp.num,]$Species == "Clavigralla tomentosicollis Burkina Faso") { TR <- 300 # must use a higher TR for this species
  } else { TR <- 293 }
  params[sp.num,]$TR <- TR
  
  
  ###################### INTRINSIC GROWTH RATE, r_m (Eq. 1a) #####################
  # Nonlinear regression
  error <- tryCatch({
    r.nls <- nls(r ~ ifelse(T_K <= Toptr, rMax*exp(-(T_K-Toptr)^2/(2*sr^2)), rMax*(1 - ((T_K-Toptr)/(Toptr-Tmaxr))^2)),
           data=sp.data, start=list(rMax=params[sp.num,]$rMax, Toptr=params[sp.num,]$Toptr, Tmaxr=params[sp.num,]$Tmaxr, sr=params[sp.num,]$sr))
    print(summary(r.nls))
    
    # Function nls does not return an error
    error <- FALSE
  }, error = function(err) { error <- TRUE; return(error) })
  
  # If there are no errors, then assign parameters
  if(error == FALSE) {
    params[sp.num,]$rMax <- round(coef(r.nls)[1], 3)
    params[sp.num,]$Toptr <- round(coef(r.nls)[2], 1)
    params[sp.num,]$Tmaxr <- round(coef(r.nls)[3], 1)
    params[sp.num,]$sr <- round(coef(r.nls)[4], 2)
  } else { # If nls returns an error, then try alternative method of splitting data and fitting each part of the piecewise function separately
    # # Split data at maximum value (which goes to both piecewise parts)
    # for(i in nrow(sp.data):1) {
    #     if(sp.data[i,"r"] == max(sp.data[,"r"])) { break }
    # }
    # data1 <- sp.data[1:i,]
    # data2 <- sp.data[i:nrow(sp.data),]
    # 
    # # Estimate first part of piecewise function
    # r.nls1 <- nls(r ~ rMax*exp(-(T_K-Toptr)^2/(2*sr^2)), data=data1,
    #              start=list(rMax=params[sp.num,]$rMax, Toptr=params[sp.num,]$Toptr, sr=params[sp.num,]$sr))
    # summary(r.nls1)
    # 
    # # Set Topt and rMax (NOTE: Topt cannot equal Tmax in nls)
    # #rMax <- max(sp.data[,"r"])
    # #Toptr <- sp.data[sp.data$r==rMax,"T_K"]
    # 
    # # Estimate all other parameters
    # r.nls <- nls(r ~ ifelse(T_K <= Toptr, rMax*exp(-(T_K-Toptr)^2/(2*sr^2)), rMax*(1 - ((T_K-Toptr)/(Toptr-Tmaxr))^2)),
    #          data=sp.data, start=list(Tmaxr=params[sp.num,]$Tmaxr, sr=params[sp.num,]$sr))
    # summary(r.nls)
    # 
    # # Assign parameters
    # params[sp.num,]$rMax <- round(coef(r.nls)[1], 3)
    # params[sp.num,]$Toptr <- round(coef(r.nls)[2], 1)
    # params[sp.num,]$Tmaxr <- round(coef(r.nls)[3], 1)
    # params[sp.num,]$sr <- round(coef(r.nls)[4], 2)
  }
  
  # Plot model fit
  plot(sp.data$T_K, sp.data$r)
  points(seq(Tmin,Tmax,1), ifelse(seq(Tmin,Tmax,1) <= params[sp.num,]$Toptr, params[sp.num,]$rMax*exp(-(seq(Tmin,Tmax,1)-params[sp.num,]$Toptr)^2/(2*params[sp.num,]$sr^2)),
                                  params[sp.num,]$rMax*(1 - ((seq(Tmin,Tmax,1)-params[sp.num,]$Toptr)/(params[sp.num,]$Toptr-params[sp.num,]$Tmaxr))^2)), type="l", col="blue")
  
  
  ##################### NET REPRODUCTIVE RATE, R0 (Eq. 1b) #######################
  # Nonlinear regression
  R0.nls <- nls(R0 ~ R0Topt*exp(-(T_K-ToptR0)^2/(2*sR0^2)), data=sp.data,
                start=list(R0Topt=params[sp.num,]$R0Topt, ToptR0=params[sp.num,]$ToptR0, sR0=params[sp.num,]$sR0))
  summary(R0.nls)
  
  # Assign parameters
  params[sp.num,]$R0Topt <- round(coef(R0.nls)[1], 1)
  params[sp.num,]$ToptR0 <- round(coef(R0.nls)[2], 1)
  params[sp.num,]$sR0 <- round(coef(R0.nls)[3], 2)
  
  # Plot model fit
  plot(sp.data$T_K, sp.data$R0)
  points(seq(Tmin,Tmax,1), params[sp.num,]$R0Topt*exp(-(seq(Tmin,Tmax,1)-params[sp.num,]$ToptR0)^2/(2*params[sp.num,]$sR0^2)), type="l", col="blue")
  
  
  ################### PER CAPITA BIRTH RATE, b[T] (Eq. 2a) #######################
  # Nonlinear regression
  b.nls<- nls(Birth_Rate ~ bTopt*exp(-(T_K-Toptb)^2/(2*sb^2)), data=sp.data,
             start=list(bTopt=params[sp.num,]$bTopt, Toptb=params[sp.num,]$Toptb, sb=params[sp.num,]$sb))
  summary(b.nls)
  
  # Assign parameters
  params[sp.num,]$bTopt <- round(coef(b.nls)[1], 2)
  params[sp.num,]$Toptb <- round(coef(b.nls)[2], 1)
  params[sp.num,]$sb <- round(coef(b.nls)[3], 2)
  
  # Plot model fit
  plot(sp.data$T_K, sp.data$Birth_Rate)
  points(seq(Tmin,Tmax,1), params[sp.num,]$bTopt*exp(-(seq(Tmin,Tmax,1)-params[sp.num,]$Toptb)^2/(2*params[sp.num,]$sb^2)), type="l", col="blue")
  
  
  ################ STAGE-SPECIFIC MORTALITY RATE, d_i[T] (Eq. 2b) ################
  # Juvenile per capita mortality rate
  # Nonlinear regression
  error <- tryCatch({
    dJ.nls <- nls(Juv_Mortality ~ dJTR*exp(AdJ*(1/TR-1/T_K)), data=sp.data,
                  start=list(dJTR=params[sp.num,]$dJTR, AdJ=params[sp.num,]$AdJ))
    print(summary(dJ.nls))
  
    # Function nls does not return an error
    error <- FALSE
    
    # If mortality estimate has high p-value (> 0.6), then return error
    if(dJ.nls[["coefficients"]][1,4] > 0.6) { error <- TRUE }
  }, error = function(err) { error <- TRUE; return(error) })
  
  # If there are no errors, then assign parameters
  if(error == FALSE) {
    params[sp.num,]$dJTR <- round(coef(dJ.nls)[1], 4)
    params[sp.num,]$AdJ <- round(coef(dJ.nls)[2], 0)
  } else { # If nls returns an error, then try alternative method of setting dJTR to the data at TR and fitting AdJ
    dJTR <- sp.data[sp.data$T_K==TR,"Juv_Mortality"]
  
    # Nonlinear regression
    dJ.nls <- nls(Juv_Mortality ~ dJTR*exp(AdJ*(1/TR-1/T_K)), data=sp.data,
                  start=list(AdJ=params[sp.num,]$AdJ))
    print(summary(dJ.nls))
    
    # Assign parameters
    params[sp.num,]$dJTR <- round(dJTR, 4)
    params[sp.num,]$AdJ <- round(coef(dJ.nls)[1], 0)
  }
  
  # Plot model fit
  plot(sp.data$T_K, sp.data$Juv_Mortality)
  points(seq(Tmin,Tmax,1), params[sp.num,]$dJTR*exp(params[sp.num,]$AdJ*(1/TR-1/seq(Tmin,Tmax,1))), type="l", col="blue")
  
  
  # Adult per capita mortality rate
  # Nonlinear regression
  error <- tryCatch({
    dA.nls <- nls(Adult_Mortality ~ dATR*exp(AdA*(1/TR-1/T_K)), data=sp.data,
                  start=list(dATR=params[sp.num,]$dATR, AdA=params[sp.num,]$AdA))
    print(summary(dA.nls))
    
    # Function nls does not return an error
    error <- FALSE
  }, error = function(err) { error <- TRUE; return(error) })
  
  # If there are no errors, then assign parameters
  if(error == FALSE) {
    params[sp.num,]$dATR <- round(coef(dA.nls)[1], 4)
    params[sp.num,]$AdA <- round(coef(dA.nls)[2], 0)
  } else { # If nls returns an error, then try alternative method of setting dATR to the data at TR and fitting AdA
    dATR <- sp.data[sp.data$T_K==TR,"Adult_Mortality"]
    
    # Nonlinear regression
    dA.nls <- nls(Adult_Mortality ~ dATR*exp(AdA*(1/TR-1/T_K)), data=sp.data,
                  start=list(AdA=params[sp.num,]$AdA))
    print(summary(dA.nls))
    
    # Assign parameters
    params[sp.num,]$dATR <- round(dATR, 4)
    params[sp.num,]$AdA <- round(coef(dA.nls)[1], 0)
  }
  
  # Plot model fit
  plot(sp.data$T_K, sp.data$Adult_Mortality)
  points(seq(Tmin,Tmax,1), params[sp.num,]$dATR*exp(params[sp.num,]$AdA*(1/TR-1/seq(Tmin,Tmax,1))), type="l", col="blue")
  
  
  ####################### DEVELOPMENT RATE, g[T] (Eq. 2c) ########################
  # Quantify gMax, Toptg, Tmaxg directly from data
  # gMax
  gMax <- max(sp.data[,"Development"])
  # Toptg
  for(i in nrow(sp.data):1) {
       if(sp.data[i,"Development"] == gMax) { break }
    }
  Toptg <- sp.data[i,"T_K"]
  # Tmaxg
  Tmaxg <- sp.data[nrow(sp.data),"T_K"]
  for(j in nrow(sp.data):i) {
    if(sp.data[j,"Development"] == 0) { Tmaxg <- sp.data[j,"T_K"] }
  }
  
  # Assign Toptg and Tmaxg parameters
  params[sp.num,]$Toptg <- round(Toptg, 1)
  params[sp.num,]$Tmaxg <- round(Tmaxg, 1)
  
  # Try to fit remaining parameters to data (excluding data beyond Toptg)
  error1 <- tryCatch({
    g.nls <- nls(Development ~ gTR*(T_K/TR)*exp(Ag*(1/TR-1/T_K))/(1+exp(AL*(1/TL-1/T_K))), data=sp.data[1:i,],
                              start=list(gTR=params[sp.num,]$gTR, Ag=params[sp.num,]$Ag, AL=params[sp.num,]$AL, TL=params[sp.num,]$TL))
    print(summary(g.nls))
    
    # Function nls does not return an error
    error1 <- FALSE
  }, error = function(err) { error1 <- TRUE; return(error1) })
  
  # If there are no errors, then assign parameters
  if(error1 == FALSE) {
    params[sp.num,]$gTR <- round(coef(g.nls)[1], 3)
    params[sp.num,]$Ag <- round(coef(g.nls)[2], 0)
    params[sp.num,]$AL <- round(coef(g.nls)[3], 0)
    params[sp.num,]$TL <- round(coef(g.nls)[4], 1)
    params[sp.num,]$gMax <- round(params[sp.num,]$gTR*(params[sp.num,]$Toptg/TR)*exp(params[sp.num,]$Ag*(1/TR-1/params[sp.num,]$Toptg))/
                                    (1+exp(params[sp.num,]$AL*(1/params[sp.num,]$TL-1/params[sp.num,]$Toptg))), 3)
  } else { # If nls returns an error, then try alternative method of first fitting gTR and Ag and then fitting AL and TL afterwards
    error2 <- tryCatch({
      # Nonlinear regression for gTR and Ag
      g.nls1 <- nls(Development ~ gTR*(T_K/TR)*exp(Ag*(1/TR-1/T_K)), data=sp.data[1:i,],
                   start=list(gTR=params[sp.num,]$gTR, Ag=params[sp.num,]$Ag))
      print(summary(g.nls1))
      
      # Assign parameters for next nonlinear regression
      gTR.test <- coef(g.nls1)[1]
      Ag.test <- coef(g.nls1)[2]
      
      # Nonlinear regression for AL and TL
      g.nls2 <- nls(Development ~ gTR.test*(T_K/TR)*exp(Ag.test*(1/TR-1/T_K))/(1+exp(AL*(1/TL-1/T_K))), data=sp.data[1:i,],
                   start=list(AL=params[sp.num,]$AL, TL=params[sp.num,]$TL))
      print(summary(g.nls2))
      
      # Function nls does not return an error
      error2 <- FALSE
    }, error = function(err) { error2 <- TRUE; return(error2) })
    
    # If there are no errors in the alternative method above, then assign parameters
    if(error2 == FALSE) {  
      # Assign parameters
      params[sp.num,]$gTR <- round(gTR.test, 3)
      params[sp.num,]$Ag <- round(Ag.test, 0)
      params[sp.num,]$AL <- round(coef(g.nls2)[1], 0)
      params[sp.num,]$TL <- round(coef(g.nls2)[2], 1)
      params[sp.num,]$gMax <- round(params[sp.num,]$gTR*(params[sp.num,]$Toptg/TR)*exp(params[sp.num,]$Ag*(1/TR-1/params[sp.num,]$Toptg))/
                                      (1+exp(params[sp.num,]$AL*(1/params[sp.num,]$TL-1/params[sp.num,]$Toptg))), 3)
    } else { # If nls still returns an error, then try alternative method of first fitting gTR and Ag and then fitting AL and then TL
      # Nonlinear regression for gTR and Ag
      g.nls1 <- nls(Development ~ gTR*(T_K/TR)*exp(Ag*(1/TR-1/T_K)), data=sp.data[1:i,],
                    start=list(gTR=params[sp.num,]$gTR, Ag=params[sp.num,]$Ag))
      print(summary(g.nls1))
      
      # Assign parameters for next nonlinear regression
      gTR.test <- coef(g.nls1)[1]
      Ag.test <- coef(g.nls1)[2]
      
      # Set TL to just below the coldest laboratory temperature and nonlinear regression for AL
      TL.test <- min(sp.data[,"T_K"]) - 1
      g.nls2 <- nls(Development ~ gTR.test*(T_K/TR)*exp(Ag.test*(1/TR-1/T_K))/(1+exp(AL*(1/TL.test-1/T_K))), data=sp.data[1:i,],
                    start=list(AL=params[sp.num,]$AL))
      print(summary(g.nls2))
      
      # Assign AL and nonlinear regression for TL
      AL.test <- coef(g.nls2)[1]
      g.nls3 <- nls(Development ~ gTR.test*(T_K/TR)*exp(Ag.test*(1/TR-1/T_K))/(1+exp(AL.test*(1/TL-1/T_K))), data=sp.data[1:i,],
                    start=list(TL=params[sp.num,]$TL))
      print(summary(g.nls3))
      
      # Assign parameters
      params[sp.num,]$gTR <- round(gTR.test, 3)
      params[sp.num,]$Ag <- round(Ag.test, 0)
      params[sp.num,]$AL <- round(AL.test, 0)
      params[sp.num,]$TL <- round(coef(g.nls3)[1], 1)
      params[sp.num,]$gMax <- round(params[sp.num,]$gTR*(params[sp.num,]$Toptg/TR)*exp(params[sp.num,]$Ag*(1/TR-1/params[sp.num,]$Toptg))/
                                      (1+exp(params[sp.num,]$AL*(1/params[sp.num,]$TL-1/params[sp.num,]$Toptg))), 3)
    }
  }

  # plot model fit
  plot(sp.data$T_K, sp.data$Development, xlim=c(Tmin,Tmax))
  g.funct <- function(T) { ifelse(T <= params[sp.num,]$Toptg, params[sp.num,]$gTR*(T/TR)*exp(params[sp.num,]$Ag*(1/TR-1/T))/
                                    (1+exp(params[sp.num,]$AL*(1/params[sp.num,]$TL-1/T))),
                                  ifelse(T <= params[sp.num,]$Tmaxg, params[sp.num,]$gMax, 0)) }
  points(seq(Tmin,Tmax,0.1), g.funct(seq(Tmin,Tmax,0.1)), type="l", col="blue")
  
  
  # Minimum developmental temperature (for overwintering in DDE model)
  Tmin.nls <- nls(Development ~ m*T_K+b, data=sp.data[1:i,], start=list(m=params[sp.num,]$gTR, b=0))
  summary(Tmin.nls)
  
  # Plot model fit
  plot(sp.data$T_K, sp.data$Development, xlim=c(Tmin,Tmax))
  points(seq(Tmin,Tmax,1), coef(Tmin.nls)[1]*seq(Tmin,Tmax,1)+coef(Tmin.nls)[2], type="l", col="blue")
  
  # Calculate minimum development temperature (Tmin)
  (Tmin <- (-coef(Tmin.nls)[2]/coef(Tmin.nls)[1])[[1]])
  params[sp.num,]$Tmin <- round(Tmin, 1)
  
  # Break for loop (line 23) if analyses are run for a specified species (all <- FALSE in line 14)
  if(all == FALSE) { break }
}

# Save model parameters to "Temperature response parameters.csv" (if desired)
if(save == TRUE) { write.csv(params, "Model parameters/Temperature response parameters.csv", row.names = FALSE) }

# Print model fit
if(all == FALSE) { params[sp.num,]
} else { params }