################################################################################
########## This R script fit functions to temperature response data ###########
################################################################################

# Load packages
library(tidyverse)
library(stringr)

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('..')


##### NOTE: HAVING TROUBLE GETTING PARAMETER-HEAVY DEVELOPMENT RATE TPC TO CONVERGE


# USER: enter species name and location OR set "all" to TRUE to run analysis for all species
name <- "Clavigralla shadabi"
location <- "Benin"
name <- paste(name, location)
all <- FALSE

# USER: Save model fit?
save <- FALSE

# Read data and model parameters file
data <- as.data.frame(read_csv("Biological data/Temperature response data.csv"))
params <- as.data.frame(read_csv("Model parameters/Temperature response parameters V2.csv"))
  
# Run analyses for all species
x <- 1
for(s in x:2) { #nrow(params)) {

  # Find model parameters for selected species or move iteratively through rows of “Temperature response parameters.csv”
  if(all == FALSE) {
    sp.num <- -1
    for(i in 1:nrow(params)) {
      if(params[i,]$Species == name) {
        sp.num <- i
        break }
    }
    if(sp.num == -1) {
      print("Location not found in Habitat temperature parameters.csv, please update species name")
    }
  } else { sp.num <- s }
  
  # Assign species (need to set location for Apolygus lucorum and Adelphocoris saturalis to "China" as in "Temperature response data.csv")
  if(word(params[sp.num,2],1) != "China") { sp.data <- data[data$Species == params[sp.num,1],] # find population in "Temperature response data.csv"
  } else { sp.data <- data[data$Species == word(params[sp.num,1], 1,3),] }
  
  # Remove columns that do not contain temperature data
  sp.data <- sp.data[-c(1:8,12,14,16,18,20,21,23,24,26,27,29,31,32,34,36,37,38)]
  
  # Set minimum and maximum for x-axes, and reference temperature (TR) for nls
  if(params[sp.num,]$Habitat == "Tropical") {
    Xmin <- 285
    Xmax <- 315
  } else {
    Xmin <- 275
    Xmax <- 305
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
  } else { # If nls returns an error, then set rMax and Toptr to data and use nls to estimate Tmaxr and sr
    rMax.test <- max(na.omit(sp.data[,"r"]))
    Toptr.test <- sp.data[which.max(sp.data$r), "T_K"]
    r.nls1 <- nls(r ~ ifelse(T_K <= Toptr.test, rMax.test*exp(-(T_K-Toptr.test)^2/(2*sr^2)), rMax.test*(1 - ((T_K-Toptr.test)/(Toptr.test-Tmaxr))^2)),
                 data=sp.data, start=list(Tmaxr=params[sp.num,]$Tmaxr, sr=params[sp.num,]$sr))
    summary(r.nls1)
    
    # Assign parameters
    params[sp.num,]$rMax <- round(rMax.test, 3)
    params[sp.num,]$Toptr <- round(Toptr.test, 1)
    params[sp.num,]$Tmaxr <- round(coef(r.nls1)[1], 1)
    params[sp.num,]$sr <- round(coef(r.nls1)[2], 2)
  }
  
  # Plot model fit
  plot(sp.data$T_K, sp.data$r)
  points(seq(Xmin,Xmax,1), ifelse(seq(Xmin,Xmax,1) <= params[sp.num,]$Toptr, params[sp.num,]$rMax*exp(-(seq(Xmin,Xmax,1)-params[sp.num,]$Toptr)^2/(2*params[sp.num,]$sr^2)),
                                  params[sp.num,]$rMax*(1 - ((seq(Xmin,Xmax,1)-params[sp.num,]$Toptr)/(params[sp.num,]$Toptr-params[sp.num,]$Tmaxr))^2)), type="l", col="blue")
  
  
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
  points(seq(Xmin,Xmax,1), params[sp.num,]$R0Topt*exp(-(seq(Xmin,Xmax,1)-params[sp.num,]$ToptR0)^2/(2*params[sp.num,]$sR0^2)), type="l", col="blue")
  
  
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
  points(seq(Xmin,Xmax,1), params[sp.num,]$bTopt*exp(-(seq(Xmin,Xmax,1)-params[sp.num,]$Toptb)^2/(2*params[sp.num,]$sb^2)), type="l", col="blue")
  
  
  ################ STAGE-SPECIFIC MORTALITY RATE, d_i[T] (Eq. 2b) ################
  # Juvenile per capita mortality rate
  # Nonlinear regression
  if(params[sp.num,2] != "Australia Acton") { 
    dJ.nls <- nls(Juv_Mortality ~ dJTR*exp(AdJ*(1/TR-1/T_K)), data=sp.data, start=list(dJTR=params[sp.num,]$dJTR, AdJ=params[sp.num,]$AdJ))
    print(summary(dJ.nls))
    
    # Assign parameters
    params[sp.num,]$dJTR <- round(coef(dJ.nls)[1], 4)
    params[sp.num,]$AdJ <- round(coef(dJ.nls)[2], 0)
  } else { dJTR.test <- round(sp.data[1,"Juv_Mortality"],4) # For population in Australia Acton, must set dJTR to data and fit AdJ via nls
    dJ.nls <- nls(Juv_Mortality ~ dJTR.test*exp(AdJ*(1/TR-1/T_K)), data=sp.data, start=list(AdJ=params[sp.num,]$AdJ))
    print(summary(dJ.nls))
  
    # Assign parameters
    params[sp.num,]$dJTR <- round(dJTR.test, 4)
    params[sp.num,]$AdJ <- round(coef(dJ.nls)[1], 0)
  }
  
  # Plot model fit
  plot(sp.data$T_K, sp.data$Juv_Mortality)
  points(seq(Xmin,Xmax,1), params[sp.num,]$dJTR*exp(params[sp.num,]$AdJ*(1/TR-1/seq(Xmin,Xmax,1))), type="l", col="blue")
  
  
  # Adult per capita mortality rate
  # Nonlinear regression
  dA.nls <- nls(Adult_Mortality ~ dATR*exp(AdA*(1/TR-1/T_K)), data=sp.data, start=list(dATR=params[sp.num,]$dATR, AdA=params[sp.num,]$AdA))
  print(summary(dA.nls))
  
  # Assign parameters
  params[sp.num,]$dATR <- round(coef(dA.nls)[1], 4)
  params[sp.num,]$AdA <- round(coef(dA.nls)[2], 0)
  
  # Plot model fit
  plot(sp.data$T_K, sp.data$Adult_Mortality)
  points(seq(Xmin,Xmax,1), params[sp.num,]$dATR*exp(params[sp.num,]$AdA*(1/TR-1/seq(Xmin,Xmax,1))), type="l", col="blue")
  
  
  ####################### DEVELOPMENT RATE, g[T] (Eq. 2c) ########################
  # Quantify Toptg directly from data (descend down temperature treatments until the maximum development rate is reached, then exit the for loop)
  for(i in nrow(sp.data):1) { if(sp.data[i,"Development"] == max(na.omit(sp.data[,"Development"]))) { break } }
  if(word(params[sp.num,2],1) == "Brazil") { Toptg <- sp.data[nrow(sp.data),"T_K"] # For Brazil populations, Topt is set to maximum laboratory temperature
  } else { Toptg <- sp.data[i,"T_K"] }
  
  # Quantify Tmaxg directly from data (set Tmaxg above the highest temperature treatment, then descend until the development rate is above zero and exit the for loop)
  Tmaxg <- sp.data[nrow(sp.data),"T_K"] + 2 # set Tmaxg to 2 degrees above highest temperature as default 
  for(j in nrow(sp.data):i) { if(sp.data[j,"Development"] == 0) { Tmaxg <- sp.data[j,"T_K"] } }

  # For Brazil populations and Apolygus lucorum, set i index to maximum laboratory temperature for nls below
  if(word(params[sp.num,2],1) == "Brazil" || word(params[sp.num,1],1,2) == "Apolygus lucorum") {
    for(i in nrow(sp.data):1) { if(sp.data[i,"Development"] != 0) { break } } # NOTE: for Brazil populations, i = 3 for parameters used in the manuscript
  }

  # Assign Toptg and Tmaxg parameters
  params[sp.num,]$Toptg <- round(Toptg, 1)
  params[sp.num,]$Tmaxg <- round(Tmaxg, 1)
  
  # NOTE: It is currently not possible to fit remaining 6 parameters in development rate TPC simultaneously;
  #       therefore, gTR and Ag are first fit to exponentially increasing portion of TPC and then other parameters are fit
  # Nonlinear regression for gTR and Ag
  g.nls <- nls(Development ~ gTR*(T_K/TR)*exp(Ag*(1/TR-1/T_K)), data=sp.data[1:i,],
                start=list(gTR=params[sp.num,]$gTR, Ag=params[sp.num,]$Ag))
  print(summary(g.nls))
  
  # Assign parameters for next nonlinear regression and update params
  gTR.test <- coef(g.nls)[1]
  Ag.test <- coef(g.nls)[2]
  params[sp.num,]$gTR <- round(gTR.test, 3)
  params[sp.num,]$Ag <- round(Ag.test, 0)
  
  # Try to fit remaining parameters to data
  error1 <- tryCatch({
    g.nls1 <- nls(Development ~ gTR.test*(T_K/TR)*exp(Ag.test*(1/TR-1/T_K))/(1+exp(AL*(1/TL-1/T_K))+exp(AH*(1/TH-1/T_K))), data=sp.data,
                              start=list(AL=params[sp.num,]$AL, TL=params[sp.num,]$TL, AH=params[sp.num,]$AH, TH=params[sp.num,]$TH))
    print(summary(g.nls1))
    
    # Function nls does not return an error
    error1 <- FALSE
  }, error = function(err) { error1 <- TRUE; return(error1) })
  
  # If there are no errors, then assign parameters
  if(error1 == FALSE) {
    params[sp.num,]$AL <- round(coef(g.nls1)[1], 0)
    params[sp.num,]$TL <- round(coef(g.nls1)[2], 1)
    params[sp.num,]$AH <- round(coef(g.nls1)[3], 0)
    params[sp.num,]$TH <- round(coef(g.nls1)[4], 1)
    params[sp.num,]$gMax <- round(params[sp.num,]$gTR*(params[sp.num,]$Toptg/TR)*exp(params[sp.num,]$Ag*(1/TR-1/params[sp.num,]$Toptg))/
                                    (1+exp(params[sp.num,]$AL*(1/params[sp.num,]$TL-1/params[sp.num,]$Toptg))), 3)
  } else { # If nls returns an error, then try setting TL and TH, fitting AL and AH, and then fitting TL and AH given fitted estimates of AL and AH
    # Set starting values for TL and TH
    TL.test <- params[sp.num,]$TL - 1 # NOTE: using a 1C "buffer" below start values helps nls convergence
    TH.test <- params[sp.num,]$TH
    
    # Nonlinear regression for AL and AH
    g.nls2 <- nls(Development ~ gTR.test*(T_K/TR)*exp(Ag.test*(1/TR-1/T_K))/(1+exp(AL*(1/TL.test-1/T_K))+exp(AH*(1/TH.test-1/T_K))), data=sp.data,
                  start=list(AL=params[sp.num,]$AL, AH=params[sp.num,]$AH))
    print(summary(g.nls2))
    
    # Assign parameters for next nonlinear regression
    AL.test <- coef(g.nls2)[1]
    AH.test <- coef(g.nls2)[2]
    params[sp.num,]$AL <- round(coef(g.nls2)[1], 0)
    params[sp.num,]$AH <- round(coef(g.nls2)[2], 0)
    
    error2 <- tryCatch({
      # Nonlinear regression for TL and TH
      g.nls3 <- nls(Development ~ gTR.test*(T_K/TR)*exp(Ag.test*(1/TR-1/T_K))/(1+exp(AL.test*(1/TL-1/T_K))+exp(AH.test*(1/TH-1/T_K))), data=sp.data,
                   #start=list(TL=TL.test, TH=TH.test))
                   start=list(TL=params[sp.num,]$TL, TH=params[sp.num,]$TH))
      print(summary(g.nls3))
      
      # Function nls does not return an error
      error2 <- FALSE
    }, error = function(err) { error2 <- TRUE; return(error2) })
    
    # If there are no errors in the alternative method above, then assign parameters
    if(error2 == FALSE) {  
      # Assign parameters
      params[sp.num,]$TL <- round(coef(g.nls3)[1], 1)
      params[sp.num,]$TH <- round(coef(g.nls3)[2], 1)
      params[sp.num,]$gMax <- round(params[sp.num,]$gTR*(params[sp.num,]$Toptg/TR)*exp(params[sp.num,]$Ag*(1/TR-1/params[sp.num,]$Toptg))/
                                      (1+exp(params[sp.num,]$AL*(1/params[sp.num,]$TL-1/params[sp.num,]$Toptg))), 3)
    } else { # Fit TL and TH separately
      # Nonlinear regression for TL
      TH.test <- params[sp.num,]$TH
      g.nls3 <- nls(Development ~ gTR.test*(T_K/TR)*exp(Ag.test*(1/TR-1/T_K))/(1+exp(AL.test*(1/TL-1/T_K))+exp(AH.test*(1/TH.test-1/T_K))), data=sp.data,
                    start=list(TL=params[sp.num,]$TL))
      print(summary(g.nls3))
      
      # Assign TL
      params[sp.num,]$TL <- round(coef(g.nls3)[1], 1)
      TL.test <- round(coef(g.nls3)[1], 1)
      
      # Nonlinear regression for TH
      g.nls4 <- nls(Development ~ gTR.test*(T_K/TR)*exp(Ag.test*(1/TR-1/T_K))/(1+exp(AL.test*(1/TL.test-1/T_K))+exp(AH.test*(1/TH-1/T_K))), data=sp.data,
                    start=list(TH=params[sp.num,]$TH))
      print(summary(g.nls4))
      
      # Assign AH and gMax
      params[sp.num,]$AH <- round(coef(g.nls4)[1], 1)
      params[sp.num,]$gMax <- round(params[sp.num,]$gTR*(params[sp.num,]$Toptg/TR)*exp(params[sp.num,]$Ag*(1/TR-1/params[sp.num,]$Toptg))/
                                      (1+exp(params[sp.num,]$AL*(1/params[sp.num,]$TL-1/params[sp.num,]$Toptg))), 3)
      
    }
  }

  # plot model fit
  plot(sp.data$T_K, sp.data$Development, xlim=c(Xmin,Tmaxg + 2), ylim=c(0, 1.5*round(params[sp.num,]$gMax, 2)))
  g.funct <- function(T) { ifelse(T <= params[sp.num,]$Toptg, params[sp.num,]$gTR*(T/TR)*exp(params[sp.num,]$Ag*(1/TR-1/T))/
                                    (1+exp(params[sp.num,]$AL*(1/params[sp.num,]$TL-1/T))),
                                  ifelse(T <= params[sp.num,]$Tmaxg, params[sp.num,]$gMax, 0)) }
  points(seq(Xmin,Tmaxg + 2,0.1), g.funct(seq(Xmin,Tmaxg + 2,0.1)), type="l", col="blue")
  
  
  # Minimum developmental temperature (for overwintering in DDE model)
  if(word(params[sp.num,2],1) == "Brazil") { i <- 3 } # For the Brazil populations, exclude high temperatures at which development rate declines

  # Linear regression of linear portion of development rate TPC
  Tmin.nls <- nls(Development ~ m*T_K+b, data=sp.data[1:i,], start=list(m=params[sp.num,]$gTR, b=0))
  summary(Tmin.nls)
  
  # Plot model fit
  plot(sp.data$T_K, sp.data$Development, xlim=c(Xmin,Xmax))
  points(seq(Xmin,Xmax,1), coef(Tmin.nls)[1]*seq(Xmin,Xmax,1)+coef(Tmin.nls)[2], type="l", col="blue")
  
  # Calculate minimum development temperature (Tmin)
  (Tmin <- (-coef(Tmin.nls)[2]/coef(Tmin.nls)[1])[[1]])
  if(word(params[sp.num,1],1,2) == "Apolygus lucorum") { Tmin <- 283.1 } # for this species, Tmin was reported to be 283.1 in Lu et al. 2010
  if(params[sp.num,2] == "UK Sand Hutton") { Tmin <- 283.1 } # for this species, Tmin is set to 10C (otherwise, species will not overwinter)
  params[sp.num,]$Tmin <- round(Tmin, 1)
  
  # Break for loop if analyses are run for a specified species (all <- FALSE)
  if(all == FALSE) { break }
}

# Save model parameters to "Temperature response parameters.csv" (if desired)
if(save == TRUE) { write.csv(params, "Model parameters/Temperature response parameters.csv", row.names = FALSE) }

# Print model fit
if(all == FALSE) { params[sp.num,]
} else { params[x:s,] }