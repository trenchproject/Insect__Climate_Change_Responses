################################################################################
########## This R script fit functions to temperature response data ###########
################################################################################

# Load packages
library(tidyverse)
library(stringr)

# Set working directory (if necessary)
#setwd() # enter working directory of main downloaded file (containing R project file)

# USER: enter species name and location (used in "temperature response data.csv") OR set "all" to TRUE to run analysis for all populations
species <- "Clavigralla shadabi"
location <- "Benin"
all <- FALSE

# USER: Save model fit?
save <- FALSE

# Read data and model parameters file
name <- paste(species, location)
data <- as.data.frame(read_csv("Biological data/Temperature response data.csv"))
params <- as.data.frame(read_csv("Model parameters/Temperature response parameters.csv"))

# Run analyses for all populations
for(s in 1:nrow(params)) {

  # Find model parameters for selected populations or move iteratively through rows of “Temperature response parameters.csv”
  if(all == FALSE) {
    found <- FALSE
    for(i in 1:nrow(params)) {
      if(params[i,]$Population == name) {
        s <- i
        found <- TRUE
        break }
    }
    if(found == FALSE) {
      print("Species name or location not found in Habitat temperature parameters.csv, please update species name or location")
      break
    }
  } else(found <- TRUE)
  
  # Obtain data for selected population (NOTE: location for "Apolygus lucorum China Dafeng" set to "China Langfang" as in "Temperature response data.csv")
  if(params[s,1] == "Apolygus lucorum China Dafeng") { sp.data <- data[data$Population == "Apolygus lucorum China Langfang",]
  } else { sp.data <- data[data$Population == params[s,1],] } # find population in "Temperature response data.csv"
  
  # Remove columns that do not contain temperature data
  sp.data <- sp.data[-c(1:8,11,13,15,17,19,21,23,25)]
  
  # Set minimum and maximum values for x-axes and reference temperature (TR) for nls
  if(params[s,]$Habitat == "Tropical") {
    Xmin <- 285
    Xmax <- 315
  } else {
    Xmin <- 275
    Xmax <- 305
  }
  if(params[s,]$Population == "Clavigralla tomentosicollis Burkina Faso") { TR <- 300 # must use a higher TR for this populations
  } else if(params[s,]$Population == "Aulacorthum solani Brazil") { TR <- 292 # must use a different TR for this populations b/c mortality rate must be set to data
  } else if(params[s,]$Population == "Myzus persicae Canada Chatham") { TR <- 283 # must use a different TR for this populations b/c mortality rate must be set to data
  } else { TR <- 293 }
  params[s,]$TR <- TR
  
  
  ###################### INTRINSIC GROWTH RATE, r_m (Eq. 1a) #####################
  # Nonlinear regression
  if(params[s,1] != "Clavigralla tomentosicollis Burkina Faso") {
    r.nls <- nls(r_m ~ ifelse(T_K <= Toptr, rMax*exp(-(T_K-Toptr)^2/(2*sr^2)), rMax*(1 - ((T_K-Toptr)/(Toptr-Tmaxr))^2)),
                 data=sp.data, start=list(rMax=params[s,]$rMax, Toptr=params[s,]$Toptr, Tmaxr=params[s,]$Tmaxr, sr=params[s,]$sr))
    print(summary(r.nls))
  
    # Assign parameters
    params[s,]$rMax <- round(coef(r.nls)[1], 3)
    params[s,]$Toptr <- round(coef(r.nls)[2], 1)
    params[s,]$Tmaxr <- round(coef(r.nls)[3], 1)
    params[s,]$sr <- round(coef(r.nls)[4], 2)
  # For Clavigralla tomentosicollis Burkina Faso, must set rMax and Toptr to data and use nls to estimate Tmaxr and sr
  } else {
    rMax.test <- max(na.omit(sp.data[,"r_m"]))
    Toptr.test <- sp.data[which.max(sp.data$r), "T_K"]
    r.nls1 <- nls(r_m ~ ifelse(T_K <= Toptr.test, rMax.test*exp(-(T_K-Toptr.test)^2/(2*sr^2)), rMax.test*(1 - ((T_K-Toptr.test)/(Toptr.test-Tmaxr))^2)),
                 data=sp.data, start=list(Tmaxr=params[s,]$Tmaxr, sr=params[s,]$sr))
    summary(r.nls1)

    # Assign parameters
    params[s,]$rMax <- round(rMax.test, 3)
    params[s,]$Toptr <- round(Toptr.test, 1)
    params[s,]$Tmaxr <- round(coef(r.nls1)[1], 1)
    params[s,]$sr <- round(coef(r.nls1)[2], 2)
  }

  # Plot model fit
  plot(sp.data$T_K, sp.data$r)
  points(seq(Xmin,Xmax,1), ifelse(seq(Xmin,Xmax,1) <= params[s,]$Toptr, params[s,]$rMax*exp(-(seq(Xmin,Xmax,1)-params[s,]$Toptr)^2/(2*params[s,]$sr^2)),
                                  params[s,]$rMax*(1 - ((seq(Xmin,Xmax,1)-params[s,]$Toptr)/(params[s,]$Toptr-params[s,]$Tmaxr))^2)), type="l", col="blue")


  ##################### NET REPRODUCTIVE RATE, R0 (Eq. 1b) #######################
  # Nonlinear regression
  R0.nls <- nls(R0 ~ R0Topt*exp(-(T_K-ToptR0)^2/(2*sR0^2)), data=sp.data,
                start=list(R0Topt=params[s,]$R0Topt, ToptR0=params[s,]$ToptR0, sR0=params[s,]$sR0))
  summary(R0.nls)

  # Assign parameters
  params[s,]$R0Topt <- round(coef(R0.nls)[1], 1)
  params[s,]$ToptR0 <- round(coef(R0.nls)[2], 1)
  params[s,]$sR0 <- abs(round(coef(R0.nls)[3], 2))

  # Plot model fit
  plot(sp.data$T_K, sp.data$R0, ylim=c(0,max(sp.data[,"R0"], na.rm=TRUE)+1))
  points(seq(Xmin,Xmax,1), params[s,]$R0Topt*exp(-(seq(Xmin,Xmax,1)-params[s,]$ToptR0)^2/(2*params[s,]$sR0^2)), type="l", col="blue")


  ################### PER CAPITA BIRTH RATE, b[T] (Eq. 2a) #######################
  # Nonlinear regression
  b.nls<- nls(Birth_Rate ~ bTopt*exp(-(T_K-Toptb)^2/(2*sb^2)), data=sp.data,
             start=list(bTopt=params[s,]$bTopt, Toptb=params[s,]$Toptb, sb=params[s,]$sb))
  summary(b.nls)

  # Assign parameters
  params[s,]$bTopt <- round(coef(b.nls)[1], 2)
  params[s,]$Toptb <- round(coef(b.nls)[2], 1)
  params[s,]$sb <- abs(round(coef(b.nls)[3], 2))

  # Plot model fit
  plot(sp.data$T_K, sp.data$Birth_Rate, ylim=c(0,max(sp.data[,"Birth_Rate"], na.rm=TRUE)+0.01))
  points(seq(Xmin,Xmax,1), params[s,]$bTopt*exp(-(seq(Xmin,Xmax,1)-params[s,]$Toptb)^2/(2*params[s,]$sb^2)), type="l", col="blue")


  ################ STAGE-SPECIFIC MORTALITY RATE, d_i[T] (Eq. 2b) ################
  # Juvenile per capita mortality rate
  # Nonlinear regression
  if(params[s,1] != "Myzus persicae Canada Chatham" && params[s,1] != "Aulacorthum solani US Ithaca") {
    dJ.nls <- nls(Juv_Mortality ~ dJTR*exp(AdJ*(1/TR-1/T_K)), data=sp.data, start=list(dJTR=params[s,]$dJTR, AdJ=params[s,]$AdJ))
    print(summary(dJ.nls))

    # Assign parameters
    if(coef(dJ.nls)[1] > 0.001) { params[s,]$dJTR <- round(coef(dJ.nls)[1], 3)
    } else { params[s,]$dJTR <- 0.0001*trunc(10000*coef(dJ.nls)[1]) } # For populations in which dJTR is really low, simply truncate at the 4th decimal point
    params[s,]$AdJ <- round(coef(dJ.nls)[2], 0)
  # For Myzus persicae Canada Chatham and Aulacorthum solani US Ithaca, must set dJTR to data and then fit AdJ via nls (otherwise, curve poorly fits data for most habitat temperatures)
  } else {
    if(params[s,1] == "Myzus persicae Canada Chatham") { dJTR.test <- sp.data[2,"Juv_Mortality"] }
    if(params[s,1] == "Aulacorthum solani US Ithaca") { dJTR.test <- sp.data[3,"Juv_Mortality"] }
    dJ.nls <- nls(Juv_Mortality ~ dJTR.test*exp(AdJ*(1/TR-1/T_K)), data=sp.data, start=list(AdJ=params[s,]$AdJ))
    print(summary(dJ.nls))
    
    # Assign parameters
    params[s,]$dJTR <- round(dJTR.test, 4)
    params[s,]$AdJ <- round(coef(dJ.nls)[1], 0)
  }
  
  # Plot model fit
  plot(sp.data$T_K, sp.data$Juv_Mortality, ylim=c(0,max(sp.data[,"Juv_Mortality"], na.rm=TRUE)+0.01))
  points(seq(Xmin,Xmax,1), params[s,]$dJTR*exp(params[s,]$AdJ*(1/TR-1/seq(Xmin,Xmax,1))), type="l", col="blue")
  
  
  # Adult per capita mortality rate
  # Nonlinear regression
  if(params[s,1] != "Aulacorthum solani Brazil") {
    dA.nls <- nls(Adult_Mortality ~ dATR*exp(AdA*(1/TR-1/T_K)), data=sp.data, start=list(dATR=params[s,]$dATR, AdA=params[s,]$AdA))
    print(summary(dA.nls))

    # Assign parameters
    params[s,]$dATR <- round(coef(dA.nls)[1], 3)
    params[s,]$AdA <- round(coef(dA.nls)[2], 0)
  # For Aulacorthum solani in Brazil, must set dATR to data and then fit AdA via nls (otherwise, curve poorly fits data for most habitat temperatures)
  } else {
    dATR.test <- sp.data[2,"Adult_Mortality"]
    dA.nls <- nls(Adult_Mortality ~ dATR.test*exp(AdA*(1/TR-1/T_K)), data=sp.data, start=list(AdA=params[s,]$AdA))
    print(summary(dA.nls))

    # Assign parameters
    params[s,]$dATR <- round(dATR.test, 3)
    params[s,]$AdA <- round(coef(dA.nls)[1], 0)
  }

  # Plot model fit
  plot(sp.data$T_K, sp.data$Adult_Mortality, ylim=c(0,max(sp.data[,"Adult_Mortality"], na.rm=TRUE)+0.01))
  points(seq(Xmin,Xmax,1), params[s,]$dATR*exp(params[s,]$AdA*(1/TR-1/seq(Xmin,Xmax,1))), type="l", col="blue")


  ####################### DEVELOPMENT RATE, g[T] (Eq. 2c) ########################
  # Quantify Toptg directly from data (descend down temperature treatments until the maximum development rate is reached, then exit the for loop)
  for(i in nrow(sp.data):1) { if(sp.data[i,"Dev_Rate"] == max(na.omit(sp.data[,"Dev_Rate"]))) { break } }
  Toptg <- sp.data[i,"T_K"]

  # Quantify Tmaxg directly from data (set Tmaxg above the highest temperature treatment, then descend until the development rate is above zero and exit the for loop)
  if(sp.data[nrow(sp.data),"Dev_Rate"] != 0) { Tmaxg <- sp.data[nrow(sp.data),"T_K"] # set Tmaxg to highest laboratory temperature if development rate is not zero at this temperature
  } else {
    for(j in nrow(sp.data):i) { if(sp.data[j,"Dev_Rate"] != 0) { break } } # find the highest laboratory temperature at which the development rate is not zero
    Tmaxg <- (sp.data[j,"T_K"] + sp.data[j+1,"T_K"])/2 # set Tmaxg to the average between this temperature and the next highest laboratory temperature
  }

  # Update i index (if necessary)
  if(word(params[s,2],1) == "Brazil") { i <- 3 } # For Brazil populations, set i index to minimum number of temperature treatments for nls below
  if(word(params[s,1],1,2) == "Apolygus lucorum") { # For Apolygus lucorum, set i index to maximum laboratory temperature for nls below
    for(i in nrow(sp.data):1) { if(sp.data[i,"Dev_Rate"] != 0) { break } }
  }

  # Assign Toptg and Tmaxg parameters
  params[s,]$Toptg <- round(Toptg, 1)
  params[s,]$Tmaxg <- round(Tmaxg, 1)

  # Try to fit remaining parameters to data (excluding data at high temperatures at which the development rate declines towards zero)
  error1 <- tryCatch({
    g.nls <- nls(Dev_Rate ~ gTR*(T_K/TR)*exp(Ag*(1/TR-1/T_K))/(1+exp(AL*(1/TL-1/T_K))), data=sp.data[1:i,],
                              start=list(gTR=params[s,]$gTR, Ag=params[s,]$Ag, AL=params[s,]$AL, TL=params[s,]$TL))
    print(summary(g.nls))

    # Function nls does not return an error
    error1 <- FALSE
  }, error = function(err) { error1 <- TRUE; return(error1) })

  # If there are no errors, then assign parameters
  if(error1 == FALSE) {
    params[s,]$gTR <- round(coef(g.nls)[1], 3)
    params[s,]$Ag <- round(coef(g.nls)[2], 0)
    params[s,]$AL <- round(coef(g.nls)[3], 0)
    params[s,]$TL <- round(coef(g.nls)[4], 1)
    params[s,]$gMax <- round(params[s,]$gTR*(params[s,]$Toptg/TR)*exp(params[s,]$Ag*(1/TR-1/params[s,]$Toptg))/
                                    (1+exp(params[s,]$AL*(1/params[s,]$TL-1/params[s,]$Toptg))), 3)
  } else { # If nls returns an error, then try alternative method of first fitting gTR and Ag and then fitting AL and TL afterwards
    # Nonlinear regression for gTR and Ag
    g.nls1 <- nls(Dev_Rate ~ gTR*(T_K/TR)*exp(Ag*(1/TR-1/T_K)), data=sp.data[1:i,],
                  start=list(gTR=params[s,]$gTR, Ag=params[s,]$Ag))
    print(summary(g.nls1))

    # Assign parameters for next nonlinear regression and update parameters
    gTR.test <- coef(g.nls1)[1]
    Ag.test <- coef(g.nls1)[2]
    params[s,]$gTR <- round(gTR.test, 3)
    params[s,]$Ag <- round(Ag.test, 0)

    # Nonlinear regression for AL and TL
    error2 <- tryCatch({
      g.nls2 <- nls(Dev_Rate ~ gTR.test*(T_K/TR)*exp(Ag.test*(1/TR-1/T_K))/(1+exp(AL*(1/TL-1/T_K))), data=sp.data[1:i,],
                   start=list(AL=params[s,]$AL, TL=params[s,]$TL))
      print(summary(g.nls2))

      # Function nls does not return an error
      error2 <- FALSE
    }, error = function(err) { error2 <- TRUE; return(error2) })

    # If there are no errors in the alternative method above, then assign parameters
    if(error2 == FALSE) {
      # Assign parameters
      params[s,]$AL <- round(coef(g.nls2)[1], 0)
      params[s,]$TL <- round(coef(g.nls2)[2], 1)
      params[s,]$gMax <- round(params[s,]$gTR*(params[s,]$Toptg/TR)*exp(params[s,]$Ag*(1/TR-1/params[s,]$Toptg))/
                                      (1+exp(params[s,]$AL*(1/params[s,]$TL-1/params[s,]$Toptg))), 3)
    } else {
      # If nls still returns an error, then try alternative method of first fitting AL and then TL
      TL.test <- min(sp.data[,"T_K"]) + 2
      g.nls2 <- nls(Dev_Rate ~ gTR.test*(T_K/TR)*exp(Ag.test*(1/TR-1/T_K))/(1+exp(AL*(1/TL.test-1/T_K))), data=sp.data[1:i,],
                    start=list(AL=params[s,]$AL))
      print(summary(g.nls2))

      # Nonlinear regression for TL
      AL.test <- coef(g.nls2)[1]
      g.nls3 <- nls(Dev_Rate ~ gTR.test*(T_K/TR)*exp(Ag.test*(1/TR-1/T_K))/(1+exp(AL.test*(1/TL-1/T_K))), data=sp.data[1:i,],
                    start=list(TL=params[s,]$TL))
      print(summary(g.nls3))

      # Assign parameters
      params[s,]$AL <- round(AL.test, 0)
      params[s,]$TL <- round(coef(g.nls3)[1], 1)
      params[s,]$gMax <- round(params[s,]$gTR*(params[s,]$Toptg/TR)*exp(params[s,]$Ag*(1/TR-1/params[s,]$Toptg))/
                                      (1+exp(params[s,]$AL*(1/params[s,]$TL-1/params[s,]$Toptg))), 3)
    }
  }

  # plot model fit
  plot(sp.data$T_K, sp.data$Dev_Rate, xlim=c(Xmin,Tmaxg + 2), ylim=c(0, 1.5*round(params[s,]$gMax, 2)))
  g.funct <- function(T) { ifelse(T <= params[s,]$Toptg, params[s,]$gTR*(T/TR)*exp(params[s,]$Ag*(1/TR-1/T))/
                                    (1+exp(params[s,]$AL*(1/params[s,]$TL-1/T))),
                                  ifelse(T <= params[s,]$Tmaxg, params[s,]$gMax, 0)) }
  points(seq(Xmin,Tmaxg + 2,0.1), g.funct(seq(Xmin,Tmaxg + 2,0.1)), type="l", col="blue")


  # Minimum developmental temperature (for overwintering in DDE model)
  # Linear regression of linear portion of development rate TPC
  Tmin.nls <- nls(Dev_Rate ~ m*T_K+b, data=sp.data[1:i,], start=list(m=params[s,]$gTR, b=0))
  summary(Tmin.nls)
  
  # Plot model fit
  plot(sp.data$T_K, sp.data$Dev_Rate, xlim=c(Xmin,Xmax))
  points(seq(Xmin,Xmax,1), coef(Tmin.nls)[1]*seq(Xmin,Xmax,1)+coef(Tmin.nls)[2], type="l", col="blue")

  # Calculate minimum development temperature (Tmin)
  (Tmin <- (-coef(Tmin.nls)[2]/coef(Tmin.nls)[1])[[1]])
  if(word(params[s,1],1,2) == "Apolygus lucorum" || # for Apolygus lucorum, Tmin was reported to be 283.1 in Lu et al. 2010
     params[s,2] == "UK Sand Hutton") { params[s,]$Tmin <- 278.1 # for Acyrthosiphon pisum, Tmin set to 5C (rough habitat temperature when the population emerges in the field); otherwise, population doesn't overwinter
  } else { params[s,]$Tmin <- round(Tmin, 1) }
  
  # Break for loop (line 26) if analyses are run for a specified population (all <- FALSE)
  if(all == FALSE) { break }
}

# Save model parameters to "Temperature response parameters.csv" (if desired)
if(save == TRUE && found == TRUE) { write.csv(params, "Model parameters/Temperature response parameters.csv", row.names = FALSE) }

# Print model fit (if the population was found in "Temperature response parameters.csv")
if(found == TRUE) {
  if(all == FALSE) { params[s,]
  } else { params[1:s,] }
}