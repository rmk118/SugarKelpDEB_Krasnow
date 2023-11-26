### Intro ################################################################################################################# 
#Site names here begin with names other than those used in the manuscript
#Sled = Pt Judith Pond N
#Dredge = Pt Judith Pond S

#File originally created by Celeste Venolia in March 2018-December 2019 for Venolia et al., (2020)
# https://doi.org/10.1016/j.ecolmodel.2020.109151 

# Modified by Ruby Krasnow in November-December 2023
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Import libraries
library(deSolve)
library(tidyverse)
library(lubridate)
library(gridExtra)
library(gdata)
library(Metrics)
library(pracma)
library(patchwork)

#Required for model runs
source("SolveR_R.R")
source("KelpDEB_model.R")
#Required for Calibration Code
source("N_uptake_Calibration.R")
source("Photosynthesis_Calibration.R")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### Minerals and Organics Section #####
#Conversion coefficients, organics (n = matrix of chemical indices)
# "food N" "food C" Structure "N reserves" "C reserves" products
#     X_N   X_C      V    E_N    E_C    P
n_O <- matrix(
    + c(0.00, 1.00, 1.00, 0.00, 1.00, 1.00,  #C/C, equals 1 by definition
      + 0.00, 0.50, 1.33, 0.00, 2.00, 1.80,  #H/C, these values show that we consider dry-mass
      + 3.00, 2.50, 1.00, 2.50, 1.00, 0.50,  #O/C
      + 1.00, 0.00, 0.04, 1.00, 0.00, 0.04), nrow=4, ncol=6, byrow = TRUE) #N/C
#V is the C-mol structure of alginate (Alginic acid: (C6H8O6)n)
#E_N is N03- and N02- averaged
#E_C is glucose C6H12O6 (Laminarin: c18h32o16 and mannitol c6h14o6)
#We aren't using the X_N, X_C, or P columns here

#Molecular weights
#t() is a matrix transpose function
#organics structure matrix multiplied by the atomic masses (mass in grams of one mole of an element) of C H O N
w_O_step <- t(n_O)*matrix(c(12, 1, 16, 14), nrow=6, ncol=4, byrow= TRUE) #g/mol, molecular weights for organics
w_O <- rowSums(w_O_step) #this provides g/mol of each of the six "pockets of mass" (i.e. X_N, X_C)

#define molecular weights
w_V <- w_O[3]  # g/mol       #molecular weight of structure
w_EN <- w_O[4]  # g/mol      #molecular weight of N reserve
w_EC <- w_O[5]  #g/mol       #molecular weight of C reserve
w_O2 <- 32 #g/mol
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                  #Ruby's Note
#Note that * performs element-wise multiplication, NOT matrix multiplication, which is denoted by %*%. 
# This means the value in n_0[1,1] is multiplied by the value in the [1,1] position in the 
# matrix of molecular weights: matrix(c(12, 1, 16, 14), nrow=6, ncol=4, byrow= TRUE).
# I found it more intuitive to achieve the same result by doing: w_O <- t(n_O) %*% c(12,1,16,14)
# Then you can define the molecular weights as w_V <- w_O[3], w_EN <- w_O[4], etc.

# Not entirely sure where the 0.04 N/C comes from in the structure if V is assumed to have the C-mol structure of alginate

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##### Parameters compiled #####
params_Lo <- c(#maximum volume-specific assimilation rate of N before temperature correction
               JENAM = 1.5e-4, #mol N / molV / h
               #half saturation constant of N uptake
               K_N = 2.5e-6, #molNO3 and NO2/L
               #max volume-specific carbon dioxide assimilation rate
               JCO2M = 0.0075, #mol DIC/molV/h
               #half saturation constant of C uptake
               K_C = 4e-7, #mol DIC/L
               #maximum volume-specific carbon assimilation rate
               JECAM = 0.282, #molC/molV/h
               #Photosynthetic unit density
               rho_PSU = 0.5, #mol PSU/ mol V
               #binding probability of photons to a Light SU
               b_I = 0.5, #dimensionless
               #Specific photon arrival cross section
               alpha = 1, #m^2 mol PSU–1
               #dissociation rate
               k_I = 0.075, #molγ molPS–1 h–1
               #Yield factor of C reserve to photon
               y_I_C = 10, #mol γ mol C-1
               #Yield factor of C reserve to DIC
               y_CO2_C = 1, #mol DIC mol C-1
               #Yield factor of photon to O2
               y_LO2 = 0.125, #molO2 molγ –1
               #reserve turnover
               kE_C = 0.02, #0.05, #1/h
               kE_N = 0.04, #0.01, #1/h
               #fraction of rejection flux from growth SU incorporated back into i-reserve
               kappa_Ei = 0.9, #dimensionless
               #yield of structure on N reserve (percent of N in structure)
               y_EN_V = 0.04, #mol N/mol V
               #yield of structure on C reserve (percent of C in structure)
               y_EC_V = 1, #mol C/mol V
               #specific maintenance costs requiring N before temp correction
               JENM = 4*10^-6, #4e-6, #mol N/molM_V/h
               #specific maintenance costs requiring C before temp correction
               JECM = 1*10^-6, #1e-6, #mol C/molM_V/h
               #Arrhenius temperature
               T_A = 6314.3, # K
               #Upper boundary of temperature tolerance
               T_H = 13.386 + 273.15, # K
               #Lower boundary of temperature tolerance
               T_L = 273.15, # K
               #Arrhenius temperature outside T_H
               T_AH = 18702, #K
               #Arrhenius temperature outside T_L
               T_AL = 4391.9, #K
               #temperature at which rate parameters are given
               T_0 = 20 + 273.15) # K
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####### Initial conditions ############
#Initial conditions of state variables
#these values are not coming from any field data or literature information, estimated
state_Lo <- c(m_EC = 0.002, #0.1, #mol C/molM_V  #Reserve density of C reserve (initial mass of C reserve per initial mass of structure)
              m_EN = 0.01, #mol N/molM_V #Reserve density of N reserve (initial mass of N reserve per initial mass of structure)
              M_V = 0.05/(w_V+0.01*w_EN+0.002*w_EC)) #molM_V #initial mass of structure

state_Johansson <- c(m_EC = 0.3, #mol C/molM_V  #Reserve density of C reserve (initial mass of C reserve per initial mass of structure)
                     m_EN = 0.01, #mol N/molM_V #Reserve density of N reserve (initial mass of N reserve per initial mass of structure)
                     M_V = 0.005/(w_V+0.01*w_EN+0.3*w_EC)) #molM_V #initial mass of structure
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####### Time steps #######
#(First number of time step, last number of time step, interval to step)
times_Lo_Sled1 <- seq(0, 4008, 1) #167 days stepped hourly
times_Lo_Sled2 <- seq(0, 3336, 1) #139 days stepped hourly
times_Lo_Dredge1 <- seq(0, 4128, 1) #172 days stepped hourly
times_Lo_Dredge2 <- seq(0, 3456, 1) #144 days stepped hourly

#### Irradiance forcing set-up ####
#NOAA irradiance data set-up: NOAASurfaceIrradiance
NOAA_Irradiance <- read.csv("NOAASurfaceIrradiance.csv", header = TRUE, fileEncoding="UTF-8-BOM")
NOAA_Irradiance$DateTime <- dmy_hms(NOAA_Irradiance$DateTime, tz = "UTC") #NOAA data in UTC (5 hours ahead)
NOAA_Irradiance <- with_tz(NOAA_Irradiance, "America/New_York") #Convert from UTC to EST
NOAA_Irradiance$DownMinusUp <- NOAA_Irradiance$dswrf-NOAA_Irradiance$uswrf #net shortwave radiation at the surface (W/m^2) is obtained by subtracting the upward short wave flux (uswrf) from the downward flux (dswrf)
#PAR = NSW*PAR_frac*C*exp(-k*z)*3600
#NSW=dswrf-uswrf
#PAR_frac is the fraction of the incident flux that is useable for photosynthesis
#C is a conversion factor = 4.56 umol photons/s/W
#k is the extinction coefficient
#3600 converts from s^-1 to h^-1
#1e-6 converts from micomoles to moles
NOAA_Irradiance$PAR <- NOAA_Irradiance$DownMinusUp*0.43*4.56*exp(-0.46*1)*3600*1e-6

##### N forcing set-up ####
WSA2_Y1 <- read.csv("WaterSampleAnalysis2Y1.csv", header = TRUE, fileEncoding="UTF-8-BOM") #Import water Q data
WSA2_Y1$Date <- mdy(WSA2_Y1$Date) #convert dates
names(WSA2_Y1)[1] <- "Site" #only necessary for some computers running this code

Sled_WSA <- filter(WSA2_Y1, Site == "Sled") #filter for Sled site
Dredge_WSA <- filter(WSA2_Y1, Site == "Dredge") #filter for Dredge site


##### DIC forcing set-up ###########
DIC <- read.csv("Ninigret_EPA_DIC.csv", header = TRUE, fileEncoding="UTF-8-BOM") #Import Ninigret DIC data
CO_2 <- mean(DIC$DIC.uMkg.mean) #micromole DIC/kg (Jason said it was okay to assume that 1kg of seawater is 1L of seawater (actual conversion requires density calc from salinity and T))
#need units to match K_C (molDIC/L)
CO_2 <- CO_2/1000000

##### Temp forcing set-up ###########
Sled_Y1_hobotemp_orig <- read.csv("Sled_Y1_TempLogger2.csv", header = TRUE, fileEncoding="UTF-8-BOM") #import
Sled_Y1_hobotemp_orig$DateTime <- mdy_hms(Sled_Y1_hobotemp_orig$Date_Time) #convert time field
Sled_Y1_hobotemp_orig$Temp_K <- Sled_Y1_hobotemp_orig$Temp_C+273.15 #create column with temp in K

Dredge_Y1_hobo_orig <- read.csv("Dredge_Y1_hobo.csv", header = TRUE, fileEncoding="UTF-8-BOM") #importing neighboring temp file to replace corrupted section
Dredge_Y1_hobo_orig$DateTime <- mdy_hms(Dredge_Y1_hobo_orig$Date_Time) #convert time field
Dredge_Y1_hobo_orig$Temp_K <- Dredge_Y1_hobo_orig$Temp_C+273.15 #create column with temp in K

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### FIELD DATA MODEL RUNS ####

### Judith N (sled) line 1 ####
W <- 0.05 #inital biomass for conversions (cannot put in initial conditions)

# N forcing set-up Judith N 1
N <- Sled_WSA[c("Date","NitrateNitrite_uM")] #subset
N$NitrateNitrite_uM <- N$NitrateNitrite_uM/1000000
#Converted to hourly by multiply by 24
N_field <- approxfun(x = c(161*24, 139*24, 105*24, 0, 28*24, 84*24, 172*24), y = N$NitrateNitrite_uM, method = "linear", rule = 2) #N forcing function

##### Irradiance set-up Judith N 1 ####
NOAA_Irradiance_Sledy1 <-  NOAA_Irradiance$PAR[2438:3774] # subset based on as_datetime("2017-11-1 12:00:00"), as_datetime("2018-04-17 12:00:00")
I_field <- approxfun(x = seq(from = 0, to = 4008, by = 3), y = NOAA_Irradiance_Sledy1, method = "linear", rule = 2) #irradiance forcing

##### Temp set-up Judith N 1 #############
Sled_Y1_hobotemp <- Sled_Y1_hobotemp_orig[14:16049,] #subset
SledT_hourly <- ceiling_date(Sled_Y1_hobotemp$DateTime, unit = "hour") #set the values to aggregate around
AvgTempKbyhr <- aggregate(Sled_Y1_hobotemp$Temp_K, by=list(SledT_hourly), mean) #calculate average hourly temp
AvgTempKbyhr_part1 <- AvgTempKbyhr$x[0:334] #subset

Dredge_Y1_hobo <- Dredge_Y1_hobo_orig[3:16531,] #subset
DredgeT_hourly <- ceiling_date(Dredge_Y1_hobo$DateTime, unit = "hour") #set the values to aggregate around
AvgTempKbyhr4FD <- aggregate(Dredge_Y1_hobo$Temp_K, by=list(DredgeT_hourly), mean) #calculate average hourly temp
AvgTempKbyhr4FD <- AvgTempKbyhr4FD[4:4132, ] #subset
fd <- AvgTempKbyhr4FD$x[335:859] #526 data points needed from dredge to replace a weird glitch in the sled temp data
AvgTempKbyhr_part2 <- AvgTempKbyhr$x[860:4009] #later part of original temp file
T_field <- approxfun(x = c(0:4008), y = c(AvgTempKbyhr_part1, fd, AvgTempKbyhr_part2), method = "linear", rule = 2) #the temp forcing function
T_Sled1_Y1 <- T_field(0:4008) #saving the forcing this way for ease of later visualization

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Model run (the differential equation solver)
sol_Sled1 <- ode(y = state_Lo, t = times_Lo_Sled1, func = rates_Lo, parms = params_Lo)

W <- 0.05
sol_high_Sled1 <- ode(y = state_Lo, t = times_Lo_Sled1, func = rates_Lo, parms = high_params_for_model)
W <- 0.05
sol_med_Sled1 <- ode(y = state_Lo, t = times_Lo_Sled1, func = rates_Lo, parms = med_params_for_model)
W <- 0.05
sol_low_Sled1 <- ode(y = state_Lo, t = times_Lo_Sled1, func = rates_Lo, parms = low_params_for_model)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setting up the forcing functions with field data for
### Judith N (sled) line 2 ####
W <- 0.05 #inital biomass for conversions (cannot put in initial conditions)

# N forcing set-up Judith N 2 #
N <- Sled_WSA[c("Date","NitrateNitrite_uM")] #create a dataframe with just the relevant collums
N$NitrateNitrite_uM <- N$NitrateNitrite_uM/1000000 #convert from micromoles/L to moles/L
#multiplying by 24 to set as hourly
N_field <- approxfun(x = c(133*24, 111*24, 77*24, -28*24, 0, 56*24, 144*24), y = N$NitrateNitrite_uM, method = "linear", rule = 2) #N forcing function
###### DIC forcing set-up ###########

###### NOAA Irradiance forcing set-up Judith N 2 ####
NOAA_Irradiance_Sledy1_L2 <-  NOAA_Irradiance$PAR[2662:3774] #subset by seq(as_datetime("2017-11-29 12:00:00"), as_datetime("2018-04-17 12:00:00"), by="hour")
I_field <- approxfun(x = seq(from = 0, to = 3336, by = 3), y = NOAA_Irradiance_Sledy1_L2, method = "linear", rule = 2) #irradiance forcing function

###### Temp forcing set-Up Judith N 2 #############
Sled_Y1_hobotemp <- Sled_Y1_hobotemp_orig[6:16051,] #subset
SledT_hourly <- ceiling_date(Sled_Y1_hobotemp$DateTime, unit = "hour") #set values to aggregate around
AvgTempKbyhr <- aggregate(Sled_Y1_hobotemp$Temp_K, by=list(SledT_hourly), mean) #calculate average hourly temp
AvgTempKbyhr <- AvgTempKbyhr[677:4011,] #subset

Dredge_Y1_hobo <- Dredge_Y1_hobo_orig[3:16531,] #subset
DredgeT_hourly <- ceiling_date(Dredge_Y1_hobo$DateTime, unit = "hour") #set values to aggregate around
AvgTempKbyhr4FD <- aggregate(Dredge_Y1_hobo$Temp_K, by=list(DredgeT_hourly), mean) #calculate average hourly temp
AvgTempKbyhr4FD <- AvgTempKbyhr4FD[4:4132, ] #subset
fd <- AvgTempKbyhr4FD$x[858:859] #526 data points needed from dredge to replace a weird glitch in the sled temp data
T_field <- approxfun(x = c(0:3336), y = c(fd, AvgTempKbyhr$x), method = "linear", rule = 2) #the temp forcing function
T_Sled2_Y1 <- T_field(0:3336) #for later ease in plotting the forcing

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Model run (the differential equation solver)
sol_Sled2 <- ode(y = state_Lo, t = times_Lo_Sled2, func = rates_Lo, parms = params_Lo)

W <- 0.05
sol_high_Sled2 <- ode(y = state_Lo, t = times_Lo_Sled2, func = rates_Lo, parms = high_params_for_model)
W <- 0.05
sol_med_Sled2 <- ode(y = state_Lo, t = times_Lo_Sled2, func = rates_Lo, parms = med_params_for_model)
W <- 0.05
sol_low_Sled2 <- ode(y = state_Lo, t = times_Lo_Sled2, func = rates_Lo, parms = low_params_for_model)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setting up the forcing functions with field data for
### Judith S (dredge) line 1 ####
W <- 0.05 #inital biomass for conversions (cannot put in initial conditions)

###### N forcing set-up Judith S 1 ##############
N <- Dredge_WSA[c("Date","NitrateNitrite_uM")] #new dataframe with the relevant collumns
N$NitrateNitrite_uM <- N$NitrateNitrite_uM/1000000 #convert from micromoles/L to moles/L
#multipling by 24 to take from daily to hourly
N_field <- approxfun(x = c(139*24, 161*24, 84*24, 0, 105*24), y = N$NitrateNitrite_uM, method = "linear", rule = 2) #N forcing function

###### NOAA Irradiance forcing set-up Judith S 1 ####
NOAA_Irradiance_Dredgey1 <-  NOAA_Irradiance$PAR[2438:3814] #subset by 11/1/17 to 2018-04-22 12:00:00
I_field <- approxfun(x = seq(from = 0, to = 4128, by = 3), y = NOAA_Irradiance_Dredgey1, method = "linear", rule = 2) #irradiance forcing function

###### Temp forcing set-Up Judith S 1 #############
Dredge_Y1_hobo <- Dredge_Y1_hobo_orig[3:16531,] #subset
DredgeT_hourly <- ceiling_date(Dredge_Y1_hobo$DateTime, unit = "hour") #set values to aggregate around
AvgTempKbyhr <- aggregate(Dredge_Y1_hobo$Temp_K, by=list(DredgeT_hourly), mean) #calculate average hourly temp
AvgTempKbyhr <- AvgTempKbyhr[4:4132, ] #subset
T_field <- approxfun(x = c(0:4128), y = AvgTempKbyhr$x, method = "linear", rule = 2) #the temp forcing function
T_Dredge1_Y1 <- T_field(0:4128) #for ease in later plotting of the forcing

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Model run (the differential equation solver)
sol_Dredge1 <- ode(y = state_Lo, t = times_Lo_Dredge1, func = rates_Lo, parms = params_Lo)

W <- 0.05
sol_high_Dredge1 <- ode(y = state_Lo, t = times_Lo_Dredge1, func = rates_Lo, parms = high_params_for_model)
W <- 0.05
sol_med_Dredge1 <- ode(y = state_Lo, t = times_Lo_Dredge1, func = rates_Lo, parms = med_params_for_model)
W <- 0.05
sol_low_Dredge1 <- ode(y = state_Lo, t = times_Lo_Dredge1, func = rates_Lo, parms = low_params_for_model)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setting up the forcing functions with field data for
### Judith S line 2 ####
W <- 0.05 #inital biomass for conversions (cannot put in initial conditions)

# N forcing set-up Judith S 2
N <- Dredge_WSA[c("Date","NitrateNitrite_uM")] #create new dataframe with just the relevant collumns
N$NitrateNitrite_uM <- N$NitrateNitrite_uM/1000000 #convert from micromoles/L to moles/L
#multiplied by 24 to convert from daily to hourly
N_field <- approxfun(x = c(111*24, 133*24, 56*24, -28*24, 77*24), y = N$NitrateNitrite_uM, method = "linear", rule = 2) #N forcing function

###### NOAA Irradiance forcing set-up Judith S 2 ####
NOAA_Irradiance_Dredgey1_L2 <-  NOAA_Irradiance$PAR[2662:3814] #subset by seq(as_datetime("2017-11-29 12:00:00"), as_datetime("2018-04-22 12:00:00"), by="hour")
I_field <- approxfun(x = seq(from = 0, to = 3456, by = 3), y = NOAA_Irradiance_Dredgey1_L2, method = "linear", rule = 2) #irradiance forcing function

###### Temp forcing set-up Judith S 2 #############
Dredge_Y1_hobo <- Dredge_Y1_hobo_orig[3:16531,] #cut 2 points in beginning, logger not yet in water
DredgeT_hourly <- ceiling_date(Dredge_Y1_hobo$DateTime, unit = "hour") #determine dates to aggregate around
AvgTempKbyhr <- aggregate(Dredge_Y1_hobo$Temp_K, by=list(DredgeT_hourly), mean) #calculate average hourly temp
AvgTempKbyhr_sub <- AvgTempKbyhr[676:4132,] #subset
T_field <- approxfun(x = c(0:3456), y = c(AvgTempKbyhr_sub$x), method = "linear", rule = 2) #the temp forcing function
T_Dredge2_Y1 <- T_field(0:3456) #for ease of later plotting the temperature forcing

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Model run (the differential equation solver)
sol_Dredge2 <- ode(y = state_Lo, t = times_Lo_Dredge2, func = rates_Lo, parms = params_Lo)

W <- 0.05
sol_high_Dredge2 <- ode(y = state_Lo, t = times_Lo_Dredge2, func = rates_Lo, parms = high_params_for_model)
W <- 0.05
sol_med_Dredge2 <- ode(y = state_Lo, t = times_Lo_Dredge2, func = rates_Lo, parms = med_params_for_model)
W <- 0.05
sol_low_Dredge2 <- ode(y = state_Lo, t = times_Lo_Dredge2, func = rates_Lo, parms = low_params_for_model)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### END FIELD DATA, START LITERATURE DATA FOR CALIBRATION ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setting up the forcing functions for Espinoza and Chapman (1983) nitrogen uptake #9C (282.15 K)

Nmax <- 73.1221719/1000000 #M
T_dat <- 9 #C (conversion in Nuptake function to K)

#Model run (the differential equation solver)
sol_EspinozaChapman1983_N_9 <- Nuptake(params_Lo, T_dat, Nmax, w_EN) #function from N_uptake_Calibration.R code
sol_EspinozaChapman1983_N_9 <- as.data.frame(sol_EspinozaChapman1983_N_9) #conversion to dataframe for later use
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setting up the forcing functions for Espinoza and Chapman (1983) nitrogen uptake #18C (282.15 K)

Nmax <- 76.9543147/1000000 #M
T_dat <- 18 #C (conversion in Nuptake function to K)

#Model run (the differential equation solver)
sol_EspinozaChapman1983_N_18 <- Nuptake(params_Lo, T_dat, Nmax, w_EN) #function from N_uptake_Calibration.R code
sol_EspinozaChapman1983_N_18 <- as.data.frame(sol_EspinozaChapman1983_N_18) #conversion to dataframe for later use
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Photosynthesis model calibration
T_dat <- 14 #C (maintained for entire experiment)
I_max <- 3233205 #micromol photons m-2 h-1
sol_Johansson2002 <- Photosynthesis(params_Lo, state_Johansson, w_V, w_EN, w_EC, w_O2, T_dat, I_max) #function from Photosynthesis_Calibration.R
sol_Johansson2002 <- as.data.frame(sol_Johansson2002) #conversion to dataframe for later use
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###### Convert DeSolve solutions into data frame for broader plotting use ####
#conversions to dataframes
sol_Sled1 <- as.data.frame(sol_Sled1)
sol_Sled2 <- as.data.frame(sol_Sled2)
sol_Dredge1 <- as.data.frame(sol_Dredge1)
sol_Dredge2 <- as.data.frame(sol_Dredge2)

sol_high_Sled1 <- as.data.frame(sol_high_Sled1)
sol_high_Sled2 <- as.data.frame(sol_high_Sled2)
sol_high_Dredge1 <- as.data.frame(sol_high_Dredge1)
sol_high_Dredge2 <- as.data.frame(sol_high_Dredge2)

sol_med_Sled1 <- as.data.frame(sol_med_Sled1)
sol_med_Sled2 <- as.data.frame(sol_med_Sled2)
sol_med_Dredge1 <- as.data.frame(sol_med_Dredge1)
sol_med_Dredge2 <- as.data.frame(sol_med_Dredge2)

sol_low_Sled1 <- as.data.frame(sol_low_Sled1)
sol_low_Sled2 <- as.data.frame(sol_Sled2)
sol_low_Dredge1 <- as.data.frame(sol_low_Dredge1)
sol_low_Dredge2 <- as.data.frame(sol_low_Dredge2)

#addition of a date variable
sol_Sled1$Date <- seq(as_datetime("2017-11-1 12:00:00"), as_datetime("2018-04-17 12:00:00"), by="hour")
sol_high_Sled1$Date <- seq(as_datetime("2017-11-1 12:00:00"), as_datetime("2018-04-17 12:00:00"), by="hour")
sol_med_Sled1$Date <- seq(as_datetime("2017-11-1 12:00:00"), as_datetime("2018-04-17 12:00:00"), by="hour")
sol_low_Sled1$Date <- seq(as_datetime("2017-11-1 12:00:00"), as_datetime("2018-04-17 12:00:00"), by="hour")


sol_Sled2$Date <- seq(as_datetime("2017-11-29 12:00:00"), as_datetime("2018-04-17 12:00:00"), by="hour")
sol_Dredge1$Date <- seq(as_datetime("2017-11-1 12:00:00"), as_datetime("2018-04-22 12:00:00"), by="hour")
sol_Dredge2$Date <-seq(as_datetime("2017-11-29 12:00:00"), as_datetime("2018-04-22 12:00:00"), by="hour")

#conversion back to Celsius from Kelvin
sol_Sled1$Temp_C <- T_Sled1_Y1 - 273.15

sol_high_Sled1$Temp_C <- T_Sled1_Y1 - 273.15
sol_med_Sled1$Temp_C <- T_Sled1_Y1 - 273.15
sol_low_Sled1$Temp_C <- T_Sled1_Y1 - 273.15

sol_Sled2$Temp_C <- T_Sled2_Y1 - 273.15
sol_Dredge1$Temp_C <- T_Dredge1_Y1 - 273.15
sol_Dredge2$Temp_C <- T_Dredge2_Y1 - 273.15

#create source column to prepare for binding all these dataframes together
sol_Sled1$source <- "Point Judith Pond N 1"

sol_high_Sled1$source <- "Point Judith Pond N 1"
sol_med_Sled1$source <- "Point Judith Pond N 1"
sol_low_Sled1$source <- "Point Judith Pond N 1"

sol_Sled2$source <- "Point Judith Pond N 2"
sol_Dredge1$source <- "Point Judith Pond S 1"
sol_Dredge2$source <- "Point Judith Pond S 2"

#combine all Y1 field data into one dataframe
sol_all <- rbind(sol_Dredge1, sol_Dredge2, sol_Sled1, sol_Sled2)

##### Model Plots (Fig 3, 6, 8, 9) #####
#Figure 3: combining all irradiance forcings
plot_I_Y1 <- ggplot() + 
  geom_line(data = sol_all[sol_all$source == "Point Judith Pond N 1",], aes(Date, I), color = "gray0") +
  theme_bw() +
  xlim(as.POSIXct(c("2017-10-30 23:00:00", "2018-06-01 23:00:00"))) +
  ylim(0,4) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  labs(x= "Date (2017-2018)", y = bquote('Irradiance (mol γ m'^"-2"*' h'^"-1)"))
plot_I_Y1

#Figure 6
#Temperature y1 plot
plot_T_Y1 <- ggplot(data = sol_all, aes(Date, Temp_C, color = source)) + 
  geom_line() +
  scale_color_manual(values = c("coral", "darkgoldenrod1", "firebrick", "black")) +
  ylim(-5, 20) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.title = element_blank()) +
 # theme(legend.position="none") + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  labs(x= "Date (2017-2018)", y = "Temperature (°C)")
plot_T_Y1

#N forcing y1 plot
Sled_WSA$Date <- as.POSIXct(c("2018-04-11 23:0:0", "2018-03-20 23:0:0", "2018-02-14 23:0:0", "2017-11-01 23:0:0", "2017-11-29 23:0:0", "2018-01-24 23:0:0", "2018-04-22 23:0:0"))
Dredge_WSA$Date <- as.POSIXct(c("2018-03-20 23:0:0", "2018-04-11 23:0:0", "2018-01-24 23:0:0", "2017-11-01 23:0:0", "2018-02-14 23:0:0"))

dredge_wsa_test <- Dredge_WSA %>% select(Date, NitrateNitrite_uM) %>% 
  add_row(Date=as_datetime("2017-11-29 23:0:0"), NitrateNitrite_uM= sol_all %>% filter(Date=="2017-11-29 23:00:00" & source=="Point Judith Pond S 1") %>% pull(N)*10^6) %>% 
  add_row(Date=as_datetime("2018-04-22 23:00:00"), NitrateNitrite_uM= sol_all %>% filter(Date=="2018-04-12 3:00:00" & source=="Point Judith Pond S 2") %>% pull(N)*10^6)

sled_wsa_test <- Sled_WSA %>% select(Date, NitrateNitrite_uM) %>% filter(Date < "2018-04-20")%>% 
  add_row(Date=as_datetime("2018-04-17 12:00:00"), NitrateNitrite_uM= sol_all %>% filter(Date==as.POSIXct("2018-04-17 12:00:00", tz="GMT") & source=="Point Judith Pond N 2") %>% pull(N)*10^6)

plot_N <- ggplot() + 
  geom_point(data = sled_wsa_test, aes(Date, NitrateNitrite_uM)) +
  geom_point(data = dredge_wsa_test, aes(Date, NitrateNitrite_uM)) +
  geom_line(data = sol_all, aes(Date, N*1000000, color = source), size = 1) +
  scale_color_manual(values = c("coral","darkgoldenrod1","black","firebrick")) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  scale_x_datetime(date_breaks = "1 month", date_labels = numform::f_month)+
  ylim(0, 10) +
  labs(x= "Date (2017-2018)", y = bquote('N concentration mol' ~NO[3]^{"-"}~ 'and' ~NO[2]^{"-"}~ 'L'^"-1"), color=NULL)
plot_N

#Figure 8
plot_J_EC_R_PJ <- ggplot() +
  geom_line(data = sol_all[sol_all$source == "Point Judith Pond S 1",], aes(Date, J_EC_R, color = source)) +
  geom_line(data = sol_all[sol_all$source == "Point Judith Pond N 1",], aes(Date, J_EC_R, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2017-10-30 23:00:00", "2018-06-01 23:00:00"))) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2017-2018)", y = bquote('Rejected C (mol C mol V'^"-1"*' h'^"-1"*')')) +
  ggtitle("A)")


plot_J_EN_R_PJ <- ggplot() +
  geom_line(data = sol_all[sol_all$source == "Point Judith Pond S 1",], aes(Date, J_EN_R, color = source)) +
  geom_line(data = sol_all[sol_all$source == "Point Judith Pond N 1",], aes(Date, J_EN_R, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2017-10-30 23:00:00", "2018-06-01 23:00:00"))) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2017-2018)", y = bquote('Rejected N (mol N mol V'^"-1"*' h'^"-1"*')')) +
  ggtitle("C)")

plot_J_EC_R_PJ/plot_J_EN_R_PJ

##### Kelp Field Data Comparison plot (Figure 7) ####
#import field data
KelpY1 <- read.csv("Year1kelpdata.csv", header = TRUE, fileEncoding="UTF-8-BOM")
names(KelpY1)[2] <- "Site"
KelpY1 <- filter(KelpY1, Site != "Fox Island")
KelpY1$Date <- mdy(KelpY1$SamplingDate)
KelpY1$SiteLine <- paste(KelpY1$Site, KelpY1$Line)
KelpY1 <- filter(KelpY1, SiteLine != "Narragansett Bay N 2")


PJS1_meandat <- KelpY1[KelpY1$SiteLine == "Point Judith Pond S 1",] %>%
  group_by(Date) %>%
  summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length, na.rm = TRUE))
PJS1_meandat$Date <- as.POSIXct(PJS1_meandat$Date)
PJS1_meandat_sub <- PJS1_meandat[2:6,]
erPJS1 <- merge(PJS1_meandat_sub, sol_all[sol_all$source == "Point Judith Pond S 1",], all.x = TRUE)
PJS1_rmse <- rmse(erPJS1$mean_length, erPJS1$L_allometric)

PJS2_meandat <- KelpY1[KelpY1$SiteLine == "Point Judith Pond S 2",] %>%
  group_by(Date) %>%
  summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length, na.rm = TRUE))
PJS2_meandat$Date <- as.POSIXct(PJS2_meandat$Date)
PJS2_meandat_sub <- PJS2_meandat[2:6,]
erPJS2 <- merge(PJS2_meandat_sub, sol_all[sol_all$source == "Point Judith Pond S 2",], all.x = TRUE)
PJS2_rmse <- rmse(erPJS2$mean_length, erPJS2$L_allometric)

PJS1_2 <- ggplot() + 
  geom_point(data = PJS1_meandat, aes(Date, mean_length), shape = 'diamond', size = 3) +
  geom_errorbar(PJS1_meandat, mapping = aes(x = Date, ymin = mean_length-sd_length, ymax = mean_length+sd_length), width = 1) +
  geom_smooth(data = sol_all[sol_all$source == "Point Judith Pond S 1",], aes(Date, L_allometric, color = source)) +
  geom_label(aes(as.POSIXct("2018-04-01 01:00:00"), 250, hjust = 1), label = sprintf("RMSE: %f", PJS1_rmse)) +
  geom_point(data = PJS2_meandat, aes(Date, mean_length), color = "gray50", size = 3) +
  geom_errorbar(PJS2_meandat, mapping = aes(x = Date, ymin = mean_length-sd_length, ymax = mean_length+sd_length), width = 1, color = "gray50") +
  geom_smooth(data = sol_all[sol_all$source == "Point Judith Pond S 2",], aes(Date, L_allometric, color = source)) +
  geom_label(aes(as.POSIXct("2018-04-01 01:00:00"), 230, hjust = 1), label = sprintf("RMSE: %f", PJS2_rmse), color = "gray50") +
  ylim(0, 258) +
  xlim(as.POSIXct(c("2017-10-30 23:00:00", "2018-04-24 23:00:00"))) +
  scale_color_manual(values = c("gray0", "gray50")) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  labs(x= "Date (2017-2018)", y = "Blade length (cm)") +
  ggtitle("Point Judith Pond S") +
  theme(legend.position="none")

PJN1_meandat <- KelpY1[KelpY1$SiteLine == "Point Judith Pond N 1",] %>%
  group_by(Date) %>%
  summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length, na.rm = TRUE))
PJN1_meandat$Date <- as.POSIXct(PJN1_meandat$Date)
PJN1_meandat_sub <- PJN1_meandat[2:6,]
erPJN1 <- merge(PJN1_meandat_sub, sol_all[sol_all$source == "Point Judith Pond N 1",], all.x = TRUE)
PJN1_rmse <- rmse(erPJN1$mean_length, erPJN1$L_allometric)

PJN2_meandat <- KelpY1[KelpY1$SiteLine == "Point Judith Pond N 2",] %>%
  group_by(Date) %>%
  summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length, na.rm = TRUE))
PJN2_meandat$Date <- as.POSIXct(PJN2_meandat$Date)
PJN2_meandat_sub <- PJN2_meandat[2:6,]
erPJN2 <- merge(PJN2_meandat_sub, sol_all[sol_all$source == "Point Judith Pond N 2",], all.x = TRUE)
PJN2_rmse <- rmse(erPJN2$mean_length, erPJN2$L_allometric)


PJN1_2 <- ggplot() + 
  geom_point(data = PJN1_meandat, aes(Date, mean_length), shape = 'diamond', size = 3) +
  geom_errorbar(PJN1_meandat, mapping = aes(x = Date, ymin = mean_length-sd_length, ymax = mean_length+sd_length), width = 1) +
  geom_smooth(data = sol_all[sol_all$source == "Point Judith Pond N 1",], aes(Date, L_allometric, color = source)) +
  geom_label(aes(as.POSIXct("2018-04-24 01:00:00"), 250, hjust = 1), label = sprintf("RMSE: %f", PJN1_rmse)) +
  geom_point(data = PJN2_meandat, aes(Date, mean_length), color ="gray50", size = 3) +
  geom_errorbar(PJN2_meandat, mapping = aes(x = Date, ymin = mean_length-sd_length, ymax = mean_length+sd_length), width = 1, color = "gray50") +
  geom_smooth(data = sol_all[sol_all$source == "Point Judith Pond N 2",], aes(Date, L_allometric, color = source)) +
  geom_label(aes(as.POSIXct("2018-04-24 01:00:00"), 230, hjust = 1), label = sprintf("RMSE: %f", PJN2_rmse), color = "gray50") +
  ylim(0, 258) +
  xlim(as.POSIXct(c("2017-10-30 23:00:00", "2018-04-24 23:00:00"))) +
  scale_color_manual(values = c("gray0", "gray50")) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  labs(x= "Date (2017-2018)", y = "Blade length (cm)") +
  ggtitle("Point Judith Pond N") +
  theme(legend.position="none")

PJS1_2+PJN1_2

#### Literature data for comparison/Calibration ####
###### Nitrate uptake #######
#Espinoza and Chapman (1983) and Ahn et al. (1998)
EC1983_9C_Nuptake_StM <- read.csv("EspinozaChapman1983_Nuptake_9C_StMargaretsBay.csv", header = TRUE, fileEncoding="UTF-8-BOM")
EC1983_18C_Nuptake_StM <- read.csv("EspinozaChapman1983_Nuptake_18C_StMargaretsBay.csv", header = TRUE, fileEncoding="UTF-8-BOM")

#conversions 9C
EC1983_9C_Nuptake_StM$N <- EC1983_9C_Nuptake_StM$ResidualNitrateConcentration
EC1983_9C_Nuptake_StM$N <- round(EC1983_9C_Nuptake_StM$N, digits = 2)
EC1983_9C_Nuptake_StM$N <- EC1983_9C_Nuptake_StM$N/1000000 #microM to M
EC1983_9C_Nuptake_StM$NuptakeRate <- EC1983_9C_Nuptake_StM$NuptakeRate/1000000/w_EN #convert micro g N gDW–1 h–1 to mol N gDW–1 h–1

#conversions 18C
EC1983_18C_Nuptake_StM$N <- EC1983_18C_Nuptake_StM$ResidualNitrateConcentration
EC1983_18C_Nuptake_StM$N <- round(EC1983_18C_Nuptake_StM$N, digits = 2)
EC1983_18C_Nuptake_StM$N <- EC1983_18C_Nuptake_StM$N/1000000 #microM to M
EC1983_18C_Nuptake_StM$NuptakeRate <- EC1983_18C_Nuptake_StM$NuptakeRate/1000000/w_EN

#testing rounding
sol_EspinozaChapman1983_N_9$N <- round(sol_EspinozaChapman1983_N_9$N*1000000, digits = 3)/1000000
sol_EspinozaChapman1983_N_18$N <- round(sol_EspinozaChapman1983_N_18$N*1000000, digits = 3)/1000000

N_calibration <- ggplot() +
  geom_line(data = sol_EspinozaChapman1983_N_9, mapping = aes(x = N*1000000, y = J_EN_A*1000000, color = "Model of Espinoza and Chapman (1983) at 9°C")) +
  geom_line(data = sol_EspinozaChapman1983_N_18, mapping = aes(x = N*1000000, y = J_EN_A*1000000, color = "Model of Espinoza and Chapman (1983) at 18°C")) +
  geom_point(data = EC1983_9C_Nuptake_StM, mapping = aes(x = N*1000000, y = NuptakeRate*1000000, color="Espinoza and Chapman (1983), St. Margaret's Bay, 9°C"), size=3) +
  geom_point(data = EC1983_18C_Nuptake_StM, mapping = aes(x = N*1000000, y = NuptakeRate*1000000, color="Espinoza and Chapman (1983), St. Margaret's Bay, 18°C"), shape = 23, fill = 'grey', size=3) +
  xlim(0, 80) +
  scale_color_manual(values = c("gray60", "gray0", "gray60", "gray0")) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  labs(x= bquote('N Concentration (μmol' ~NO[3]^{"-"}~ 'L'^"-1"*')'), y = bquote('N uptake (μmol' ~NO[3]^{"-"}~ 'g DW'^"-1"*' h'^"-1"*')')) +
  ggtitle('a)')

#Error calculations
er9 <- merge(EC1983_9C_Nuptake_StM, sol_EspinozaChapman1983_N_9, all.x = TRUE)
rmse(er9$NuptakeRate, er9$J_EN_A) #3.683799e-07

er18 <- merge(EC1983_18C_Nuptake_StM, sol_EspinozaChapman1983_N_18, all.x = TRUE)
rmse(er18$NuptakeRate, er18$J_EN_A) #2.606024e-07

######## Photosynthesis related ####
#Johansson2002
Johansson2002 <- read.csv("Johansson2002.csv", header = TRUE, fileEncoding="UTF-8-BOM")
#conversions
Johansson2002$Irradiance <- Johansson2002$Irradiance*3600*1e-6 #micromol photons m-2 s-1 to mol photons m-2 h-1
Johansson2002$O2production <- Johansson2002$O2production/1e+6*32/1000*3600 #micromol O2 kg DW-1 s-1 to g O2/g/h
Johansson2002$O2productionSHIFT <- Johansson2002$O2production + 0.001720976 #from net to gross

Photosynthesis_calibration <- ggplot(data = Johansson2002) +
  geom_line(data = sol_Johansson2002, mapping = aes(x = I, y = J_O*1000)) +
  geom_point(mapping = aes(x = Irradiance, y = O2productionSHIFT*1000), size = 3) +
  scale_color_grey() +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  labs(x= bquote('Irradiance (mol γ m'^"-2"*' h'^"-1)"), y = bquote('Oxygen production (mg' ~O[2]~ 'g DW'^"-1"*' h'^"-1"*')')) +
  ggtitle('b)')

#error calculations
Johansson2002$I <- round(Johansson2002$Irradiance, digits = 6)
sol_Johansson2002$I <- round(sol_Johansson2002$I, digits = 6)
erPhoto <- merge(Johansson2002, sol_Johansson2002, all.x = TRUE)
rmse(erPhoto$O2productionSHIFT, erPhoto$J_O)

######## Combine calibration plot (Figure 5) #######
#Figure 5
grid.arrange(N_calibration, Photosynthesis_calibration, ncol=2)
