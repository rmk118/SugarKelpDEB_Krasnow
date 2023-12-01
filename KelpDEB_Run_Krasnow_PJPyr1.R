### Point Judith Point Year 1 ################################################################################################################# 

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
library(Metrics)
library(patchwork)

#Required for model runs
source("SolveR_R.R")
source("KelpDEB_model.R")
source("./outdoorExpt/outdoorHOBO/outdoor_HOBO.R")
source("./outdoorExpt/outdoor_expt.R")

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

# New (refit from literature)
new_params_for_model <- params_Lo
new_params_for_model[c("T_A", "T_H", "T_AH")]<-c(new_T_A, new_T_H, new_T_AH)

# Means
high_params_for_model <- params_Lo
med_params_for_model <- params_Lo
low_params_for_model <- params_Lo

high_params_for_model[c("T_A", "T_H", "T_AH")]<-params %>% filter(type=="ctrl", level=="high", res=="means") %>% select(T_A, T_H, T_AH)
med_params_for_model[c("T_A", "T_H", "T_AH")]<-params %>% filter(type=="ctrl", level=="med", res=="means") %>% select(T_A, T_H, T_AH)
low_params_for_model[c("T_A", "T_H", "T_AH")]<-params %>% filter(type=="ctrl", level=="low", res=="means") %>% select(T_A, T_H, T_AH)

# Crosses
high_params_for_model_cross <- params_Lo
med_params_for_model_cross <- params_Lo
low_params_for_model_cross <- params_Lo

high_params_for_model_cross[c("T_A", "T_H", "T_AH")]<-params %>% filter(type=="ctrl", level=="high", res=="all") %>% select(T_A, T_H, T_AH)
med_params_for_model_cross[c("T_A", "T_H", "T_AH")]<-params %>% filter(type=="ctrl", level=="med", res=="all") %>% select(T_A, T_H, T_AH)
low_params_for_model_cross[c("T_A", "T_H", "T_AH")]<-params %>% filter(type=="ctrl", level=="low", res=="all") %>% select(T_A, T_H, T_AH)

# Replicates
high_params_for_model_rep <- params_Lo
med_params_for_model_rep <- params_Lo
low_params_for_model_rep <- params_Lo

high_params_for_model_rep[c("T_A", "T_H", "T_AH")]<-params %>% filter(type=="ctrl", level=="high", res=="all") %>% select(T_A, T_H, T_AH)
med_params_for_model_rep[c("T_A", "T_H", "T_AH")]<-params %>% filter(type=="ctrl", level=="med", res=="all") %>% select(T_A, T_H, T_AH)
low_params_for_model_rep[c("T_A", "T_H", "T_AH")]<-params %>% filter(type=="ctrl", level=="low", res=="all") %>% select(T_A, T_H, T_AH)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####### Initial conditions ############
#Initial conditions of state variables
#these values are not coming from any field data or literature information, estimated
state_Lo <- c(m_EC = 0.002, #0.1, #mol C/molM_V  #Reserve density of C reserve (initial mass of C reserve per initial mass of structure)
              m_EN = 0.01, #mol N/molM_V #Reserve density of N reserve (initial mass of N reserve per initial mass of structure)
              M_V = 0.05/(w_V+0.01*w_EN+0.002*w_EC)) #molM_V #initial mass of structure

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

# Judith N
N_sled <- Sled_WSA[c("Date","NitrateNitrite_uM")] #subset
N_sled$NitrateNitrite_uM <- N_sled$NitrateNitrite_uM/1000000

# Judith S
N_dredge <- Dredge_WSA[c("Date","NitrateNitrite_uM")] #new dataframe with the relevant columns
N_dredge$NitrateNitrite_uM <- N_dredge$NitrateNitrite_uM/1000000 #convert from micromoles/L to moles/L

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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### FIELD DATA MODEL RUNS ####

### Judith N (sled) line 1 ####
W <- 0.05 #inital biomass for conversions (cannot put in initial conditions)

#Converted to hourly by multiply by 24
N_field <- approxfun(x = c(161*24, 139*24, 105*24, 0, 28*24, 84*24, 172*24), y = N_sled$NitrateNitrite_uM, method = "linear", rule = 2) #N forcing function

###### Irradiance set-up Judith N 1 #######
NOAA_Irradiance_Sledy1 <-  NOAA_Irradiance$PAR[2438:3774] # subset based on as_datetime("2017-11-1 12:00:00"), as_datetime("2018-04-17 12:00:00")
I_field <- approxfun(x = seq(from = 0, to = 4008, by = 3), y = NOAA_Irradiance_Sledy1, method = "linear", rule = 2) #irradiance forcing

######  Temp set-up Judith N 1 #############
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

###### Model runs standard ######
# (the differential equation solver)
#T_field <- approxfun(x = c(0:4008), y = c(AvgTempKbyhr_part1+2, fd, AvgTempKbyhr_part2+2), method = "linear", rule = 2)

sol_Sled1 <- ode(y = state_Lo, t = times_Lo_Sled1, func = rates_Lo, parms = params_Lo)

W <- 0.05
sol_new_Sled1 <- ode(y = state_Lo, t = times_Lo_Sled1, func = rates_Lo, parms = new_params_for_model)
W <- 0.05
sol_high_Sled1 <- ode(y = state_Lo, t = times_Lo_Sled1, func = rates_Lo, parms = high_params_for_model)
W <- 0.05
sol_med_Sled1 <- ode(y = state_Lo, t = times_Lo_Sled1, func = rates_Lo, parms = med_params_for_model)
W <- 0.05
sol_low_Sled1 <- ode(y = state_Lo, t = times_Lo_Sled1, func = rates_Lo, parms = low_params_for_model)

W <- 0.05
sol_high_Sled1_cross <- ode(y = state_Lo, t = times_Lo_Sled1, func = rates_Lo, parms = high_params_for_model_cross)
W <- 0.05
sol_med_Sled1_cross <- ode(y = state_Lo, t = times_Lo_Sled1, func = rates_Lo, parms = med_params_for_model_cross)
W <- 0.05
sol_low_Sled1_cross <- ode(y = state_Lo, t = times_Lo_Sled1, func = rates_Lo, parms = low_params_for_model_cross)

W <- 0.05
sol_high_Sled1_rep <- ode(y = state_Lo, t = times_Lo_Sled1, func = rates_Lo, parms = high_params_for_model_rep)
W <- 0.05
sol_med_Sled1_rep <- ode(y = state_Lo, t = times_Lo_Sled1, func = rates_Lo, parms = med_params_for_model_rep)
W <- 0.05
sol_low_Sled1_rep <- ode(y = state_Lo, t = times_Lo_Sled1, func = rates_Lo, parms = low_params_for_model_rep)


###### Model runs +2°C ######
T_field <- approxfun(x = c(0:4008), y = c(AvgTempKbyhr_part1+2, fd+2, AvgTempKbyhr_part2+2), method = "linear", rule = 2)
T_Sled1_plus2 <- T_field(0:4008) #saving the forcing this way for ease of later visualization

W <- 0.05
sol_Sled1_plus2 <- ode(y = state_Lo, t = times_Lo_Sled1, func = rates_Lo, parms = params_Lo)
W <- 0.05
sol_new_Sled1_plus2 <- ode(y = state_Lo, t = times_Lo_Sled1, func = rates_Lo, parms = new_params_for_model)
W <- 0.05
sol_high_Sled1_plus2 <- ode(y = state_Lo, t = times_Lo_Sled1, func = rates_Lo, parms = high_params_for_model)
W <- 0.05
sol_med_Sled1_plus2 <- ode(y = state_Lo, t = times_Lo_Sled1, func = rates_Lo, parms = med_params_for_model)
W <- 0.05
sol_low_Sled1_plus2 <- ode(y = state_Lo, t = times_Lo_Sled1, func = rates_Lo, parms = low_params_for_model)

W <- 0.05
sol_high_Sled1_plus2_cross <- ode(y = state_Lo, t = times_Lo_Sled1, func = rates_Lo, parms = high_params_for_model_cross)
W <- 0.05
sol_med_Sled1_plus2_cross <- ode(y = state_Lo, t = times_Lo_Sled1, func = rates_Lo, parms = med_params_for_model_cross)
W <- 0.05
sol_low_Sled1_plus2_cross <- ode(y = state_Lo, t = times_Lo_Sled1, func = rates_Lo, parms = low_params_for_model_cross)

W <- 0.05
sol_high_Sled1_plus2_rep <- ode(y = state_Lo, t = times_Lo_Sled1, func = rates_Lo, parms = high_params_for_model_rep)
W <- 0.05
sol_med_Sled1_plus2_rep <- ode(y = state_Lo, t = times_Lo_Sled1, func = rates_Lo, parms = med_params_for_model_rep)
W <- 0.05
sol_low_Sled1_plus2_rep <- ode(y = state_Lo, t = times_Lo_Sled1, func = rates_Lo, parms = low_params_for_model_rep)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setting up the forcing functions with field data for
### Judith N (sled) line 2 ####
W <- 0.05 #inital biomass for conversions (cannot put in initial conditions)

# N forcing set-up Judith N 2 #
N_field <- approxfun(x = c(133*24, 111*24, 77*24, -28*24, 0, 56*24, 144*24), y = N_sled$NitrateNitrite_uM, method = "linear", rule = 2) #N forcing function

###### Irradiance set-up Judith N 2 ####
NOAA_Irradiance_Sledy1_L2 <-  NOAA_Irradiance$PAR[2662:3774] #subset by seq(as_datetime("2017-11-29 12:00:00"), as_datetime("2018-04-17 12:00:00"), by="hour")
I_field <- approxfun(x = seq(from = 0, to = 3336, by = 3), y = NOAA_Irradiance_Sledy1_L2, method = "linear", rule = 2) #irradiance forcing function

###### Temp set-Up Judith N 2 #############
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

##### Model runs ####
sol_Sled2 <- ode(y = state_Lo, t = times_Lo_Sled2, func = rates_Lo, parms = params_Lo)
W <- 0.05
sol_new_Sled2 <- ode(y = state_Lo, t = times_Lo_Sled2, func = rates_Lo, parms = new_params_for_model)
W <- 0.05
sol_high_Sled2 <- ode(y = state_Lo, t = times_Lo_Sled2, func = rates_Lo, parms = high_params_for_model)
W <- 0.05
sol_med_Sled2 <- ode(y = state_Lo, t = times_Lo_Sled2, func = rates_Lo, parms = med_params_for_model)
W <- 0.05
sol_low_Sled2 <- ode(y = state_Lo, t = times_Lo_Sled2, func = rates_Lo, parms = low_params_for_model)

W <- 0.05
sol_high_Sled2_cross <- ode(y = state_Lo, t = times_Lo_Sled2, func = rates_Lo, parms = high_params_for_model_cross)
W <- 0.05
sol_med_Sled2_cross <- ode(y = state_Lo, t = times_Lo_Sled2, func = rates_Lo, parms = med_params_for_model_cross)
W <- 0.05
sol_low_Sled2_cross <- ode(y = state_Lo, t = times_Lo_Sled2, func = rates_Lo, parms = low_params_for_model_cross)

W <- 0.05
sol_high_Sled2_rep <- ode(y = state_Lo, t = times_Lo_Sled2, func = rates_Lo, parms = high_params_for_model_rep)
W <- 0.05
sol_med_Sled2_rep <- ode(y = state_Lo, t = times_Lo_Sled2, func = rates_Lo, parms = med_params_for_model_rep)
W <- 0.05
sol_low_Sled2_rep <- ode(y = state_Lo, t = times_Lo_Sled2, func = rates_Lo, parms = low_params_for_model_rep)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setting up the forcing functions with field data for
### Judith S (dredge) line 1 ####

# N forcing set-up Judith S 1 #
N_field <- approxfun(x = c(139*24, 161*24, 84*24, 0, 105*24), y = N_dredge$NitrateNitrite_uM, method = "linear", rule = 2)

###### Irradiance set-up Judith S 1 ####
NOAA_Irradiance_Dredgey1 <-  NOAA_Irradiance$PAR[2438:3814] #subset by 11/1/17 to 2018-04-22 12:00:00
I_field <- approxfun(x = seq(from = 0, to = 4128, by = 3), y = NOAA_Irradiance_Dredgey1, method = "linear", rule = 2) #irradiance forcing function

###### Temp set-Up Judith S 1 #############
Dredge_Y1_hobo <- Dredge_Y1_hobo_orig[3:16531,] #subset
DredgeT_hourly <- ceiling_date(Dredge_Y1_hobo$DateTime, unit = "hour") #set values to aggregate around
AvgTempKbyhr <- aggregate(Dredge_Y1_hobo$Temp_K, by=list(DredgeT_hourly), mean) #calculate average hourly temp
AvgTempKbyhr <- AvgTempKbyhr[4:4132, ] #subset
T_field <- approxfun(x = c(0:4128), y = AvgTempKbyhr$x, method = "linear", rule = 2) #the temp forcing function
T_Dredge1_Y1 <- T_field(0:4128) #for ease in later plotting of the forcing

#### Model runs ####

W <- 0.05 #initial biomass for conversions (cannot put in initial conditions)
sol_Dredge1 <- ode(y = state_Lo, t = times_Lo_Dredge1, func = rates_Lo, parms = params_Lo)
W <- 0.05
sol_new_Dredge1 <- ode(y = state_Lo, t = times_Lo_Dredge1, func = rates_Lo, parms = new_params_for_model)
W <- 0.05
sol_high_Dredge1 <- ode(y = state_Lo, t = times_Lo_Dredge1, func = rates_Lo, parms = high_params_for_model)
W <- 0.05
sol_med_Dredge1 <- ode(y = state_Lo, t = times_Lo_Dredge1, func = rates_Lo, parms = med_params_for_model)
W <- 0.05
sol_low_Dredge1 <- ode(y = state_Lo, t = times_Lo_Dredge1, func = rates_Lo, parms = low_params_for_model)

W <- 0.05
sol_high_Dredge1_cross <- ode(y = state_Lo, t = times_Lo_Dredge1, func = rates_Lo, parms = high_params_for_model_cross)
W <- 0.05
sol_med_Dredge1_cross <- ode(y = state_Lo, t = times_Lo_Dredge1, func = rates_Lo, parms = med_params_for_model_cross)
W <- 0.05
sol_low_Dredge1_cross <- ode(y = state_Lo, t = times_Lo_Dredge1, func = rates_Lo, parms = low_params_for_model_cross)

W <- 0.05
sol_high_Dredge1_rep <- ode(y = state_Lo, t = times_Lo_Dredge1, func = rates_Lo, parms = high_params_for_model_rep)
W <- 0.05
sol_med_Dredge1_rep <- ode(y = state_Lo, t = times_Lo_Dredge1, func = rates_Lo, parms = med_params_for_model_rep)
W <- 0.05
sol_low_Dredge1_rep <- ode(y = state_Lo, t = times_Lo_Dredge1, func = rates_Lo, parms = low_params_for_model_rep)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setting up the forcing functions with field data for
### Judith S line 2 ####

# N forcing set-up Judith S 2
N_field <- approxfun(x = c(111*24, 133*24, 56*24, -28*24, 77*24), y = N_dredge$NitrateNitrite_uM, method = "linear", rule = 2) #N forcing function

###### Irradiance set-up Judith S 2 ####
NOAA_Irradiance_Dredgey1_L2 <-  NOAA_Irradiance$PAR[2662:3814] #subset by seq(as_datetime("2017-11-29 12:00:00"), as_datetime("2018-04-22 12:00:00"), by="hour")
I_field <- approxfun(x = seq(from = 0, to = 3456, by = 3), y = NOAA_Irradiance_Dredgey1_L2, method = "linear", rule = 2) #irradiance forcing function

###### Temp set-up Judith S 2 #############
Dredge_Y1_hobo <- Dredge_Y1_hobo_orig[3:16531,] #cut 2 points in beginning, logger not yet in water
DredgeT_hourly <- ceiling_date(Dredge_Y1_hobo$DateTime, unit = "hour") #determine dates to aggregate around
AvgTempKbyhr <- aggregate(Dredge_Y1_hobo$Temp_K, by=list(DredgeT_hourly), mean) #calculate average hourly temp
AvgTempKbyhr_sub <- AvgTempKbyhr[676:4132,] #subset
T_field <- approxfun(x = c(0:3456), y = c(AvgTempKbyhr_sub$x), method = "linear", rule = 2) #the temp forcing function
T_Dredge2_Y1 <- T_field(0:3456) #for ease of later plotting the temperature forcing

#### Model runs ####

W <- 0.05 #initial biomass for conversions (cannot put in initial conditions)
sol_Dredge2 <- ode(y = state_Lo, t = times_Lo_Dredge2, func = rates_Lo, parms = params_Lo)
W <- 0.05
sol_new_Dredge2 <- ode(y = state_Lo, t = times_Lo_Dredge2, func = rates_Lo, parms = new_params_for_model)
W <- 0.05
sol_high_Dredge2 <- ode(y = state_Lo, t = times_Lo_Dredge2, func = rates_Lo, parms = high_params_for_model)
W <- 0.05
sol_med_Dredge2 <- ode(y = state_Lo, t = times_Lo_Dredge2, func = rates_Lo, parms = med_params_for_model)
W <- 0.05
sol_low_Dredge2 <- ode(y = state_Lo, t = times_Lo_Dredge2, func = rates_Lo, parms = low_params_for_model)

W <- 0.05
sol_high_Dredge2_cross <- ode(y = state_Lo, t = times_Lo_Dredge2, func = rates_Lo, parms = high_params_for_model_cross)
W <- 0.05
sol_med_Dredge2_cross <- ode(y = state_Lo, t = times_Lo_Dredge2, func = rates_Lo, parms = med_params_for_model_cross)
W <- 0.05
sol_low_Dredge2_cross <- ode(y = state_Lo, t = times_Lo_Dredge2, func = rates_Lo, parms = low_params_for_model_cross)

W <- 0.05
sol_high_Dredge2_rep <- ode(y = state_Lo, t = times_Lo_Dredge2, func = rates_Lo, parms = high_params_for_model_rep)
W <- 0.05
sol_med_Dredge2_rep <- ode(y = state_Lo, t = times_Lo_Dredge2, func = rates_Lo, parms = med_params_for_model_rep)
W <- 0.05
sol_low_Dredge2_rep <- ode(y = state_Lo, t = times_Lo_Dredge2, func = rates_Lo, parms = low_params_for_model_rep)


###### Convert DeSolve solutions into data frame for broader plotting use ####

clean_ode_sol <- function(sol, date_seq, temp_seq, source){
  df <- as.data.frame(sol) #convert to data frame
  df <- df %>% mutate(Temp_C = temp_seq-273.15, #conversion back to Celsius from Kelvin
                      Date=date_seq,
                      source=source) #addition of a date variable
  return(df)
}

sled1_date_seq <- seq(as_datetime("2017-11-1 12:00:00"), as_datetime("2018-04-17 12:00:00"), by="hour")
sled2_date_seq <- seq(as_datetime("2017-11-29 12:00:00"), as_datetime("2018-04-17 12:00:00"), by="hour")
dredge1_date_seq<- seq(as_datetime("2017-11-1 12:00:00"), as_datetime("2018-04-22 12:00:00"), by="hour")
dredge2_date_seq<-seq(as_datetime("2017-11-29 12:00:00"), as_datetime("2018-04-22 12:00:00"), by="hour")

#Example applied to one deSolve solution
#sol_Sled1 <- clean_ode_sol(sol_Sled1, sled1_date_seq, T_Sled1_Y1)

sled1_sols <- list(orig=sol_Sled1, new=sol_new_Sled1, high=sol_high_Sled1, med=sol_med_Sled1, low=sol_low_Sled1, high_cross=sol_high_Sled1_cross, med_cross=sol_med_Sled1_cross, low_cross=sol_low_Sled1_cross, high_rep=sol_high_Sled1_rep, med_rep=sol_med_Sled1_rep, low_rep=sol_low_Sled1_rep)

sled2_sols <- list(orig=sol_Sled2, new=sol_new_Sled2, high=sol_high_Sled2, med=sol_med_Sled2, low=sol_low_Sled2, high_cross=sol_high_Sled2_cross, med_cross=sol_med_Sled2_cross, low_cross=sol_low_Sled2_cross, high_rep=sol_high_Sled2_rep, med_rep=sol_med_Sled2_rep, low_rep=sol_low_Sled2_rep)

dredge1_sols <- list(orig=sol_Dredge1, new=sol_new_Dredge1, high=sol_high_Dredge1, med=sol_med_Dredge1, low=sol_low_Dredge1, high_cross=sol_high_Dredge1_cross, med_cross=sol_med_Dredge1_cross, low_cross=sol_low_Dredge1_cross, high_rep=sol_high_Dredge1_rep, med_rep=sol_med_Dredge1_rep, low_rep=sol_low_Dredge1_rep)

dredge2_sols <- list(orig=sol_Dredge2, new=sol_new_Dredge2, high=sol_high_Dredge2, med=sol_med_Dredge2, low=sol_low_Dredge2, high_cross=sol_high_Dredge2_cross, med_cross=sol_med_Dredge2_cross, low_cross=sol_low_Dredge2_cross, high_rep=sol_high_Dredge2_rep, med_rep=sol_med_Dredge2_rep, low_rep=sol_low_Dredge2_rep)

sled1_sols_plus2 <- list(orig=sol_Sled1_plus2, new=sol_new_Sled1_plus2, high=sol_high_Sled1_plus2, med=sol_med_Sled1_plus2, low=sol_low_Sled1_plus2, high_rep=sol_high_Sled1_plus2_rep, med_rep=sol_med_Sled1_plus2_rep, low_rep=sol_low_Sled1_plus2_rep)

sled1_sols<- map(sled1_sols, clean_ode_sol, sled1_date_seq, T_Sled1_Y1, "Point Judith Pond N 1") %>% bind_rows(.id = "params")
sled2_sols<- map(sled2_sols, clean_ode_sol, sled2_date_seq, T_Sled2_Y1, "Point Judith Pond N 2") %>% bind_rows(.id = "params")
dredge1_sols<- map(dredge1_sols, clean_ode_sol, dredge1_date_seq, T_Dredge1_Y1, "Point Judith Pond S 1") %>% bind_rows(.id = "params")
dredge2_sols<- map(dredge2_sols, clean_ode_sol, dredge2_date_seq, T_Dredge2_Y1, "Point Judith Pond S 2") %>% bind_rows(.id = "params")

sled1_sols_plus2<- map(sled1_sols_plus2, clean_ode_sol, sled1_date_seq, T_Sled1_plus2, "Point Judith Pond N 1") %>% bind_rows(.id = "params")

#combine all original field data into one dataframe
sol_all_orig <- rbind(dredge1_sols %>% filter(params=="orig"),
                 dredge2_sols %>% filter(params=="orig"),
                 sled1_sols %>% filter(params=="orig"),
                 sled2_sols %>% filter(params=="orig"))

##### Model Plots (Fig 3, 6, 8, 9) #####
#Figure 3: combining all irradiance forcings
plot_I_Y1 <- ggplot() + 
  geom_line(data = sol_all_orig[sol_all_orig$source == "Point Judith Pond N 1",], aes(Date, I), color = "gray0") +
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
plot_T_Y1 <- ggplot(data = sol_all_orig, aes(Date, Temp_C, color = source)) + 
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
  add_row(Date=as_datetime("2017-11-29 23:0:0"), NitrateNitrite_uM= sol_all_orig %>% filter(Date=="2017-11-29 23:00:00" & source=="Point Judith Pond S 1") %>% pull(N)*10^6) %>% 
  add_row(Date=as_datetime("2018-04-22 23:00:00"), NitrateNitrite_uM= sol_all_orig %>% filter(Date=="2018-04-12 3:00:00" & source=="Point Judith Pond S 2") %>% pull(N)*10^6)

sled_wsa_test <- Sled_WSA %>% select(Date, NitrateNitrite_uM) %>% filter(Date < "2018-04-20")%>% 
  add_row(Date=as_datetime("2018-04-17 12:00:00"), NitrateNitrite_uM= sol_all_orig %>% filter(Date==as.POSIXct("2018-04-17 12:00:00", tz="GMT") & source=="Point Judith Pond N 2") %>% pull(N)*10^6)

plot_N <- ggplot() + 
  geom_point(data = sled_wsa_test, aes(Date, NitrateNitrite_uM)) +
  geom_point(data = dredge_wsa_test, aes(Date, NitrateNitrite_uM)) +
  geom_line(data = sol_all_orig, aes(Date, N*1000000, color = source), size = 1) +
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
  geom_line(data = sol_all_orig[sol_all_orig$source == "Point Judith Pond S 1",], aes(Date, J_EC_R, color = source)) +
  geom_line(data = sol_all_orig[sol_all_orig$source == "Point Judith Pond N 1",], aes(Date, J_EC_R, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2017-10-30 23:00:00", "2018-06-01 23:00:00"))) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2017-2018)", y = bquote('Rejected C (mol C mol V'^"-1"*' h'^"-1"*')')) +
  ggtitle("A)")


plot_J_EN_R_PJ <- ggplot() +
  geom_line(data = sol_all_orig[sol_all_orig$source == "Point Judith Pond S 1",], aes(Date, J_EN_R, color = source)) +
  geom_line(data = sol_all_orig[sol_all_orig$source == "Point Judith Pond N 1",], aes(Date, J_EN_R, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2017-10-30 23:00:00", "2018-06-01 23:00:00"))) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2017-2018)", y = bquote('Rejected N (mol N mol V'^"-1"*' h'^"-1"*')')) +
  ggtitle("C)")

plot_J_EC_R_PJ/plot_J_EN_R_PJ

### Kelp Field Data Comparison plot (Figure 7) ####

#import field data
KelpY1 <- read.csv("Year1kelpdata.csv", header = TRUE, fileEncoding="UTF-8-BOM")
names(KelpY1)[2] <- "Site"
KelpY1 <- filter(KelpY1, Site != "Fox Island")
KelpY1$Date <- mdy(KelpY1$SamplingDate)
KelpY1$SiteLine <- paste(KelpY1$Site, KelpY1$Line)
KelpY1 <- filter(KelpY1, SiteLine != "Narragansett Bay N 2")

###### Combine with model data #####

#Point Judith Point N (sled) Line 1
PJN1_meandat <- KelpY1[KelpY1$SiteLine == "Point Judith Pond N 1",] %>%
  group_by(Date) %>%
  summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length, na.rm = TRUE))
PJN1_meandat$Date <- as.POSIXct(PJN1_meandat$Date)
PJN1_meandat_sub <- PJN1_meandat[2:6,]
erPJN1 <- merge(PJN1_meandat_sub, sol_all_orig[sol_all_orig$source == "Point Judith Pond N 1",], all.x = TRUE)
PJN1_rmse <- rmse(erPJN1$mean_length, erPJN1$L_allometric)

#Point Judith Point N (sled) Line 2
PJN2_meandat <- KelpY1[KelpY1$SiteLine == "Point Judith Pond N 2",] %>%
  group_by(Date) %>%
  summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length, na.rm = TRUE))
PJN2_meandat$Date <- as.POSIXct(PJN2_meandat$Date)
PJN2_meandat_sub <- PJN2_meandat[2:6,]
erPJN2 <- merge(PJN2_meandat_sub, sol_all_orig[sol_all_orig$source == "Point Judith Pond N 2",], all.x = TRUE)
PJN2_rmse <- rmse(erPJN2$mean_length, erPJN2$L_allometric)

#Point Judith Point S (dredge) Line 1
PJS1_meandat <- KelpY1[KelpY1$SiteLine == "Point Judith Pond S 1",] %>%
  group_by(Date) %>%
  summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length, na.rm = TRUE))
PJS1_meandat$Date <- as.POSIXct(PJS1_meandat$Date)
PJS1_meandat_sub <- PJS1_meandat[2:6,]
erPJS1 <- merge(PJS1_meandat_sub, sol_all_orig[sol_all_orig$source == "Point Judith Pond S 1",], all.x = TRUE)
# erPJS1 <- merge(PJS1_meandat_sub, dredge1_sols %>% filter(source=="Point Judith Pond S 1", params=="low"), all.x = TRUE)
PJS1_rmse <- rmse(erPJS1$mean_length, erPJS1$L_allometric)

#Point Judith Point S (dredge) Line 2
PJS2_meandat <- KelpY1[KelpY1$SiteLine == "Point Judith Pond S 2",] %>%
  group_by(Date) %>%
  summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length, na.rm = TRUE))
PJS2_meandat$Date <- as.POSIXct(PJS2_meandat$Date)
PJS2_meandat_sub <- PJS2_meandat[2:6,]
erPJS2 <- merge(PJS2_meandat_sub, sol_all_orig[sol_all_orig$source == "Point Judith Pond S 2",], all.x = TRUE)
erPJS2 <- merge(PJS2_meandat_sub, dredge2_sols %>% filter(source=="Point Judith Pond S 2", params=="low"), all.x = TRUE)
PJS2_rmse <- rmse(erPJS2$mean_length, erPJS2$L_allometric)

###### Plots #####
PJS1_2 <- ggplot() + 
  geom_point(data = PJS1_meandat, aes(Date, mean_length), shape = 'diamond', size = 3) +
  geom_errorbar(PJS1_meandat, mapping = aes(x = Date, ymin = mean_length-sd_length, ymax = mean_length+sd_length), width = 1) +
  geom_smooth(data = sol_all_orig[sol_all_orig$source == "Point Judith Pond S 1",], aes(Date, L_allometric, color = source)) +
  geom_label(aes(as.POSIXct("2018-04-01 01:00:00"), 250, hjust = 1), label = sprintf("RMSE: %f", PJS1_rmse)) +
  geom_point(data = PJS2_meandat, aes(Date, mean_length), color = "gray50", size = 3) +
  geom_errorbar(PJS2_meandat, mapping = aes(x = Date, ymin = mean_length-sd_length, ymax = mean_length+sd_length), width = 1, color = "gray50") +
  geom_smooth(data = sol_all_orig[sol_all_orig$source == "Point Judith Pond S 2",], aes(Date, L_allometric, color = source)) +
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

PJN1_2 <- ggplot() + 
  geom_point(data = PJN1_meandat, aes(Date, mean_length), shape = 'diamond', size = 3) +
  geom_errorbar(PJN1_meandat, mapping = aes(x = Date, ymin = mean_length-sd_length, ymax = mean_length+sd_length), width = 1) +
  geom_smooth(data = sol_all_orig[sol_all_orig$source == "Point Judith Pond N 1",], aes(Date, L_allometric, color = source)) +
  geom_label(aes(as.POSIXct("2018-04-24 01:00:00"), 250, hjust = 1), label = sprintf("RMSE: %f", PJN1_rmse)) +
  geom_point(data = PJN2_meandat, aes(Date, mean_length), color ="gray50", size = 3) +
  geom_errorbar(PJN2_meandat, mapping = aes(x = Date, ymin = mean_length-sd_length, ymax = mean_length+sd_length), width = 1, color = "gray50") +
  geom_smooth(data = sol_all_orig[sol_all_orig$source == "Point Judith Pond N 2",], aes(Date, L_allometric, color = source)) +
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

#combine all original field data into one dataframe
sol_all <- rbind(dredge1_sols,dredge2_sols,sled1_sols, sled2_sols)

ggplot() +
  geom_smooth(data = sol_all, aes(Date, L_allometric, color = params), se=FALSE)+
  labs(x= "Date", y = "Blade length (cm)")+
  facet_wrap(~source)
