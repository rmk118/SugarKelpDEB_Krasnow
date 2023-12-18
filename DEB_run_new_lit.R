### ALL SITES AND YEARS ################################################################################################################# 

#Site names here begin with names other than those used in the manuscript
#Sled = Pt Judith Pond N
#Dredge = Pt Judith Pond S
#Wickford = Narragansett Bay N
#Rome Point = Narragansett Bay S

#Some material originally created by Celeste Venolia in March 2018-December 2019 for Venolia et al., (2020)
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
library(furrr) #parallel processing version of purrr, to speed up model runs
library(tidytext)
library(scales)
library(clipr)
library(rstatix)

#Required for model runs
source("SolveR_R.R")
source("KelpDEB_model.R")
source("./new_lit.R")

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

params_nested <- params_new_lit %>% nest(data = c(T_A, T_H, T_AH))
### Year 1 ################################################################################################################# 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###### Initial conditions year 1 ############
#Initial conditions of state variables
#these values are not coming from any field data or literature information, estimated
state_Lo <- c(m_EC = 0.002, #0.1, #mol C/molM_V  #Reserve density of C reserve (initial mass of C reserve per initial mass of structure)
              m_EN = 0.01, #mol N/molM_V #Reserve density of N reserve (initial mass of N reserve per initial mass of structure)
              M_V = 0.05/(w_V+0.01*w_EN+0.002*w_EC)) #molM_V #initial mass of structure

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####### Time steps year 1 #######
#(First number of time step, last number of time step, interval to step)
times_Lo_Sled1 <- seq(0, 4008, 1) #167 days stepped hourly - Pt Judith Pond N 1
times_Lo_Sled2 <- seq(0, 3336, 1) #139 days stepped hourly - Pt Judith Pond N 2
times_Lo_Dredge1 <- seq(0, 4128, 1) #172 days stepped hourly - Pt Judith Pond S 1
times_Lo_Dredge2 <- seq(0, 3456, 1) #144 days stepped hourly - Pt Judith Pond S 2
times_Y1_W <- seq(0, 3312, 1) #138 days stepped hourly - Narragansett Bay N
times_Y1_R1 <- seq(0, 4104, 1) #171 days stepped hourly - Narragansett Bay S 1
times_Y1_R2 <- seq(0, 3264, 1) #136 days stepped hourly - Narragansett Bay S 2

sled1_date_seq <- seq(as_datetime("2017-11-1 12:00:00"), as_datetime("2018-04-17 12:00:00"), by="hour")
sled2_date_seq <- seq(as_datetime("2017-11-29 12:00:00"), as_datetime("2018-04-17 12:00:00"), by="hour")
dredge1_date_seq<- seq(as_datetime("2017-11-1 12:00:00"), as_datetime("2018-04-22 12:00:00"), by="hour")
dredge2_date_seq<-seq(as_datetime("2017-11-29 12:00:00"), as_datetime("2018-04-22 12:00:00"), by="hour")
W_date_seq_Y1 <- seq(as_datetime("2017-12-4 12:00:00"), as_datetime("2018-04-21 12:00:00"), by="hour")
R1_date_seq_Y1 <- seq(as_datetime("2017-11-1 12:00:00"), as_datetime("2018-04-21 12:00:00"), by="hour")
R2_date_seq_Y1 <-seq(as_datetime("2017-12-6 12:00:00"), as_datetime("2018-04-21 12:00:00"), by="hour")

#### Irradiance set-up year 1 ####
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

##### PJP N set-up year 1 ####
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

##### PJP DIC set-up year 1###########
DIC <- read.csv("Ninigret_EPA_DIC.csv", header = TRUE, fileEncoding="UTF-8-BOM") #Import Ninigret DIC data
CO_2 <- mean(DIC$DIC.uMkg.mean) #micromole DIC/kg (Jason said it was okay to assume that 1kg of seawater is 1L of seawater (actual conversion requires density calc from salinity and T))
#need units to match K_C (molDIC/L)
CO_2 <- CO_2/1000000

##### PJP temp set-up year 1 ###########
Sled_Y1_hobotemp_orig <- read.csv("Sled_Y1_TempLogger2.csv", header = TRUE, fileEncoding="UTF-8-BOM") #import
Sled_Y1_hobotemp_orig$DateTime <- mdy_hms(Sled_Y1_hobotemp_orig$Date_Time) #convert time field
Sled_Y1_hobotemp_orig$Temp_K <- Sled_Y1_hobotemp_orig$Temp_C+273.15 #create column with temp in K

Dredge_Y1_hobo_orig <- read.csv("Dredge_Y1_hobo.csv", header = TRUE, fileEncoding="UTF-8-BOM") #importing neighboring temp file to replace corrupted section
Dredge_Y1_hobo_orig$DateTime <- mdy_hms(Dredge_Y1_hobo_orig$Date_Time) #convert time field
Dredge_Y1_hobo_orig$Temp_K <- Dredge_Y1_hobo_orig$Temp_C+273.15 #create column with temp in K


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Point Judith Year 1 ####

#### Sled line 1 Year 1####
#initial biomass for conversions (cannot put in initial conditions)
W=0.05
#Converted to hourly by multiply by 24
N_field <- approxfun(x = c(161*24, 139*24, 105*24, 0, 28*24, 84*24, 172*24), y = N_sled$NitrateNitrite_uM, method = "linear", rule = 2) #N forcing function

###### Irradiance set-up Judith N 1 Yr 1
NOAA_Irradiance_Sledy1 <-  NOAA_Irradiance$PAR[2438:3774] # subset based on as_datetime("2017-11-1 12:00:00"), as_datetime("2018-04-17 12:00:00")
I_field <- approxfun(x = seq(from = 0, to = 4008, by = 3), y = NOAA_Irradiance_Sledy1, method = "linear", rule = 2) #irradiance forcing

######  Temp set-up Judith N 1 Yr 1
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

###### Model runs ######
# (the differential equation solver)

plan(multisession, workers=availableCores())

output_sled1_yr1 <- params_nested %>% mutate(std_L = future_map(data, function(df) {
  temp_params <- params_Lo
  temp_params[c("T_A", "T_H", "T_AH")] <- c(df$T_A, df$T_H, df$T_AH)
  ode_output<- as.data.frame(ode(y = state_Lo, t = times_Lo_Sled1, func = rates_Lo, parms = temp_params)) %>% select(time, W, L_allometric)
  ode_output
})) %>% select(-data)

output_sled1_yr1_clean <- output_sled1_yr1 %>%
  unnest(cols=std_L) %>% 
  ungroup() %>% 
  #group_by(type, level, res) %>% 
  group_by(level) %>% 
  mutate(Temp_C = T_Sled1_Y1-273.15, #conversion back to Celsius from Kelvin
         Date=sled1_date_seq,
         source="Point Judith Pond N 1")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setting up the forcing functions with field data for
#### Sled line 2 Year 1####
#inital biomass for conversions (cannot put in initial conditions)

# N forcing set-up Judith N 2 #
N_field <- approxfun(x = c(133*24, 111*24, 77*24, -28*24, 0, 56*24, 144*24), y = N_sled$NitrateNitrite_uM, method = "linear", rule = 2) #N forcing function

###### Irradiance set-up Judith N 2 Year 1
NOAA_Irradiance_Sledy1_L2 <-  NOAA_Irradiance$PAR[2662:3774] #subset by seq(as_datetime("2017-11-29 12:00:00"), as_datetime("2018-04-17 12:00:00"), by="hour")
I_field <- approxfun(x = seq(from = 0, to = 3336, by = 3), y = NOAA_Irradiance_Sledy1_L2, method = "linear", rule = 2) #irradiance forcing function

###### Temp set-Up Judith N 2 Year 1
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

###### Model runs ####

output_sled2_yr1 <- params_nested %>% mutate(std_L = future_map(data, function(df) {
  temp_params <- params_Lo
  temp_params[c("T_A", "T_H", "T_AH")] <- c(df$T_A, df$T_H, df$T_AH)
  ode_output<- as.data.frame(ode(y = state_Lo, t = times_Lo_Sled2, func = rates_Lo, parms = temp_params)) %>% select(time, W, L_allometric)
  ode_output
})) %>% select(-data)

output_sled2_yr1_clean <- output_sled2_yr1 %>%
  unnest(cols=std_L) %>% 
  ungroup() %>% 
 #group_by(type, level, res) %>%
  group_by(level) %>%
  mutate(Temp_C = T_Sled2_Y1-273.15, #conversion back to Celsius from Kelvin
         Date=sled2_date_seq,
         source="Point Judith Pond N 2")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Dredge line 1 Year 1####

# N forcing set-up Judith S 1 #
N_field <- approxfun(x = c(139*24, 161*24, 84*24, 0, 105*24), y = N_dredge$NitrateNitrite_uM, method = "linear", rule = 2)

###### Irradiance set-up Judith S 1 Year 1
NOAA_Irradiance_Dredgey1 <-  NOAA_Irradiance$PAR[2438:3814] #subset by 11/1/17 to 2018-04-22 12:00:00
I_field <- approxfun(x = seq(from = 0, to = 4128, by = 3), y = NOAA_Irradiance_Dredgey1, method = "linear", rule = 2) #irradiance forcing function

###### Temp set-Up Judith S 1 Year 1
Dredge_Y1_hobo <- Dredge_Y1_hobo_orig[3:16531,] #subset
DredgeT_hourly <- ceiling_date(Dredge_Y1_hobo$DateTime, unit = "hour") #set values to aggregate around
AvgTempKbyhr <- aggregate(Dredge_Y1_hobo$Temp_K, by=list(DredgeT_hourly), mean) #calculate average hourly temp
AvgTempKbyhr <- AvgTempKbyhr[4:4132, ] #subset
T_field <- approxfun(x = c(0:4128), y = AvgTempKbyhr$x, method = "linear", rule = 2) #the temp forcing function
T_Dredge1_Y1 <- T_field(0:4128) #for ease in later plotting of the forcing

##### Model runs ####
output_dredge1_yr1 <- params_nested %>% mutate(std_L = future_map(data, function(df) {
  temp_params <- params_Lo
  temp_params[c("T_A", "T_H", "T_AH")] <- c(df$T_A, df$T_H, df$T_AH)
  ode_output<- as.data.frame(ode(y = state_Lo, t = times_Lo_Dredge1, func = rates_Lo, parms = temp_params)) %>% select(time, W, L_allometric)
  ode_output
})) %>% select(-data)

output_dredge1_yr1_clean <- output_dredge1_yr1 %>%
  ungroup() %>% 
  unnest(cols=std_L) %>% 
  #group_by(type, level, res) %>%
  group_by(level) %>%
  mutate(Temp_C = T_Dredge1_Y1-273.15, #conversion back to Celsius from Kelvin
         Date=dredge1_date_seq,
         source="Point Judith Pond S 1")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setting up the forcing functions with field data for
#### Dredge line 2 Year 1####

# N forcing set-up Judith S 2
N_field <- approxfun(x = c(111*24, 133*24, 56*24, -28*24, 77*24), y = N_dredge$NitrateNitrite_uM, method = "linear", rule = 2) #N forcing function

###### Irradiance set-up Judith S 2 Year 1
NOAA_Irradiance_Dredgey1_L2 <-  NOAA_Irradiance$PAR[2662:3814] #subset by seq(as_datetime("2017-11-29 12:00:00"), as_datetime("2018-04-22 12:00:00"), by="hour")
I_field <- approxfun(x = seq(from = 0, to = 3456, by = 3), y = NOAA_Irradiance_Dredgey1_L2, method = "linear", rule = 2) #irradiance forcing function

###### Temp set-up Judith S 2 Year 1
Dredge_Y1_hobo <- Dredge_Y1_hobo_orig[3:16531,] #cut 2 points in beginning, logger not yet in water
DredgeT_hourly <- ceiling_date(Dredge_Y1_hobo$DateTime, unit = "hour") #determine dates to aggregate around
AvgTempKbyhr <- aggregate(Dredge_Y1_hobo$Temp_K, by=list(DredgeT_hourly), mean) #calculate average hourly temp
AvgTempKbyhr_sub <- AvgTempKbyhr[676:4132,] #subset
T_field <- approxfun(x = c(0:3456), y = c(AvgTempKbyhr_sub$x), method = "linear", rule = 2) #the temp forcing function
T_Dredge2_Y1 <- T_field(0:3456) #for ease of later plotting the temperature forcing

##### Model runs ####

output_dredge2_yr1 <- params_nested %>% mutate(std_L = future_map(data, function(df) {
  temp_params <- params_Lo
  temp_params[c("T_A", "T_H", "T_AH")] <- c(df$T_A, df$T_H, df$T_AH)
  ode_output<- as.data.frame(ode(y = state_Lo, t = times_Lo_Dredge2, func = rates_Lo, parms = temp_params)) %>% select(time, W, L_allometric)
  ode_output
})) %>% select(-data)

output_dredge2_yr1_clean <- output_dredge2_yr1 %>%
  unnest(cols=std_L) %>% 
  ungroup() %>% 
  #group_by(type, level, res) %>%
  group_by(level) %>%
  mutate(Temp_C = T_Dredge2_Y1-273.15, #conversion back to Celsius from Kelvin
         Date=dredge2_date_seq,
         source="Point Judith Pond S 2")

#combine all output into one dataframe
all_output_yr1 <- rbind(output_sled1_yr1_clean, output_sled2_yr1_clean, output_dredge1_yr1_clean, output_dredge2_yr1_clean) %>% 
  ungroup() #%>% 
#   mutate(params = case_when(
#     # res=="orig" ~ "orig",
#     # res=="lit" ~ "new",
#     # level=="high" & res=="means" ~ "high",
#     # level=="med" & res=="means" ~ "med",
#     # level=="low" & res=="means" ~ "low",
#     # level=="high" & res=="cross" ~ "high_cross",
#     # level=="med" & res=="cross" ~ "med_cross",
#     # level=="low" & res=="cross" ~ "low_cross",
#     # level=="high" & res=="all" ~ "high_rep",
#     # level=="med" & res=="all" ~ "med_rep",
#     # level=="low" & res=="all" ~ "low_rep"
#     level=="high"  ~ "high",
#     level=="low" ~ "low",
#     level=="orig" ~ "orig",
#     level=="lit" ~ "lit"
#   )
# ) 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Narragansett Bay Year 1 ####
### NB N set-up year 1 ####

#Wickford
WSA2_Y1 <- read.csv("WaterSampleAnalysis2Y1.csv", header = TRUE, fileEncoding="UTF-8-BOM") #Import water quality data
WSA2_Y1$Date <- mdy(WSA2_Y1$Date) #convert dates
names(WSA2_Y1)[1] <- "Site" #only necessary to run this code on some computers
Wickford_WSA <- filter(WSA2_Y1, Site == "Wickford") #filter by site

#Rome Pt
RomePt_WSA <- filter(WSA2_Y1, Site == "Rome Point") #filter by site

### NB DIC set-up year 1 ###########
Segarra2002Carbon <- read.csv("BrentonPoint_Segarra2002CarbonData.csv", header = TRUE, fileEncoding="UTF-8-BOM") #Import lit TCO2 data
names(Segarra2002Carbon)[1] <- "Date" #Only necessary for running the code on some computer
Segarra2002Carbon$Date <- mdy(Segarra2002Carbon$Date) #time field conversion
CO_2 <- mean(Segarra2002Carbon$TCO2_micromolPERkg)/1000000 #(mol CO2/L)

# Using same data frame for irradiance, just different subset

### NB temp set-up year 1 #############
# NB N (Wickford)
Wickford_Y1_hobo_orig <- read.csv("Wickford_Y1_hobo.csv", header = TRUE, fileEncoding="UTF-8-BOM") #Import Wickford Hobo data
Wickford_Y1_hobo_orig$DateTime <- mdy_hms(Wickford_Y1_hobo_orig$DateTime) #convert time field
Wickford_Y1_hobo_orig$Temp_K <- Wickford_Y1_hobo_orig$Temp_C+273.15 #create column with temp in K

# NB S (Rome)
RomePoint_Y1_hobotemp_orig <- read.csv("RomePoint_Y1_hobotemp.csv", header = TRUE, fileEncoding="UTF-8-BOM") #import
RomePoint_Y1_hobotemp_orig$DateTime <- mdy_hms(RomePoint_Y1_hobotemp_orig$DateTime) #convert time field
RomePoint_Y1_hobotemp_orig$Temp_K <- RomePoint_Y1_hobotemp_orig$Temp_C+273.15 #create column with temp in K


### Wickford Year 1 ####
N <- Wickford_WSA[c("Date","NitrateNitrite_uM")] #new data frame with the relevant columns
N$NitrateNitrite_uM <- N$NitrateNitrite_uM/1000000 #convert from micromoles/L to moles/L

N_field <- approxfun(x = c(115*24, 0, 81*24, 59*24, 38*24, 138*24), y = c(N$NitrateNitrite_uM[1], N$NitrateNitrite_uM[3:7]), method = "linear", rule = 2) #N forcing function

###### Irradiance forcing set-up #
NOAA_Irradiance_Wickfordy1 <-  NOAA_Irradiance$PAR[2702:3806] #subset by seq(as_datetime("2017-12-4 12:00:00"), as_datetime("2018-04-21 12:00:00"), by="hour")
I_field <- approxfun(x = seq(from = 0, to = 3312, by = 3), y = NOAA_Irradiance_Wickfordy1, method = "linear", rule = 2) #irradiance forcing function

###### Temp forcing set-up
Wickford_Y1_hobo <- Wickford_Y1_hobo_orig[607:13839,] #subset based on 2017-12-04 15:30:00 to 2018-04-21 11:30:00 
WickfordT_hourly <- ceiling_date(Wickford_Y1_hobo$DateTime, unit = "hour") #determine values to aggregate around
AvgTempKbyhr <- aggregate(Wickford_Y1_hobo$Temp_K, by=list(WickfordT_hourly), mean) #calculate average hourly temp
fd <- AvgTempKbyhr[1:4,] #a few replacement data points at the front of the forcing
T_field <- approxfun(x = c(0:3312), y = c(fd$x, AvgTempKbyhr$x), method = "linear", rule = 2) #the temp forcing function
T_W_Y1 <- T_field(0:3312) #For ease of plotting the temp forcing

###### Model runs ######

output_W_yr1 <- params_nested %>% mutate(std_L = future_map(data, function(df) {
  temp_params <- params_Lo
  temp_params[c("T_A", "T_H", "T_AH")] <- c(df$T_A, df$T_H, df$T_AH)
  ode_output<- as.data.frame(ode(y = state_Lo, t = times_Y1_W, func = rates_Lo, parms = temp_params)) %>% select(time, W, L_allometric)
  return(ode_output)
})) %>% select(-data)

output_W_yr1_clean <- output_W_yr1 %>%
  unnest(cols = std_L) %>% 
  ungroup() %>% 
  #group_by(type, level, res) %>% 
  group_by(level) %>%
  mutate(Temp_C = T_W_Y1-273.15, #conversion back to Celsius from Kelvin
         Date=W_date_seq_Y1,
         source="Narragansett Bay N")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Rome Pt line 1 Year 1####

###### N forcing set-up #
N <- RomePt_WSA[c("Date","NitrateNitrite_uM")] #new dataframe with relevant collumns
N$NitrateNitrite_uM <- N$NitrateNitrite_uM/1000000 #convert from micromoles/L to moles/Lmean
#multiplied by 24 to take from daily to hourly
N_field <- approxfun(x = c(114*24, 148*24, 0, 71*24, 92*24, 171*24), y = c(N$NitrateNitrite_uM[1:3], N$NitrateNitrite_uM[5:7]), method = "linear", rule = 2) #N forcing function

###### Irradiance forcing set-up #
NOAA_Irradiance_RomePty1 <-  NOAA_Irradiance$PAR[2438:3806] #subset by seq(as_datetime("2017-11-1 12:00:00"), as_datetime("2018-04-21 12:00:00"), by="hour")
I_field <- approxfun(x = seq(from = 0, to = 4104, by = 3), y = NOAA_Irradiance_RomePty1, method = "linear", rule = 2) #irradiance forcing function

RomePoint_Y1_hobotemp <- RomePoint_Y1_hobotemp_orig[6:16425,] #subset based on 2017-11-01 13:15:00 start and 2018-04-21 14:00:00 end
RomePointT_hourly <- ceiling_date(RomePoint_Y1_hobotemp$DateTime, unit = "hour") #determine dates to aggregate around
AvgTempKbyhr <- aggregate(RomePoint_Y1_hobotemp$Temp_K, by=list(RomePointT_hourly), mean) #calculate average hourly temp
AvgTempKbyhr <- AvgTempKbyhr[1:4103,] #subset
fd <- AvgTempKbyhr[1:2,] #two points of simulated data
T_field <- approxfun(x = c(0:4104), y = c(fd$x, AvgTempKbyhr$x), method = "linear", rule = 2) #the temp forcing function
T_R1_Y1 <- T_field(0:4104) #for ease in plotting the temperature forcing

#### Model runs ####
output_R1_yr1 <- params_nested %>% mutate(std_L = future_map(data, function(df) {
  temp_params <- params_Lo
  temp_params[c("T_A", "T_H", "T_AH")] <- c(df$T_A, df$T_H, df$T_AH)
  ode_output<- as.data.frame(ode(y = state_Lo, t = times_Y1_R1, func = rates_Lo, parms = temp_params)) %>% select(time, W, L_allometric)
  ode_output
})) %>% select(-data)

output_R1_yr1_clean <- output_R1_yr1 %>%
  unnest(cols = std_L) %>% 
  ungroup() %>% 
  #group_by(type, level, res) %>% 
  group_by(level) %>%
  mutate(Temp_C = T_R1_Y1-273.15, #conversion back to Celsius from Kelvin
         Date=R1_date_seq_Y1,
         source="Narragansett Bay S 1")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Rome Pt line 2 Year 1 ####
N <- RomePt_WSA[c("Date","NitrateNitrite_uM")] #new dataframe with relevant collumns
N$NitrateNitrite_uM <- N$NitrateNitrite_uM/1000000 #convert from micromoles/L to moles/L
N_field <- approxfun(x = c(79*24, 113*24, -25*24, 57*24, 136*24), y = c(N$NitrateNitrite_uM[1:2], N$NitrateNitrite_uM[5:7]), method = "linear", rule = 2) #N forcing function

NOAA_Irradiance_RomePty1_L2 <-  NOAA_Irradiance$PAR[2718:3806] #subset by seq(as_datetime("2017-12-6 12:00:00"), as_datetime("2018-04-21 12:00:00"), by="hour")
I_field <- approxfun(x = seq(from = 0, to = 3264, by = 3), y = NOAA_Irradiance_RomePty1_L2, method = "linear", rule = 2) #irradiance forcing function

RomePoint_Y1_hobotemp <- RomePoint_Y1_hobotemp_orig[3313:16425,] #subset based on 2017-12-06 00:00:00 - 2018-04-21 14:00:00
RomePointT_hourly <- ceiling_date(RomePoint_Y1_hobotemp$DateTime, unit = "hour") #determine dates to aggregate around
AvgTempKbyhr <- aggregate(RomePoint_Y1_hobotemp$Temp_K, by=list(RomePointT_hourly), mean) #calculate average hourly temp
AvgTempKbyhr <- AvgTempKbyhr[13:3277, ] #subset
T_field <- approxfun(x = c(0:3264), y = c(AvgTempKbyhr$x), method = "linear", rule = 2) #the temp forcing function
T_R2_Y1 <- T_field(0:3264) #for ease in later plotting

#### Model runs ####
output_R2_yr1 <- params_nested %>% mutate(std_L = future_map(data, function(df) {
  temp_params <- params_Lo
  temp_params[c("T_A", "T_H", "T_AH")] <- c(df$T_A, df$T_H, df$T_AH)
  ode_output<- as.data.frame(ode(y = state_Lo, t = times_Y1_R2, func = rates_Lo, parms = temp_params)) %>% select(time, W, L_allometric)
  ode_output
})) %>% select(-data)

output_R2_yr1_clean <- output_R2_yr1 %>%
  unnest(cols = std_L) %>% 
  ungroup() %>% 
  #group_by(type, level, res) %>%
  group_by(level) %>%
  mutate(Temp_C = T_R2_Y1-273.15, #conversion back to Celsius from Kelvin
         Date=R2_date_seq_Y1,
         source="Narragansett Bay S 2")

#combine all output into one dataframe
all_output_NB_yr1 <- rbind(output_W_yr1_clean, output_R1_yr1_clean, output_R2_yr1_clean) %>% ungroup() %>% mutate(
  params = case_when(
    # res=="orig" ~ "orig",
    # res=="lit" ~ "new",
    # level=="high" & res=="means" ~ "high",
    # level=="med" & res=="means" ~ "med",
    # level=="low" & res=="means" ~ "low",
    # level=="high" & res=="cross" ~ "high_cross",
    # level=="med" & res=="cross" ~ "med_cross",
    # level=="low" & res=="cross" ~ "low_cross",
    # level=="high" & res=="all" ~ "high_rep",
    # level=="med" & res=="all" ~ "med_rep",
    # level=="low" & res=="all" ~ "low_rep"
    level=="high"  ~ "high",
    level=="low" ~ "low",
    level=="orig" ~ "orig",
    level=="lit" ~ "lit"
  )
) 


### Import field data ####
KelpY1 <- read.csv("Year1kelpdata.csv", header = TRUE, fileEncoding="UTF-8-BOM")
names(KelpY1)[2] <- "Site"
KelpY1 <- filter(KelpY1, Site != "Fox Island")
KelpY1$Date <- mdy(KelpY1$SamplingDate)
KelpY1$SiteLine <- paste(KelpY1$Site, KelpY1$Line)
KelpY1 <- filter(KelpY1, SiteLine != "Narragansett Bay N 2")

field_data_fun <- function(df) {
  df_out <- df %>%
    group_by(Date) %>%
    summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length, na.rm = TRUE)) %>%
    mutate(Date= as.POSIXct(Date)) %>%
    na.omit()
  df_out
}


PJN1_meandat <- KelpY1[KelpY1$SiteLine == "Point Judith Pond N 1",] %>% field_data_fun()
PJN2_meandat <- KelpY1[KelpY1$SiteLine == "Point Judith Pond N 2",] %>% field_data_fun()
PJS1_meandat <- KelpY1[KelpY1$SiteLine == "Point Judith Pond S 1",] %>% field_data_fun()
PJS2_meandat <- KelpY1[KelpY1$SiteLine == "Point Judith Pond S 2",] %>% field_data_fun()
NBN1_Y1_meandat <- KelpY1[KelpY1$SiteLine == "Narragansett Bay N 1",] %>% field_data_fun()
NBS1_Y1_meandat <- KelpY1[KelpY1$SiteLine == "Narragansett Bay S 1",] %>% field_data_fun()
NBS2_Y1_meandat <- KelpY1[KelpY1$SiteLine == "Narragansett Bay S 2",] %>% field_data_fun()

### Combine with model data #####

field_data <- bind_rows(list("Point Judith Pond N 1" = PJN1_meandat, "Point Judith Pond N 2" = PJN2_meandat, "Point Judith Pond S 1" = PJS1_meandat, "Point Judith Pond S 2" =PJS2_meandat), .id="source")

field_data_NB_Y1 <- bind_rows(list("Narragansett Bay N" = NBN1_Y1_meandat, "Narragansett Bay S 1" = NBS1_Y1_meandat, "Narragansett Bay S 2" =NBS2_Y1_meandat), .id="source")



### Year 2 ################################################################################################################# 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Initial conditions year 2 ############
state_LoY2 <- c(m_EC = 0.01, #0.9 #mol C/molM_V  #Reserve density of C reserve (initial mass of C reserve per initial mass of structure)
                m_EN = 0.09, #mol N/molM_V #Reserve density of N reserve (initial mass of N reserve per initial mass of structure)
                M_V = 0.05/(w_V+0.09*w_EN+0.01*w_EC)) #molM_V #initial mass of structure

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####### Time steps year 2 #######
#(First number of time step, last number of time step, interval to step)
times_Y2_Sled1 <- seq(0, 3408, 1) #142 days stepped hourly
times_Y2_Sled2 <- seq(0, 2064, 1) #86 days stepped hourly
times_Y2_Dredge1 <- seq(0, 3408, 1) #142 days stepped hourly
times_Y2_Dredge2 <- seq(0, 2064, 1) #86 days stepped hourly
times_Y2_W <- seq(0, 3720, 1) #155 days stepped hourly
times_Y2_R1 <- seq(0, 3720, 1) #155 days stepped hourly
times_Y2_R2 <- seq(0, 2208, 1) #92 days stepped hourly

L1_date_seq_Y2 <- seq(as_datetime("2018-12-12 12:00:00"), as_datetime("2019-05-03 12:00:00"), by="hour") #for Sled 1 and Dredge 1
L2_date_seq_Y2 <- seq(as_datetime("2019-02-06 12:00:00"), as_datetime("2019-05-03 12:00:00"), by="hour") #for Sled 2 and Dredge 2
W_date_seq_Y2 <- seq(as_datetime("2018-12-19 12:00:00"), as_datetime("2019-05-23 12:00:00"), by="hour")
R1_date_seq_Y2 <- seq(as_datetime("2018-12-20 12:00:00"), as_datetime("2019-05-24 12:00:00"), by="hour")
R2_date_seq_Y2 <-seq(as_datetime("2019-2-21 12:00:00"), as_datetime("2019-05-24 12:00:00"), by="hour")

#### PJP N set-up year 2 ####
WSA_Y2 <- read.csv("WaterSamplesY2.csv", header = TRUE, fileEncoding="UTF-8-BOM") #Import water Q data
WSA_Y2$Date <- mdy(WSA_Y2$Date) #convert dates
names(WSA_Y2)[1] <- "Site" #only necessary for some computers running this code

Sled_WSA2 <- filter(WSA_Y2, Site == "Moonstone Sled") #filter by site
Dredge_WSA2 <- filter(WSA_Y2, Site == "Moonstone Dredge")

Sled_WSA2$NO3NO2_µM <- Sled_WSA2$NO3NO2_µM/1000000 #convert from micromoles/L to moles/L
Dredge_WSA2$NO3NO2_µM <- Dredge_WSA2$NO3NO2_µM/1000000

#### PJP DIC set-up year 2 ###########
DIC <- read.csv("Ninigret_EPA_DIC.csv", header = TRUE, fileEncoding="UTF-8-BOM") #Import Ninigret DIC data
CO_2 <- mean(DIC$DIC.uMkg.mean) #micromole DIC/kg (Jason said it was okay to assume that 1kg of seawater is 1L of seawater (actual conversion requires density calc from salinity and T))
#need units to match K_C (molDIC/L)
CO_2 <- CO_2/1000000

#Using same data frame for irradiance, just different subset

#### PJP temp set-up year 2 #############
# Point Judith Pond N (sled)
Sled_Y2_Hobo_orig <- read.csv("Sled_Y2_HoboLightTemp.csv", header = TRUE, fileEncoding="UTF-8-BOM") #import
Sled_Y2_Hobo_orig$DateTime <- mdy_hms(Sled_Y2_Hobo_orig$DateTime) #convert date time field
Sled_Y2_Hobo_orig$Temp_K <- Sled_Y2_Hobo_orig$Temp_C+273.15 #create column with temp in K
SledY2T_hourly <- ceiling_date(Sled_Y2_Hobo_orig$DateTime, unit = "hour") #determine times to aggregate around

# Point Judith Pond S (dredge)
Dredge_Y2_Hobo_orig <- read.csv("Dredge_Y2_HoboTempLight.csv", header = TRUE, fileEncoding="UTF-8-BOM") #import
Dredge_Y2_Hobo_orig$DateTime <- mdy_hms(Dredge_Y2_Hobo_orig$DateTime) #convert date time field
Dredge_Y2_Hobo_orig$Temp_K <- Dredge_Y2_Hobo_orig$Temp_C+273.15 #create column with temp in K
DredgeY2T_hourly <- ceiling_date(Dredge_Y2_Hobo_orig$DateTime, unit = "hour") #determine what times to aggregate around
AvgTempKbyhr <- aggregate(Dredge_Y2_Hobo_orig$Temp_K, by=list(DredgeY2T_hourly), mean) #calculate average hourly temp



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Point Judith Year 2 ####

### Sled line 1 Year 2 ####
N_field <- approxfun(x = c(1*24, 57*24, 93*24, 124*24, 142*24, 163*24), y = c(Sled_WSA2$NO3NO2_µM), method = "linear", rule = 2) #N forcing function

###### Irradiance forcing set-up #
NOAA_Irradiance_Sledy2 <-  NOAA_Irradiance$PAR[5686:6822] #subset by seq(as_datetime("2018-12-12 12:00:00"), as_datetime("2019-05-03 12:00:00"), by="hour")
I_field <- approxfun(x = seq(from = 0, to = 3408, by = 3), y = NOAA_Irradiance_Sledy2, method = "linear", rule = 2) #irradiance forcing function

###### Temp forcing set-up #
AvgTempKbyhr <- aggregate(Sled_Y2_Hobo_orig$Temp_K, by=list(SledY2T_hourly), mean) #calculate average hourly temp
AvgTempKbyhr <- AvgTempKbyhr[2:3385,] #subset
fd <- rep(285, 25) #small section of replacement
T_field <- approxfun(x = c(0:3408), y = c(AvgTempKbyhr$x, fd), method = "linear", rule = 2) #the temp forcing function
T_Sled1_Y2 <- T_field(0:3408) #for ease in plotting the temperature forcing

###### Model runs ######
output_sled1_yr2 <- params_nested %>% mutate(std_L = future_map(data, function(df) {
  temp_params <- params_Lo
  temp_params[c("T_A", "T_H", "T_AH")] <- c(df$T_A, df$T_H, df$T_AH)
  ode_output<- as.data.frame(ode(y = state_LoY2, t = times_Y2_Sled1, func = rates_Lo, parms = temp_params)) %>% select(time, W, L_allometric)
  return(ode_output)
})) %>% select(-data)

output_sled1_yr2_clean <- output_sled1_yr2 %>%
  unnest(cols = std_L) %>% 
  ungroup() %>% 
  #group_by(type, level, res) %>% 
  group_by(level) %>%
  mutate(Temp_C = T_Sled1_Y2-273.15, #conversion back to Celsius from Kelvin
         Date=L1_date_seq_Y2,
         source="Point Judith Pond N 1")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Sled line 2 Year 2 ####

# N forcing set-up Judith N 2 #
Sled_WSA2_sub <- Sled_WSA2[2:6,] #remove the point before the relevant time range
N_field <- approxfun(x = c(1*24, 37*24, 68*24, 86*24, 107*24), y = c(Sled_WSA2_sub$NO3NO2_µM), method = "linear", rule = 2) #N forcing function

###### Irradiance forcing set-up #
NOAA_Irradiance_Sledy2_L2 <-  NOAA_Irradiance$PAR[6134:6822] #subset by seq(as_datetime("2019-02-06 12:00:00"), as_datetime("2019-05-03 12:00:00"), by="hour")
I_field <- approxfun(x = seq(from = 0, to = 2064, by = 3), y = NOAA_Irradiance_Sledy2_L2, method = "linear", rule = 2) #irradiance forcing function

###### Temp forcing set-Up #
AvgTempKbyhr <- aggregate(Sled_Y2_Hobo_orig$Temp_K, by=list(SledY2T_hourly), mean) #calculate average hourly temp
#calculate average hourly temp
AvgTempKbyhr <- AvgTempKbyhr[1346:3385,] #subset
fd <- rep(285, 25) #small data replacement
T_field <- approxfun(x = c(0:2064), y = c(AvgTempKbyhr$x, fd), method = "linear", rule = 2) #the temp forcing function
T_Sled2_Y2 <- T_field(0:2064) #for ease in later plotting

##### Model runs ####
output_sled2_yr2 <- params_nested %>% mutate(std_L = future_map(data, function(df) {
  temp_params <- params_Lo
  temp_params[c("T_A", "T_H", "T_AH")] <- c(df$T_A, df$T_H, df$T_AH)
  ode_output<- as.data.frame(ode(y = state_LoY2, t = times_Y2_Sled2, func = rates_Lo, parms = temp_params)) %>% select(time, W, L_allometric)
  ode_output
})) %>% select(-data)

output_sled2_yr2_clean <- output_sled2_yr2 %>%
  unnest(cols = std_L) %>% 
  ungroup() %>% 
 # group_by(type, level, res) %>% 
  group_by(level) %>%
  mutate(Temp_C = T_Sled2_Y2-273.15, #conversion back to Celsius from Kelvin
         Date=L2_date_seq_Y2,
         source="Point Judith Pond N 2")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Dredge line 1 Year 2 ####

N_field <- approxfun(x = c(1*24, 93*24, 124*24, 142*24, 163*24), y = c(Dredge_WSA2$NO3NO2_µM), method = "linear", rule = 2) #N forcing function

###### Irradiance forcing set-up #
NOAA_Irradiance_Dredgey2 <-  NOAA_Irradiance$PAR[5686:6822] #subset by seq(as_datetime("2018-12-12 12:00:00"), as_datetime("2019-05-03 12:00:00"), by="hour")
I_field <- approxfun(x = seq(from = 0, to = 3408, by = 3), y = NOAA_Irradiance_Dredgey2, method = "linear", rule = 2) #irradiance forcing function

AvgTempKbyhr <- aggregate(Dredge_Y2_Hobo_orig$Temp_K, by=list(DredgeY2T_hourly), mean) #calculate average hourly temp
AvgTempKbyhr <- AvgTempKbyhr[2:3384,] #subset
fd <- rep(285, 26) #small amount of replacement data
T_field <- approxfun(x = c(0:3408), y = c(AvgTempKbyhr$x, fd), method = "linear", rule = 2) #the temp forcing function
T_Dredge1_Y2 <- T_field(0:3408) #for ease of plotting

#### Model runs ####

output_dredge1_yr2 <- params_nested %>% mutate(std_L = future_map(data, function(df) {
  temp_params <- params_Lo
  temp_params[c("T_A", "T_H", "T_AH")] <- c(df$T_A, df$T_H, df$T_AH)
  ode_output<- as.data.frame(ode(y = state_LoY2, t = times_Y2_Dredge1, func = rates_Lo, parms = temp_params)) %>% select(time, W, L_allometric)
  ode_output
})) %>% select(-data)

output_dredge1_yr2_clean <- output_dredge1_yr2 %>%
  unnest(cols = std_L) %>% 
  ungroup() %>% 
  #group_by(type, level, res) %>%
  group_by(level) %>%
  mutate(Temp_C = T_Dredge1_Y2-273.15, #conversion back to Celsius from Kelvin
         Date=L1_date_seq_Y2,
         source="Point Judith Pond S 1")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Dredge line 2 Year 2 ####

N_field <- approxfun(x = c(-25*24, 37*24, 68*24, 86*24, 107*24), y = c(Dredge_WSA2$NO3NO2_µM), method = "linear", rule = 2) #N forcing function

NOAA_Irradiance_Dredgey2_L2 <-  NOAA_Irradiance$PAR[6134:6822] #subset by seq(as_datetime("2019-02-06 12:00:00"), as_datetime("2019-05-03 12:00:00"), by="hour")
I_field <- approxfun(x = seq(from = 0, to = 2064, by = 3), y = NOAA_Irradiance_Dredgey2_L2, method = "linear", rule = 2) #irradiance forcing function

AvgTempKbyhr <- aggregate(Dredge_Y2_Hobo_orig$Temp_K, by=list(DredgeY2T_hourly), mean) #calculate average hourly temp
AvgTempKbyhr <- AvgTempKbyhr[1346:3385,] #subset
fd <- rep(285, 25) #estimation to fill in gap in data
T_field <- approxfun(x = c(0:2064), y = c(AvgTempKbyhr$x, fd), method = "linear", rule = 2) #the temp forcing function
T_Dredge2_Y2 <- T_field(0:2064) #for ease in plotting

#### Model runs ####
output_dredge2_yr2 <- params_nested %>% mutate(std_L = future_map(data, function(df) {
  temp_params <- params_Lo
  temp_params[c("T_A", "T_H", "T_AH")] <- c(df$T_A, df$T_H, df$T_AH)
  ode_output<- as.data.frame(ode(y = state_LoY2, t = times_Y2_Dredge2, func = rates_Lo, parms = temp_params)) %>% select(time, W, L_allometric)
  ode_output
})) %>% select(-data)

output_dredge2_yr2_clean <- output_dredge2_yr2 %>%
  unnest(cols = std_L) %>% 
  ungroup() %>% 
  #group_by(type, level, res) %>% 
  group_by(level) %>%
  mutate(Temp_C = T_Dredge2_Y2-273.15, #conversion back to Celsius from Kelvin
         Date=L2_date_seq_Y2,
         source="Point Judith Pond S 2")

#combine all output into one data frame
all_output_yr2 <- rbind(output_sled1_yr2_clean, output_sled2_yr2_clean, output_dredge1_yr2_clean, output_dredge2_yr2_clean) %>% 
  ungroup() #%>% 
#   mutate(
#   params = case_when(
#     # res=="orig" ~ "orig",
#     # res=="lit" ~ "new",
#     # level=="high" & res=="means" ~ "high",
#     # level=="med" & res=="means" ~ "med",
#     # level=="low" & res=="means" ~ "low",
#     # level=="high" & res=="cross" ~ "high_cross",
#     # level=="med" & res=="cross" ~ "med_cross",
#     # level=="low" & res=="cross" ~ "low_cross",
#     # level=="high" & res=="all" ~ "high_rep",
#     # level=="med" & res=="all" ~ "med_rep",
#     # level=="low" & res=="all" ~ "low_rep"
#     level=="high"  ~ "high",
#     level=="low" ~ "low",
#     level=="orig" ~ "orig",
#     level=="lit" ~ "lit"
# ) )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Narragansett Bay Year 2 ####
##### NB set-up year 2 ####

#Wickford
WSA_Y2 <- read.csv("WaterSamplesY2.csv", header = TRUE, fileEncoding="UTF-8-BOM") #Import water Q data
WSA_Y2$Date <- mdy(WSA_Y2$Date) #convert dates
names(WSA_Y2)[1] <- "Site" #only necessary for some computers running this code

Wickford_WSA2 <- filter(WSA_Y2, Site == "Wickford") #filter by site
Wickford_WSA2$NO3NO2_µM <- Wickford_WSA2$NO3NO2_µM/1000000 #convert from micromoles/L to moles/L

##### NB DIC set-up year 2 ###########
Segarra2002Carbon <- read.csv("BrentonPoint_Segarra2002CarbonData.csv", header = TRUE, fileEncoding="UTF-8-BOM") #Import lit TCO2 data
names(Segarra2002Carbon)[1] <- "Date" #Only necessary for running the code on some computer
Segarra2002Carbon$Date <- mdy(Segarra2002Carbon$Date) #time field conversion
CO_2 <- mean(Segarra2002Carbon$TCO2_micromolPERkg)/1000000 #(mol CO2/L)

# Using same data frame for irradiance, just different subset

#### NB temp set-up year 2 #############
# NB N (Wickford)
Wickford_Y2_Hobo <- read.csv("Wickford_Y2_HoboLightTemp.csv", header = TRUE, fileEncoding="UTF-8-BOM") #import
Wickford_Y2_Hobo$DateTime <- mdy_hms(Wickford_Y2_Hobo$DateTime) #convert date time field
Wickford_Y2_Hobo$Temp_K <- Wickford_Y2_Hobo$Temp_C+273.15 #create column with temp in K
WickfordY2T_hourly <- ceiling_date(Wickford_Y2_Hobo$DateTime, unit = "hour") #determine the times to aggregate around

# NB S (Rome)
RomePt_Y2_Hobo_orig <- read.csv("RomePt_Y2_HoboTempLight.csv", header = TRUE, fileEncoding="UTF-8-BOM") #import
RomePt_Y2_Hobo_orig$DateTime <- mdy_hm(RomePt_Y2_Hobo_orig$DateTime) #convert date time field
RomePt_Y2_Hobo_orig$Temp_K <- RomePt_Y2_Hobo_orig$Temp_C+273.15 #create collumn with temp in K
RomePtY2T_hourly <- ceiling_date(RomePt_Y2_Hobo_orig$DateTime, unit = "hour") #determine the times to aggregate around


### Wickford Year 2####
N_field <- approxfun(x = c(1*24, 55*24, 85*24, 156*24), y = c(Wickford_WSA2$NO3NO2_µM), method = "linear", rule = 2) #N forcing function

###### Irradiance forcing set-up #
NOAA_Irradiance_Wickfordy2 <-  NOAA_Irradiance$PAR[5742:6982] #subset by seq(as_datetime("2018-12-19 12:00:00"), as_datetime("2019-05-23 12:00:00"), by="hour")
I_field <- approxfun(x = seq(from = 0, to = 3720, by = 3), y = NOAA_Irradiance_Wickfordy2, method = "linear", rule = 2) #irradiance forcing function

###### Temp forcing set-up
AvgTempKbyhr <- aggregate(Wickford_Y2_Hobo$Temp_K, by=list(WickfordY2T_hourly), mean) #calculate average hourly temp
fd <- rep(278, 4) #replacement data
AvgTempKbyhr_sub <- AvgTempKbyhr[4:3716,] #subset
fd2 <- rep(287, 4) #second small section of replacement data
T_field <- approxfun(x = c(0:3720), y = c(fd, AvgTempKbyhr_sub$x, fd2), method = "linear", rule = 2) #the temp forcing function
T_W_Y2 <- T_field(0:3720) #for ease of plotting

###### Model runs ######

output_W_yr2 <- params_nested %>% mutate(std_L = future_map(data, function(df) {
  temp_params <- params_Lo
  temp_params[c("T_A", "T_H", "T_AH")] <- c(df$T_A, df$T_H, df$T_AH)
  ode_output<- as.data.frame(ode(y = state_LoY2, t = times_Y2_W, func = rates_Lo, parms = temp_params)) %>% select(time, W, L_allometric)
  return(ode_output)
})) %>% select(-data)

output_W_yr2_clean <- output_W_yr2 %>%
  unnest(cols = std_L) %>% 
  ungroup() %>% 
  #group_by(type, level, res) %>% 
  group_by(level) %>%
  mutate(Temp_C = T_W_Y2-273.15, #conversion back to Celsius from Kelvin
         Date=W_date_seq_Y2,
         source="Narragansett Bay N")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Rome Pt line 1 Year 2####

###### N forcing set-up #
GSO_N1 <- read.csv("T98BayNitrate.csv", header = TRUE, fileEncoding="UTF-8-BOM") #Import water quality data
GSO_N1$Date <- mdy(GSO_N1$Date) #convert dates
GSO_N1 <- GSO_N1[103:124,]
GSO_N1$NO3NO2 <- GSO_N1$NO3NO2/1000000 #convert from micromoles/L to moles/L
N_field <- approxfun(x = c(8*24, 14*24, 20*24, 29*24, 35*24, 41*24, 48*24, 56*24, 62*24, 69*24, 79*24, 83*24, 93*24, 99*24, 107*24, 111*24, 118*24, 125*24, 132*24, 139*24, 146*24, 153*24), y = c(GSO_N1$NO3NO2), method = "linear", rule = 2) #N forcing function

###### Irradiance forcing set-up #
NOAA_Irradiance_RomePty2 <-  NOAA_Irradiance$PAR[5750:6990] #subset by seq(as_datetime("2018-12-20 12:00:00"), as_datetime("2019-05-24 12:00:00"), by="hour")
I_field <- approxfun(x = seq(from = 0, to = 3720, by = 3), y = NOAA_Irradiance_RomePty2, method = "linear", rule = 2) #irradiance forcing function

AvgTempKbyhr <- aggregate(RomePt_Y2_Hobo_orig$Temp_K, by=list(RomePtY2T_hourly), mean) #calculate average hourly temp
fd <- rep(280, 4) #small bit of replacement data
AvgTempKbyhr_sub <- AvgTempKbyhr[28:2414,] #subset
#Using Wickford temp to fill in gap in the Rome Pt temp
AvgTempKbyhr_Wickford <- aggregate(Wickford_Y2_Hobo$Temp_K, by=list(WickfordY2T_hourly), mean) #calculate average hourly temp
AvgTempKbyhr_W <- AvgTempKbyhr_Wickford[2415:3716,]
fd2 <- rep(287, 28)
T_field <- approxfun(x = c(0:3720), y = c(fd, AvgTempKbyhr_sub$x, AvgTempKbyhr_W$x, fd2), method = "linear", rule = 2) #the temp forcing function
T_R1_Y2 <- T_field(0:3720) #for ease of plotting

#### Model runs ####
output_R1_yr2 <- params_nested %>% mutate(std_L = future_map(data, function(df) {
  temp_params <- params_Lo
  temp_params[c("T_A", "T_H", "T_AH")] <- c(df$T_A, df$T_H, df$T_AH)
  ode_output<- as.data.frame(ode(y = state_LoY2, t = times_Y2_R1, func = rates_Lo, parms = temp_params)) %>% select(time, W, L_allometric)
  ode_output
})) %>% select(-data)

output_R1_yr2_clean <- output_R1_yr2 %>%
  unnest(cols = std_L) %>% 
  ungroup() %>% 
  #group_by(type, level, res) %>% 
  group_by(level) %>%
  mutate(Temp_C = T_R1_Y2-273.15, #conversion back to Celsius from Kelvin
         Date=R1_date_seq_Y2,
         source="Narragansett Bay S 1")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Rome Pt line 2 Year 2 ####
GSO_N1 <- read.csv("T98BayNitrate.csv", header = TRUE, fileEncoding="UTF-8-BOM") #Import water quality data
GSO_N1$Date <- mdy(GSO_N1$Date) #convert dates
GSO_N1 <- GSO_N1[112:124,]
GSO_N1$NO3NO2 <- GSO_N1$NO3NO2/1000000 #convert from micromoles/L to moles/L
N_field <- approxfun(x = c(6*24, 16*24, 20*24, 30*24, 36*24, 44*24, 48*24, 55*24, 62*24, 69*24, 76*24, 83*24, 90*24), y = c(GSO_N1$NO3NO2), method = "linear", rule = 2) #N forcing function

NOAA_Irradiance_RomePty2_L2 <-  NOAA_Irradiance$PAR[6254:6990] #subset by seq(as_datetime("2019-2-21 12:00:00"), as_datetime("2019-05-24 12:00:00"), by="hour")
I_field <- approxfun(x = seq(from = 0, to = 2208, by = 3), y = NOAA_Irradiance_RomePty2_L2, method = "linear", rule = 2) #irradiance forcing function

AvgTempKbyhr <- aggregate(RomePt_Y2_Hobo_orig$Temp_K, by=list(RomePtY2T_hourly), mean) #calculate average hourly temp
AvgTempKbyhr_sub <- AvgTempKbyhr[1536:2414,] #subset
#Using Wickford temp to fill in th gap in the Rome Pt temp
AvgTempKbyhr_W <- AvgTempKbyhr_Wickford[2415:3716,]
fd2 <- rep(287, 28)
T_field <- approxfun(x = c(0:2208), y = c(AvgTempKbyhr_sub$x, AvgTempKbyhr_W$x, fd2), method = "linear", rule = 2) #the temp forcing function
T_R2_Y2 <- T_field(0:2208) #for ease in plotting

#### Model runs ####
output_R2_yr2 <- params_nested %>% mutate(std_L = future_map(data, function(df) {
  temp_params <- params_Lo
  temp_params[c("T_A", "T_H", "T_AH")] <- c(df$T_A, df$T_H, df$T_AH)
  ode_output<- as.data.frame(ode(y = state_LoY2, t = times_Y2_R2, func = rates_Lo, parms = temp_params)) %>% select(time, W, L_allometric)
  ode_output
})) %>% select(-data)

output_R2_yr2_clean <- output_R2_yr2 %>%
  unnest(cols = std_L) %>% 
  ungroup() %>% 
  #group_by(type, level, res) %>% 
  group_by(level) %>%
  mutate(Temp_C = T_R2_Y2-273.15, #conversion back to Celsius from Kelvin
         Date=R2_date_seq_Y2,
         source="Narragansett Bay S 2")

#combine all output into one data frame
all_output_NB_yr2 <- rbind(output_W_yr2_clean, output_R1_yr2_clean, output_R2_yr2_clean) %>% ungroup() %>% mutate(
  params = case_when(
    # res=="orig" ~ "orig",
    # res=="lit" ~ "new",
    # level=="high" & res=="means" ~ "high",
    # level=="med" & res=="means" ~ "med",
    # level=="low" & res=="means" ~ "low",
    # level=="high" & res=="cross" ~ "high_cross",
    # level=="med" & res=="cross" ~ "med_cross",
    # level=="low" & res=="cross" ~ "low_cross",
    # level=="high" & res=="all" ~ "high_rep",
    # level=="med" & res=="all" ~ "med_rep",
    # level=="low" & res=="all" ~ "low_rep"
    level=="high"  ~ "high",
    level=="low" ~ "low",
    level=="orig" ~ "orig",
    level=="lit" ~ "lit"
  ))

### Import field data ####
KelpY2 <- read.csv("Year2kelpdata.csv", header = TRUE, fileEncoding="UTF-8-BOM")
names(KelpY2)[2] <- "Site"
KelpY2 <- filter(KelpY2, Site != "Fox Island")
KelpY2$Date <- mdy(KelpY2$SamplingDate)
KelpY2$SiteLine <- paste(KelpY2$Site, KelpY2$Line)
KelpY2 <- filter(KelpY2, SiteLine != "Narragansett Bay N 2")


PJN1_Y2_meandat <- KelpY2[KelpY2$SiteLine == "Point Judith Pond N 1",] %>% field_data_fun()
PJN2_Y2_meandat <- KelpY2[KelpY2$SiteLine == "Point Judith Pond N 2",] %>% field_data_fun()
PJS1_Y2_meandat <- KelpY2[KelpY2$SiteLine == "Point Judith Pond S 1",] %>% field_data_fun()
PJS2_Y2_meandat <- KelpY2[KelpY2$SiteLine == "Point Judith Pond S 2",] %>% field_data_fun()
NBN1_Y2_meandat <- KelpY2[KelpY2$SiteLine == "Narragansett Bay N 1",] %>% field_data_fun()
NBS1_Y2_meandat <- KelpY2[KelpY2$SiteLine == "Narragansett Bay S 1",] %>% field_data_fun()
NBS2_Y2_meandat <- KelpY2[KelpY2$SiteLine == "Narragansett Bay S 2",] %>% field_data_fun()

### Combine with model data ####


field_data_Y2 <- bind_rows(list("Point Judith Pond N 1" = PJN1_Y2_meandat, "Point Judith Pond N 2" = PJN2_Y2_meandat, "Point Judith Pond S 1" = PJS1_Y2_meandat, "Point Judith Pond S 2" =PJS2_Y2_meandat), .id="source")

field_data_NB_Y2 <- bind_rows(list("Narragansett Bay N" = NBN1_Y2_meandat, "Narragansett Bay S 1" = NBS1_Y2_meandat, "Narragansett Bay S 2" =NBS2_Y2_meandat), .id="source")

all_field_data <- bind_rows(list("1"=field_data, "1"=field_data_NB_Y1, "2"=field_data_Y2, "2"=field_data_NB_Y2), .id="year")

### rmse ####

rmse_dat_new <- field_data %>% 
  group_by(source) %>% 
  left_join((all_output_yr1 %>% group_by(source))) %>%
  group_by(source, level) %>% 
  summarise(rmse = rmse(mean_length, L_allometric))

rmse_NB_Y1 <- field_data_NB_Y1 %>% group_by(source) %>% left_join((all_output_NB_yr1 %>% group_by(source))) %>% 
  group_by(source, level) %>% 
  summarise(rmse = rmse(mean_length, L_allometric))

rmse_dat_Y2 <- field_data_Y2 %>% group_by(source) %>% 
  left_join((all_output_yr2 %>% group_by(source))) %>%
  group_by(source, level) %>% 
  summarise(rmse = rmse(mean_length, L_allometric))

rmse_NB_Y2 <- field_data_NB_Y2 %>% group_by(source) %>% 
  left_join((all_output_NB_yr2 %>% group_by(source))) %>% 
  group_by(source, level) %>% 
  summarise(rmse = rmse(mean_length, L_allometric))

all_rmse <- bind_rows(rmse_dat_new %>% mutate(year=1), 
                      rmse_dat_Y2 %>% mutate(year=2),
                      rmse_NB_Y1 %>% mutate(year=1),
                      rmse_NB_Y2%>% mutate(year=2))

orig_rmse <- all_rmse %>% filter(level =="orig") %>% ungroup()

all_rmse <- all_rmse %>% ungroup() %>% 
  left_join(orig_rmse %>% 
              mutate(orig_rmse = rmse) %>% 
              select(source, year, orig_rmse), by=c("source", "year")) %>%
  mutate(improvement = orig_rmse-rmse) %>% 
  filter(level!="lit") %>% 
  mutate(level=fct_drop(level))

perc_imp <- all_rmse %>% mutate(imp = if_else(improvement>0, TRUE, FALSE)) %>% 
 group_by(level) %>% 
  summarize(num_imp = sum(imp), 
            perc_imp = sum(imp)/length(imp), 
            mean_imp = mean(if_else(improvement>0, improvement, 0)))

pjp_plot1 <- ggplot(all_output_yr1)+
  theme_classic()+
  geom_smooth(aes(x=Date, y=L_allometric, color=level))+
  geom_point(data=field_data, aes(x=Date, y=mean_length, size="obs"))+
  facet_grid(~source)+
  labs(y="Kelp length (cm)", x="", color=NULL, size=NULL)+
  scale_size_manual(values=c("obs"=2), breaks=c("obs"), labels=c("obs"="Observations"))+
  scale_color_manual(values=c("low"='#0f85a0',"high"="#dd4124","orig"="black"),
                     breaks=c("low","high", "orig"),
                     labels=c("high"="Warm", "low"="Cold" ,"orig"="Original"))+
  scale_x_datetime(limits=as.POSIXct(c("2017-10-30 23:00:00", "2018-04-24 23:00:00")),breaks="2 months", date_labels="%b")+
  theme(text = element_text(size=18),
        title = element_text(size=15),
        axis.title.y = element_text(margin = margin(t = 0, r = 9, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 9, r = 0, b = 0, l = 0)),
        legend.spacing = unit(0, "pt"))

pjp_plot2 <-ggplot(all_output_yr2)+
  theme_classic()+
  geom_smooth(aes(x=Date, y=L_allometric, color=level))+
  geom_point(data=field_data_Y2, aes(x=Date, y=mean_length, size="obs"))+
  labs(y="Kelp length (cm)", x="", color=NULL, size=NULL)+
  facet_grid(~source, scales = "free")+
  scale_size_manual(values=c("obs"=2), breaks=c("obs"), labels=c("obs"="Observations"))+
  scale_color_manual(values=c("low"='#0f85a0',"high"="#dd4124","orig"="black"),
                     breaks=c("low","high", "orig"),
                     labels=c("high"="Warm", "low"="Cold" ,"orig"="Original"))+
  scale_x_datetime(limits=as.POSIXct(c("2018-11-29 12:00:00", "2019-05-30 12:00:00")),breaks="2 months", date_labels="%b")+
  ylim(0,100)

(pjp_plot1/pjp_plot2) +
  plot_layout(guides = 'collect') &
  theme(text = element_text(size=18),
        title = element_text(size=15),
        axis.title.y = element_text(margin = margin(t = 0, r = 9, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 9, r = 0, b = 0, l = 0)))


nb_plot1 <-ggplot(all_output_NB_yr1)+
  theme_classic()+
  geom_smooth(aes(x=Date, y=L_allometric, color=level))+
  geom_point(data=field_data_NB_Y1, aes(x=Date, y=mean_length, size="obs"))+
  facet_wrap(~source)+
  labs(y="Kelp length (cm)", x="", color=NULL, size=NULL)+
  scale_size_manual(values=c("obs"=2), breaks=c("obs"), labels=c("obs"="Observations"))+
  scale_color_manual(values=c("low"='#0f85a0',"high"="#dd4124","orig"="black"),
                     breaks=c("low","high", "orig"),
                     labels=c("high"="Warm", "low"="Cold" ,"orig"="Original"))+
 scale_x_datetime(limits=as.POSIXct(c("2017-10-30 23:00:00", "2018-04-24 23:00:00")),breaks="2 months", date_labels="%b")+
  ylim(0,180)

nb_plot2 <- ggplot(all_output_NB_yr2)+
  theme_classic()+
  geom_smooth(aes(x=Date, y=L_allometric, color=level))+
  geom_point(data=field_data_NB_Y2, aes(x=Date, y=mean_length, size="obs"))+
  facet_wrap(~source)+
  labs(y="Kelp length (cm)", x=NULL, color=NULL, size=NULL)+
  scale_size_manual(values=c("obs"=2), breaks=c("obs"), labels=c("obs"="Observations"))+
  scale_color_manual(values=c("low"='#0f85a0',"high"="#dd4124","orig"="black"),
                     breaks=c("low","high", "orig"),
                     labels=c("high"="Warm", "low"="Cold" ,"orig"="Original"))+
  scale_x_datetime(limits=as.POSIXct(c("2018-11-29 12:00:00", "2019-05-30 12:00:00")),breaks="2 months", date_labels="%b")+
  ylim(0,70)

(nb_plot1/nb_plot2) +
  plot_layout(guides = 'collect') &
  theme(text = element_text(size=18),
        title = element_text(size=15),
        axis.title.y = element_text(margin = margin(t = 0, r = 9, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 9, r = 0, b = 0, l = 0)),
        legend.spacing = unit(0, "pt"))

(nb_plot1+pjp_plot1)/(nb_plot2+pjp_plot2) +
  plot_layout(guides = 'collect')

all_output <- bind_rows(all_output_yr1 %>% mutate(year=1), all_output_yr2 %>% mutate(year=2), all_output_NB_yr1 %>% mutate(year=1), all_output_NB_yr2 %>% mutate(year=2))

all_output<- all_output %>% 
  mutate(date2017=paste(sep="", month(Date), "-", day(Date), " ", hour(Date), ":00:00")) %>% mutate(date2017 = case_when(year(Date)==2017 ~ Date,
    year(Date)==2018 & month(Date) <7 ~ Date,
    year(Date)==2018 & month(Date) >=7 ~ as_datetime(paste(sep="", "2017-", as.character(date2017))), 
    year(Date)==2019 ~ as_datetime(paste(sep="", "2018-", as.character(date2017)))))

all_field_data<- all_field_data %>% 
  mutate(date2017=paste(sep="", month(Date), "-", day(Date), " ", hour(Date), ":00:00")) %>% mutate(date2017 = case_when(year(Date)==2017 ~ Date,
                                                                                                                         year(Date)==2018 & month(Date) <7 ~ Date,
                                                                                                                         year(Date)==2018 & month(Date) >=7 ~ as_datetime(paste(sep="", "2017-", as.character(date2017))), 
                                                                                                                         year(Date)==2019 ~ as_datetime(paste(sep="", "2018-", as.character(date2017)))))

ggplot(data=all_output %>% filter(year==1), aes(x=Date, y=L_allometric, color=level)) +
  geom_smooth()+
  facet_wrap(~source, scales="free")


ggplot() +
  theme_bw()+
  geom_point(data=all_field_data %>% filter(str_detect(source, "Bay")), aes(x=date2017, y=mean_length,size="obs"))+
  geom_smooth(data=all_output %>% filter(str_detect(source, "Bay"), level!="lit"), aes(x=date2017, y=L_allometric, color=level))+
  facet_grid(year~source, labeller=labeller(year = c("1"="Year 1", "2"="Year 2")))+
  
  scale_size_manual(values=c("obs"=2), breaks=c("obs"), labels=c("obs"="Observations"))+
  scale_color_manual(values=c("cold"='#0f85a0',"warm"="#dd4124", "orig"="black"),
                     breaks=c("cold","warm","orig"),
                     labels=c("warm"="Warm", "cold"="Cold" ,"orig"="Original"))+
  labs(x=NULL, y="Length (cm)", color=NULL, size=NULL)+
  theme(strip.background=element_rect(fill="white"),
        text = element_text(size=16))

ggplot() +
  theme_bw()+
  geom_point(data=all_field_data %>% filter(str_detect(source, "Point")), aes(x=date2017, y=mean_length,size="obs"))+
  geom_smooth(data=all_output %>% filter(str_detect(source, "Point"), level!="lit"), aes(x=date2017, y=L_allometric, color=level))+
  facet_grid(year~source, labeller=labeller(year = c("1"="Year 1", "2"="Year 2")))+
  scale_size_manual(values=c("obs"=2), breaks=c("obs"), labels=c("obs"="Observations"))+
  scale_color_manual(values=c("cold"='#0f85a0',"warm"="#dd4124", "orig"="black"),
                     breaks=c("cold","warm","orig"),
                     labels=c("warm"="Warm", "cold"="Cold" ,"orig"="Original"))+
  labs(x=NULL, y="Length (cm)", color=NULL, size=NULL)+
  theme(strip.background=element_rect(fill="white"),
        text = element_text(size=16))


all_rmse %>% 
  filter(level!="orig") %>% 
  ungroup() %>% 
  mutate(imp = if_else(improvement>0, TRUE, FALSE)) %>% 
  group_by(level, year) %>% 
  summarize(num_imp = sum(imp), 
            perc_imp = sum(imp)/length(imp), 
            med_imp = median(improvement),
            mean_imp = mean(improvement))

