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
library(furrr) #parallel processing version of purrr, to speed up model runs
library(tidytext)

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

params_allometric <- bind_rows(list("rates_Lo"=params))
params_nested <- params_allometric %>% nest(data = c(T_A, T_H, T_AH))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####### Initial conditions year 1 ############
#Initial conditions of state variables
#these values are not coming from any field data or literature information, estimated
state_Lo <- c(m_EC = 0.002, #0.1, #mol C/molM_V  #Reserve density of C reserve (initial mass of C reserve per initial mass of structure)
              m_EN = 0.01, #mol N/molM_V #Reserve density of N reserve (initial mass of N reserve per initial mass of structure)
              M_V = 0.05/(w_V+0.01*w_EN+0.002*w_EC)) #molM_V #initial mass of structure

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####### Time steps year 1 #######
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

##### N forcing set-up year 1 ####
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

##### Temp forcing set-up year 1 ###########
Sled_Y1_hobotemp_orig <- read.csv("Sled_Y1_TempLogger2.csv", header = TRUE, fileEncoding="UTF-8-BOM") #import
Sled_Y1_hobotemp_orig$DateTime <- mdy_hms(Sled_Y1_hobotemp_orig$Date_Time) #convert time field
Sled_Y1_hobotemp_orig$Temp_K <- Sled_Y1_hobotemp_orig$Temp_C+273.15 #create column with temp in K

Dredge_Y1_hobo_orig <- read.csv("Dredge_Y1_hobo.csv", header = TRUE, fileEncoding="UTF-8-BOM") #importing neighboring temp file to replace corrupted section
Dredge_Y1_hobo_orig$DateTime <- mdy_hms(Dredge_Y1_hobo_orig$Date_Time) #convert time field
Dredge_Y1_hobo_orig$Temp_K <- Dredge_Y1_hobo_orig$Temp_C+273.15 #create column with temp in K

sled1_date_seq <- seq(as_datetime("2017-11-1 12:00:00"), as_datetime("2018-04-17 12:00:00"), by="hour")
sled2_date_seq <- seq(as_datetime("2017-11-29 12:00:00"), as_datetime("2018-04-17 12:00:00"), by="hour")
dredge1_date_seq<- seq(as_datetime("2017-11-1 12:00:00"), as_datetime("2018-04-22 12:00:00"), by="hour")
dredge2_date_seq<-seq(as_datetime("2017-11-29 12:00:00"), as_datetime("2018-04-22 12:00:00"), by="hour")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### FIELD DATA MODEL RUNS YR 1 ####

### Judith N (sled) line 1 ####
 #initial biomass for conversions (cannot put in initial conditions)
W=0.05
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
  group_by(type, level, res) %>% 
  mutate(Temp_C = T_Sled1_Y1-273.15, #conversion back to Celsius from Kelvin
         Date=sled1_date_seq,
         source="Point Judith Pond N 1")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setting up the forcing functions with field data for
### Judith N (sled) line 2 ####
 #inital biomass for conversions (cannot put in initial conditions)

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

output_sled2_yr1 <- params_nested %>% mutate(std_L = future_map(data, function(df) {
  temp_params <- params_Lo
  temp_params[c("T_A", "T_H", "T_AH")] <- c(df$T_A, df$T_H, df$T_AH)
  ode_output<- as.data.frame(ode(y = state_Lo, t = times_Lo_Sled2, func = rates_Lo, parms = temp_params)) %>% select(time, W, L_allometric)
  ode_output
})) %>% select(-data)

output_sled2_yr1_clean <- output_sled2_yr1 %>%
  unnest(cols=std_L) %>% 
  ungroup() %>% 
  group_by(type, level, res) %>%
  mutate(Temp_C = T_Sled2_Y1-273.15, #conversion back to Celsius from Kelvin
         Date=sled2_date_seq,
         source="Point Judith Pond N 2")
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

output_dredge1_yr1 <- params_nested %>% mutate(std_L = future_map(data, function(df) {
  temp_params <- params_Lo
  temp_params[c("T_A", "T_H", "T_AH")] <- c(df$T_A, df$T_H, df$T_AH)
  ode_output<- as.data.frame(ode(y = state_Lo, t = times_Lo_Dredge1, func = rates_Lo, parms = temp_params)) %>% select(time, W, L_allometric)
  ode_output
})) %>% select(-data)

output_dredge1_yr1_clean <- output_dredge1_yr1 %>%
  ungroup() %>% 
  unnest(cols=std_L) %>% 
  group_by(type, level, res) %>%
  mutate(Temp_C = T_Dredge1_Y1-273.15, #conversion back to Celsius from Kelvin
         Date=dredge1_date_seq,
         source="Point Judith Pond S 1")
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

output_dredge2_yr1 <- params_nested %>% mutate(std_L = future_map(data, function(df) {
  temp_params <- params_Lo
  temp_params[c("T_A", "T_H", "T_AH")] <- c(df$T_A, df$T_H, df$T_AH)
  ode_output<- as.data.frame(ode(y = state_Lo, t = times_Lo_Dredge2, func = rates_Lo, parms = temp_params)) %>% select(time, W, L_allometric)
  ode_output
})) %>% select(-data)

output_dredge2_yr1_clean <- output_dredge2_yr1 %>%
  unnest(cols=std_L) %>% 
  ungroup() %>% 
  group_by(type, level, res) %>%
  mutate(Temp_C = T_Dredge2_Y1-273.15, #conversion back to Celsius from Kelvin
         Date=dredge2_date_seq,
         source="Point Judith Pond S 2")


#combine all output into one dataframe
all_output_yr1 <- rbind(output_sled1_yr1_clean, output_sled2_yr1_clean, output_dredge1_yr1_clean, output_dredge2_yr1_clean) %>% ungroup() %>% mutate(
  params = case_when(
    res=="orig" ~ "orig",
    res=="lit" ~ "new",
    level=="high" & res=="means" ~ "high",
    level=="med" & res=="means" ~ "med",
    level=="low" & res=="means" ~ "low",
    level=="high" & res=="cross" ~ "high_cross",
    level=="med" & res=="cross" ~ "med_cross",
    level=="low" & res=="cross" ~ "low_cross",
    level=="high" & res=="all" ~ "high_rep",
    level=="med" & res=="all" ~ "med_rep",
    level=="low" & res=="all" ~ "low_rep",
  )
) 


### Import field data ####
KelpY1 <- read.csv("Year1kelpdata.csv", header = TRUE, fileEncoding="UTF-8-BOM")
names(KelpY1)[2] <- "Site"
KelpY1 <- filter(KelpY1, Site != "Fox Island")
KelpY1$Date <- mdy(KelpY1$SamplingDate)
KelpY1$SiteLine <- paste(KelpY1$Site, KelpY1$Line)
KelpY1 <- filter(KelpY1, SiteLine != "Narragansett Bay N 2")

### Combine with model data #####

#Point Judith Point N (sled) Line 1
PJN1_meandat <- KelpY1[KelpY1$SiteLine == "Point Judith Pond N 1",] %>%
  group_by(Date) %>%
  summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length, na.rm = TRUE)) %>%
  mutate(Date= as.POSIXct(Date)) %>%
  na.omit()

#Point Judith Point N (sled) Line 2
PJN2_meandat <- KelpY1[KelpY1$SiteLine == "Point Judith Pond N 2",] %>%
  group_by(Date) %>%
  summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length, na.rm = TRUE))  %>%
  mutate(Date= as.POSIXct(Date)) %>%
  na.omit()

#Point Judith Point S (dredge) Line 1
PJS1_meandat <- KelpY1[KelpY1$SiteLine == "Point Judith Pond S 1",] %>%
  group_by(Date) %>%
  summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length, na.rm = TRUE))  %>%
  mutate(Date= as.POSIXct(Date)) %>%
  na.omit()

#Point Judith Point S (dredge) Line 2
PJS2_meandat <- KelpY1[KelpY1$SiteLine == "Point Judith Pond S 2",] %>%
  group_by(Date) %>%
  summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length, na.rm = TRUE))  %>%
  mutate(Date= as.POSIXct(Date)) %>%
  na.omit()

field_data <- bind_rows(list("Point Judith Pond N 1" = PJN1_meandat, "Point Judith Pond N 2" = PJN2_meandat, "Point Judith Pond S 1" = PJS1_meandat, "Point Judith Pond S 2" =PJS2_meandat), .id="source")

rmse_dat_new <- field_data %>% group_by(source) %>% left_join((all_output_yr1 %>% group_by(source))) %>% group_by(source, params, type) %>% summarise(rmse = rmse(mean_length, L_allometric))

rmse_dat_new %>% filter(params=="orig") #check to make sure it's the same as original paper


ggplot(data=rmse_dat_new %>% ungroup() %>% filter(type!="stress"), aes(x=reorder(params, rmse), y=rmse, fill=params)) +
  geom_col()+
  geom_text(aes(label = round(rmse,1), vjust = -0.2))+
  facet_wrap(~source)

ggplot(data=all_output_yr1 %>% filter(res!="cross"), aes(x=Date, y=L_allometric, color=params)) +
  geom_line()+
  facet_grid(type~source)

ggplot(data=rmse_dat_new %>% ungroup()%>% filter(str_detect(params, "cross", TRUE)), aes(x=reorder_within(params, rmse,list(type,source)), y=rmse, fill=params)) +
  geom_col()+
  geom_text(aes(label = round(rmse,1), vjust = -0.2))+
  scale_x_reordered()+
  facet_wrap(type~source, scales = "free_x")

pjp_plot1 <- ggplot(all_output_yr1 %>% filter(params %in% c("orig", "high", "high_rep"), type!="stress"))+
  geom_smooth(aes(x=Date, y=L_allometric, color=params))+
  geom_point(data=field_data, aes(x=Date, y=mean_length))+
  facet_wrap(~source, scales = "free")

pjp_plot2 <-ggplot(all_output_yr2 %>% filter(params %in% c("orig", "high", "high_rep"), type!="stress"))+
  geom_smooth(aes(x=Date, y=L_allometric, color=params))+
  geom_point(data=field_data_Y2, aes(x=Date, y=mean_length))+
  facet_wrap(~source, scales = "free")


nb_plot1 <-ggplot(all_output_NB_yr1 %>% filter(params %in% c("orig", "high", "high_rep"), type!="stress"))+
  geom_smooth(aes(x=Date, y=L_allometric, color=params))+
  geom_point(data=field_data_NB_Y1, aes(x=Date, y=mean_length))+
  facet_wrap(~source, scales = "free_x")

nb_plot2 <- ggplot(all_output_NB_yr2 %>% filter(params %in% c("orig", "high", "high_rep"), type!="stress"))+
  geom_smooth(aes(x=Date, y=L_allometric, color=params))+
  geom_point(data=field_data_NB_Y2, aes(x=Date, y=mean_length))+
  facet_wrap(~source, scales = "free_x")

(nb_plot1+pjp_plot1)/(nb_plot2+pjp_plot2) +
  plot_layout(guides = 'collect')
