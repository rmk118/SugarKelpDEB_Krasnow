### Point Judith Point Year 2 ################################################################################################################# 

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
source("KelpDEB_Run_Krasnow_PJPyr1.R")
source("SolveR_R.R")
source("KelpDEB_model.R")
source("./outdoorExpt/outdoorHOBO/outdoor_HOBO.R")
source("./outdoorExpt/outdoor_expt.R")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####### Initial conditions year 2 ############
#Initial conditions of state variables

state_LoY2 <- c(m_EC = 0.01, #0.9 #mol C/molM_V  #Reserve density of C reserve (initial mass of C reserve per intital mass of structure)
                m_EN = 0.09, #mol N/molM_V #Reserve density of N reserve (initial mass of N reserve per intital mass of structure)
                M_V = 0.05/(w_V+0.09*w_EN+0.01*w_EC)) #molM_V #initial mass of structure

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####### Time steps year 2 #######
#(First number of time step, last number of time step, interval to step)
times_Y2_Sled1 <- seq(0, 3408, 1) #142 days stepped hourly
times_Y2_Sled2 <- seq(0, 2064, 1) #86 days stepped hourly
times_Y2_Dredge1 <- seq(0, 3408, 1) #142 days stepped hourly
times_Y2_Dredge2 <- seq(0, 2064, 1) #86 days stepped hourly

##### N forcing set-up year 2 ####
WSA_Y2 <- read.csv("WaterSamplesY2.csv", header = TRUE, fileEncoding="UTF-8-BOM") #Import water Q data
WSA_Y2$Date <- mdy(WSA_Y2$Date) #convert dates
names(WSA_Y2)[1] <- "Site" #only necessary for some computers running this code

Sled_WSA2 <- filter(WSA_Y2, Site == "Moonstone Sled") #filter by site
Dredge_WSA2 <- filter(WSA_Y2, Site == "Moonstone Dredge")

Sled_WSA2$NO3NO2_µM <- Sled_WSA2$NO3NO2_µM/1000000 #convert from micromoles/L to moles/L
Dredge_WSA2$NO3NO2_µM <- Dredge_WSA2$NO3NO2_µM/1000000

#No new DIC data, year 1 variables re-used

###### Temp forcing set-up year 2 #############

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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### FIELD DATA MODEL RUNS YR 2 ####

### Judith N (sled) line 1 ####