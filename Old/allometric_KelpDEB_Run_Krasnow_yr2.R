### Point Judith Point Year 2 ################################################################################################################# 

# INCLUDES NEW ALLOMETRIC EQUATION

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
source("KelpDEB_Run_Krasnow_yr1_parallel.R")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####### Initial conditions year 2 ############
state_LoY2 <- c(m_EC = 0.01, #0.9 #mol C/molM_V  #Reserve density of C reserve (initial mass of C reserve per intital mass of structure)
                m_EN = 0.09, #mol N/molM_V #Reserve density of N reserve (initial mass of N reserve per intital mass of structure)
                M_V = 0.05/(w_V+0.09*w_EN+0.01*w_EC)) #molM_V #initial mass of structure

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

#No new DIC data, year 1 variables re-used. Using same data frame for irradiance, just different subset

###### Temp forcing set-up year 2 #############
# Point Judith Pond N (sled)
Sled_Y2_Hobo_orig <- read.csv("Sled_Y2_HoboLightTemp.csv", header = TRUE, fileEncoding="UTF-8-BOM") #import
Sled_Y2_Hobo_orig$DateTime <- mdy_hms(Sled_Y2_Hobo_orig$DateTime) #convert date time field
Sled_Y2_Hobo_orig$Temp_K <- Sled_Y2_Hobo_orig$Temp_C+273.15 #create column with temp in K
SledY2T_hourly <- ceiling_date(Sled_Y2_Hobo_orig$DateTime, unit = "hour") #determine times to aggregate around
AvgTempKbyhr <- aggregate(Sled_Y2_Hobo_orig$Temp_K, by=list(SledY2T_hourly), mean) #calculate average hourly temp

# Point Judith Pond S (dredge)
Dredge_Y2_Hobo_orig <- read.csv("Dredge_Y2_HoboTempLight.csv", header = TRUE, fileEncoding="UTF-8-BOM") #import
Dredge_Y2_Hobo_orig$DateTime <- mdy_hms(Dredge_Y2_Hobo_orig$DateTime) #convert date time field
Dredge_Y2_Hobo_orig$Temp_K <- Dredge_Y2_Hobo_orig$Temp_C+273.15 #create column with temp in K
DredgeY2T_hourly <- ceiling_date(Dredge_Y2_Hobo_orig$DateTime, unit = "hour") #determine what times to aggregate around
AvgTempKbyhr <- aggregate(Dredge_Y2_Hobo$Temp_K, by=list(DredgeY2T_hourly), mean) #calculate average hourly temp

L1_date_seq_Y2 <- seq(as_datetime("2018-12-12 12:00:00"), as_datetime("2019-05-03 12:00:00"), by="hour") 
L2_date_seq_Y2 <- seq(as_datetime("2019-02-06 12:00:00"), as_datetime("2019-05-03 12:00:00"), by="hour")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### FIELD DATA MODEL RUNS YR 2 ####

### Sled line 1 ####
N_field <- approxfun(x = c(1*24, 57*24, 93*24, 124*24, 142*24, 163*24), y = c(Sled_WSA2$NO3NO2_µM), method = "linear", rule = 2) #N forcing function

###### Irradiance forcing set-up #
NOAA_Irradiance_Sledy2 <-  NOAA_Irradiance$PAR[5686:6822] #subset by seq(as_datetime("2018-12-12 12:00:00"), as_datetime("2019-05-03 12:00:00"), by="hour")
I_field <- approxfun(x = seq(from = 0, to = 3408, by = 3), y = NOAA_Irradiance_Sledy2, method = "linear", rule = 2) #irradiance forcing function

###### Temp forcing set-up #
AvgTempKbyhr <- AvgTempKbyhr[2:3385,] #subset
fd <- rep(285, 25) #small section of replacement
T_field <- approxfun(x = c(0:3408), y = c(AvgTempKbyhr$x, fd), method = "linear", rule = 2) #the temp forcing function
T_Sled1_Y2 <- T_field(0:3408) #for ease in plotting the temperature forcing

###### Model runs ######
# (the differential equation solver)

plan(multisession, workers=availableCores())

output_sled1_yr2 <- params_nested %>% mutate(std_L = future_map(data, function(df) {
  temp_params <- params_Lo
  temp_params[c("T_A", "T_H", "T_AH")] <- c(df$T_A, df$T_H, df$T_AH)
  ode_output<- as.data.frame(ode(y = state_LoY2, t = times_Y2_Sled1, func = rates_Lo, parms = temp_params)) %>% select(time, W, L_allometric)
  return(ode_output)
})) %>% mutate(new_L = future_map(data, function(df) {
  temp_params <- params_Lo
  temp_params[c("T_A", "T_H", "T_AH")] <- c(df$T_A, df$T_H, df$T_AH)
  ode_output<- as.data.frame(ode(y = state_LoY2, t = times_Y2_Sled1, func = rates_L_new, parms = temp_params)) %>% select(time, L_allometric) %>% mutate(L_allometric_new = L_allometric,.keep="unused")
  return(ode_output)
})) %>% select(-data)

output_sled1_yr2 <- output_sled1_yr2 %>%
  mutate(data = future_map2(std_L, new_L, ~ cbind(.x,.y))) %>%
  select(-c(std_L, new_L, scaling)) %>% 
  distinct() %>% 
  rowwise() %>% 
  mutate(data = list(data[,c(1,2,3,5)])) %>% 
  unnest(cols=data) %>%
  rename(L_allometric_old = L_allometric) %>%
  pivot_longer(cols=c(L_allometric_old, L_allometric_new), names_to="L_formula", values_to = "L_allometric", names_prefix = "L_allometric_")

output_sled1_yr2_clean <- output_sled1_yr2 %>%
  ungroup() %>% 
  group_by(type, level, res, L_formula) %>% 
  mutate(Temp_C = T_Sled1_Y2-273.15, #conversion back to Celsius from Kelvin
         Date=L1_date_seq_Y2,
         source="Point Judith Pond N 1")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setting up the forcing functions with field data for
### Sled line 2 ####

# N forcing set-up Judith N 2 #
Sled_WSA2_sub <- Sled_WSA2[2:6,] #remove the point before the relevant time range
N_field <- approxfun(x = c(1*24, 37*24, 68*24, 86*24, 107*24), y = c(Sled_WSA2_sub$NO3NO2_µM), method = "linear", rule = 2) #N forcing function

###### Irradiance forcing set-up #
NOAA_Irradiance_Sledy2_L2 <-  NOAA_Irradiance$PAR[6134:6822] #subset by seq(as_datetime("2019-02-06 12:00:00"), as_datetime("2019-05-03 12:00:00"), by="hour")
I_field <- approxfun(x = seq(from = 0, to = 2064, by = 3), y = NOAA_Irradiance_Sledy2_L2, method = "linear", rule = 2) #irradiance forcing function

###### Temp forcing set-Up #
AvgTempKbyhr <- aggregate(Dredge_Y2_Hobo$Temp_K, by=list(DredgeY2T_hourly), mean) #calculate average hourly temp
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
})) %>% mutate(new_L = future_map(data, function(df) {
  temp_params <- params_Lo
  temp_params[c("T_A", "T_H", "T_AH")] <- c(df$T_A, df$T_H, df$T_AH)
  ode_output<- as.data.frame(ode(y = state_LoY2, t = times_Y2_Sled2, func = rates_L_new, parms = temp_params)) %>% select(time, L_allometric) %>% mutate(L_allometric_new = L_allometric,.keep="unused")
  ode_output
})) %>% select(-data)

output_sled2_yr2 <- output_sled2_yr2 %>%
  mutate(data = future_map2(std_L, new_L, ~ cbind(.x,.y))) %>%
  select(-c(std_L, new_L, scaling)) %>% 
  distinct() %>% 
  rowwise() %>% 
  mutate(data = list(data[,c(1,2,3,5)])) %>% 
  unnest(cols=data) %>%
  rename(L_allometric_old = L_allometric) %>%
  pivot_longer(cols=c(L_allometric_old, L_allometric_new), names_to="L_formula", values_to = "L_allometric", names_prefix = "L_allometric_")

output_sled2_yr2_clean <- output_sled2_yr2 %>%
  ungroup() %>% 
  group_by(type, level, res, L_formula) %>%
  mutate(Temp_C = T_Sled2_Y2-273.15, #conversion back to Celsius from Kelvin
         Date=L2_date_seq_Y2,
         source="Point Judith Pond N 2")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setting up the forcing functions with field data for
### Dredge line 1 ####

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
})) %>% mutate(new_L = future_map(data, function(df) {
  temp_params <- params_Lo
  temp_params[c("T_A", "T_H", "T_AH")] <- c(df$T_A, df$T_H, df$T_AH)
  ode_output<- as.data.frame(ode(y = state_LoY2, t = times_Y2_Dredge1, func = rates_L_new, parms = temp_params)) %>% select(time, L_allometric) %>% mutate(L_allometric_new = L_allometric,.keep="unused")
  ode_output
})) %>% select(-data)

output_dredge1_yr2 <- output_dredge1_yr2 %>%
  mutate(data = future_map2(std_L, new_L, ~ cbind(.x,.y))) %>%
  select(-c(std_L, new_L, scaling)) %>% 
  distinct() %>% 
  rowwise() %>% 
  mutate(data = list(data[,c(1,2,3,5)])) %>% 
  unnest(cols=data) %>%
  rename(L_allometric_old = L_allometric) %>%
  pivot_longer(cols=c(L_allometric_old, L_allometric_new), names_to="L_formula", values_to = "L_allometric", names_prefix = "L_allometric_")

output_dredge1_yr2_clean <- output_dredge1_yr2 %>%
  ungroup() %>%
  group_by(type, level, res, L_formula) %>%
  mutate(Temp_C = T_Dredge1_Y2-273.15, #conversion back to Celsius from Kelvin
         Date=L1_date_seq_Y2,
         source="Point Judith Pond S 1")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setting up the forcing functions with field data for
### Dredge line 2 ####

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
})) %>% mutate(new_L = future_map(data, function(df) {
  temp_params <- params_Lo
  temp_params[c("T_A", "T_H", "T_AH")] <- c(df$T_A, df$T_H, df$T_AH)
  ode_output<- as.data.frame(ode(y = state_LoY2, t = times_Y2_Dredge2, func = rates_L_new, parms = temp_params)) %>% select(time, L_allometric) %>% mutate(L_allometric_new = L_allometric,.keep="unused")
  ode_output
})) %>% select(-data)

output_dredge2_yr2 <- output_dredge2_yr2 %>%
  mutate(data = future_map2(std_L, new_L, ~ cbind(.x,.y))) %>%
  select(-c(std_L, new_L, scaling)) %>% 
  distinct() %>% 
  rowwise() %>% 
  mutate(data = list(data[,c(1,2,3,5)])) %>% 
  unnest(cols=data) %>%
  rename(L_allometric_old = L_allometric) %>%
  pivot_longer(cols=c(L_allometric_old, L_allometric_new), names_to="L_formula", values_to = "L_allometric", names_prefix = "L_allometric_")

output_dredge2_yr2_clean <- output_dredge2_yr2 %>%
  ungroup() %>%
  group_by(type, level, res, L_formula) %>%
  mutate(Temp_C = T_Dredge2_Y2-273.15, #conversion back to Celsius from Kelvin
         Date=L2_date_seq_Y2,
         source="Point Judith Pond S 2")

#combine all output into one dataframe
all_output_yr2 <- rbind(output_sled1_yr2_clean, output_sled2_yr2_clean, output_dredge1_yr2_clean, output_dredge2_yr2_clean) %>% ungroup() %>% mutate(
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
KelpY2 <- read.csv("Year2kelpdata.csv", header = TRUE, fileEncoding="UTF-8-BOM")
names(KelpY2)[2] <- "Site"
KelpY2 <- filter(KelpY2, Site != "Fox Island")
KelpY2$Date <- mdy(KelpY2$SamplingDate)
KelpY2$SiteLine <- paste(KelpY2$Site, KelpY2$Line)
KelpY2 <- filter(KelpY2, SiteLine != "Narragansett Bay N 2")


### Combine with model data ####
#Point Judith Point N (sled) Line 1
PJN1_Y2_meandat <- KelpY2[KelpY2$SiteLine == "Point Judith Pond N 1",] %>%
  group_by(Date) %>%
  summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length, na.rm = TRUE)) %>%
  mutate(Date= as.POSIXct(Date)) %>%
  na.omit()

#Point Judith Point N (sled) Line 2
PJN2_Y2_meandat <- KelpY2[KelpY2$SiteLine == "Point Judith Pond N 2",] %>%
  group_by(Date) %>%
  summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length, na.rm = TRUE))  %>%
  mutate(Date= as.POSIXct(Date)) %>%
  na.omit()

#Point Judith Point S (dredge) Line 1
PJS1_Y2_meandat <- KelpY2[KelpY2$SiteLine == "Point Judith Pond S 1",] %>%
  group_by(Date) %>%
  summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length, na.rm = TRUE))  %>%
  mutate(Date= as.POSIXct(Date)) %>%
  na.omit()

#Point Judith Point S (dredge) Line 2
PJS2_Y2_meandat <- KelpY2[KelpY2$SiteLine == "Point Judith Pond S 2",] %>%
  group_by(Date) %>%
  summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length, na.rm = TRUE))  %>%
  mutate(Date= as.POSIXct(Date)) %>%
  na.omit()

field_data_Y2 <- bind_rows(list("Point Judith Pond N 1" = PJN1_Y2_meandat, "Point Judith Pond N 2" = PJN2_Y2_meandat, "Point Judith Pond S 1" = PJS1_Y2_meandat, "Point Judith Pond S 2" =PJS2_Y2_meandat), .id="source")

rmse_dat_Y2 <- field_data_Y2 %>% group_by(source) %>% left_join((all_output_yr2 %>% group_by(source))) %>% group_by(source, params, L_formula, type) %>% summarise(rmse = rmse(mean_length, L_allometric))

rmse_dat_Y2 %>% filter(params=="orig", L_formula=="old") #check to make sure it's the same as original paper

ggplot(data=rmse_dat_Y2 %>% ungroup() %>% filter(type!="stress"), aes(x=reorder_within(params, rmse, list(L_formula, source)), y=rmse, fill=params)) +
  geom_col()+
  geom_text(aes(label = round(rmse,1), vjust = -0.2))+
  facet_wrap(L_formula~source, scales = "free_x")+
  scale_x_reordered()
