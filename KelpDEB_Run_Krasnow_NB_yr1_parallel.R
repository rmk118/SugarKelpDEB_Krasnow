### Narragansett Bay Year 1 ################################################################################################################# 

#Site names here begin with names other than those used in the manuscript
#Wickford = Narragansett Bay N
#Rome Point = Narragansett Bay S

#File originally created by Celeste Venolia in March 2018-December 2019 for Venolia et al., (2020)
# https://doi.org/10.1016/j.ecolmodel.2020.109151 

# Modified by Ruby Krasnow in December 2023
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
#source("KelpDEB_Run_Krasnow_yr1_parallel.R")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####### Initial conditions year 1 ############
#Initial conditions of state variables
#these values are not coming from any field data or literature information, estimated
state_Lo <- c(m_EC = 0.002, #0.1, #mol C/molM_V  #Reserve density of C reserve (initial mass of C reserve per initial mass of structure)
              m_EN = 0.01, #mol N/molM_V #Reserve density of N reserve (initial mass of N reserve per initial mass of structure)
              M_V = 0.05/(w_V+0.01*w_EN+0.002*w_EC)) #molM_V #initial mass of structure

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####### Time steps year 1 #######
#(First number of time step, last number of time step, interval to step)
times_Y1_W <- seq(0, 3312, 1) #138 days stepped hourly
times_Y1_R1 <- seq(0, 4104, 1) #171 days stepped hourly
times_Y1_R2 <- seq(0, 3264, 1) #136 days stepped hourly

##### N forcing set-up year 1 ####

#Wickford
WSA2_Y1 <- read.csv("WaterSampleAnalysis2Y1.csv", header = TRUE, fileEncoding="UTF-8-BOM") #Import water quality data
WSA2_Y1$Date <- mdy(WSA2_Y1$Date) #convert dates
names(WSA2_Y1)[1] <- "Site" #only necessary to run this code on some computers
Wickford_WSA <- filter(WSA2_Y1, Site == "Wickford") #filter by site

#Rome Pt
RomePt_WSA <- filter(WSA2_Y1, Site == "Rome Point") #filter by site

###### DIC forcing set-up year 1 ###########
Segarra2002Carbon <- read.csv("BrentonPoint_Segarra2002CarbonData.csv", header = TRUE, fileEncoding="UTF-8-BOM") #Import lit TCO2 data
names(Segarra2002Carbon)[1] <- "Date" #Only necessary for running the code on some computer
Segarra2002Carbon$Date <- mdy(Segarra2002Carbon$Date) #time field conversion
CO_2 <- mean(Segarra2002Carbon$TCO2_micromolPERkg)/1000000 #(mol CO2/L)

# Using same data frame for irradiance, just different subset

###### Temp forcing set-up year 1 #############
# NB N (Wickford)
Wickford_Y1_hobo_orig <- read.csv("Wickford_Y1_hobo.csv", header = TRUE, fileEncoding="UTF-8-BOM") #Import Wickford Hobo data
Wickford_Y1_hobo_orig$DateTime <- mdy_hms(Wickford_Y1_hobo_orig$DateTime) #convert time field
Wickford_Y1_hobo_orig$Temp_K <- Wickford_Y1_hobo_orig$Temp_C+273.15 #create column with temp in K

# NB S (Rome)
RomePoint_Y1_hobotemp_orig <- read.csv("RomePoint_Y1_hobotemp.csv", header = TRUE, fileEncoding="UTF-8-BOM") #import
RomePoint_Y1_hobotemp_orig$DateTime <- mdy_hms(RomePoint_Y1_hobotemp_orig$DateTime) #convert time field
RomePoint_Y1_hobotemp_orig$Temp_K <- RomePoint_Y1_hobotemp_orig$Temp_C+273.15 #create column with temp in K

W_date_seq_Y1 <- seq(as_datetime("2017-12-4 12:00:00"), as_datetime("2018-04-21 12:00:00"), by="hour")
R1_date_seq_Y1 <- seq(as_datetime("2017-11-1 12:00:00"), as_datetime("2018-04-21 12:00:00"), by="hour")
R2_date_seq_Y1 <-seq(as_datetime("2017-12-6 12:00:00"), as_datetime("2018-04-21 12:00:00"), by="hour")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### FIELD DATA MODEL RUNS YR 1 ####

### Wickford ####
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
# (the differential equation solver)

#plan(multisession, workers=availableCores())

output_W_yr1 <- params_nested %>% mutate(std_L = future_map(data, function(df) {
  temp_params <- params_Lo
  temp_params[c("T_A", "T_H", "T_AH")] <- c(df$T_A, df$T_H, df$T_AH)
  ode_output<- as.data.frame(ode(y = state_Lo, t = times_Y1_W, func = rates_Lo, parms = temp_params)) %>% select(time, W, L_allometric)
  return(ode_output)
})) %>% select(-data)

output_W_yr1_clean <- output_W_yr1 %>%
  unnest(cols = std_L) %>% 
  ungroup() %>% 
  group_by(type, level, res) %>% 
  mutate(Temp_C = T_W_Y1-273.15, #conversion back to Celsius from Kelvin
         Date=W_date_seq_Y1,
         source="Narragansett Bay N")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setting up the forcing functions with field data for
### Rome Pt line 1 ####

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
  group_by(type, level, res) %>% 
  mutate(Temp_C = T_R1_Y1-273.15, #conversion back to Celsius from Kelvin
         Date=R1_date_seq_Y1,
         source="Narragansett Bay S 1")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setting up the forcing functions with field data for
### Rome Pt line 2 ####
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
  group_by(type, level, res) %>% 
  mutate(Temp_C = T_R2_Y1-273.15, #conversion back to Celsius from Kelvin
         Date=R2_date_seq_Y1,
         source="Narragansett Bay S 2")

#combine all output into one dataframe
all_output_NB_yr1 <- rbind(output_W_yr1_clean, output_R1_yr1_clean, output_R2_yr1_clean) %>% ungroup() %>% mutate(
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

field_data_fun <- function(df) {
  df_out <- df %>%
    group_by(Date) %>%
    summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length, na.rm = TRUE)) %>%
    mutate(Date= as.POSIXct(Date)) %>%
    na.omit()
  df_out
}

NBN1_Y1_meandat <- KelpY1[KelpY1$SiteLine == "Narragansett Bay N 1",] %>% field_data_fun()

NBS1_Y1_meandat <- KelpY1[KelpY1$SiteLine == "Narragansett Bay S 1",] %>% field_data_fun()

NBS2_Y1_meandat <- KelpY1[KelpY1$SiteLine == "Narragansett Bay S 2",] %>% field_data_fun()

field_data_NB_Y1 <- bind_rows(list("Narragansett Bay N" = NBN1_Y1_meandat, "Narragansett Bay S 1" = NBS1_Y1_meandat, "Narragansett Bay S 2" =NBS2_Y1_meandat), .id="source")

### Combine with model data ####
rmse_NB_Y1 <- field_data_NB_Y1 %>% group_by(source) %>% left_join((all_output_NB_yr1 %>% group_by(source))) %>% group_by(source, params,type) %>% summarise(rmse = rmse(mean_length, L_allometric))

rmse_NB_Y1 %>% filter(params=="orig") #check to make sure it's the same as original paper


ggplot(data=rmse_NB_Y1 %>% ungroup() %>% filter(type!="stress"), aes(x=reorder_within(params, rmse, source), y=rmse, fill=params)) +
  geom_col()+
  geom_text(aes(label = round(rmse,1), vjust = -0.2))+
  facet_wrap(~source, scales = "free_x")+
  scale_x_reordered()

ggplot(data=all_output_NB_yr1 %>% filter(params=="high"|params=="orig"), aes(x=Date, y=L_allometric, color=params, linetype=type)) +
  geom_smooth()+
  facet_wrap(~source, scales="free_y")
