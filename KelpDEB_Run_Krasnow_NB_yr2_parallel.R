### Narragansett Bay Year 2 #################################################################################################################

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
####### Initial conditions year 2 ############
state_LoY2 <- c(m_EC = 0.01, #0.9 #mol C/molM_V  #Reserve density of C reserve (initial mass of C reserve per intital mass of structure)
                m_EN = 0.09, #mol N/molM_V #Reserve density of N reserve (initial mass of N reserve per intital mass of structure)
                M_V = 0.05/(w_V+0.09*w_EN+0.01*w_EC)) #molM_V #initial mass of structure

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####### Time steps year 2 #######
#(First number of time step, last number of time step, interval to step)
times_Y2_W <- seq(0, 3720, 1) #155 days stepped hourly
times_Y2_R1 <- seq(0, 3720, 1) #155 days stepped hourly
times_Y2_R2 <- seq(0, 2208, 1) #92 days stepped hourly

##### N forcing set-up year 2 ####

#Wickford
WSA_Y2 <- read.csv("WaterSamplesY2.csv", header = TRUE, fileEncoding="UTF-8-BOM") #Import water Q data
WSA_Y2$Date <- mdy(WSA_Y2$Date) #convert dates
names(WSA_Y2)[1] <- "Site" #only necessary for some computers running this code

Wickford_WSA2 <- filter(WSA_Y2, Site == "Wickford") #filter by site
Wickford_WSA2$NO3NO2_µM <- Wickford_WSA2$NO3NO2_µM/1000000 #convert from micromoles/L to moles/L

###### DIC forcing set-up year 2 ###########
Segarra2002Carbon <- read.csv("BrentonPoint_Segarra2002CarbonData.csv", header = TRUE, fileEncoding="UTF-8-BOM") #Import lit TCO2 data
names(Segarra2002Carbon)[1] <- "Date" #Only necessary for running the code on some computer
Segarra2002Carbon$Date <- mdy(Segarra2002Carbon$Date) #time field conversion
CO_2 <- mean(Segarra2002Carbon$TCO2_micromolPERkg)/1000000 #(mol CO2/L)

# Using same data frame for irradiance, just different subset

###### Temp forcing set-up year 2 #############
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

W_date_seq_Y2 <- seq(as_datetime("2018-12-19 12:00:00"), as_datetime("2019-05-23 12:00:00"), by="hour")
R1_date_seq_Y2 <- seq(as_datetime("2018-12-20 12:00:00"), as_datetime("2019-05-24 12:00:00"), by="hour")
R2_date_seq_Y2 <-seq(as_datetime("2019-2-21 12:00:00"), as_datetime("2019-05-24 12:00:00"), by="hour")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### FIELD DATA MODEL RUNS YR 2 ####

### Wickford ####
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
# (the differential equation solver)

plan(multisession, workers=availableCores())

output_W_yr2 <- params_nested %>% mutate(std_L = future_map(data, function(df) {
  temp_params <- params_Lo
  temp_params[c("T_A", "T_H", "T_AH")] <- c(df$T_A, df$T_H, df$T_AH)
  ode_output<- as.data.frame(ode(y = state_LoY2, t = times_Y2_W, func = rates_Lo, parms = temp_params)) %>% select(time, W, L_allometric)
  return(ode_output)
})) %>% select(-data)

output_W_yr2_clean <- output_W_yr2 %>%
  unnest(cols = std_L) %>% 
  ungroup() %>% 
  group_by(type, level, res) %>% 
  mutate(Temp_C = T_W_Y2-273.15, #conversion back to Celsius from Kelvin
         Date=W_date_seq_Y2,
         source="Narragansett Bay N")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setting up the forcing functions with field data for
### Rome Pt line 1 ####

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
  group_by(type, level, res) %>% 
  mutate(Temp_C = T_R1_Y2-273.15, #conversion back to Celsius from Kelvin
         Date=R1_date_seq_Y2,
         source="Narragansett Bay S 1")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setting up the forcing functions with field data for
### Rome Pt line 2 ####
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
  group_by(type, level, res) %>% 
  mutate(Temp_C = T_R2_Y2-273.15, #conversion back to Celsius from Kelvin
         Date=R2_date_seq_Y2,
         source="Narragansett Bay S 2")

#combine all output into one dataframe
all_output_NB_yr2 <- rbind(output_W_yr2_clean, output_R1_yr2_clean, output_R2_yr2_clean) %>% ungroup() %>% mutate(
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
NBN1_Y2_meandat <- KelpY2[KelpY2$SiteLine == "Narragansett Bay N 1",] %>%
  group_by(Date) %>%
  summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length, na.rm = TRUE)) %>%
  mutate(Date= as.POSIXct(Date)) %>%
  na.omit()

NBS1_Y2_meandat <- KelpY2[KelpY2$SiteLine == "Narragansett Bay S 1",] %>%
  group_by(Date) %>%
  summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length, na.rm = TRUE)) %>%
  mutate(Date= as.POSIXct(Date)) %>%
  na.omit()

NBS2_Y2_meandat <- KelpY2[KelpY2$SiteLine == "Narragansett Bay S 2",] %>%
  group_by(Date) %>%
  summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length, na.rm = TRUE)) %>%
  mutate(Date= as.POSIXct(Date)) %>%
  na.omit()


field_data_NB_Y2 <- bind_rows(list("Narragansett Bay N" = NBN1_Y2_meandat, "Narragansett Bay S 1" = NBS1_Y2_meandat, "Narragansett Bay S 2" =NBS2_Y2_meandat), .id="source")

### Combine with model data ####
rmse_NB_Y2 <- field_data_NB_Y2 %>% group_by(source) %>% 
  left_join((all_output_NB_yr2 %>% group_by(source))) %>% 
  group_by(source, params, type) %>% summarise(rmse = rmse(mean_length, L_allometric))

rmse_NB_Y2 %>% filter(params=="orig") #check to make sure it's the same as original paper

ggplot(data=rmse_NB_Y2 %>% ungroup() %>% filter(type!="stress"), aes(x=reorder_within(params, rmse, source), y=rmse, fill=params)) +
  geom_col()+
  geom_text(aes(label = round(rmse,1), vjust = -0.2))+
  facet_wrap(~source, scales = "free_x")+
  scale_x_reordered()

ggplot(data=rmse_NB_Y2 %>% ungroup() %>% filter(type!="ctrl"), aes(x=reorder_within(params, rmse, source), y=rmse, fill=params)) +
  geom_col()+
  geom_text(aes(label = round(rmse,1), vjust = -0.2))+
  facet_wrap(~source, scales = "free_x")+
  scale_x_reordered()+ggtitle("Stress Y2")+labs(y="RMSE",x=NULL)

ggplot(all_output_NB_yr1 %>% filter(params %in% c("orig", "high", "high_rep", "high_cross")))+
  geom_smooth(aes(x=Date, y=L_allometric, color=params, linetype=type))+
  geom_point(data=field_data_NB_Y1, aes(x=Date, y=mean_length))+
  facet_wrap(~source, scales = "free_x")
  
View(rmse_NB_Y2 %>% ungroup() %>% filter(params %in% c("orig", "high", "high_rep", "high_cross")))


all_rmse <- bind_rows(rmse_dat_new %>% mutate(year=1), rmse_dat_Y2 %>% mutate(year=2), rmse_NB_Y1 %>% mutate(year=1), rmse_NB_Y2%>% mutate(year=2))

orig_rmse <- all_rmse %>% filter(params =="orig") %>% ungroup()

all_rmse <- all_rmse %>% bind_rows(orig_rmse %>% mutate(type="stress")) %>% bind_rows(orig_rmse %>% mutate(type="ctrl"))

all_rmse <- all_rmse %>% ungroup() %>% left_join(orig_rmse %>% mutate(orig_rmse = rmse) %>% select(source, year, orig_rmse), by=c("source", "year")) %>% mutate(improvement = orig_rmse-rmse)

perc_imp <- all_rmse %>% mutate(imp = if_else(improvement>0, TRUE, FALSE)) %>% group_by(params, type) %>% summarize(num_imp = sum(imp), perc_imp = sum(imp)/14, mean_imp = mean(if_else(improvement>0, improvement, 0)))
