#Comparing the effects of temperature correction parameters
#Ruby Krasnow
#Last modified: Nov 30, 2023

ggplot() +
  geom_smooth(data = sled1_sols, aes(Date, L_allometric, color = params), se=FALSE)+
  labs(x= "Date", y = "Blade length (cm)")

ggplot() +
  geom_smooth(data = sled2_sols, aes(Date, L_allometric, color = params), se=FALSE)+
  labs(x= "Date", y = "Blade length (cm)")

ggplot() +
  #geom_smooth(data = sled1_sols_plus2, aes(Temp_C, J_EC_A, color = params), se=FALSE)
  geom_smooth(data = sled1_sols_plus2, aes(Date, L_allometric, color = params), se=FALSE)+
  labs(x= "Date", y = "Blade length (cm)")

#Pt. Judith Pond N 1
ggplot() +
  geom_bar(data = sled1_sols  %>% filter(Date==as.POSIXct("2018-04-17 12:00:00", tz="GMT")),
           aes(x=reorder(params, -L_allometric), y= L_allometric), stat="identity")+
  geom_hline(data=PJN1_meandat %>% 
               mutate(plus=mean_length+sd_length, minus=mean_length-sd_length) %>% 
               select(Date, mean_length, plus,minus) %>% 
               pivot_longer(2:4) %>% filter(Date==as_datetime("2018-04-17")), aes(yintercept=value))+ggtitle("Pt. Judith Pond N 1")

#Pt. Judith Pond N 2
ggplot() +
  geom_bar(data = sled2_sols  %>% filter(Date==as.POSIXct("2018-04-17 12:00:00", tz="GMT")),
           aes(x=reorder(params, -L_allometric), y= L_allometric), stat="identity")+
  geom_hline(data=PJN2_meandat %>% 
               mutate(plus=mean_length+sd_length, minus=mean_length-sd_length) %>% 
               select(Date, mean_length, plus,minus) %>% 
               pivot_longer(2:4) %>% filter(Date==as_datetime("2018-04-17")), aes(yintercept=value))+ggtitle("Pt. Judith Pond N 2")

#Pt. Judith Pond S 1
ggplot() +
  geom_bar(data = dredge1_sols  %>% filter(Date==as.POSIXct("2018-04-22 12:00:00", tz="GMT")),
           aes(x=reorder(params, -L_allometric), y= L_allometric), stat="identity")+
  geom_hline(data=PJS1_meandat %>% 
               mutate(plus=mean_length+sd_length, minus=mean_length-sd_length) %>% 
               select(Date, mean_length, plus,minus) %>% 
               pivot_longer(2:4) %>% filter(Date==as_datetime("2018-04-22")), aes(yintercept=value))+ggtitle("Pt. Judith Pond S 1")

#Pt. Judith Pond S 2
ggplot() +
  geom_bar(data = dredge2_sols  %>% filter(Date==as.POSIXct("2018-04-22 12:00:00", tz="GMT")),
           aes(x=reorder(params, -L_allometric), y= L_allometric), stat="identity")+
  geom_hline(data=PJS2_meandat %>% 
               mutate(plus=mean_length+sd_length, minus=mean_length-sd_length) %>% 
               select(Date, mean_length, plus,minus) %>% 
               pivot_longer(2:4) %>% filter(Date==as_datetime("2018-04-22")), aes(yintercept=value))+ggtitle("Pt. Judith Pond S 2")


# Allometric --------------------------------------------------------------
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

params_allometric <- bind_rows(list("rates_Lo"=params, "rates_L_new"=params), .id="scaling")
params_nested <- params_allometric %>% nest(data = c(T_A, T_H, T_AH))

plan(multisession, workers=availableCores())

output_sled1_yr1 <- params_nested %>% mutate(std_L = future_map(data, function(df) {
  temp_params <- params_Lo
  temp_params[c("T_A", "T_H", "T_AH")] <- c(df$T_A, df$T_H, df$T_AH)
  ode_output<- as.data.frame(ode(y = state_Lo, t = times_Lo_Sled1, func = rates_Lo, parms = temp_params)) %>% select(time, W, L_allometric)
  ode_output
})) %>% mutate(new_L = future_map(data, function(df) {
  temp_params <- params_Lo
  temp_params[c("T_A", "T_H", "T_AH")] <- c(df$T_A, df$T_H, df$T_AH)
  ode_output<- as.data.frame(ode(y = state_Lo, t = times_Lo_Sled1, func = rates_L_new, parms = temp_params)) %>% select(time, L_allometric) %>% mutate(L_allometric_new = L_allometric,.keep="unused")
  ode_output
})) %>% select(-data)

output_sled1_yr1 <- output_sled1_yr1 %>%
  mutate(data = future_map2(std_L, new_L, ~ cbind(.x,.y))) %>%
  select(-c(std_L, new_L))

output_sled1_yr1 <- output_sled1_yr1 %>% unnest(cols=data, names_repair = "unique") %>%
  select(-time...8) %>%
  rename(time=time...5, L_allometric_old = L_allometric) %>%
  pivot_longer(cols=c(L_allometric_old, L_allometric_new), names_to="L_formula", values_to = "L_allometric", names_prefix = "L_allometric_")

output_sled1_yr1_clean <- output_sled1_yr1 %>%
  ungroup() %>% select(-scaling) %>%
  group_by(type, level, res, L_formula) %>% distinct() %>%
  mutate(Temp_C = T_Sled1_Y1-273.15, #conversion back to Celsius from Kelvin
             Date=sled1_date_seq,
             source="Point Judith Pond N 1")

for_comp <- output_sled1_yr1_clean %>% filter(type!="ctrl") %>% ungroup() %>% mutate(
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


comp_rmse <- field_data %>% group_by(source) %>% left_join((for_comp %>% group_by(source)%>% filter(L_formula=="old"))) %>% group_by(params) %>% summarise(rmse = rmse(mean_length, L_allometric)) %>% na.omit() #%>% filter(params=="orig") %>% pull(rmse)

comp_rmse_newL <- field_data %>% group_by(source) %>% left_join((for_comp %>% group_by(source)%>% filter(L_formula=="new"))) %>% group_by(params) %>% summarise(rmse = rmse(mean_length, L_allometric)) %>% na.omit() #%>% filter(params=="orig") %>% pull(rmse)

(ggplot(data=comp_rmse %>% ungroup(), aes(x=reorder(params, rmse), y=rmse, fill=params)) + geom_col()+
  scale_x_reordered()+
  geom_text(aes(label = round(rmse,2), vjust = -0.2))+ggtitle("old L"))/
  
  (ggplot(data=comp_rmse_newL %>% ungroup(), aes(x=reorder(params, rmse), y=rmse, fill=params)) + geom_col()+
     scale_x_reordered()+
     geom_text(aes(label = round(rmse,2), vjust = -0.2))+ggtitle("new L"))

ggplot(data=rmse_dat %>% ungroup() %>% filter(source=="Point Judith Pond N 1"), aes(x=reorder(params, rmse), y=rmse, fill=params)) + 
  geom_col()+
  geom_text(aes(label = round(rmse,1), vjust = -0.2))+
  scale_x_reordered()+ggtitle("first way")
