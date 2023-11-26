#Comparing the effects of temperature correction parameters
#Ruby Krasnow
#Last modified: Nov 25, 2023

#Required for model runs
source("KelpDEB_Run_Krasnow.R")
source("./outdoorExpt/outdoorHOBO/outdoor_HOBO.R")
source("./outdoorExpt/outdoor_expt.R")

high_params_for_model <- params_Lo
med_params_for_model <- params_Lo
low_params_for_model <- params_Lo

high_params_for_model[c("T_A", "T_H", "T_AH")]<-c(high_T_A, high_T_H,high_T_AH)
med_params_for_model[c("T_A", "T_H", "T_AH")]<-c(med_T_A, med_T_H,med_T_AH)
low_params_for_model[c("T_A", "T_H", "T_AH")]<-c(low_T_A, low_T_H,low_T_AH)

ggplot() + 
  geom_smooth(data = sol_all, aes(Date, L_allometric, color = source))

ggplot() + 
  geom_smooth(data = sol_all %>% filter(source=="Point Judith Pond N 1"), aes(Date, L_allometric, color="orig"))+
  geom_smooth(data=sol_low_Sled1, aes(Date, L_allometric, color="low"))+
geom_smooth(data=sol_med_Sled1, aes(Date, L_allometric, color="med"))+
geom_smooth(data=sol_high_Sled1, aes(Date, L_allometric, color="high"))


sol_high_Sled1 %>% filter(Date==as.POSIXct("2018-04-17 12:00:00", tz="GMT")) %>% pull(W)
sol_med_Sled1 %>% filter(Date==as.POSIXct("2018-04-17 12:00:00", tz="GMT")) %>% pull(W)
sol_low_Sled1 %>% filter(Date==as.POSIXct("2018-04-17 12:00:00", tz="GMT")) %>% pull(W)
sol_Sled1 %>% filter(Date==as.POSIXct("2018-04-17 12:00:00", tz="GMT")) %>% pull(W)

sol_high_Sled1 %>% filter(Date==as.POSIXct("2018-04-17 12:00:00", tz="GMT")) %>% pull(L_allometric)
sol_med_Sled1 %>% filter(Date==as.POSIXct("2018-04-17 12:00:00", tz="GMT")) %>% pull(L_allometric)
sol_low_Sled1 %>% filter(Date==as.POSIXct("2018-04-17 12:00:00", tz="GMT")) %>% pull(L_allometric)
sol_Sled1 %>% filter(Date==as.POSIXct("2018-04-17 12:00:00", tz="GMT")) %>% pull(L_allometric)