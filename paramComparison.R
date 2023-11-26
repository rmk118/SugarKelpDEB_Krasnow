#Comparing the effects of temperature correction parameters
#Ruby Krasnow
#Last modified: Nov 25, 2023

#Required for model runs

ggplot() +
  geom_smooth(data = sled2_sols, aes(Date, M_V, color = params))
# 
# ggplot() + 
#   geom_smooth(data = sol_all %>% filter(source=="Point Judith Pond N 1"), aes(Date, L_allometric, color="orig"))+
#   geom_smooth(data=sol_low_Sled1, aes(Date, L_allometric, color="low"))+
# geom_smooth(data=sol_med_Sled1, aes(Date, L_allometric, color="med"))+
# geom_smooth(data=sol_high_Sled1, aes(Date, L_allometric, color="high"))
# 
# 
# sol_high_Sled1 %>% filter(Date==as.POSIXct("2018-04-17 12:00:00", tz="GMT")) %>% pull(W)
# sol_med_Sled1 %>% filter(Date==as.POSIXct("2018-04-17 12:00:00", tz="GMT")) %>% pull(W)
# sol_low_Sled1 %>% filter(Date==as.POSIXct("2018-04-17 12:00:00", tz="GMT")) %>% pull(W)
# sol_Sled1 %>% filter(Date==as.POSIXct("2018-04-17 12:00:00", tz="GMT")) %>% pull(W)
# 
# sol_high_Sled1 %>% filter(Date==as.POSIXct("2018-04-17 12:00:00", tz="GMT")) %>% pull(L_allometric)
# sol_med_Sled1 %>% filter(Date==as.POSIXct("2018-04-17 12:00:00", tz="GMT")) %>% pull(L_allometric)
# sol_low_Sled1 %>% filter(Date==as.POSIXct("2018-04-17 12:00:00", tz="GMT")) %>% pull(L_allometric)
# sol_Sled1 %>% filter(Date==as.POSIXct("2018-04-17 12:00:00", tz="GMT")) %>% pull(L_allometric)