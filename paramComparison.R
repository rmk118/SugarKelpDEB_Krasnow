#Comparing the effects of temperature correction parameters
#Ruby Krasnow
#Last modified: Nov 27, 2023

#Required for model runs

ggplot() +
  #geom_smooth(data = sled1_sols_plus2, aes(Temp_C, J_EC_A, color = params), se=FALSE)
geom_line(data = sled1_sols_plus2, aes(Date, L_allometric, color = params), se=FALSE)


ggplot() +
  geom_smooth(data = dredge1_sols, aes(Date, L_allometric, color = params), se=FALSE)
# 
ggplot() +
  #geom_smooth(data = sol_all %>% filter(source=="Point Judith Pond N 1"), aes(Date, L_allometric, color="orig"))+
  geom_smooth(data=clean_ode_sol(sol_low_Sled1,sled1_date_seq, T_Sled1_Y1,"Point Judith Pond N 1"), aes(Date, L_allometric, color="low"))+
geom_smooth(data=clean_ode_sol(sol_med_Sled1,sled1_date_seq, T_Sled1_Y1,"Point Judith Pond N 1"), aes(Date, L_allometric, color="med"))+
geom_smooth(data=clean_ode_sol(sol_high_Sled1,sled1_date_seq, T_Sled1_Y1,"Point Judith Pond N 1"), aes(Date, L_allometric, color="high"))+
  geom_smooth(data=clean_ode_sol(sol_Sled1,sled1_date_seq, T_Sled1_Y1,"Point Judith Pond N 1"), aes(Date, L_allometric, color="orig"))
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