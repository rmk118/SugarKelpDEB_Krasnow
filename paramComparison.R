#Comparing the effects of temperature correction parameters
#Ruby Krasnow
#Last modified: Nov 27, 2023

#Required for model runs

ggplot() +
  #geom_smooth(data = sled1_sols_plus2, aes(Temp_C, J_EC_A, color = params), se=FALSE)
geom_smooth(data = sled1_sols_plus2, aes(Date, L_allometric, color = params), se=FALSE)


ggplot() +
  geom_smooth(data = sled1_sols, aes(Date, L_allometric, color = params), se=FALSE)

ggplot() +
  geom_smooth(data = sled1_sols_plus2, aes(Date, L_allometric, color = params), se=FALSE)

# sol_high_Sled1 %>% filter(Date==as.POSIXct("2018-04-17 12:00:00", tz="GMT")) %>% pull(W)
# sol_med_Sled1 %>% filter(Date==as.POSIXct("2018-04-17 12:00:00", tz="GMT")) %>% pull(W)
# sol_low_Sled1 %>% filter(Date==as.POSIXct("2018-04-17 12:00:00", tz="GMT")) %>% pull(W)
# sol_Sled1 %>% filter(Date==as.POSIXct("2018-04-17 12:00:00", tz="GMT")) %>% pull(W)
# 
sled1_sols %>% filter(Date==as.POSIXct("2018-04-17 12:00:00", tz="GMT")) %>% select(L_allometric, params) %>% arrange(L_allometric)
