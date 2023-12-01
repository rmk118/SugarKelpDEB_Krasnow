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

