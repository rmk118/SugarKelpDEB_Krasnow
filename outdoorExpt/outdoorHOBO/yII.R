#### PAM data setup ####
pam_data <- read.csv("~/Downloads/MBL_SES/outdoorExpt/PAM.csv")
pam_data <- pam_data %>% filter(Time!="175 was measured twice") %>% 
  filter((!is.na(YII..meristem.))) %>% 
  select(!starts_with("X"))

time_fixed <- pam_data  %>% 
  mutate(date = mdy(Date),
         time_period = hm(Time, quiet = TRUE), 
         time_duration = as.duration(time_period), 
         .keep="unused", .before = Sampling.week) %>%
  mutate(heat_test_f = na_if(Heat.tested.Female.GP,"#N/A"),
         heat_test_m = na_if(Heat.tested.Male.GP,"#N/A"), .keep="unused") %>% 
  arrange(orig.order) %>% 
  mutate(num = row_number(), .before = 1)

#average time 14:31 (~2:30pm) for all sampling days    
avg_time_all <- as.duration(mean(time_fixed$time_duration, na.rm=TRUE))
seconds_to_period(avg_time_all)

#no times were recorded for 6/30/23 or 7/07/23, so we'll assume 14:31
#for the other dates, we'll use linear interpolation to fill in the measurement times

time_filled <- time_fixed %>% 
  mutate(time_duration = if_else(date == "2023-06-30" | date == "2023-07-07", 
                                 avg_time_all, time_duration))

time_filled <- time_filled %>% 
  mutate(time_interp = as.duration(approx(x=time_duration, xout=num)[["y"]]), .before=2)

#Check that we didn't overwrite any points with times originally recorded
any(time_filled$time_interp!=time_filled$time_duration, na.rm = TRUE)

time_done <- time_filled %>% 
  mutate(time = seconds_to_period(time_interp),
         meristem = YII..meristem., 
         cross=Cross,
         trt= case_when(
           Treatment=="A" ~ "control",
           Treatment=="B" ~ "highN",
           Treatment=="C" ~ "lowN"),
         above = YII..above.tipward.of.10cm.holepunch.,
         .keep="unused", .before=2) %>% 
  select(-c(time_duration, time_period))

#### PAM data summary ####
time_done <- time_done %>%
  mutate(ctrl_group = case_when(
    as.character(cross) %in% levels(high_ctrl$cross) ~ "high",
    as.character(cross) %in% levels(low_ctrl$cross) ~ "low",
    .default = "med"),
    stress_group = case_when(
      as.character(cross) %in% levels(high_stress$cross) ~ "high",
      as.character(cross) %in% levels(low_stress$cross) ~ "low",
      .default = "med")) %>% 
  mutate(ctrl_group = as_factor(ctrl_group),
         stress_group = as_factor(stress_group)) %>% 
  mutate(ctrl_group = fct_relevel(ctrl_group, c("high", "med", "low")),
         stress_group = fct_relevel(stress_group, c("high", "med", "low"))) %>% 
  ungroup() %>% 
  mutate(trt=fct_drop(trt)) %>% 
  mutate(week = case_when(
    date %within% week1 ~ 1,
    date %within% week2 ~ 2,
    date %within% week3 ~ 3,
    date %within% week4 ~ 4,
    date %within% week5 ~ 5
  ))

ggplot(data = time_done %>% group_by(trt, date) %>% 
         summarise(avg = mean(meristem)) %>% ungroup())+ 
  labs(y="YII meristem", x=NULL)+
  geom_line(aes(x=date, y=avg, color=trt))



weekly_PAM <- time_done %>% group_by(stress_group, week, trt) %>% 
  summarise(mean_yII_meristem = mean(meristem, na.rm=TRUE),
            mean_yII_blade = mean(above, na.rm=TRUE)) %>% ungroup()


weekly_means_2on<-full_join(weekly_means_growth, weekly_PAM, 
                            by=c("stress_group", "week", "trt"))  %>% 
  mutate( stress_group = fct_relevel(stress_group, c("high", "med", "low"))) %>% 
  filter(week != 1) %>% ungroup()

#growth rates under stress not driven by yII
(ggplot(data = weekly_means_2on)+
    geom_boxplot(aes(y=mean_yII_meristem, x=stress_group, fill=stress_group))+facet_wrap(~trt)) /
  (ggplot(data = weekly_means_2on)+
     geom_boxplot(aes(y=mean_yII_blade, x=stress_group, fill=stress_group))+
     facet_wrap(~trt))

#ggbetweenstats(data=weekly_means_2on, x=stress_group, y=mean_yII_meristem)
# ggbetweenstats(data=time_done %>% ungroup() %>% filter(week!=1), x=stress_group, y=above)
# ggbetweenstats(data=time_done %>% ungroup() %>% filter(week!=1& week !=2), x=stress_group, y=meristem)
# ggbetweenstats(data=time_done %>% ungroup()%>% filter(week!=1& week !=2), x=ctrl_group, y=meristem)
# ggbetweenstats(data=time_done %>% ungroup()%>% filter(week!=1), x=ctrl_group, y=above)
# ggbetweenstats(data=time_done %>% ungroup()%>% filter(week!=1), x=ctrl_group, y=meristem)
