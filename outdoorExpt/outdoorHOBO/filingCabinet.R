#library(parameters)

##### temperature #####

# weeks_3to5_hourly <- weeks_3to5 %>% 
#   group_by(trt) %>% 
#   mutate(hourly_mean_temp = rollmean(temp,k=12, fill=NA)) %>% na.omit()

# weeks_3to5_daily <- weeks_3to5 %>% 
#   group_by(trt) %>% 
#   mutate(daily_mean_temp = rollmean(temp,k=24, fill=NA)) %>% na.omit()
# 
# daily_no_roll <- weeks_3to5 %>% 
#   mutate(trunc_date = round_date(datetime, "day")) %>% 
#   group_by(trt, trunc_date) %>% 
#   summarize(daily_no_roll = mean(temp, na.rm=TRUE))

# hourly_roll_mean<- ggplot(data=weeks_3to5)+
#   geom_line(aes(x=datetime, y=hourly_mean_temp, color=trt))
# 
# daily_roll_mean<-ggplot(data=weeks_3to5_daily)+
#   geom_line(aes(x=datetime, y=daily_mean_temp, color=trt))

#daily_no_roll_mean<-ggplot(data=daily_no_roll)+
#  geom_line(aes(x=trunc_date, y=daily_no_roll, color=trt))

#hourly_roll_mean/daily_roll_mean/hourly_no_roll_mean

weeks_3to5 <- dat_full %>% 
  filter(datetime %within% interval(ymd("2023-06-26"), ymd_hms(("2023-07-14 13:00:00"))))

hourly_no_roll <- weeks_3to5 %>% 
  mutate(trunc_date = round_date(datetime, "hour")) %>% 
  group_by(trt, trunc_date) %>% 
  summarize(hourly_no_roll = mean(temp, na.rm=TRUE))

hourly_no_roll_mean<-ggplot(data=hourly_no_roll)+
  geom_line(aes(x=trunc_date, y=hourly_no_roll, color=trt))+
  scale_x_datetime(breaks="4 days", date_labels = "%b %d")+
  labs(x=NULL, y="Temperature (°C)", color="Treatment")+
  theme_light()+
  theme(text = element_text(size=14),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title = element_text(size=13), axis.text = element_text(size=14))+
  scale_color_hue(labels = c("Control + N", "Heat + N", "Heat"))
hourly_no_roll_mean

dat_full %>% filter(temp<20.5 & temp>19.5, trt!="control") %>% pull(datetime) %>% round_date("day") %>% as.character() %>% fct_count() %>% print(n=29)

weeks_1to2 <- dat_full %>% 
  filter(datetime %within% interval(ymd("2023-06-12"), ymd("2023-06-25")))
weeks_1to2 %>% filter(trt!="control") %>% summarise(mean(temp))

dat_full %>% 
  filter(datetime %within% interval(ymd_hm("2023-06-12 14:00"), ymd_hm("2023-06-30 14:00"))) %>% filter(trt=="lowN") %>% summarise(mean(temp)) #20.0

dat_full %>% 
  filter(datetime %within% interval(ymd_hm("2023-06-12 14:00"), ymd_hm("2023-06-30 14:00"))) %>% filter(trt!="control") %>% summarise(mean(temp)) #19.9

#week3_yII
dat_full %>% 
  filter(datetime %within% interval(ymd_hm("2023-06-29 14:30"), ymd_hm("2023-06-30 14:30"))) %>% filter(trt=="highN") %>% summarise(mean(temp))

#week3_full A
dat_full %>% 
  filter(datetime %within% interval(ymd("2023-06-20"), ymd("2023-06-28"))) %>% filter(trt!="control") %>% summarise(mean(temp)) #20.1

dat_full %>% 
  filter(datetime %within% interval(ymd("2023-06-19"), ymd("2023-06_30"))) %>% 
  filter(!(datetime %within% interval(ymd_hm("2023-06-20 09:05"), ymd_hm("2023-06-20 09:25")))) %>% 
  filter(trt!="control") %>% summarise(mean(temp)) #20.2

dat_full %>% 
  filter(datetime %within% interval(ymd("2023-06-19"), ymd("2023-06-30"))) %>% 
  filter(trt=="highN") %>% summarise(mean(temp)) #20.2

#week2_yII
dat_full %>% 
  filter(datetime %within% interval(ymd("2023-06-20"), ymd("2023-06-21"))) %>% 
  filter(trt!="control") %>% summarise(mean(temp)) #19.0

#week1_yII
dat_full %>% 
  filter(datetime %within% interval(ymd("2023-06-12"), ymd("2023-06-13"))) %>% 
  filter(trt!="control") %>% summarise(mean(temp)) #18.8

dat_full %>% 
  filter(datetime %within% interval(ymd("2023-06-14"), ymd("2023-06-25"))) %>% 
  filter(trt!="control") %>% summarise(mean(temp)) #19.3

##### growth #####

# View(net_growth %>% filter(cross==40) %>% group_by(cross, trt, rep) %>% mutate(growth_rate = (blade_len - lag(blade_len))/lag(blade_len)) %>% filter(growth_rate!=-1) %>% na.omit())

growth_data %>% summarise(num = n_distinct(id))
growth_data %>% filter(blade_len==0) %>% group_by(id) %>% summarise(num = n_distinct(date)) %>% filter(num>1)

growth_rates_zero %>% ungroup() %>% group_by(date) %>% summarise(tot=mean(growth_rate_tot),
                                                                 hp = mean(growth_rate_hp))

mean_growth_by_cross <- growth_rates_zero %>% ungroup() %>% group_by(cross) %>% 
  summarise(g_tot=mean(growth_rate_tot)) %>% arrange(-g_tot)

high_performing <- mean_growth_by_cross %>% slice_max(g_tot, n=10)
low_performing <- mean_growth_by_cross %>% slice_min(g_tot, n=10)
mid_performing <- anti_join(mean_growth_by_cross, high_performing) %>% anti_join(low_performing)

model <- lm(growth_rate_tot ~ date, data = growth_rates_zero)
model_parameters(model)

ggplot(data=growth_rates %>% filter(blade_len!=0))+
  geom_point(aes(x=blade_len, y=growth_rate_tot))+
  geom_smooth(aes(x=blade_len, y=growth_rate_tot), method="gam")

growth_no_zero <- growth_rates %>% 
  filter(blade_len!=0) %>% 
  group_by(id) %>% 
  mutate(rel_growth_tot = (blade_len - lag(blade_len))/lag(blade_len), 
         rel_growth_hp = if_else(is.na(len),NA,(len-lag(len))/lag(len)),
         time_elapsed = as.integer(date-lag(date)),
         growth_rate_tot = (blade_len-lag(blade_len))/time_elapsed,
         growth_rate_hp = if_else(is.na(len),NA,(len-lag(len))/time_elapsed)) %>% 
  filter(growth_rate_tot!=-1) %>% na.omit()

ggplot(data=growth_no_zero, aes(x=blade_len, y=growth_rate_tot, label=id, color=as.factor(date)))+
  geom_point()+
  geom_text()
#geom_smooth()

ggplot(data=growth_no_zero, aes(x=blade_len, y=rel_growth_hp, color=trt))+
  geom_point()+
  geom_smooth(method = "lm", se=FALSE)

ggplot(data=growth_no_zero, aes(x=len, y=rel_growth_hp, color=trt))+
  geom_point()+
  geom_smooth(method = "lm", se=FALSE)

ggplot(data=growth_no_zero, aes(x=trt, y=growth_rate_hp))+
  geom_boxplot()

ggplot(data=growth_no_zero, aes(x=trt, y=len))+facet_wrap(~date)+
  geom_boxplot()

growth_no_zero %>% group_by(trt) %>% summarise(num=n())

#growth_rates_no_flags %>% ungroup() %>% 
#summarize(mean(diff), sd(diff))
# summarize(mean(growth_hp), sd(growth_hp))
# ggplot()+geom_density(aes(x=growth_hp, fill=group), alpha=0.6)+facet_wrap(~trt)

#library(kSamples)
# ad.test(growth_hp ~ trt, data=growth_rates_no_flags)
# ad.test(growth_hp ~ cross, data=growth_rates_no_flags)
# ad.test(growth_hp ~ group, data=growth_rates_no_flags)
# ad.test(growth_hp ~ group, data=growth_rates_no_flags %>% filter(trt=="A"))
# ad.test(gr_hp ~ group, data=mean_growth_by_cross)

##### ggstatsplot #####
#ggbetweenstats(data = growth_rates_no_flags, x=group, y=growth_hp)

##### NLS #####

#temp correction arrhenius function - basic shape
#C_T <- exp((T_A/T_0)-(T_A/T_field(t))) * 
#(1+exp((T_AL/T_0)-(T_AL/T_L))+ exp((T_AH/T_H)-(T_AH/T_0))) * 
#((1+exp((T_AL/T_field(t))-(T_AL/T_L))+exp((T_AH/T_H)-(T_AH/T_field(t))))^-1)

# y_smooth_minpack1=exp((7629/293.15)-(7629/temp_smooth)) *
#   (1+exp((196979/293.15)-(196979/272)) + exp((20035/288)-(20035/293.15))) * 
#   ((1+exp((196979/temp_smooth)-(196979/272))+exp((20035/288)-(20035/temp_smooth)))^-1)
# ggplot()+geom_line(aes(x=temp_smooth, y=y_smooth_minpack1))

nls(std_rate ~ exp((T_a/T_0)-(T_a/temp_K)) *
      (1+exp((T_al/T_0)-(T_al/T_l)) + exp((T_ah/T_h)-(T_ah/T_0))) *
      ((1+exp((T_al/temp_K)-(T_al/T_l))+exp((T_ah/T_h)-(T_ah/temp_K)))^-1),
    start=list(T_a=T_a, T_h=T_h, T_ah=T_ah, T_al=T_al, T_l=T_l),
    data = lit_data)

library(nlsr)

try(nlxb(formula=std_rate ~ exp((T_a/T_0)-(T_a/temp_K)) *
           (1+exp((T_al/T_0)-(T_al/T_l)) + exp((T_ah/T_h)-(T_ah/T_0))) *
           ((1+exp((T_al/temp_K)-(T_al/T_l))+exp((T_ah/T_h)-(T_ah/temp_K)))^-1),
         start=list(T_a=T_a, T_h=T_h, T_ah=T_ah, T_al=T_al),
         data = lit_data))

try(nlxb(formula=std_rate~exp((T_a/T_0)-(T_a/temp_K)) *
           (1+exp((T_AL/T_0)-(T_AL/T_L)) + exp((T_ah/T_h)-(T_ah/T_0))) *
           ((1+exp((T_AL/temp)-(T_AL/T_L))+exp((T_ah/T_h)-(T_ah/temp_K)))^-1),
         start=list(T_a=T_a, T_h=T_h, T_ah=T_ah),
         data = lit_data))


try(nlsLM(formula=std_rate ~ exp((T_a/T_0)-(T_a/temp_K)) *
            (1+exp((T_al/T_0)-(T_al/T_l)) + exp((T_ah/T_h)-(T_ah/T_0))) *
            ((1+exp((T_al/temp_K)-(T_al/T_l))+exp((T_ah/T_h)-(T_ah/temp_K)))^-1),
          start=list(T_a=T_a, T_h=T_h, T_ah=T_ah, T_al=T_al),
          data = lit_data))