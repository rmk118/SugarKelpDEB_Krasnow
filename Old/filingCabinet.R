#library(parameters)
#library(gdata)
#library(Metrics)
#library(gridExtra)

# Add column to data frame that identifies which group each cross is in (control temps)
mean_growth_by_cross_ctrl <- mean_growth_by_cross_ctrl %>%
  mutate(group = case_when(
    as.character(cross) %in% levels(high_stress$cross) ~ "high",
    as.character(cross) %in% levels(low_stress$cross) ~ "low",
    .default = "med"
  ), group = as.factor(group))

# Add column to data frame that identifies which group each cross is in (stress temps)
mean_growth_by_cross_stress <- mean_growth_by_cross_stress %>%
  mutate(group = case_when(
    as.character(cross) %in% levels(high_stress$cross) ~ "high",
    as.character(cross) %in% levels(low_stress$cross) ~ "low",
    .default = "med"
  ), group = as.factor(group))

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

# growth_no_zero <- growth_rates %>% 
#   filter(blade_len!=0) %>% 
#   group_by(id) %>% 
#   mutate(rel_growth_tot = (blade_len - lag(blade_len))/lag(blade_len), 
#          rel_growth_hp = if_else(is.na(len),NA,(len-lag(len))/lag(len)),
#          time_elapsed = as.integer(date-lag(date)),
#          growth_rate_tot = (blade_len-lag(blade_len))/time_elapsed,
#          growth_rate_hp = if_else(is.na(len),NA,(len-lag(len))/time_elapsed)) %>% 
#   filter(growth_rate_tot!=-1) %>% na.omit()

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

ggstatsplot::ggbetweenstats(unh_stress2,x=ctrl_group, y=std_rgr)

mean(growth_rates$blade_len) #in cm
mean(growth_rates$hp) #in cm
mean(growth_rates$net_perc_change) #as % inc.
mean(growth_rates$perc_change_hp) #as % inc.
mean(growth_rates$growth_rate_tot) #in cm/day
mean(growth_rates$growth_rate_hp) #in cm/day
mean(growth_rates$rgr) #in %/day

mean(growth_rates_no_flags$blade_len) #in cm
mean(growth_rates_no_flags$hp) #in cm
mean(growth_rates_no_flags$net_perc_change) #as % inc.
mean(growth_rates_no_flags$perc_change_hp) #as % inc.
mean(growth_rates_no_flags$growth_rate_tot) #in cm/day
mean(growth_rates_no_flags$growth_rate_hp) #in cm/day
mean(growth_rates_no_flags$rgr) #in %/day

### CONTROL PLOT MEANS ###
ggplot()+
  #geom_line(data=stress_groups, aes(x=temp, y=std_rate, color=level, linetype=type),linewidth=1.1)+
  geom_line(data=ctrl_groups, aes(x=temp, y=std_rate, color=level),linewidth=1.1)+
  theme_classic()+
  scale_x_continuous(breaks = seq(-5,35,5))+
  #scale_color_manual(values = c("#dd4124", "#edd746",'#0f85a0'))+
  labs(x="Temperature (°C)", y="Standardized rate", shape=NULL,color=NULL, linetype=NULL, fill=NULL)+
  guides(color = guide_legend(
    override.aes=list(shape = "-")))+
  theme(text = element_text(size=25, color="black"),
        plot.margin = margin(0.4,0.7,0.4,0.4, "cm"),
        axis.title.y = element_text(margin = margin(t = 0, r = 9, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 9, r = 0, b = 0, l = 0)),
        axis.title = element_text(size=18), axis.text = element_text(size=18, color="black"),
        #legend theming
        legend.margin = margin(-0.7,0, -0.6, 0, "lines"),
        legend.position = c(.85, .8),
        legend.text = element_text(size=9, face="bold"),
        legend.box.background = element_rect(color="black"),
        legend.box.margin = margin(0.1, 0.3, 0.8, 0.2, "lines"),
        legend.key.size = unit(0.9, "lines"))

#### Plot: means control points/curves ####
ctrl_params_plot_means<-ggplot()+
  geom_point(data=lit_data, aes(x=temp, y=std_rate, shape=paper), size=4)+
  geom_point(data=unh_ctrl %>% mutate(ctrl_group=factor(ctrl_group, levels=c("low", "med", "high"))), aes(x=temp, y=std_rgr,color=ctrl_group, shape="unh_ctrl"), size=3)+
  geom_line(data=ctrl_groups %>% filter(type=="means"), aes(x=temp, y=std_rate, color=group),linewidth=1.1)+
  theme_classic()+
  scale_x_continuous(breaks = seq(-5,35,5))+
  scale_shape_manual(values=c(15,18,8,16,17), 
                     labels=c("Bolton and Lüning (1982)",
                              "Fortes and Lüning (1980)",
                              "Davison and Davison (1987)",
                              "Davison (1987)",
                              "This study"))+
  scale_color_manual(values = c("#dd4124", "#edd746",'#0f85a0'))+
  labs(x="Temperature (°C)", y="Standardized rate", color=NULL, shape=NULL, fill=NULL)+
  guides(color = guide_legend( 
    override.aes=list(shape = "-")))+
  theme(text = element_text(size=15, color="black"),
        plot.margin = margin(0.4,0.7,0.4,0.4, "cm"),
        axis.title.y = element_text(margin = margin(t = 0, r = 9, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 9, r = 0, b = 0, l = 0)),
        axis.title = element_text(size=18), axis.text = element_text(size=18, color="black"))+
        #legend theming
        #legend.margin = margin(-0.7,0, -0.6, 0, "lines"),
        #legend.position = c(.85, .8),
        #legend.position = "none")+
        #legend.text = element_text(size=9, face="bold"),
        #legend.box.background = element_rect(color="black"),
        #legend.box.margin = margin(0.1, 0.3, 0.8, 0.2, "lines"),
        #legend.key.size = unit(0.9, "lines"))+
  ggtitle("Control - means")

(ctrl_params_plot_means+new_params_plot_means)/(control_all_plot+stress_all_plot)

nls_fun(df=lit_data_plus, type="stress", level="high")$std_rate==high_smooth
#IT WORKS


# #### Plot: means stress points/curves ####
# new_params_plot_means <-ggplot()+
#   geom_point(data=lit_data, aes(x=temp, y=std_rate, shape=paper), size=4)+
#   geom_point(data=unh_stress %>% mutate(stress_group=factor(stress_group, levels=c("low", "med", "high"))), aes(x=temp, y=std_rgr,color=stress_group, shape="unh_stress"), size=3)+
#   geom_line(data=stress_groups %>% filter(type=="means"), aes(x=temp, y=std_rate, color=group),linewidth=1.1)+
#   theme_classic()+
#   #coord_cartesian(xlim=c(-5, 35), ylim=c(0,2.5), expand=FALSE)+
#   scale_x_continuous(breaks = seq(-5,35,5))+
#   #scale_y_continuous(breaks = seq(0,2.5,0.5))+
#   scale_shape_manual(values=c(15,18,8,16,17), 
#                      labels=c("Bolton and Lüning (1982)",
#                               "Fortes and Lüning (1980)",
#                               "Davison and Davison (1987)",
#                               "Davison (1987)",
#                               "This study"))+
#   scale_color_manual(values = c("#dd4124", "#edd746",'#0f85a0'), labels=c("Low tolerance", "Medium tolerance", "High tolerance"))+
#   labs(x="Temperature (°C)", y="Standardized rate", color=NULL, shape=NULL, fill=NULL)+
#   guides(color = guide_legend( 
#     override.aes=list(shape = "-")))+
#   theme(text = element_text(size=15, color="black"),
#         plot.margin = margin(0.4,0.7,0.4,0.4, "cm"),
#         axis.title.y = element_text(margin = margin(t = 0, r = 9, b = 0, l = 0)),
#         axis.title.x = element_text(margin = margin(t = 9, r = 0, b = 0, l = 0)),
#         axis.title = element_text(size=18), axis.text = element_text(size=18, color="black"),
#         #legend theming
#         #legend.margin = margin(-0.7,0, -0.6, 0, "lines"),
#         #legend.position = c(.85, .8),
#         legend.position = "none")+
#        # legend.text = element_text(size=9, face="bold"),
#        # legend.box.background = element_rect(color="black"),
#        # legend.box.margin = margin(0.1, 0.3, 0.8, 0.2, "lines"),
#         #legend.key.size = unit(0.9, "lines"))+
#   ggtitle("Stress - means")
# new_params_plot_means
# 
# ggsave(
#   filename="./figures/threeCurves.png",
#   plot=new_params_plot_means, 
#   device="png",
#   width = 855, height = 750, units = "px",scale=2.6
# )
# 
# #### Consolidated refit parameters ###
# params<-data.frame(curve=c("Venolia", "Refit", "High stress means", "Med stress means", "Low stress means", "High stress all", "Med stress all", "Low stress all", "High control all", "Med control all", "Low control all"),
#            T_A=c(T_A, new_T_A, high_T_A, med_T_A, low_T_A, high_T_A_all, med_T_A_all, low_T_A_all, high_T_A_all_ctrl, med_T_A_all_ctrl,low_T_A_all_ctrl),
#            T_H=c(T_H, new_T_H, high_T_H, med_T_H, low_T_H, high_T_H_all, med_T_H_all, low_T_H_all, high_T_H_all_ctrl, med_T_H_all_ctrl,low_T_H_all_ctrl),
#            T_AH = c(T_AH, new_T_AH, high_T_AH, med_T_AH, low_T_AH, high_T_AH_all, med_T_AH_all, low_T_AH_all,high_T_AH_all_ctrl, med_T_AH_all_ctrl,low_T_AH_all_ctrl)) %>% 
#   mutate(T_L=T_L, T_AL=T_AL, T_0=T_0)
# 
# 
# #### Plot: all stress points/curves ####
# stress_all_plot<- ggplot()+
#   geom_line(data=stress_groups %>% filter(type=="all"), aes(x=temp, y=std_rate, color=level),linewidth=1.1)+
#   geom_line(aes(x=temp_smooth-273.15, y=y_smooth_venolia, color="Venolia"))+
#   theme_classic()+
#   scale_shape_manual(values=c(15,18,8,16,20), 
#                      labels=c("Bolton and Lüning (1982)",
#                               "Fortes and Lüning (1980)",
#                               "Davison and Davison (1987)",
#                               "Davison (1987)",
#                               "This study"))+
#   geom_point(data=lit_data, aes(x=temp, y=std_rate, shape=paper), size=4)+
#   geom_point(data=unh_stress_all %>% mutate(stress_group=factor(stress_group, levels=c("low", "med", "high"))), aes(x=temp, y=std_rgr,color=stress_group, shape="unh_stress"), size=3)+
#   scale_x_continuous(breaks = seq(-5,35,5))+
#   scale_color_manual(values = c("#dd4124", "#edd746",'#0f85a0', "black"), labels=c("Low tolerance", "Medium tolerance", "High tolerance", "Venolia"))+
#   labs(x="Temperature (°C)", y="Standardized rate", shape=NULL,color=NULL, linetype=NULL, fill=NULL)+
#   guides(color = guide_legend( 
#     override.aes=list(shape = "-")))+
#   theme(text = element_text(size=15, color="black"),
#         plot.margin = margin(0.4,0.7,0.4,0.4, "cm"),
#         axis.title.y = element_text(margin = margin(t = 0, r = 9, b = 0, l = 0)),
#         axis.title.x = element_text(margin = margin(t = 9, r = 0, b = 0, l = 0)),
#         axis.title = element_text(size=18), axis.text = element_text(size=18, color="black"),
#         #legend theming
#         legend.margin = margin(-0.7,0, -0.6, 0, "lines"),
#         #legend.position = c(.85, .8),
#         legend.position = "right",
#         legend.text = element_text(size=9, face="bold"),
#         legend.box.background = element_rect(color="black"),
#         legend.box.margin = margin(0.1, 0.3, 0.8, 0.2, "lines"),
#         legend.key.size = unit(0.9, "lines"))+ggtitle("Stress - all")
# 
# #### Plot: Control stress points/curves ####
# control_all_plot<- ggplot()+
#   geom_line(data=ctrl_groups %>% filter(type=="all"), aes(x=temp, y=std_rate, color=level),linewidth=1.1)+
#   theme_classic()+
#   scale_shape_manual(values=c(15,18,8,16,20), 
#                      labels=c("Bolton and Lüning (1982)",
#                               "Fortes and Lüning (1980)",
#                               "Davison and Davison (1987)",
#                               "Davison (1987)",
#                               "This study"))+
#   geom_point(data=lit_data, aes(x=temp, y=std_rate, shape=paper), size=4)+
#   geom_point(data=unh_ctrl_all %>% mutate(ctrl_group=factor(ctrl_group, levels=c("low", "med", "high"))), aes(x=temp, y=std_rgr,color=ctrl_group, shape="unh_ctrl"), size=3)+
#   scale_x_continuous(breaks = seq(-5,35,5))+
#   scale_color_manual(values = c("#dd4124", "#edd746",'#0f85a0'), labels=c("Low tolerance", "Medium tolerance", "High tolerance"))+
#   labs(x="Temperature (°C)", y="Standardized rate", shape=NULL,color=NULL, linetype=NULL, fill=NULL)+
#   guides(color = guide_legend( 
#     override.aes=list(shape = "-")))+
#   theme(text = element_text(size=15, color="black"),
#         plot.margin = margin(0.4,0.7,0.4,0.4, "cm"),
#         axis.title.y = element_text(margin = margin(t = 0, r = 9, b = 0, l = 0)),
#         axis.title.x = element_text(margin = margin(t = 9, r = 0, b = 0, l = 0)),
#         axis.title = element_text(size=18), axis.text = element_text(size=18, color="black"))+
#         #legend theming
#         #legend.margin = margin(-0.7,0, -0.6, 0, "lines"),
#         #legend.position = c(.85, .8),
#         #legend.position = "none")+
#         #legend.text = element_text(size=9, face="bold"),
#         #legend.box.background = element_rect(color="black"),
#         #legend.box.margin = margin(0.1, 0.3, 0.8, 0.2, "lines"),
#         #legend.key.size = unit(0.9, "lines"))+
#   ggtitle("Control - all")

#### old way of doing crosses ####
cross144 <- c(growth_rates_no_flags %>% filter(cross==144, trt!="control") %>% mutate(rgr= if_else(rgr<0 | is.na(rgr), 0, rgr)) %>% summarise(mean(rgr)))

ref_cross <- growth_rates_no_flags %>% group_by(cross, date, trt) %>%
  filter(growth_hp > 0) %>% 
  summarise(mean_growth = mean(growth_hp, na.rm=TRUE),
            mean_rgr = mean(growth_rate_hp, na.rm=TRUE),
            mean_rgr = mean(rgr, na.rm=TRUE),
            max_rgr = max(rgr, na.rm=TRUE)) %>% add_week_fun() %>% 
  left_join(weekly_means_degC) %>% 
  filter(!(week %in% c(1,5)))  %>% 
  group_by(cross) %>% 
  filter(round(mean_temp,0)==20) %>% 
  #summarize(ref_rgr = mean(max_rgr)) %>% 
  summarize(ref_rgr = max(max_rgr)) %>% 
  select(cross, ref_rgr)  %>%
  add_row(cross=as.factor(144), ref_rgr=cross144[[1]])

unh_cross<-growth_rates_no_flags %>% 
  group_by(cross, date, trt, mean_temp) %>%
  summarize(rgr = mean(rgr, na.rm=TRUE)) %>% 
  left_join(ref_cross) %>% 
  mutate(std_rate = rgr/ref_rgr) %>% 
  add_ctrl_fun() %>% add_stress_fun() %>% 
  mutate(ctrl_group = fct_relevel(ctrl_group, c("high", "med", "low")),stress_group = fct_relevel(stress_group, c("high", "med", "low"))) %>% 
  ungroup() %>% 
  mutate(paper="unh_cross2023", paper_full = paste(paper, stress_group, sep="_"), temp=mean_temp, rate=rgr, extra_info=paste(trt, cross,sep="_"), temp_K= temp + 273.15) %>%  #adding columns with identifying info
  select(paper, paper_full, temp, rate,std_rate, extra_info, temp_K, ctrl_group, stress_group)

# stress_cross_df <- map(.x=c(high="high", med="med", low="low"), .f=nls_fun, df=unh_cross, type="stress") %>% bind_rows(.id="level") %>% mutate(res="cross")

# ctrl_cross_df <- map(.x=c(high="high", med="med", low="low"), .f=nls_fun, df=unh_cross, type="ctrl") %>% bind_rows(.id="level") %>% mutate(res="cross")

#start_Wd <- 0.00387*33^(1.469) 
# M_V = W/(w_V+M_EN*w_EN+M_EC*w_EC)
# W = (w_V+M_EN*w_EN+M_EC*w_EC)* M_V
#(w_V+0.01*w_EN+0.1*w_EC)* 5

#state_Lo_UNH <- c(m_EC = 0.08, m_EN = 0.008, M_V = 1/(w_V+0.008*w_EN+0.08*w_EC))

# Create zoo objects
# zT <- zoo(temp_forcing, order.by = temp_forcing$date)[,2:4]
# zPAR <- zoo(PAR_forcing$PAR, PAR_forcing$date) # high freq
# zN <- zoo(nitrate$nitrate, as.POSIXct(nitrate$date))  # low freq
# 
# # Merge series into one object
# z <- merge(zT, zPAR,zN)
# # Interpolate calibration data
# z$zN <- na.approx(z$zN, rule=2)
# z$zPAR <- na.approx(z$zPAR, rule=2)
# # Only keep index values from sample data
# z <- z[index(zT),]
# Z <- as_tibble(z) %>% mutate(PAR=as.double(zPAR), nitrate=as.double(zN), .keep="unused", date=temp_forcing$date) %>% mutate(across(where(is.character), ~as.double(.x)))

View(growth_data_orig %>% select(cross, Female.GP, Male.GP) %>% distinct() %>% na.omit() %>% left_join(growth_rates_no_flags, by="cross") %>% select(cross, Female.GP, Male.GP, ctrl_group, stress_group) %>% distinct())

##### OLD HALVES #########

#testing halves instead of thirds

top_half_ctrl <- mean_growth_by_cross_ctrl %>% slice_max(gr_hp, n=15) %>% mutate(cross = fct_drop(cross))
bottom_half_ctrl <- mean_growth_by_cross_ctrl %>% slice_min(gr_hp, n=15) %>% mutate(cross = fct_drop(cross))

top_half_stress <- mean_growth_by_cross_stress %>% slice_max(gr_hp, n=15) %>% mutate(cross = fct_drop(cross))
bottom_half_stress <- mean_growth_by_cross_stress %>% slice_min(gr_hp, n=15) %>% mutate(cross = fct_drop(cross))

growth_rates_halves <- growth_rates_no_flags %>%
  mutate(ctrl_group = case_when(
    as.character(cross) %in% levels(top_half_ctrl$cross) ~ "top",
    as.character(cross) %in% levels(bottom_half_ctrl$cross) ~ "bottom"), 
    ctrl_group = as.factor(ctrl_group),
    stress_group = case_when(
      as.character(cross) %in% levels(top_half_stress$cross) ~ "top",
      as.character(cross) %in% levels(bottom_half_stress$cross) ~ "bottom"), 
    stress_group = as.factor(stress_group)) %>%
  ungroup()

growth_no_zero <- growth_data %>% 
  group_by(id) %>% 
  mutate(time_elapsed = as.integer(date-lag(date)), #days between sampling events
         rgr=((log(hp)-log(lag(hp)))/time_elapsed)*100) %>% #relative growth rate, in % per day
  mutate(rgr=if_else(rgr<0,0,rgr)) %>% na.omit() %>% 
  add_trt_fun() #%>% filter(trt!="lowN")

crosses<-growth_no_zero %>% group_by(cross, date, trt) %>%
  summarise(mean_rgr = mean(rgr, na.rm=TRUE),
            max_rgr = max(rgr, na.rm=TRUE)) %>% add_week_fun() %>% 
  left_join(weekly_means_degC %>% mutate(mean_temp=round(mean_temp,0))) %>% 
  group_by(cross) %>% 
  filter(mean_temp==20) %>% 
  summarize(ref_rgr = max(max_rgr)) %>% 
  select(cross, ref_rgr) %>% 
  filter(ref_rgr>0.2)

top_half_ctrl <- mean_growth_by_cross_ctrl %>% filter(cross %in% c(fct_drop(crosses$cross))) %>% 
  slice_max(gr_hp, n=14) %>% mutate(cross = fct_drop(cross))
bottom_half_ctrl <- mean_growth_by_cross_ctrl %>% filter(cross %in% c(fct_drop(crosses$cross))) %>% slice_min(gr_hp, n=14) %>% mutate(cross = fct_drop(cross))

top_half_stress <- mean_growth_by_cross_stress %>% filter(cross %in% c(fct_drop(crosses$cross))) %>% slice_max(gr_hp, n=14) %>% mutate(cross = fct_drop(cross))
bottom_half_stress <- mean_growth_by_cross_stress %>% filter(cross %in% c(fct_drop(crosses$cross))) %>% slice_min(gr_hp, n=14) %>% mutate(cross = fct_drop(cross))

crosses<- crosses %>%
  mutate(ctrl_group = case_when(
    as.character(cross) %in% levels(top_half_ctrl$cross) ~ "top",
    as.character(cross) %in% levels(bottom_half_ctrl$cross) ~ "bottom"), 
    ctrl_group = as.factor(ctrl_group),
    stress_group = case_when(
      as.character(cross) %in% levels(top_half_stress$cross) ~ "top",
      as.character(cross) %in% levels(bottom_half_stress$cross) ~ "bottom"), 
    stress_group = as.factor(stress_group)) %>%
  ungroup()


half_test_stress<-growth_no_zero %>%
  group_by(cross, trt, date) %>%
  summarise(mean_rgr = mean(rgr, na.rm=TRUE)) %>% 
  add_week_fun() %>% #add week labels
  filter(!(week %in% c(1,5))) %>% 
  left_join(crosses) %>% 
  left_join(weekly_means_degC) %>% 
  na.omit() %>% 
  mutate(std_rate=mean_rgr/ref_rgr) %>%
  mutate(ctrl_group = case_when(
    as.character(cross) %in% levels(top_half_ctrl$cross) ~ "top",
    as.character(cross) %in% levels(bottom_half_ctrl$cross) ~ "bottom"), 
    ctrl_group = as.factor(ctrl_group),
    stress_group = case_when(
      as.character(cross) %in% levels(top_half_stress$cross) ~ "top",
      as.character(cross) %in% levels(bottom_half_stress$cross) ~ "bottom"), 
    stress_group = as.factor(stress_group)) %>%
  ungroup() %>% 
  mutate(paper="unh_cross_stress2023", paper_full = paste(paper, stress_group, sep="_"), temp=mean_temp, rate=mean_rgr, extra_info=paste(trt, cross,sep="_"), temp_K= temp + 273.15) %>%  #adding columns with identifying info
  select(paper, paper_full, temp, rate,std_rate, extra_info, temp_K, ctrl_group, stress_group)

half_test_ctrl<-growth_no_zero %>%
  group_by(cross, trt, date) %>%
  summarise(mean_rgr = mean(rgr, na.rm=TRUE)) %>% 
  add_week_fun() %>% #add week labels
  filter(!(week %in% c(1,5))) %>% 
  left_join(crosses) %>% 
  left_join(weekly_means_degC) %>% 
  na.omit() %>% 
  mutate(std_rate=mean_rgr/ref_rgr) %>%
  mutate(ctrl_group = case_when(
    as.character(cross) %in% levels(top_half_ctrl$cross) ~ "top",
    as.character(cross) %in% levels(bottom_half_ctrl$cross) ~ "bottom"), 
    ctrl_group = as.factor(ctrl_group),
    stress_group = case_when(
      as.character(cross) %in% levels(top_half_stress$cross) ~ "top",
      as.character(cross) %in% levels(bottom_half_stress$cross) ~ "bottom"), 
    stress_group = as.factor(stress_group)) %>%
  ungroup() %>% 
  mutate(paper="unh_cross_ctrl2023", paper_full = paste(paper, ctrl_group, sep="_"), temp=mean_temp, rate=mean_rgr, extra_info=paste(trt, cross,sep="_"), temp_K= temp + 273.15) %>%  #adding columns with identifying info
  select(paper, paper_full, temp, rate,std_rate, extra_info, temp_K, ctrl_group, stress_group)




# calculate weekly mean and maximum relative growth rates for each half, where grouping was done based on growth under non-control temperatures
weekly_halves_stress <- growth_rates_halves %>% 
  group_by(stress_group, trt, date, mean_temp) %>%
  summarise(mean_growth = mean(growth_hp, na.rm=TRUE),
            mean_rgr = mean(growth_rate_hp, na.rm=TRUE),
            mean_rgr = mean(rgr, na.rm=TRUE),
            max_rgr = max(rgr, na.rm=TRUE)) %>% 
  add_week_fun() %>% #add week labels
  filter(!(week %in% c(1,5))) #remove weeks without growth rates

# calculate weekly mean and maximum relative growth rates for each of the 3 groups, where grouping was done based on growth under control temperatures
weekly_halves_ctrl <- growth_rates_halves %>% 
  group_by(ctrl_group, trt, date, mean_temp) %>% 
  summarise(mean_growth = mean(growth_hp, na.rm=TRUE),
            mean_rgr = mean(growth_rate_hp, na.rm=TRUE),
            mean_rgr = mean(rgr, na.rm=TRUE),
            max_rgr = max(rgr, na.rm=TRUE)) %>% add_week_fun() %>%
  filter(!(week %in% c(1,5)))


ref_stress_half <- weekly_halves_stress %>% 
  group_by(stress_group) %>% 
  filter(round(mean_temp,1)==20) %>% 
  mutate(ref_rgr = max_rgr) %>% 
  select(stress_group, ref_rgr)

ref_ctrl_half <- weekly_halves_ctrl %>% 
  group_by(ctrl_group) %>% 
  filter(round(mean_temp,1)==20) %>% 
  mutate(ref_rgr = max_rgr) %>% 
  select(ctrl_group, ref_rgr)

# Finding standardized RGR at the group level
unh_stress_half <- weekly_halves_stress %>% 
  left_join(ref_stress_half) %>%  #add a column to the weekly means (stress grouping) with reference RGRs
  mutate(across(c(mean_rgr, max_rgr),
                ~if_else(.x<0 | is.na(.x), 0, .x))) %>% # Replace negative or missing RGRs with 0
  mutate(std_rgr_max = max_rgr/ref_rgr, # standardized max RGR (max RGR divided by ref temp RGR)
         std_rgr_mean = mean_rgr/ref_rgr, # standardized mean RGR (mean RGR divided by ref temp RGR)
         std_rgr = case_when(week==2 & trt!="control"~ std_rgr_max, #use max for week 2 non-control
                             .default= std_rgr_mean)) %>% 
  ungroup() %>%
  mutate(paper="unh_stress_all2023", paper_full=paste(paper, stress_group, sep="_"), temp=mean_temp, rate=mean_rgr, extra_info=NA, temp_K=temp+273.15) %>% #adding columns with identifying info
  select(paper, paper_full, temp, rate, std_rgr, extra_info, temp_K, stress_group)

unh_ctrl_half <- weekly_halves_ctrl %>% 
  left_join(ref_ctrl_half) %>%  #add a column to the weekly means (ctrl grouping) with reference RGRs
  mutate(across(c(mean_rgr, max_rgr),
                ~if_else(.x<0 | is.na(.x), 0, .x))) %>% # Replace negative or missing RGRs with 0
  mutate(std_rgr_max = max_rgr/ref_rgr, # standardized max RGR (max RGR divided by ref temp RGR)
         std_rgr_mean = mean_rgr/ref_rgr, # standardized mean RGR (mean RGR divided by ref temp RGR)
         std_rgr = case_when(week==2 & trt!="control"~ std_rgr_max, #use max for week 2 non-control
                             .default= std_rgr_mean)) %>% 
  ungroup() %>%
  mutate(paper="unh_ctrl_all2023", paper_full=paste(paper, ctrl_group, sep="_"), temp=mean_temp, rate=mean_rgr, extra_info=NA, temp_K=temp+273.15) %>% #adding columns with identifying info
  select(paper, paper_full, temp, rate, std_rgr, extra_info, temp_K, ctrl_group)


##### Cross level ####################################################################

###### Stress #####
unh_cross_stress_half <-growth_rates_halves %>% 
  group_by(cross, date, trt, mean_temp) %>% #group by cross, week, and treatment
  summarize(rgr = mean(rgr, na.rm=TRUE)) %>%  #find mean RGR for each cross under each trt in each week
  mutate(stress_group = case_when(
    as.character(cross) %in% levels(top_half_stress$cross) ~ "top",
    as.character(cross) %in% levels(bottom_half_stress$cross) ~ "bottom"), 
    stress_group = as.factor(stress_group)) %>% # add stress grouping labels
  left_join(ref_stress_half) %>% # add reference RGRs for stress groupings
  mutate(std_rate = rgr/ref_rgr) %>% # find standardized rate by dividing RGR at non-20°C temps by the reference rate, which is the stress group RGR at 20°C
  mutate(std_rate = if_else(std_rate<0 | is.na(std_rate), 0, std_rate)) %>% #replace negative rates with 0
  mutate(ctrl_group = case_when(
    as.character(cross) %in% levels(top_half_ctrl$cross) ~ "top",
    as.character(cross) %in% levels(bottom_half_ctrl$cross) ~ "bottom"), 
    ctrl_group = as.factor(ctrl_group)) %>%  #add control grouping labels
  ungroup() %>% 
  mutate(paper="unh_cross_stress2023", paper_full = paste(paper, stress_group, sep="_"), temp=mean_temp, rate=rgr, extra_info=paste(trt, cross,sep="_"), temp_K= temp + 273.15) %>%  #adding columns with identifying info
  select(paper, paper_full, temp, rate,std_rate, extra_info, temp_K, ctrl_group, stress_group)

###### Control #####
unh_cross_ctrl_half <-growth_rates_halves %>% 
  group_by(cross, date, trt, mean_temp) %>% #group by cross, week, and treatment
  summarize(rgr = mean(rgr, na.rm=TRUE)) %>%  #find mean RGR for each cross under each trt in each week
  mutate(ctrl_group = case_when(
    as.character(cross) %in% levels(top_half_ctrl$cross) ~ "top",
    as.character(cross) %in% levels(bottom_half_ctrl$cross) ~ "bottom"), 
    ctrl_group = as.factor(ctrl_group)) %>% # add ctrl grouping labels
  left_join(ref_ctrl_half) %>% # add reference RGRs for ctrl groupings
  mutate(std_rate = rgr/ref_rgr) %>% # find standardized rate by dividing RGR at non-20°C temps by the reference rate, which is the ctrl group RGR at 20°C
  mutate(std_rate = if_else(std_rate<0 | is.na(std_rate), 0, std_rate)) %>% #replace negative rates with 0
  mutate(stress_group = case_when(
    as.character(cross) %in% levels(top_half_stress$cross) ~ "top",
    as.character(cross) %in% levels(bottom_half_stress$cross) ~ "bottom"), 
    stress_group = as.factor(stress_group)) %>%  #add stress grouping labels
  ungroup() %>% 
  mutate(paper="unh_cross_ctrl2023", paper_full = paste(paper, ctrl_group, sep="_"), temp=mean_temp, rate=rgr, extra_info=paste(trt, cross,sep="_"), temp_K= temp + 273.15) %>%  #adding columns with identifying info
  select(paper, paper_full, temp, rate,std_rate, extra_info, temp_K, ctrl_group, stress_group)


#### Combine new and lit data #### 
lit_data_plus_half <- bind_rows(lit_data, 
                                unh_stress_half %>% mutate(std_rate=std_rgr,ctrl_group = NA,.keep="unused")) %>% 
  mutate(type="stress", res="means")

lit_data_plus_ctrl_half <- bind_rows(lit_data,
                                     unh_ctrl_half %>% mutate(std_rate=std_rgr, stress_group = NA,.keep="unused")) %>%
  mutate(type="ctrl", res="means")

lit_data_plus_cross_half <- bind_rows(lit_data,
                                      #                                unh_cross_stress_half) %>%  mutate(type="stress", res="cross")
                                      half_test_stress) %>%  mutate(type="stress", res="cross")

lit_data_plus_cross_ctrl_half <- bind_rows(lit_data,
                                           #                                 unh_cross_ctrl_half) %>% mutate(type="ctrl", res="cross")
                                           half_test_ctrl) %>%  mutate(type="ctrl", res="cross")


all_lit_data_half <- rbind(lit_data_plus_half, lit_data_plus_cross_half, lit_data_plus_ctrl_half, lit_data_plus_cross_ctrl_half) %>% distinct() %>% # removes duplicate rows (i.e., the original literature data)
  mutate(level = case_when(str_ends(paper_full, "top") ~ "top",
                           str_ends(paper_full, "bottom") ~ "bottom",
                           .default = "lit"),
         level = as_factor(level),
         level = fct_relevel(level, c("top", "bottom", "lit")))

ggplot()+geom_point(data=all_lit_data_half, aes(x=temp, y=std_rate, color=level))+facet_grid(type~res)


#### Stress means
stress_means_half_df <- map(.x=c(top="top", bottom="bottom"), .f=nls_fun, df=lit_data_plus_half, type="stress") %>% bind_rows(.id="level") %>% mutate(res="means")

#### Stress crosses
stress_cross_half_df <- map(.x=c(top="top", bottom="bottom"), .f=nls_fun, df=lit_data_plus_cross_half, type="stress") %>% bind_rows(.id="level") %>% mutate(res="cross")

##### Control means
ctrl_means_half_df <- map(.x=c(top="top", bottom="bottom"), .f=nls_fun, df=lit_data_plus_ctrl_half, type="ctrl") %>% bind_rows(.id="level") %>% mutate(res="means")

#### Control crosses
ctrl_cross_half_df <- map(.x=c(top="top", bottom="bottom"), .f=nls_fun, df=lit_data_plus_cross_ctrl_half, type="ctrl") %>% bind_rows(.id="level") %>% mutate(res="cross")

halves_calibrations <- rbind(stress_means_half_df, ctrl_means_half_df, stress_cross_half_df, ctrl_cross_half_df) %>% 
  mutate(level = as_factor(level),
         level = fct_expand(level, "lit"),
         level = fct_relevel(level, c("top", "bottom", "lit")))

params_halves<-halves_calibrations %>% select(T_A, T_H, T_AH,type,level,res) %>% distinct() %>% 
  add_row(T_A=T_A, T_H=T_H,T_AH=T_AH,type="orig", level="orig", res="orig")%>% 
  add_row(T_A=new_T_A, T_H=new_T_H,T_AH=new_T_AH,type="lit", level="lit", res="lit")

ggplot()+
  geom_line(data=halves_calibrations %>% ungroup(), aes(x=temp_smooth-273.15, y=std_rate, color=level), linewidth=2)+
  theme_bw()+
  facet_grid(type~res)+
  annotate(geom='line', x=temp_smooth-273.15,y=y_smooth_venolia)+
  labs(x="Temperature (°C)", y="Standardized rate", color=NULL)+
  scale_color_manual(values=c("top"="#dd4124", "bottom"='#0f85a0'),
                     breaks=c("top", "bottom"),
                     labels=c("top"="Top half", "bottom"="Bottom half"))


##### END OLD HALVES #########

new_lit %>%  left_join(new_lit20)%>% 
  mutate(std_rate2 = round(rate/rate20,2)) %>% mutate(check=std_rate!=std_rate2) %>% pull(check) %>% sum()

new_lit20 <- new_lit %>% 
  filter(temp==20) %>% 
  select(paper_full, rate) %>% 
  rename(rate20=rate) 

new_lit <-new_lit%>%  left_join(new_lit20)%>% 
  mutate(std_rate2 = round(rate/rate20,2)) %>% mutate(check=std_rate==std_rate2) #%>% pull(check) %>% sum()

new_lit <- new_lit %>%  left_join(new_lit20)%>% 
  mutate(std_rate = round(rate/rate20,2))


ggplot(data=all_rmse %>% ungroup() %>% filter(year==2), aes(x=reorder_within(params, rmse, source), y=rmse, fill=params)) +
  geom_col()+
  geom_text(aes(label = round(rmse,1), vjust = -0.2))+
  facet_wrap(~source, scales="free_x")+
  scale_x_reordered()




# Sig tests ---------------------------------------------------------------

all_rmse %>% 
  filter(level!="orig") %>% 
  ungroup() %>% 
  mutate(imp = if_else(improvement>0, TRUE, FALSE)) %>% 
  group_by(level, year) %>% 
  summarize(num_imp = sum(imp), 
            perc_imp = sum(imp)/length(imp), 
            #mean_imp = mean(if_else(improvement>0, improvement, 0))) %>% left_join(
            med_imp = median(improvement),
            mean_imp = mean(improvement)) %>% left_join(
              (all_rmse %>% filter(level!="orig") %>% group_by(level, year) %>% 
                 wilcox_test(improvement ~ 0, alternative="greater") %>%
                 adjust_pvalue() %>% 
                 select(level, year, p.adj))) %>%
  arrange(year) #%>% write_clip()

all_rmse %>% 
  filter(level!="orig") %>% 
  ungroup() %>% 
  mutate(imp = if_else(improvement>0, TRUE, FALSE)) %>% 
  group_by(level, year) %>% 
  summarize(num_imp = sum(imp), 
            perc_imp = sum(imp)/length(imp),
            med_imp = median(improvement),
            mean_imp = mean(improvement)) %>% 
  left_join((all_rmse %>% filter(level!="orig") %>% group_by(level, year) %>% 
               wilcox_test(improvement ~ 0))) %>%
  #adjust_pvalue() %>% 
  #select(level, year, p.adj))) %>%
  arrange(year)

all_rmse %>% 
  select(source, year, level, rmse) %>% 
  ungroup() %>% 
  levene_test(rmse~level)

all_rmse %>% 
  select(source, year, level, rmse) %>% 
  ungroup() %>% 
  group_by(level) %>% 
  shapiro_test(rmse)

all_rmse %>% 
  #filter(level!="orig") %>% 
  select(source, year, level, rmse) %>% 
  #group_by(year) %>% 
  group_by(year) %>% 
  kruskal_test(rmse~level)