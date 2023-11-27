#Analysis of outdoor tank sporophyte heat stress growth expt
#Ruby Krasnow
#Last updated: Nov 27, 2023

library(tidyverse)
library(patchwork)
library(lubridate)
library(clipr)
library(tseries)
library(zoo)
library(ggstatsplot)
library(minpack.lm)

source("./outdoorExpt/outdoorHOBO/outdoor_HOBO.R")

#### Nitrate data ####
nitrate_data<- read.csv("~/Downloads/MBL_SES/outdoorExpt/nitrate.csv")

nitrate_data <- nitrate_data %>% mutate(date = mdy(Date),
                                        tank = Tank,
                                        ppm = Level..ppm.,
                                        umol = Level..umol., .keep="none") %>% na.omit()

nitrate <- nitrate_data %>% 
  group_by(date) %>% 
  summarise(nitrate = mean(umol)) %>%
  mutate(tank = "Ambient")

nitrate_enriched <- nitrate_data %>% 
  group_by(date) %>% 
  summarise(nitrate = mean(umol)+7) %>%
  mutate(tank = "Enriched")

nitrate <- rbind(nitrate, nitrate_enriched)
rm(nitrate_enriched, nitrate_data)

ggplot(data=nitrate)+
  geom_line(aes(x=date, y=nitrate, color=tank))+
  labs(x=NULL, y="Nitrate (µmol)", color="Treatment")+
  theme_light()

#### Growth data ####
growth_data <- read.csv("~/Downloads/MBL_SES/outdoorExpt/growth.csv") %>% 
  mutate(date = mdy(Date),
         heat_test_f = na_if(Heat.testesd.Female.GP.ID,"#N/A"),
         heat_test_m = na_if(Heat.tested.Male.GP.ID,"#N/A"),
         blade_len = Total.length.of.blade,
         hp = as.double(na_if(Distance.base.of.blade.to.10cm.hole.punch,"--")),
         cross = as.factor(Cross),
         trt = as.factor(Treatment),
         rep = as.factor(Rep),
         id = paste(cross, trt, rep, sep="."),
         .before = 3, .keep="unused")

growth_data <- growth_data %>% 
  select(c(date, blade_len,hp, cross, trt, rep, id, comments)) %>% #select relevant columns
  complete(date, id) %>% #make implicit NAs explicit
  filter(hp>9) #the hole punch was done 10cm from the base, this gets rid of a few instances of potential measuring error or degradation while allowing slight room for error
  

growth_rates <- growth_data %>% 
  group_by(id) %>% 
  mutate(
    net_growth = (blade_len - lag(blade_len)), #difference in total blade length between sampling dates
    growth_hp = if_else(is.na(hp),NA,(hp-lag(hp))), #difference in distance between hole punch and stipe between sampling dates (provides estimate of growth without effect of erosion from the end of the blade)
    net_perc_change = (blade_len - lag(blade_len))/lag(blade_len), #percent change in total blade length between sampling dates
    perc_change_hp = if_else(is.na(hp),NA,(hp-lag(hp))/lag(hp)), #percent change in stipe to hole punch distance between sampling dates
    time_elapsed = as.integer(date-lag(date)), #days between sampling events
    growth_rate_tot = (blade_len-lag(blade_len))/time_elapsed, #growth rate of overall blade, in cm/day
    growth_rate_hp = if_else(is.na(hp),NA,(hp-lag(hp))/time_elapsed), #growth rate near meristem, in cm/day
     rgr=((log(hp)-log(lag(hp)))/time_elapsed)*100) %>% #relative growth rate, in % per day
    na.omit() # %>% filter(growth_rate_tot!=-1)


growth_rates %>% mutate(diff = net_growth-growth_hp) %>% 
  ungroup() %>% summarize(mean(diff, na.rm=TRUE), sd(diff, na.rm=TRUE))

growth_rates <- growth_rates %>% mutate(diff = net_growth-growth_hp,
                                               flag_diff = (diff > 2),
                                        big_loss = (growth_hp < -2))

mean(growth_rates$blade_len) #in cm
mean(growth_rates$hp) #in cm
mean(growth_rates$net_perc_change) #as % inc.
mean(growth_rates$perc_change_hp) #as % inc.
mean(growth_rates$growth_rate_tot) #in cm/day
mean(growth_rates$growth_rate_hp) #in cm/day
mean(growth_rates$rgr) #in %/day

sum(growth_rates$flag_diff)
sum(growth_rates$big_loss)

growth_rates_no_flags <- growth_rates %>% filter(big_loss==FALSE, flag_diff==FALSE) %>% 
  mutate(trt = fct_drop(trt), cross=fct_drop(cross))

mean_growth_by_cross_ctrl <- growth_rates_no_flags %>% 
  filter(trt=="A") %>%
  ungroup() %>% group_by(cross) %>% 
  summarise(gr_hp=mean(growth_rate_hp)) %>% 
  arrange(-gr_hp) %>% 
  mutate(cross = fct_drop(cross))

mean_growth_by_cross_stress <- growth_rates_no_flags %>% 
  filter(trt!="A") %>%
  ungroup() %>% group_by(cross) %>% 
  summarise(gr_hp=mean(perc_change_hp)) %>% 
  arrange(-gr_hp) %>% 
  mutate(cross = fct_drop(cross))

high_ctrl <- mean_growth_by_cross_ctrl %>% slice_max(gr_hp, n=10) %>% mutate(cross = fct_drop(cross))
low_ctrl <- mean_growth_by_cross_ctrl %>% slice_min(gr_hp, n=10) %>% mutate(cross = fct_drop(cross))
mid_ctrl <- anti_join(mean_growth_by_cross_ctrl, high_ctrl) %>% anti_join(low_ctrl) %>% mutate(cross = fct_drop(cross))

high_stress <- mean_growth_by_cross_stress %>% slice_max(gr_hp, n=10) %>% mutate(cross = fct_drop(cross))
low_stress <- mean_growth_by_cross_stress %>% slice_min(gr_hp, n=10) %>% mutate(cross = fct_drop(cross))
mid_stress <- anti_join(mean_growth_by_cross_stress, high_stress) %>% anti_join(low_stress) %>% mutate(cross = fct_drop(cross))

mean_growth_by_cross_ctrl <- mean_growth_by_cross_ctrl %>%
  mutate(group = case_when(
    as.character(cross) %in% levels(high_ctrl$cross) ~ "high",
    as.character(cross) %in% levels(low_ctrl$cross) ~ "low",
    .default = "med"
  ), group = as.factor(group))

mean_growth_by_cross_stress <- mean_growth_by_cross_stress %>%
  mutate(group = case_when(
    as.character(cross) %in% levels(high_stress$cross) ~ "high",
    as.character(cross) %in% levels(low_stress$cross) ~ "low",
    .default = "med"
  ), group = as.factor(group))


growth_rates_no_flags <- growth_rates_no_flags %>%
  mutate(ctrl_group = case_when(
    as.character(cross) %in% levels(high_ctrl$cross) ~ "high",
    as.character(cross) %in% levels(low_ctrl$cross) ~ "low",
    .default = "med"),
    stress_group = case_when(
      as.character(cross) %in% levels(high_stress$cross) ~ "high",
      as.character(cross) %in% levels(low_stress$cross) ~ "low",
      .default = "med") ) %>%
  mutate(ctrl_group = as_factor(ctrl_group),
         stress_group = as_factor(stress_group)) %>% 
  mutate(ctrl_group = fct_relevel(ctrl_group, c("high", "med", "low")),stress_group = fct_relevel(stress_group, c("high", "med", "low"))) %>% 
ungroup()

sum(growth_rates_no_flags$ctrl_group==growth_rates_no_flags$stress_group)/length(growth_rates_no_flags$stress_group)

ggplot(data=growth_rates_no_flags)+
  geom_boxplot(aes(y=growth_hp, x=stress_group, fill=stress_group))+
  facet_wrap(~trt)+
  scale_fill_viridis_d()+
  theme_classic()

# library(ggstatsplot)
# ggbetweenstats(data = growth_rates_no_flags %>% filter(trt=="control") %>% ungroup(), x=stress_group, y=growth_rate_hp)

weekly_growth <- growth_rates_no_flags %>% group_by(stress_group, date, trt) %>% 
  summarise(mean_growth = mean(growth_hp, na.rm=TRUE),
            mean_growth_rate = mean(growth_rate_hp, na.rm=TRUE),
            max_growth_rate = max(growth_rate_hp, na.rm=TRUE)) %>% ungroup() %>% 
  mutate(week = case_when(
    date=="2023-06-28" ~ 2,
    date=="2023-07-05" ~ 3,
    date=="2023-07-13" ~ 4), 
    trt = case_when(
      trt=="A" ~ "control",
      trt=="B" ~ "highN",
      trt=="C" ~ "lowN"), .keep="unused")

growth_rates_no_flags <- growth_rates_no_flags %>% 
  mutate(week = case_when(
    date=="2023-06-28" ~ 2,
    date=="2023-07-05" ~ 3,
    date=="2023-07-13" ~ 4), 
    trt = case_when(
      trt=="A" ~ "control",
      trt=="B" ~ "highN",
      trt=="C" ~ "lowN"), .keep="unused") %>% 
  left_join(weekly_means)

weekly_means_growth<-full_join(weekly_means_all, weekly_growth) %>% 
  filter(!(week %in% c(1,5)))

ggplot(data=growth_rates_no_flags %>% na.omit())+
  geom_point(aes(x=mean, y=growth_rate_hp,color=stress_group, fill=stress_group))+
  geom_smooth(aes(x=mean, y=growth_rate_hp,color=stress_group, fill=stress_group), se=FALSE)

#### Literature data ####
lit_data <- read.csv("~/Downloads/MBL_SES/arrhenius_lit_data.csv") %>% 
  select(c(paper, paper_full, temp, rate, std_rate, extra_info)) %>% 
  filter(temp!=23) %>% 
  mutate(extra_info = if_else(paper!="bl182" & extra_info=="Germany", NA, extra_info),
         paper_full=as_factor(paper_full)) %>% 
 mutate(paper_full = fct_relevel(paper_full, c("bl1982_France", "bl1982_Norway", "bl1982_Germany", "bl1982_UK", "fl1980","dd1987",  "d1987_0", "d1987_5", "d1987_10", "d1987_15", "d1987_20")),
        temp_K = temp+273.15)

# Venolia et al., 2020 values
#Arrhenius temperature
T_A <- 6314.3 # K
#Upper boundary of temperature tolerance
T_H <- 13.386 + 273.15 # K, 286.54
#Lower boundary of temperature tolerance
T_L <- 273.15 # K
#Arrhenius temperature outside T_H
T_AH <- 18702 #K
#Arrhenius temperature outside T_L
T_AL <- 4391.9 #K
#temperature at which rate parameters are given
T_0 <- 20 + 273.15 # K, 293.15

# Initial values for least-squares fitting
#Arrhenius temperature
T_a <- 6000 # K
#Upper boundary of temperature tolerance
T_h <- 286 # K
#Arrhenius temperature outside T_H
T_ah <- 18000 #K

lit_data_fit <- lit_data %>% mutate(y=exp((T_A/T_0)-(T_A/temp_K)) *
            (1+exp((T_AL/T_0)-(T_AL/T_L)) + exp((T_AH/T_H)-(T_AH/T_0))) * 
            ((1+exp((T_AL/temp_K)-(T_AL/T_L))+exp((T_AH/T_H)-(T_AH/temp_K)))^-1))

temp_smooth <- seq(-5,35,by=0.1)+273.15

y_smooth_venolia=exp((T_A/T_0)-(T_A/temp_smooth)) *
  (1+exp((T_AL/T_0)-(T_AL/T_L)) + exp((T_AH/T_H)-(T_AH/T_0))) * 
  ((1+exp((T_AL/temp_smooth)-(T_AL/T_L))+exp((T_AH/T_H)-(T_AH/temp_smooth)))^-1)

#### Venolia plot #####
venolia_plot<- ggplot()+
  geom_point(data=lit_data, aes(x=temp, y=std_rate, group=paper_full,shape=paper_full, color=paper_full), size=4)+
  geom_line(aes(x=temp_smooth-273.15, y=y_smooth_venolia, linetype="Arrhenius function"))+
  theme_classic()+
  coord_cartesian(xlim=c(-5, 35), ylim=c(0,2.5), expand=FALSE)+
  scale_x_continuous(breaks = seq(-5,35,5))+
  scale_y_continuous(breaks = seq(0,2.5,0.5))+
  scale_shape_manual(values=c(15,15,15,15,18,8,16,16,16,16,16), 
                     labels=c("Bolton and Lüning (1982), Growth, France", 
                              "Bolton and Lüning (1982), Growth, Norway",
                              "Bolton and Lüning (1982), Growth, Germany",
                              "Bolton and Lüning (1982), Growth, UK",
                              "Fortes and Lüning (1980), Growth",
                              "Davison and Davison (1987), Growth",
                              "Davison (1987), Photosynthesis (grown at 0°C)",
                              "Davison (1987), Photosynthesis (grown at 5°C)",
                              "Davison (1987), Photosynthesis (grown at 10°C)",
                              "Davison (1987), Photosynthesis (grown at 15°C)",
                              "Davison (1987), Photosynthesis (grown at 20°C)"))+
  scale_color_manual(values=c('orange',"yellow", "maroon3", "green3","deepskyblue", "red", "blue","orange", "yellow","maroon3","green3"),
                     labels=c("Bolton and Lüning (1982), Growth, France", 
                              "Bolton and Lüning (1982), Growth, Norway",
                              "Bolton and Lüning (1982), Growth, Germany",
                              "Bolton and Lüning (1982), Growth, UK",
                              "Fortes and Lüning (1980), Growth",
                              "Davison and Davison (1987), Growth",
                              "Davison (1987), Photosynthesis (grown at 0°C)",
                              "Davison (1987), Photosynthesis (grown at 5°C)",
                              "Davison (1987), Photosynthesis (grown at 10°C)",
                              "Davison (1987), Photosynthesis (grown at 15°C)",
                              "Davison (1987), Photosynthesis (grown at 20°C)"))+
  labs(x="Temperature (°C)", y="Standardized rate", linetype=NULL,color=NULL, shape=NULL)+
  theme(text = element_text(size=25, color="black"),
        plot.margin = margin(0.4,0.7,0.4,0.4, "cm"),
        axis.title.y = element_text(margin = margin(t = 0, r = 9, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 9, r = 0, b = 0, l = 0)),
        axis.title = element_text(size=18), axis.text = element_text(size=18, color="black"),
        #legend theming
        legend.margin = margin(-0.7,0, -0.7, 0, "lines"),
        legend.position = c(.8, .8),
        legend.text = element_text(size=8, face="bold"),
        legend.box.background = element_rect(color="black"),
        legend.box.margin = margin(0.1, 0.3, 0.8, 0.2, "lines"),
        legend.key.size = unit(0.9, "lines"))

ggsave(
  filename="./figures/venolia.png",
  plot=venolia_plot, 
  device="png",
  width = 855, height = 750, units = "px",scale=2.6
)

#### NLS to literature data only #####
params_orig <- nlsLM(formula=std_rate ~ exp((T_a/T_0)-(T_a/temp_K)) *
            (1+exp((T_AL/T_0)-(T_AL/T_L)) + exp((T_ah/T_h)-(T_ah/T_0))) *
            ((1+exp((T_AL/temp_K)-(T_AL/T_L))+exp((T_ah/T_h)-(T_ah/temp_K)))^-1),
          start=list(T_a=T_a, T_h=T_h, T_ah=T_ah),
          data = lit_data)
summary(params_orig)

new_T_A <- params_orig$m$getPars()[[1]]
new_T_H <- params_orig$m$getPars()[[2]]
new_T_AH <- params_orig$m$getPars()[[3]]

new_smooth <- exp((new_T_A/T_0)-(new_T_A/temp_smooth))*(1+exp((T_AL/T_0)-(T_AL/T_L))+exp((new_T_AH/new_T_H)-(new_T_AH/T_0))) * ((1+exp((T_AL/temp_smooth)-(T_AL/T_L))+exp((new_T_AH/new_T_H)-(new_T_AH/temp_smooth)))^-1)

t <- weekly_means_growth %>% group_by(stress_group) %>% 
  filter(round(mean_temp,1)==20) %>% mutate(ref_rate = max_growth_rate) %>% 
  select(stress_group, ref_rate)

unh <- weekly_means_growth %>% left_join(t) %>% 
  mutate(across(c(mean_growth_rate, max_growth_rate),
                ~if_else(.x<0 | is.na(.x), 0, .x))) %>%  
  mutate(std_rate_max = max_growth_rate/ref_rate,
         std_rate_mean = mean_growth_rate/ref_rate,
         std_rate = case_when(week==2 & trt!="control"~ std_rate_max,
                              .default= std_rate_mean)) %>% 
  ungroup() %>% 
  mutate(paper="unh2023", paper_full=paste(paper, stress_group, sep="_"), temp=mean_temp, rate=mean_growth_rate, extra_info=trt, temp_K=temp+273.15) %>% 
  filter(trt!="control") %>% 
  select(paper, paper_full, temp, rate, std_rate, extra_info, temp_K, stress_group)

lit_data_plus <- bind_rows(lit_data %>% mutate(stress_group="lit", std_rate), unh)

ggplot()+
  geom_point(data=lit_data_plus, aes(x=temp, y=std_rate, color=stress_group))

ggstatsplot::ggbetweenstats(unh,x=stress_group, y=std_rate)

#### NLS for lit data + most heat-tolerant strains  #####
params_high <- nlsLM(formula=std_rate ~ exp((T_a/T_0)-(T_a/temp_K)) *
                       (1+exp((T_AL/T_0)-(T_AL/T_L)) + exp((T_ah/T_h)-(T_ah/T_0))) *
                       ((1+exp((T_AL/temp_K)-(T_AL/T_L))+exp((T_ah/T_h)-(T_ah/temp_K)))^-1),
                     start=list(T_a=T_a, T_h=T_h, T_ah=T_ah),
                     data = lit_data_plus %>% filter(stress_group=="lit" | stress_group=="high"))
summary(params_high)

high_T_A <- params_high$m$getPars()[[1]]
high_T_H <- params_high$m$getPars()[[2]]
high_T_AH <- params_high$m$getPars()[[3]]

high_smooth <- exp((high_T_A/T_0)-(high_T_A/temp_smooth))*(1+exp((T_AL/T_0)-(T_AL/T_L))+exp((high_T_AH/high_T_H)-(high_T_AH/T_0))) * ((1+exp((T_AL/temp_smooth)-(T_AL/T_L))+exp((high_T_AH/high_T_H)-(high_T_AH/temp_smooth)))^-1)

#### NLS for lit data + medium heat-tolerant strains  #####
params_med <- nlsLM(formula=std_rate ~ exp((T_a/T_0)-(T_a/temp_K)) *
                       (1+exp((T_AL/T_0)-(T_AL/T_L)) + exp((T_ah/T_h)-(T_ah/T_0))) *
                       ((1+exp((T_AL/temp_K)-(T_AL/T_L))+exp((T_ah/T_h)-(T_ah/temp_K)))^-1),
                     start=list(T_a=T_a, T_h=T_h, T_ah=T_ah),
                     data = lit_data_plus %>% 
                      filter(stress_group=="lit" | stress_group=="med"))
summary(params_med)

med_T_A <- params_med$m$getPars()[[1]]
med_T_H <- params_med$m$getPars()[[2]]
med_T_AH <- params_med$m$getPars()[[3]]

med_smooth <- exp((med_T_A/T_0)-(med_T_A/temp_smooth))*(1+exp((T_AL/T_0)-(T_AL/T_L))+exp((med_T_AH/med_T_H)-(med_T_AH/T_0))) * ((1+exp((T_AL/temp_smooth)-(T_AL/T_L))+exp((med_T_AH/med_T_H)-(med_T_AH/temp_smooth)))^-1)

#### NLS for lit data + low heat-tolerant strains  #####
params_low <- nlsLM(formula=std_rate ~ exp((T_a/T_0)-(T_a/temp_K)) *
                      (1+exp((T_AL/T_0)-(T_AL/T_L)) + exp((T_ah/T_h)-(T_ah/T_0))) *
                      ((1+exp((T_AL/temp_K)-(T_AL/T_L))+exp((T_ah/T_h)-(T_ah/temp_K)))^-1),
                    start=list(T_a=T_a, T_h=T_h, T_ah=T_ah),
                    data = lit_data_plus %>% 
                      filter(stress_group=="lit" | stress_group=="low"))
summary(params_low)

low_T_A <- params_low$m$getPars()[[1]]
low_T_H <- params_low$m$getPars()[[2]]
low_T_AH <- params_low$m$getPars()[[3]]

low_smooth <- exp((low_T_A/T_0)-(low_T_A/temp_smooth))*(1+exp((T_AL/T_0)-(T_AL/T_L))+exp((low_T_AH/low_T_H)-(low_T_AH/T_0))) * ((1+exp((T_AL/temp_smooth)-(T_AL/T_L))+exp((low_T_AH/low_T_H)-(low_T_AH/temp_smooth)))^-1)

stress_groups <- data.frame(temp = temp_smooth-273.15, high=high_smooth, med=med_smooth, low=low_smooth) %>% pivot_longer(cols=c(2:4), names_to = "group", values_to = "std_rate") %>% mutate(group=factor(group, levels=c("low", "med", "high")))

ggplot()+
  #geom_line(aes(x=temp_smooth, y=y_smooth_venolia, color="venolia"))+
  #geom_line(aes(x=temp_smooth, y=new_smooth, color="new"))+
  #geom_line(aes(x=lit_data_fit$temp_K, y=lit_data_fit$y, color="discrete"))+
  geom_line(aes(x=temp_smooth, y=high_smooth, color="high"))+
  geom_line(aes(x=temp_smooth, y=med_smooth, color="med"))+
  geom_line(aes(x=temp_smooth, y=low_smooth, color="low"))


#### Plot w/ 3 curves ####
new_params_plot<-ggplot()+
  geom_point(data=lit_data, aes(x=temp, y=std_rate, shape=paper), size=4)+
  geom_point(data=unh %>% mutate(stress_group=factor(stress_group, levels=c("low", "med", "high"))), aes(x=temp, y=std_rate,color=stress_group, shape="unh"), size=4)+
  geom_line(data=stress_groups, aes(x=temp, y=std_rate, color=group),linewidth=1.1)+
  theme_classic()+
  coord_cartesian(xlim=c(-5, 35), ylim=c(0,2.5), expand=FALSE)+
  scale_x_continuous(breaks = seq(-5,35,5))+
  scale_y_continuous(breaks = seq(0,2.5,0.5))+
  scale_shape_manual(values=c(15,18,8,16,17), 
                     labels=c("Bolton and Lüning (1982)",
                              "Fortes and Lüning (1980)",
                              "Davison and Davison (1987)",
                              "Davison (1987)",
                              "This study"))+
  scale_color_manual(values = c("#dd4124", "#edd746",'#0f85a0'), labels=c("Low tolerance", "Medium tolerance", "High tolerance"))+
  labs(x="Temperature (°C)", y="Standardized rate", color=NULL, shape=NULL, fill=NULL)+
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

ggsave(
  filename="./figures/threeCurves.png",
  plot=new_params_plot, 
  device="png",
  width = 855, height = 750, units = "px",scale=2.6
)


params<-data.frame(curve=c("Venolia", "Refit", "High", "Med", "Low"),
           T_A=c(T_A, new_T_A, high_T_A, med_T_A, low_T_A),
           T_H=c(T_H, new_T_H, high_T_H, med_T_H, low_T_H),
           T_AH = c(T_AH, new_T_AH, high_T_AH, med_T_AH, low_T_AH)) %>% 
  mutate(T_L=T_L, T_AL=T_AL, T_0=T_0)
