#Analysis of outdoor tank sporophyte heat stress growth expt
#Ruby Krasnow
#Last updated: Nov 28, 2023

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

# Here we check the average and SD of the difference between total (entire blade) and hole punch growth measurements.
growth_rates %>% mutate(diff = net_growth-growth_hp) %>% 
  ungroup() %>% summarize(mean(diff, na.rm=TRUE), sd(diff, na.rm=TRUE))

#Adding new columns to dataframe
growth_rates <- growth_rates %>% mutate(diff = net_growth-growth_hp,
                                        flag_diff = (diff > 2),
                                        big_loss = (growth_hp < -2))

sum(growth_rates$flag_diff)
sum(growth_rates$big_loss)

#We remove instances where the entire blade growth was substantially more than the hole punch growth or the hole punch length shrunk substantially between weeks
growth_rates_no_flags <- growth_rates %>% filter(big_loss==FALSE, flag_diff==FALSE) %>% 
  mutate(trt = fct_drop(trt), cross=fct_drop(cross))

# Note: overall, we are using the three treatments as a way to get growth rates under different temperatures. There were no substantial differences between the high N and low N tanks, so we are considering treatment A to be "control" conditions and everything else to be "stress" conditions

# This summarizes the mean growth of each cross under control temps (~9-12°C)
mean_growth_by_cross_ctrl <- growth_rates_no_flags %>% 
  filter(trt=="A") %>%
  ungroup() %>% group_by(cross) %>% 
  summarise(gr_hp=mean(rgr)) %>% 
  arrange(-gr_hp) %>% 
  mutate(cross = fct_drop(cross))

# This summarizes the mean growth of each cross under stress/ambient temps, which increased throughout the summer
mean_growth_by_cross_stress <- growth_rates_no_flags %>% 
  filter(trt!="A") %>%
  ungroup() %>% group_by(cross) %>% 
  summarise(gr_hp=mean(rgr)) %>% 
  arrange(-gr_hp) %>% 
  mutate(cross = fct_drop(cross))

##### Grouping strains by growth ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #### 

# Create vectors splitting the 30 crosses into 3 groups of 10 based on their relative growth rates under the control temps
high_ctrl <- mean_growth_by_cross_ctrl %>% slice_max(gr_hp, n=10) %>% mutate(cross = fct_drop(cross))
low_ctrl <- mean_growth_by_cross_ctrl %>% slice_min(gr_hp, n=10) %>% mutate(cross = fct_drop(cross))
mid_ctrl <- anti_join(mean_growth_by_cross_ctrl, high_ctrl) %>% anti_join(low_ctrl) %>% mutate(cross = fct_drop(cross))

# Create vectors splitting the 30 crosses into 3 groups of 10 based on their relative growth rates under the stress temps
high_stress <- mean_growth_by_cross_stress %>% slice_max(gr_hp, n=10) %>% mutate(cross = fct_drop(cross))
low_stress <- mean_growth_by_cross_stress %>% slice_min(gr_hp, n=10) %>% mutate(cross = fct_drop(cross))
mid_stress <- anti_join(mean_growth_by_cross_stress, high_stress) %>% anti_join(low_stress) %>% mutate(cross = fct_drop(cross))

# Add column to data frame that identifies which group each cross is in (control temps)
mean_growth_by_cross_ctrl <- mean_growth_by_cross_ctrl %>%
  mutate(group = case_when(
    as.character(cross) %in% levels(high_ctrl$cross) ~ "high",
    as.character(cross) %in% levels(low_ctrl$cross) ~ "low",
    .default = "med"
  ), group = as.factor(group))

# Add column to data frame that identifies which group each cross is in (stress temps)
mean_growth_by_cross_stress <- mean_growth_by_cross_stress %>%
  mutate(group = case_when(
    as.character(cross) %in% levels(high_stress$cross) ~ "high",
    as.character(cross) %in% levels(low_stress$cross) ~ "low",
    .default = "med"
  ), group = as.factor(group))

# Add columns to data frame that identifies which group each cross is in (both control and stress temps)
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

sum(growth_rates_no_flags$ctrl_group==growth_rates_no_flags$stress_group)/length(growth_rates_no_flags$stress_group) #what percentage of crosses are in the same group when we consider growth under control temps vs growth under stress?

ggplot(data=growth_rates_no_flags)+
  geom_boxplot(aes(y=rgr, x=stress_group, fill=stress_group))+
  facet_wrap(~trt)+
  scale_fill_viridis_d()+
  theme_classic()

# Function to add week and treatment to dataframes
add_week_trt_fun <- function (df){
  df_out <- df %>% ungroup() %>% mutate(week = case_when(
    date=="2023-06-28" ~ 2, #add week column
    date=="2023-07-05" ~ 3,
    date=="2023-07-13" ~ 4), 
    trt = case_when( #add treatment column
      trt=="A" ~ "control",
      trt=="B" ~ "highN",
      trt=="C" ~ "lowN"), .keep="unused")
  df_out
}

weekly_growth_stress <- growth_rates_no_flags %>% group_by(stress_group, date, trt) %>% 
  summarise(mean_growth = mean(growth_hp, na.rm=TRUE),
            mean_rgr = mean(growth_rate_hp, na.rm=TRUE),
            mean_rgr = mean(rgr, na.rm=TRUE),
            max_rgr = max(rgr, na.rm=TRUE)) %>% add_week_trt_fun()

weekly_growth_ctrl <- growth_rates_no_flags %>% group_by(ctrl_group, date, trt) %>% 
  summarise(mean_growth = mean(growth_hp, na.rm=TRUE),
            mean_rgr = mean(growth_rate_hp, na.rm=TRUE),
            mean_rgr = mean(rgr, na.rm=TRUE),
            max_rgr = max(rgr, na.rm=TRUE)) %>% add_week_trt_fun()

growth_rates_no_flags <- growth_rates_no_flags %>% add_week_trt_fun() %>% 
  left_join(weekly_means) #join with temperature data

weekly_means_growth_stress<-full_join(weekly_means_all, weekly_growth_stress) %>% 
  filter(!(week %in% c(1,5)))

weekly_means_growth_ctrl<-full_join(weekly_means, weekly_growth_ctrl) %>% 
  filter(!(week %in% c(1,5)))

ggplot(data=growth_rates_no_flags %>% na.omit())+
  geom_point(aes(x=mean_temp, y=rgr,color=stress_group, fill=stress_group))+
  geom_smooth(aes(x=mean_temp, y=rgr,color=stress_group, fill=stress_group), se=FALSE)

#### Literature data ####
lit_data <- read.csv("~/Downloads/MBL_SES/arrhenius_lit_data.csv") %>% 
  select(c(paper, paper_full, temp, rate, std_rate, extra_info)) %>% 
  filter(temp!=23) %>% 
  mutate(extra_info = if_else(paper!="bl182" & extra_info=="Germany", NA, extra_info),
         paper_full=as_factor(paper_full)) %>% 
 mutate(paper_full = fct_relevel(paper_full, c("bl1982_France", "bl1982_Norway", "bl1982_Germany", "bl1982_UK", "fl1980","dd1987",  "d1987_0", "d1987_5", "d1987_10", "d1987_15", "d1987_20")),
        temp_K = temp+273.15,
        stress_group="lit",
        ctrl_group="lit") 

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

new_T_A <- params_orig$m$getPars()[[1]]
new_T_H <- params_orig$m$getPars()[[2]]
new_T_AH <- params_orig$m$getPars()[[3]]

new_smooth <- exp((new_T_A/T_0)-(new_T_A/temp_smooth))*(1+exp((T_AL/T_0)-(T_AL/T_L))+exp((new_T_AH/new_T_H)-(new_T_AH/T_0))) * ((1+exp((T_AL/temp_smooth)-(T_AL/T_L))+exp((new_T_AH/new_T_H)-(new_T_AH/temp_smooth)))^-1)


#### Set-up for new data NLS #####
# 20°C was the average temp in the low N experimental tank during week 2, so we use the RGR of each group after week 2 as the growth rate at the Arrhenius reference temperature. However, the experimental kelp blades were still acclimating to their new environments and were generally stressed during the first couple weeks of the experiment— this is also reflected in the PAM measurements. Many of the growth rates were negative or very low during this period, so we're using the max RGR instead of mean RGR when considering week 2 data.

# find the RGR at the reference temperature (20°C) when we split the 30 strains into groups of 3 based on their growth rates in the HEAT STRESS treatments
t_stress <- weekly_means_growth_stress %>% 
  group_by(stress_group) %>% 
  filter(round(mean_temp,1)==20) %>% 
  mutate(ref_rgr = max_rgr) %>% 
  select(stress_group, ref_rgr)

# find the RGR at the reference temperature (20°C) for the three groups when we split the 30 strains into groups of 3 based on their growth rates in the CONTROL TEMPERATURE treatment
t_ctrl <- weekly_means_growth_ctrl %>% ungroup() %>% 
  group_by(ctrl_group) %>% 
  filter(round(mean_temp,1)==20) %>% 
  mutate(ref_rgr = max_rgr) %>% 
  select(ctrl_group, ref_rgr)

#### Finding RGR at the grouping level
unh_stress <- weekly_means_growth_stress %>% 
  left_join(t_stress) %>%  #add a column to the weekly means (stress grouping) with reference RGRs
  mutate(across(c(mean_rgr, max_rgr),
                ~if_else(.x<0 | is.na(.x), 0, .x))) %>% # Replace any negative or missing RGRs with 0
  mutate(std_rgr_max = max_rgr/ref_rgr, # standardized max RGR (max RGR divided by ref temp RGR)
         std_rgr_mean = mean_rgr/ref_rgr, # standardized mean RGR (mean RGR divided by ref temp RGR)
         std_rgr = case_when(week==2 & trt!="control"~ std_rgr_max, #use max for week 2 non-control
                              .default= std_rgr_mean)) %>% ungroup() %>%
  mutate(paper="unh_stress_all023", paper_full=paste(paper, stress_group, sep="_"), temp=mean_temp, rate=mean_rgr, extra_info=trt, temp_K=temp+273.15) %>% #adding columns with identifying info
  select(paper, paper_full, temp, rate, std_rgr, extra_info, temp_K, stress_group)

unh_ctrl <- weekly_means_growth_ctrl %>%
  left_join(t_ctrl) %>% #add a column to the weekly means (stress grouping) with reference RGRs
  mutate(across(c(mean_rgr, max_rgr),
                ~if_else(.x<0 | is.na(.x), 0, .x))) %>% # Replace any negative or missing RGRs with 0
  mutate(std_rgr_max = max_rgr/ref_rgr, # standardized max RGR (max RGR divided by ref temp RGR)
         std_rgr = mean_rgr/ref_rgr) %>% # standardized mean RGR (mean RGR divided by ref temp RGR)
  ungroup() %>%
  mutate(paper="unh_ctrl_all023", paper_full=paste(paper, ctrl_group, sep="_"), temp=mean_temp, rate=mean_rgr, extra_info=trt, temp_K=temp+273.15) %>% #adding columns with identifying info
  select(paper, paper_full, temp, rate, std_rgr, extra_info, temp_K, ctrl_group)

unh_stress_all <- growth_rates_no_flags %>% 
    left_join(t_stress) %>% 
    mutate(rgr= if_else(rgr<0 | is.na(rgr), 0, rgr),
           std_rgr = rgr/ref_rgr) %>% 
    ungroup() %>%
    mutate(paper="unh_stress_all023", paper_full = paste(paper, stress_group, sep="_"), temp=mean_temp, rate=rgr, extra_info=trt, temp_K= temp + 273.15) %>% #adding columns with identifying info
    select(paper, paper_full, temp, rate, std_rgr, extra_info, temp_K, ctrl_group, stress_group)
  
unh_ctrl_all <- growth_rates_no_flags %>% 
    left_join(t_ctrl) %>% 
    mutate(rgr= if_else(rgr<0 | is.na(rgr), 0, rgr),
           std_rgr = rgr/ref_rgr) %>%
    ungroup() %>%
    mutate(paper="unh_ctrl_all023", paper_full=paste(paper, ctrl_group, sep="_"), temp=mean_temp, rate=rgr, extra_info=trt, temp_K=temp+273.15) %>% #adding columns with identifying info
    select(paper, paper_full, temp, rate, std_rgr, extra_info, temp_K, ctrl_group,stress_group)


lit_data_plus <- bind_rows(lit_data, 
              unh_stress %>% mutate(std_rate=std_rgr,ctrl_group = NA,.keep="unused")) %>% 
  mutate(type="stress", res="means")

lit_data_plus_ext <- bind_rows(lit_data,
              unh_stress_all %>% mutate(std_rate=std_rgr,.keep="unused")) %>% 
  mutate(type="stress", res="all")

lit_data_plus_ctrl <- bind_rows(lit_data,
              unh_ctrl %>% mutate(std_rate=std_rgr, stress_group = NA,.keep="unused")) %>%
  mutate(type="ctrl", res="means")

lit_data_plus_ext_ctrl <- bind_rows(lit_data,
              unh_ctrl_all %>% mutate(std_rate=std_rgr,.keep="unused")) %>%
  mutate(type="ctrl", res="all")

all_lit_data <- rbind(lit_data_plus, lit_data_plus_ext, lit_data_plus_ctrl, lit_data_plus_ext_ctrl) %>% distinct() %>% # removes duplicate rows (i.e., the original literature data)
  mutate(level = case_when(str_ends(paper_full, "high") ~ "high",
                           str_ends(paper_full, "med") ~ "med",
                           str_ends(paper_full, "low") ~ "low",
                           .default = "lit"),
         level = as_factor(level),
         level = fct_relevel(level, c("high", "med", "low", "lit")))

(ggplot()+
  geom_point(data=lit_data_plus, aes(x=temp, y=std_rate, color=stress_group)))+
(ggplot()+
  geom_point(data=lit_data_plus_ctrl, aes(x=temp, y=std_rate, color=ctrl_group)))+
(ggplot()+
  geom_point(data=lit_data_plus_ext, aes(x=temp, y=std_rate, color=stress_group)))+
(ggplot()+
  geom_point(data=lit_data_plus_ext, aes(x=temp, y=std_rate, color=ctrl_group)))

ggplot()+geom_point(data=all_lit_data, aes(x=temp, y=std_rate, color=level))+facet_grid(type~res)


### NLS incorporating new data ####

nls_fun <- function(df, type, level) {
  #type can be "ctrl" or "stress"
  #level can be "high", "med", or "low"
  #df should be a data frame that includes the original literature values (in order to fit the lower portion of the curve), along with the new points you want to include in the nls fitting
  
  relevant_grouping <- paste(type, "group", sep="_")
  data <- df %>% filter(.data[[relevant_grouping]]=="lit" | .data[[relevant_grouping]]==level)
  
  nls_res <- nlsLM(formula=std_rate ~ exp((T_a/T_0)-(T_a/temp_K)) *
                     (1+exp((T_AL/T_0)-(T_AL/T_L)) + exp((T_ah/T_h)-(T_ah/T_0))) *
                     ((1+exp((T_AL/temp_K)-(T_AL/T_L))+exp((T_ah/T_h)-(T_ah/temp_K)))^-1),
                   start=list(T_a=T_a, T_h=T_h, T_ah=T_ah),
                   data = data)
  
  temp_T_A <- nls_res$m$getPars()[[1]]
  temp_T_H <- nls_res$m$getPars()[[2]]
  temp_T_AH <- nls_res$m$getPars()[[3]]
  
  vector_smooth <- exp((temp_T_A/T_0)-(temp_T_A/temp_smooth))*(1+exp((T_AL/T_0)-(T_AL/T_L))+exp((temp_T_AH/temp_T_H)-(temp_T_AH/T_0))) * ((1+exp((T_AL/temp_smooth)-(T_AL/T_L))+exp((temp_T_AH/temp_T_H)-(temp_T_AH/temp_smooth)))^-1)
  
  df_out <- data.frame(T_A=temp_T_A, T_H=temp_T_H, T_AH=temp_T_AH,type=type,level=level, temp_smooth=temp_smooth, std_rate=vector_smooth)
  df_out
}

##### Stress means #### 
stress_means_df <- map(.x=c(high="high", med="med", low="low"), .f=nls_fun, df=lit_data_plus, type="stress") %>% bind_rows(.id="level") %>% mutate(res="means")

ggplot()+
  geom_line(data=stress_means_df %>% group_by(level), aes(x=temp_smooth, y=std_rate, color=level))

##### Stress all #### 
stress_all_df <- map(.x=c(high="high", med="med", low="low"), .f=nls_fun, df=lit_data_plus_ext, type="stress") %>% bind_rows(.id="level") %>% mutate(res="all")

ggplot()+
  geom_line(data=stress_all_df %>% group_by(level), aes(x=temp_smooth, y=std_rate, color=level))

##### Control means #### 
ctrl_means_df <- map(.x=c(high="high", med="med", low="low"), .f=nls_fun, df=lit_data_plus_ctrl, type="ctrl") %>% bind_rows(.id="level") %>% mutate(res="means")

ggplot()+
  geom_line(data=ctrl_means_df %>% group_by(level), aes(x=temp_smooth, y=std_rate, color=level))

##### Control all #### 
ctrl_all_df <- map(.x=c(high="high", med="med", low="low"), .f=nls_fun, df=lit_data_plus_ext_ctrl, type="ctrl") %>% bind_rows(.id="level") %>% mutate(res="all")

ggplot()+
  geom_line(data=ctrl_all_df %>% group_by(level), aes(x=temp_smooth, y=std_rate, color=level))

all_calibrations <- rbind(stress_means_df, stress_all_df, ctrl_means_df, ctrl_all_df) %>% 
  mutate(level = as_factor(level),
         level = fct_expand(level, "lit"),
         level = fct_relevel(level, c("high", "med", "low", "lit")))


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


ggplot()+
  geom_line(data=all_calibrations %>% group_by(type, level, res), aes(x=temp_smooth-273.15, y=std_rate, color=level), linewidth=2)+
  geom_point(data=all_lit_data, aes(x=temp, y=std_rate, color=level))+
  theme_bw()+
  facet_grid(type~res)+
  #ylim(0,5)+
  labs(x="Temperature (°C)", y="Standardized rate", color=NULL)

ggplot()+
  geom_point(data=all_lit_data %>% filter(res=="means"), aes(x=temp, y=std_rate, color=level))+
  geom_line(data=all_calibrations %>% group_by(type, level, res) %>% filter(res=="means"), aes(x=temp_smooth-273.15, y=std_rate, color=level), linewidth=2)+
  
  theme_bw()+
  facet_wrap(~type, labeller = as_labeller(c(ctrl="Grouped by overall growth", stress='Grouped by heat tolerance')))+
 scale_color_manual(values=c("high"="#dd4124", "med"="#edd746","low"='#0f85a0', "lit"="gray"),
                    breaks=c("high", "med", "low", "lit"),
                    labels=c("high"="High", "med"="Medium", "low"="Low", "lit"="Literature"))+
  labs(x="Temperature (°C)", y="Standardized rate", color=NULL)+
  theme(text = element_text(size=15),
        axis.title.y = element_text(margin = margin(t = 0, r = 9, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 9, r = 0, b = 0, l = 0)))

ggplot()+
  geom_line(data=all_calibrations %>% group_by(type, level, res) %>% filter(res=="all"), aes(x=temp_smooth-273.15, y=std_rate, color=level), linewidth=1)+
  geom_point(data=all_lit_data %>% filter(res=="all"), aes(x=temp, y=std_rate, color=level))+
  theme_bw()+
  facet_wrap(~type, labeller = as_labeller(c(ctrl="Grouped by overall growth", stress='Grouped by heat tolerance')))+
  #ylim(0,5)+
  labs(x="Temperature (°C)", y="Standardized rate", color=NULL)+
  theme(text = element_text(size=15),
  axis.title.y = element_text(margin = margin(t = 0, r = 9, b = 0, l = 0)),
  axis.title.x = element_text(margin = margin(t = 9, r = 0, b = 0, l = 0)))
