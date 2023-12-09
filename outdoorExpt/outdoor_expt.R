#Analysis of outdoor tank sporophyte heat stress growth expt
#Ruby Krasnow
#Last updated: Dec 6, 2023

library(tidyverse)
library(patchwork)
library(lubridate)
library(clipr)
library(tseries)
library(ggstatsplot)
library(minpack.lm)

source("./outdoorExpt/outdoorHOBO/outdoor_HOBO.R")

#### Nitrate data ####################################################################
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

nitrate2 <- rbind(nitrate, nitrate_enriched)

ggplot(data=nitrate2)+
  geom_line(aes(x=date, y=nitrate, color=tank))+
  labs(x=NULL, y="Nitrate (µmol/L)", color="Treatment")+
  theme_light()

#### Growth data ####################################################################
growth_data_orig <- read.csv("~/Downloads/MBL_SES/outdoorExpt/growth.csv") %>% #import data
  mutate(date = mdy(Date), #fix date column
  #Shorten variable names and change to factors as necessary
         heat_test_f = na_if(Heat.testesd.Female.GP.ID,"#N/A"),
         heat_test_m = na_if(Heat.tested.Male.GP.ID,"#N/A"),
         blade_len = Total.length.of.blade,
         hp = as.double(na_if(Distance.base.of.blade.to.10cm.hole.punch,"--")),
         cross = as.factor(Cross),
         trt = as.factor(Treatment),
         rep = as.factor(Rep),
         id = paste(cross, trt, rep, sep="."), #create ID variable as cross.treatment.replicate
         .before = 3, .keep="unused")

growth_data <- growth_data_orig %>% 
  select(c(date, blade_len,hp, cross, trt, rep, id, comments)) %>% #select relevant columns
  complete(date, id) %>%
  filter(hp>9) #the hole punch was done 10cm from the base, this gets rid of a few instances of potential measuring error or degradation while allowing slight room for error
  

growth_rates <- growth_data %>% 
  group_by(id) %>% 
  mutate(
    net_growth = (blade_len - lag(blade_len)), #find difference in total blade length between sampling dates
    growth_hp = if_else(is.na(hp),NA,(hp-lag(hp))), #find difference in distance between hole punch and stipe between sampling dates (provides estimate of growth without effect of erosion from the end of the blade)
    net_perc_change = (blade_len - lag(blade_len))/lag(blade_len), #percent change in total blade length between sampling dates
    perc_change_hp = if_else(is.na(hp),NA,(hp-lag(hp))/lag(hp)), #percent change in stipe to hole punch distance between sampling dates
    time_elapsed = as.integer(date-lag(date)), #days between sampling events
    growth_rate_tot = (blade_len-lag(blade_len))/time_elapsed, #growth rate of overall blade, in cm/day
    growth_rate_hp = if_else(is.na(hp),NA,(hp-lag(hp))/time_elapsed), #growth rate near meristem, in cm/day
     rgr=((log(hp)-log(lag(hp)))/time_elapsed)*100) %>% #relative growth rate, in % per day
    na.omit()

# Here we check the average and SD of the difference between total (entire blade) and hole punch growth measurements.
growth_rates %>% mutate(diff = net_growth-growth_hp) %>% 
  ungroup() %>% summarize(mean(diff, na.rm=TRUE), sd(diff, na.rm=TRUE))

#Adding new columns to data frame to remove most dramatic instances of measurement error
growth_rates <- growth_rates %>% mutate(diff = net_growth-growth_hp,
                                        flag_diff = (diff > 2),
                                        big_loss = (growth_hp < -2))

sum(growth_rates$flag_diff)
sum(growth_rates$big_loss)

#We remove instances where the entire blade growth was substantially more than the hole punch growth or the hole punch length shrunk substantially between weeks
growth_rates_no_flags <- growth_rates %>% filter(big_loss==FALSE, flag_diff==FALSE) %>%
  mutate(trt = fct_drop(trt), cross=fct_drop(cross))

# Note: overall, we are using the three treatments as a way to get growth rates under different temperatures. There were no substantial differences between the high N and low N tanks, so we are considering treatment A to be "control" conditions and everything else to be (thermal) "stress" conditions


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

##### Grouping strains by growth  ################################################################

# Create vectors splitting the 30 crosses into 3 groups of 10 based on their relative growth rates under the control temps
high_ctrl <- mean_growth_by_cross_ctrl %>% slice_max(gr_hp, n=10) %>% mutate(cross = fct_drop(cross))
low_ctrl <- mean_growth_by_cross_ctrl %>% slice_min(gr_hp, n=10) %>% mutate(cross = fct_drop(cross))
mid_ctrl <- anti_join(mean_growth_by_cross_ctrl, high_ctrl) %>% anti_join(low_ctrl) %>% mutate(cross = fct_drop(cross))

# Create vectors splitting the 30 crosses into 3 groups of 10 based on their relative growth rates under the stress temps
high_stress <- mean_growth_by_cross_stress %>% slice_max(gr_hp, n=10) %>% mutate(cross = fct_drop(cross))
low_stress <- mean_growth_by_cross_stress %>% slice_min(gr_hp, n=10) %>% mutate(cross = fct_drop(cross))
mid_stress <- anti_join(mean_growth_by_cross_stress, high_stress) %>% anti_join(low_stress) %>% mutate(cross = fct_drop(cross))

#function to add control grouping tag to dfs with cross in them
add_ctrl_fun <- function(df) { 
  df_out<- df %>% mutate(ctrl_group = case_when(
    as.character(cross) %in% levels(high_ctrl$cross) ~ "high",
    as.character(cross) %in% levels(low_ctrl$cross) ~ "low",
    .default = "med"
  ), ctrl_group = as.factor(ctrl_group),
  ctrl_group = fct_relevel(ctrl_group, c("high", "med", "low")))
  df_out
}

#function to add stress grouping tag to dfs with cross in them
add_stress_fun <- function(df) {
  df_out<- df %>% mutate(stress_group = case_when(
    as.character(cross) %in% levels(high_stress$cross) ~ "high",
    as.character(cross) %in% levels(low_stress$cross) ~ "low",
    .default = "med"
  ), stress_group = as.factor(stress_group),
  stress_group = fct_relevel(stress_group, c("high", "med", "low")))
  df_out
}

# Use those functions to add columns to growth rate data frame that identifies which group each cross is in (both control and stress groupings)
growth_rates_no_flags <- growth_rates_no_flags %>%
  add_ctrl_fun() %>% #add control groups
  add_stress_fun() %>%  #add stress groups
  mutate(ctrl_group = fct_relevel(ctrl_group, c("high", "med", "low")), #reordering the factors for later plotting purposes
         stress_group = fct_relevel(stress_group, c("high", "med", "low"))) %>% 
ungroup()

sum(growth_rates_no_flags$ctrl_group==growth_rates_no_flags$stress_group)/length(growth_rates_no_flags$stress_group) #Around 39% of crosses are in the same group when we consider growth under control temps vs growth under stress

ggplot(data=growth_rates_no_flags)+
  geom_boxplot(aes(y=rgr, x=stress_group, fill=stress_group))+
  facet_wrap(~trt)+
  scale_fill_viridis_d()+
  theme_classic()

# Function to add week column to data frames
add_week_fun <- function (df){
  df_out <- df %>% ungroup() %>% mutate(week = case_when(
    date=="2023-06-28" ~ 2,
    date=="2023-07-05" ~ 3,
    date=="2023-07-13" ~ 4))
  df_out
}

# Function to update treatment column so dfs match
add_trt_fun <- function (df){
  df_out <- df %>% ungroup() %>% mutate(trt = case_when(
      trt=="A" ~ "control",
      trt=="B" ~ "highN",
      trt=="C" ~ "lowN"))
  df_out
}

# calculate weekly mean and maximum relative growth rates for each of the 3 groups, where grouping was done based on growth under non-control temperatures
weekly_growth_stress <- growth_rates_no_flags %>% 
  group_by(stress_group, date, trt) %>%
  summarise(mean_growth = mean(growth_hp, na.rm=TRUE),
            mean_rgr = mean(growth_rate_hp, na.rm=TRUE),
            mean_rgr = mean(rgr, na.rm=TRUE),
            max_rgr = max(rgr, na.rm=TRUE)) %>% 
  add_week_fun() %>% #add week labels
  add_trt_fun() %>% #update treatment group column
  full_join(weekly_means_degC) %>% #join with mean temperature of each week
  filter(!(week %in% c(1,5))) #remove weeks without growth rates

# calculate weekly mean and maximum relative growth rates for each of the 3 groups, where grouping was done based on growth under control temperatures
weekly_growth_ctrl <- growth_rates_no_flags %>% group_by(ctrl_group, date, trt) %>% 
  summarise(mean_growth = mean(growth_hp, na.rm=TRUE),
            mean_rgr = mean(growth_rate_hp, na.rm=TRUE),
            mean_rgr = mean(rgr, na.rm=TRUE),
            max_rgr = max(rgr, na.rm=TRUE)) %>% add_week_fun() %>% add_trt_fun() %>% 
  full_join(weekly_means_degC) %>% 
  filter(!(week %in% c(1,5)))

growth_rates_no_flags <- growth_rates_no_flags %>% #update growth rate data frame
  add_week_fun() %>% #add week labels
  add_trt_fun() %>% #update treatment group column
  left_join(weekly_means_degC) #join with mean temperature of each week

ggplot(data=growth_rates_no_flags %>% na.omit())+
  geom_point(aes(x=mean_temp, y=rgr,color=stress_group,))+
  geom_smooth(aes(x=mean_temp, y=rgr,color=stress_group), se=FALSE)

ggplot(data=growth_rates_no_flags %>% na.omit())+
  geom_point(aes(x=mean_temp, y=rgr,color=ctrl_group))+
  geom_smooth(aes(x=mean_temp, y=rgr,color=ctrl_group), se=FALSE)

#### Set-up for new data NLS #####################################################################
# 20°C was the average temp in the low N experimental tank during week 2, so we use the RGR of each group after week 2 as the growth rate at the Arrhenius reference temperature. However, the experimental kelp blades were still acclimating to their new environments and were generally stressed during the first couple weeks of the experiment— this is also reflected in the PAM measurements. Many of the growth rates were negative or very low during this period, so we used the max RGR instead of mean RGR when considering week 2 data.

##### Groups of 10 level ####################################################################

###### Stress ######
# find the RGR at the reference temperature (20°C) when we split the 30 strains into groups of 3 based on their growth rates in the HEAT STRESS treatments
ref_stress <- weekly_growth_stress %>% 
  group_by(stress_group) %>% 
  filter(round(mean_temp,1)==20) %>% 
  mutate(ref_rgr = max_rgr) %>% 
  select(stress_group, ref_rgr)

# Finding standardized RGR at the group level
unh_stress <- weekly_growth_stress %>% 
  left_join(ref_stress) %>%  #add a column to the weekly means (stress grouping) with reference RGRs
  mutate(across(c(mean_rgr, max_rgr),
                ~if_else(.x<0 | is.na(.x), 0, .x))) %>% # Replace negative or missing RGRs with 0
  mutate(std_rgr_max = max_rgr/ref_rgr, # standardized max RGR (max RGR divided by ref temp RGR)
         std_rgr_mean = mean_rgr/ref_rgr, # standardized mean RGR (mean RGR divided by ref temp RGR)
         std_rgr = case_when(week==2 & trt!="control"~ std_rgr_max, #use max for week 2 non-control
                             .default= std_rgr_mean)) %>% 
  ungroup() %>%
  mutate(paper="unh_stress_all2023", paper_full=paste(paper, stress_group, sep="_"), temp=mean_temp, rate=mean_rgr, extra_info=trt, temp_K=temp+273.15) %>% #adding columns with identifying info
  select(paper, paper_full, temp, rate, std_rgr, extra_info, temp_K, stress_group)


###### Control ####
# find the RGR at the reference temperature (20°C) for the three groups when we split the 30 strains into groups of 3 based on their growth rates in the CONTROL TEMPERATURE treatment
ref_ctrl <- weekly_growth_ctrl %>% ungroup() %>% 
  group_by(ctrl_group) %>% 
  filter(round(mean_temp,1)==20) %>% 
  mutate(ref_rgr = max_rgr) %>% 
  select(ctrl_group, ref_rgr)

# Finding standardized RGR at the group level
unh_ctrl <- weekly_growth_ctrl %>%
  left_join(ref_ctrl) %>% #add a column to the weekly means (stress grouping) with reference RGRs
  mutate(across(c(mean_rgr, max_rgr),
                ~if_else(.x<0 | is.na(.x), 0, .x))) %>% # Replace any negative or missing RGRs with 0
  mutate(std_rgr_max = max_rgr/ref_rgr, # standardized max RGR (max RGR divided by ref temp RGR)
         std_rgr = mean_rgr/ref_rgr) %>% # standardized mean RGR (mean RGR divided by ref temp RGR)
  ungroup() %>%
  mutate(paper="unh_ctrl_all2023", paper_full=paste(paper, ctrl_group, sep="_"), temp=mean_temp, rate=mean_rgr, extra_info=trt, temp_K=temp+273.15) %>% #adding columns with identifying info
  select(paper, paper_full, temp, rate, std_rgr, extra_info, temp_K, ctrl_group)



##### Cross level ####################################################################

###### Stress #####
unh_cross_stress<-growth_rates_no_flags %>% 
  group_by(cross, date, trt, mean_temp) %>% #group by cross, week, and treatment
  summarize(rgr = mean(rgr, na.rm=TRUE)) %>%  #find mean RGR for each cross under each trt in each week
  add_stress_fun() %>% # add stress grouping labels
  left_join(ref_stress) %>% # add reference RGRs for stress groupings
  mutate(std_rate = rgr/ref_rgr) %>% # find standardized rate by dividing RGR at non-20°C temps by the reference rate, which is the stress group RGR at 20°C
  mutate(std_rate = if_else(std_rate<0 | is.na(std_rate), 0, std_rate)) %>% #replace negative rates with 0
  add_ctrl_fun() %>% #add control grouping labels
  ungroup() %>% 
  mutate(paper="unh_cross_stress2023", paper_full = paste(paper, stress_group, sep="_"), temp=mean_temp, rate=rgr, extra_info=paste(trt, cross,sep="_"), temp_K= temp + 273.15) %>%  #adding columns with identifying info
  select(paper, paper_full, temp, rate,std_rate, extra_info, temp_K, ctrl_group, stress_group)

###### Control #### 
unh_cross_ctrl<-growth_rates_no_flags %>% 
  group_by(cross, date, trt, mean_temp) %>% #group by cross, week, and treatment
  summarize(rgr = mean(rgr, na.rm=TRUE)) %>% #find mean RGR for each cross under each trt in each week
  add_ctrl_fun() %>% #add control grouping labels
  left_join(ref_ctrl) %>% # add reference RGRs for control groupings
  mutate(std_rate = rgr/ref_rgr) %>% # find standardized rate by dividing RGR at non-20°C temps by the reference rate, which is the control group RGR at 20°C
  mutate(std_rate = if_else(std_rate<0 | is.na(std_rate), 0, std_rate)) %>% #replace negative rates with 0
  add_stress_fun() %>% # add stress grouping labels
  ungroup() %>% 
  mutate(paper="unh_cross_ctrl2023", paper_full = paste(paper, ctrl_group, sep="_"), temp=mean_temp, rate=rgr, extra_info=paste(trt, cross,sep="_"), temp_K= temp + 273.15) %>%  #adding columns with identifying info
  select(paper, paper_full, temp, rate,std_rate, extra_info, temp_K, ctrl_group, stress_group)

##### Replicate level ################################################################

###### Stress #####
unh_stress_rep <- growth_rates_no_flags %>% # start with all growth rates
  left_join(ref_stress) %>% # add reference RGRs for stress groupings
  mutate(rgr= if_else(rgr<0 | is.na(rgr), 0, rgr), #replace negative rates with 0
         std_rgr = rgr/ref_rgr) %>% # find standardized rate by dividing RGR at non-20°C temps by the reference rate, which is the stress group RGR at 20°C
  ungroup() %>%
  mutate(paper="unh_stress_all2023", paper_full = paste(paper, stress_group, sep="_"), temp=mean_temp, rate=rgr, extra_info=trt, temp_K= temp + 273.15) %>% #adding columns with identifying info
  select(paper, paper_full, temp, rate, std_rgr, extra_info, temp_K, ctrl_group, stress_group)

###### Control #### 
unh_ctrl_rep <- growth_rates_no_flags %>% # start with all growth rates
  left_join(ref_ctrl) %>% # add reference RGRs for control groupings
  mutate(rgr= if_else(rgr<0 | is.na(rgr), 0, rgr), #replace negative rates with 0
         std_rgr = rgr/ref_rgr) %>% # find standardized rate by dividing RGR at non-20°C temps by the reference rate, which is the control group RGR at 20°C
  ungroup() %>%
  mutate(paper="unh_ctrl_all2023", paper_full=paste(paper, ctrl_group, sep="_"), temp=mean_temp, rate=rgr, extra_info=trt, temp_K=temp+273.15) %>% #adding columns with identifying info
  select(paper, paper_full, temp, rate, std_rgr, extra_info, temp_K, ctrl_group,stress_group)

#### Literature data ####################################################################
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

#### Combine new and lit data #### 
lit_data_plus <- bind_rows(lit_data, 
  unh_stress %>% mutate(std_rate=std_rgr,ctrl_group = NA,.keep="unused")) %>% 
  mutate(type="stress", res="means")

lit_data_plus_rep <- bind_rows(lit_data,
  unh_stress_rep %>% mutate(std_rate=std_rgr,.keep="unused")) %>% 
  mutate(type="stress", res="all")

lit_data_plus_ctrl <- bind_rows(lit_data,
  unh_ctrl %>% mutate(std_rate=std_rgr, stress_group = NA,.keep="unused")) %>%
  mutate(type="ctrl", res="means")

lit_data_plus_rep_ctrl <- bind_rows(lit_data,
  unh_ctrl_rep %>% mutate(std_rate=std_rgr,.keep="unused")) %>%
  mutate(type="ctrl", res="all")


lit_data_plus_rep_ctrl <- bind_rows(lit_data,
                                    unh_ctrl_rep %>% mutate(std_rate=std_rgr,.keep="unused")) %>%
  mutate(type="ctrl", res="all")

lit_data_plus_ctrl_cross <- bind_rows(lit_data,
                                      unh_cross_ctrl %>% mutate(stress_group = NA,.keep="unused")) %>%
  mutate(type="ctrl", res="cross")

lit_data_plus_stress_cross <- bind_rows(lit_data,
                                      unh_cross_stress %>% mutate(ctrl_group = NA,.keep="unused")) %>%
  mutate(type="stress", res="cross")

all_lit_data <- rbind(lit_data_plus, lit_data_plus_rep, lit_data_plus_ctrl, lit_data_plus_rep_ctrl, lit_data_plus_ctrl_cross, lit_data_plus_stress_cross) %>% distinct() %>% # removes duplicate rows (i.e., the original literature data)
  mutate(level = case_when(str_ends(paper_full, "high") ~ "high",
                           str_ends(paper_full, "med") ~ "med",
                           str_ends(paper_full, "low") ~ "low",
                           .default = "lit"),
         level = as_factor(level),
         level = fct_relevel(level, c("high", "med", "low", "lit")))

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
                   lower=c(6000,282.15,12000),
                   upper=c(12500,290.15,25000),
                   data = data)
  
  temp_T_A <- nls_res$m$getPars()[[1]]
  temp_T_H <- nls_res$m$getPars()[[2]]
  temp_T_AH <- nls_res$m$getPars()[[3]]
  
  vector_smooth <- exp((temp_T_A/T_0)-(temp_T_A/temp_smooth))*(1+exp((T_AL/T_0)-(T_AL/T_L))+exp((temp_T_AH/temp_T_H)-(temp_T_AH/T_0))) * ((1+exp((T_AL/temp_smooth)-(T_AL/T_L))+exp((temp_T_AH/temp_T_H)-(temp_T_AH/temp_smooth)))^-1)
  
  df_out <- data.frame(T_A=temp_T_A, T_H=temp_T_H, T_AH=temp_T_AH,type=type,level=level, temp_smooth=temp_smooth, std_rate=vector_smooth)
  df_out
}

#### Stress means
stress_means_df <- map(.x=c(high="high", med="med", low="low"), .f=nls_fun, df=lit_data_plus, type="stress") %>% bind_rows(.id="level") %>% mutate(res="means")

##### Stress reps
stress_all_df <- map(.x=c(high="high", med="med", low="low"), .f=nls_fun, df=lit_data_plus_rep, type="stress") %>% bind_rows(.id="level") %>% mutate(res="all")

#### Stress crosses
stress_cross_df <- map(.x=c(high="high", med="med", low="low"), .f=nls_fun, df=lit_data_plus_stress_cross, type="stress") %>% bind_rows(.id="level") %>% mutate(res="cross")

##### Control means
ctrl_means_df <- map(.x=c(high="high", med="med", low="low"), .f=nls_fun, df=lit_data_plus_ctrl, type="ctrl") %>% bind_rows(.id="level") %>% mutate(res="means")

#### Control crosses
ctrl_cross_df <- map(.x=c(high="high", med="med", low="low"), .f=nls_fun, df=lit_data_plus_ctrl_cross, type="ctrl") %>% bind_rows(.id="level") %>% mutate(res="cross")

##### Control reps
ctrl_all_df <- map(.x=c(high="high", med="med", low="low"), .f=nls_fun, df=lit_data_plus_rep_ctrl, type="ctrl") %>% bind_rows(.id="level") %>% mutate(res="all")

all_calibrations <- rbind(stress_means_df, stress_all_df, ctrl_means_df, ctrl_all_df, stress_cross_df, ctrl_cross_df) %>% 
  mutate(level = as_factor(level),
         level = fct_expand(level, "lit"),
         level = fct_relevel(level, c("high", "med", "low", "lit")))

params<-all_calibrations %>% select(T_A, T_H, T_AH,type,level,res) %>% distinct() %>% 
  add_row(T_A=T_A, T_H=T_H,T_AH=T_AH,type="orig", level="orig", res="orig")%>% 
  add_row(T_A=new_T_A, T_H=new_T_H,T_AH=new_T_AH,type="lit", level="lit", res="lit")

ggplot()+
  geom_line(data=all_calibrations %>% group_by(type, level, res), aes(x=temp_smooth-273.15, y=std_rate, color=level), linewidth=2)+
  geom_point(data=all_lit_data, aes(x=temp, y=std_rate, color=level))+
  theme_bw()+
  facet_grid(type~res)+
  #ylim(0,5)+
  labs(x="Temperature (°C)", y="Standardized rate", color=NULL)

means_plot<-ggplot()+
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

ind_plot<-ggplot()+
  geom_line(data=all_calibrations %>% group_by(type, level, res) %>% filter(res=="all"), aes(x=temp_smooth-273.15, y=std_rate, color=level), linewidth=1)+
  #geom_point(data=all_lit_data %>% filter(res=="all"), aes(x=temp, y=std_rate, color=level))+
  theme_bw()+
  scale_color_manual(values=c("high"="#dd4124", "med"="#edd746","low"='#0f85a0', "lit"="gray"),
                     breaks=c("high", "med", "low", "lit"),
                     labels=c("high"="High", "med"="Medium", "low"="Low", "lit"="Literature"))+
  facet_wrap(~type, labeller = as_labeller(c(ctrl="Grouped by overall growth", stress='Grouped by heat tolerance')))+
  labs(x="Temperature (°C)", y="Standardized rate", color=NULL)+
  theme(text = element_text(size=15),
  axis.title.y = element_text(margin = margin(t = 0, r = 9, b = 0, l = 0)),
  axis.title.x = element_text(margin = margin(t = 9, r = 0, b = 0, l = 0)))

means_plot/ind_plot

ggplot()+
  geom_line(data=all_calibrations %>% ungroup(), aes(x=temp_smooth-273.15, y=std_rate, color=level), linewidth=2)+
  theme_bw()+
  facet_grid(type~res)+
  annotate(geom='line', x=temp_smooth-273.15,y=y_smooth_venolia)+
  labs(x="Temperature (°C)", y="Standardized rate", color=NULL)+
  scale_color_manual(values=c("high"="#dd4124", "med"="#edd746","low"='#0f85a0'),
                     breaks=c("high", "med", "low"),
                     labels=c("high"="High", "med"="Medium", "low"="Low"))

myPlot<-ggplot()+
  geom_line(data=all_calibrations %>% group_by(type, level, res) %>% filter(type=="ctrl", res=="all" | res=="means"), aes(x=temp_smooth-273.15, y=std_rate, color=level), linewidth=2)+
  theme_bw()+
  labs(x="Temperature (°C)", y="Standardized rate", color=NULL)+
  scale_color_manual(values=c("high"="#dd4124", "med"="#edd746","low"='#0f85a0'),
                     breaks=c("high", "med", "low"),
                     labels=c("high"="High", "med"="Medium", "low"="Low"))+
  theme(legend.position = "none")+
  ylim(0,2.6)+
  facet_wrap(~res)
myPlot

