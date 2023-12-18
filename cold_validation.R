#Validation of DEB kelp model to data from Trømso, Norway
#A complete description of the Venolia et al., (2020) model is 
# available at https://doi.org/10.1016/j.ecolmodel.2020.109151

# This code was written by Ruby Krasnow between November-December 2023
# Last updated: Dec 17, 2023

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Import libraries
library(tidyverse)
library(patchwork)
library(lubridate)
library(tseries)
library(minpack.lm)
library(deSolve)

#Required for model runs
source("SolveR_R.R")
source("KelpDEB_model.R")
source("new_lit.R")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##### Minerals and Organics Section - Written by Celeste Venolia #####
#Conversion coefficients, organics (n = matrix of chemical indices)
# "food N" "food C" Structure "N reserves" "C reserves" products
#     X_N   X_C      V    E_N    E_C    P
n_O <- matrix(
  + c(0.00, 1.00, 1.00, 0.00, 1.00, 1.00,  #C/C, equals 1 by definition
      + 0.00, 0.50, 1.33, 0.00, 2.00, 1.80,  #H/C, these values show that we consider dry-mass
      + 3.00, 2.50, 1.00, 2.50, 1.00, 0.50,  #O/C
      + 1.00, 0.00, 0.04, 1.00, 0.00, 0.04), nrow=4, ncol=6, byrow = TRUE) #N/C
#V is the C-mol structure of alginate (Alginic acid: (C6H8O6)n)
#E_N is N03- and N02- averaged
#E_C is glucose C6H12O6 (Laminarin: c18h32o16 and mannitol c6h14o6)
#We aren't using the X_N, X_C, or P collumn here

#Molecular weights
#t() is a matrix transpose function
#organics structure matrix multiplied by the atomic masses (mass in grams of one mole of an element) of C H O N
w_O_step <- t(n_O)*matrix(c(12, 1, 16, 14), nrow=6, ncol=4, byrow= TRUE) #g/mol, molecular weights for organics
w_O <- rowSums(w_O_step) #this provides g/mol of each of the six "pockets of mass" (i.e. X_N, X_C)

#define molecular weights
w_V <- w_O[3]  # g/mol       #molecular weight of structure
w_EN <- w_O[4]  # g/mol      #molecular weight of N reserve
w_EC <- w_O[5]  #g/mol       #molecular weight of C reserve
w_O2 <- 32 #g/mol

##### End Minerals and Organics Section - all following code written by Ruby Krasnow ##
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Time steps #####
# Kelp was outplanted April 4 and harvested September 5 (growth period of 155 days)
hourly_seq <- seq(ymd_hms("2018-04-04 00:00:00"),ymd_hms("2018-09-05 00:00:00"), by="hour") # hourly sequence of POSIXct values
times_Matsson <- seq(0, 154*24, 1) #155 days stepped hourly

### Setting up environmental forcing data #####
matsson_temp_df <- read_csv("matsson_temp.csv", col_names=c("date", "temp"), col_types = "Dd") #from Matsson et al., 2012 in J. Appl. Phyc.
matsson_N_df <- read_csv("matssonN.csv", col_names=c("date", "N"), col_types = "Dd") #from Matsson et al., 2012
matsson_PAR_df <- read_csv("matsson_PAR.csv", col_names=c("date", "PAR"), col_types = "Dd") #from Matsson et al., 2012

#From Jones et al., 2019 (Monitoring ocean acidification in Norwegian seas in 2018)
CO_2 <- 2261/10^6 #convert from M to µM

matsson_temp <- matsson_temp_df %>% 
  group_by(date) %>% 
  summarise(temp = mean(temp,na.rm=TRUE)) %>% 
  filter(date > "2018-04-03" & date < "2018-09-07")

matsson_PAR <- matsson_PAR_df %>% 
  group_by(date) %>% 
  summarise(PAR = mean(PAR,na.rm=TRUE)) %>% 
  filter(date > "2018-04-03" & date < "2018-09-07")

matsson_N <- matsson_N_df %>% 
  filter(date > "2018-04-03" & date < "2018-09-07")

env_data <- matsson_temp %>% 
  full_join(matsson_N) %>% 
  full_join(matsson_PAR) %>% 
  mutate(date=as_datetime(date), N=N/10^6) %>% 
  full_join(data.frame(date=hourly_seq)) %>%
  arrange(date) %>% filter(date < as_datetime("2018-09-06: 00:00:00"))

env_data <- env_data %>% mutate(across(c(temp, N, PAR), ~na.approx(.x, rule=2))) %>% 
  na.omit() %>% 
  mutate(PAR=PAR*3600*1e-6,
         temp_K=temp+273.15)

# Irradiance forcing function
I_field <- approxfun(x = seq(from = 0, to = 154*24, by = 1), y = env_data$PAR, method = "linear", rule = 2) 
# Temperature forcing function
T_field <- approxfun(x = seq(from = 0, to = 154*24, by = 1), y = env_data$temp_K, method = "linear", rule = 2)
# Nitrate forcing function
N_field <- approxfun(x = seq(from = 0, to = 154*24, by = 1), y = env_data$N, method = "linear", rule = 2) 



### Initial conditions #####
state_Lo <- c(m_EC = 0.002, #Reserve density of C reserve (initial mass of C reserve per initial mass of structure)
              m_EN = 0.01,  #Reserve density of N reserve (initial mass of N reserve per initial mass of structure)
              M_V = 0.05/(w_V+0.01*w_EN+0.002*w_EC)) #molM_V #initial mass of structure
W <- 0.05 #initial biomass for conversions


### Model runs #####
output_matsson <- params_nested %>% 
  mutate(std_L = future_map(data, function(df) {
  temp_params <- params_Lo
  temp_params[c("T_A", "T_H", "T_AH")] <- c(df$T_A, df$T_H, df$T_AH)
  
  ode_output <- ode(y = state_Lo, t = times_Matsson, func = rates_Lo, parms = temp_params)
   ode_output <- as.data.frame(ode_output) %>% select(time, W, L_allometric) #convert deSolve output into data frame and keep most relevant columns
   ode_output
}))
  
output_matsson <- output_matsson %>%
  select(-data) %>% 
  unnest(cols=std_L) %>% 
  ungroup() %>% 
  group_by(level) %>% 
  mutate(Temp_C = env_data$temp, #conversion back to Celsius from Kelvin
         date=hourly_seq, #add date column
         width= L_allometric/4.5, #using mean from throughout the growing period (see supplementary fig. 1 in Matsson et al., 2021)
         area=0.289*(L_allometric*width)^1.15, #allometric relationship used by Matsson et al.
         .before=1) #putting new columns at the beginning for readability

### Import observed data #####
matsson_obs <- data.frame(date=as_datetime(c("2018-06-08", "2018-06-28", "2018-07-17", "2018-08-01", "2018-08-13", "2018-09-05")),length=c(42,50,66,73,83,87.5))

### Growth figure #####
ggplot(data=output_matsson %>% filter(level!="lit"))+
  geom_line(aes(x=date, y=L_allometric, color=level), linewidth=1)+
  geom_point(data=matsson_obs, aes(x=date, y=length, size="obs"))+
  labs(x="Date", y="Kelp frond length (cm)", color=NULL, size=NULL)+
  theme_classic()+
  scale_size_manual(values=c("obs"=2), breaks=c("obs"), labels=c("obs"="Observations"))+
  scale_color_manual(values=c("cold"='#0f85a0',"warm"="#dd4124", "orig"="black", "Observed"="black"),
                     breaks=c("cold","warm","orig", "Observed"),
                     labels=c("warm"="Warm", "cold"="Cold","orig"="Original", "Observed"="Observations"))+
  theme(text = element_text(size=18),
        axis.title.y = element_text(margin = margin(t = 0, r = 9, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 9, r = 0, b = 0, l = 0)))+
        #legend.position="bottom",
        #legend.box = "horizontal")+ 
  guides(colour = guide_legend(order = 1),
         size = guide_legend(order = 2))

### Environmental figures #####
PAR_plot_Norway<- ggplot(data=env_data)+
  geom_line(aes(x=date, y=PAR/3600*10^6), linewidth=1)+
  labs(x=NULL, y=expression(paste("PAR (μmol photons ",  m^-2, " ",s^-1, ")")), color=NULL, linetype=NULL)+
  theme_classic()+
  theme(text = element_text(size=16),
        axis.title.y = element_text(margin = margin(t = 0, r = 9, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 9, r = 0, b = 0, l = 0)))

N_plot_Norway<-ggplot(data=env_data)+
  geom_line(aes(x=date, y=N*10^6), linewidth=1)+
  labs(x=NULL, y=bquote("NO"[3]~" (µM)"), color=NULL, linetype=NULL)+
  theme_classic()+
  theme(text = element_text(size=16),
        axis.title.y = element_text(margin = margin(t = 0, r = 9, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 9, r = 0, b = 0, l = 0)))

temp_plot_Norway<-ggplot(data=env_data)+
  geom_line(aes(x=date, y=temp), linewidth=1)+
  labs(x=NULL, y="Temperature (°C)", color=NULL, linetype=NULL)+
  theme_classic()+
  theme(text = element_text(size=16),
        axis.title.y = element_text(margin = margin(t = 0, r = 9, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 9, r = 0, b = 0, l = 0)))

PAR_plot_Norway+N_plot_Norway+temp_plot_Norway

matsson_obs %>% 
  left_join(output_matsson) %>% 
  group_by(level) %>% 
  summarise(rmse = rmse(length, L_allometric),
            mae=mae(length, L_allometric)) %>% write_clip()
