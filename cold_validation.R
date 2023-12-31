#Validation of DEB kelp model to data from Trømso, Norway
#A complete description of the Venolia et al., (2020) model is 
# available at https://doi.org/10.1016/j.ecolmodel.2020.109151

# This code was written by Ruby Krasnow between November-December 2023
# Last updated: Dec 31, 2023

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Import libraries
library(tidyverse)
library(patchwork)
library(lubridate)
library(tseries)
library(minpack.lm)
library(deSolve)
library(zoo)

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

##### Parameters compiled #####
params_Lo <- c(#maximum volume-specific assimilation rate of N before temperature correction
  JENAM = 1.5e-4, #mol N / molV / h
  #half saturation constant of N uptake
  K_N = 2.5e-6, #molNO3 and NO2/L
  #max volume-specific carbon dioxide assimilation rate
  JCO2M = 0.0075, #mol DIC/molV/h
  #half saturation constant of C uptake
  K_C = 4e-7, #mol DIC/L
  #maximum volume-specific carbon assimilation rate
  JECAM = 0.282, #molC/molV/h
  #Photosynthetic unit density
  rho_PSU = 0.5, #mol PSU/ mol V
  #binding probability of photons to a Light SU
  b_I = 0.5, #dimensionless
  #Specific photon arrival cross section
  alpha = 1, #m^2 mol PSU–1
  #dissociation rate
  k_I = 0.075, #molγ molPS–1 h–1
  #Yield factor of C reserve to photon
  y_I_C = 10, #mol γ mol C-1
  #Yield factor of C reserve to DIC
  y_CO2_C = 1, #mol DIC mol C-1
  #Yield factor of photon to O2
  y_LO2 = 0.125, #molO2 molγ –1
  #reserve turnover
  kE_C = 0.02, #0.05, #1/h
  kE_N = 0.04, #0.01, #1/h
  #fraction of rejection flux from growth SU incorporated back into i-reserve
  kappa_Ei = 0.9, #dimensionless
  #yield of structure on N reserve (percent of N in structure)
  y_EN_V = 0.04, #mol N/mol V
  #yield of structure on C reserve (percent of C in structure)
  y_EC_V = 1, #mol C/mol V
  #specific maintenance costs requiring N before temp correction
  JENM = 4*10^-6, #4e-6, #mol N/molM_V/h
  #specific maintenance costs requiring C before temp correction
  JECM = 1*10^-6, #1e-6, #mol C/molM_V/h
  #Arrhenius temperature
  T_A = 6314.3, # K
  #Upper boundary of temperature tolerance
  T_H = 13.386 + 273.15, # K
  #Lower boundary of temperature tolerance
  T_L = 273.15, # K
  #Arrhenius temperature outside T_H
  T_AH = 18702, #K
  #Arrhenius temperature outside T_L
  T_AL = 4391.9, #K
  #temperature at which rate parameters are given
  T_0 = 20 + 273.15) # K

##### End Minerals and Organics Section - all following code written by Ruby Krasnow ##
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# MATSSON 2012 ------------------------------------------------------------


###### Time steps ######
# Kelp was outplanted April 4 and harvested September 5 (growth period of 155 days)
hourly_seq <- seq(ymd_hms("2018-04-04 00:00:00"),ymd_hms("2018-09-05 00:00:00"), by="hour") # hourly sequence of POSIXct values
times_Matsson <- seq(0, 154*24, 1) #155 days stepped hourly

###### Environmental forcing data #####
matsson_temp_df <- read_csv("./validation_data/matsson2021/matsson_temp.csv", col_names=c("date", "temp"), col_types = "Dd") #from Matsson et al., 2012 in J. Appl. Phyc.
matsson_N_df <- read_csv("./validation_data/matsson2021/matssonN.csv", col_names=c("date", "N"), col_types = "Dd") #from Matsson et al., 2012
matsson_PAR_df <- read_csv("./validation_data/matsson2021/matsson_PAR.csv", col_names=c("date", "PAR"), col_types = "Dd") #from Matsson et al., 2012

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

###### Initial conditions #####
state_Lo <- c(m_EC = 0.002, #Reserve density of C reserve (initial mass of C reserve per initial mass of structure)
              m_EN = 0.01,  #Reserve density of N reserve (initial mass of N reserve per initial mass of structure)
              M_V = 0.05/(w_V+0.01*w_EN+0.002*w_EC)) #molM_V #initial mass of structure
W <- 0.05 #initial biomass for conversions

###### Model runs #####
output_matsson <- params_nested %>% 
  mutate(std_L = future_map(data, function(df) {
  temp_params <- params_Lo
  temp_params[c("T_A", "T_H", "T_AH")] <- c(df$T_A, df$T_H, df$T_AH)
  
  ode_output <- ode(y = state_Lo, t = times_Matsson, func = rates_Lo, parms = temp_params)
   ode_output <- as.data.frame(ode_output) #convert deSolve output into data frame
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

###### Import observed data #####
matsson_obs <- data.frame(date=as_datetime(c("2018-06-08", "2018-06-28", "2018-07-17", "2018-08-01", "2018-08-13", "2018-09-05")),length=c(42,50,66,73,83,87.5))

###### Figures #####
# Growth figure
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
  guides(colour = guide_legend(order = 1),
         size = guide_legend(order = 2))

### Environmental figures
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
            mae=mae(length, L_allometric)) #%>% write_clip()

output_matsson %>% filter(time==3696) %>% select(date, L_allometric, level)


# JEVNE 2020 --------------------------------------------------------------
#D4: deep water low light
#S1: surface water high light

###### Time steps ######
# Experimental period was June 1-19th, 2014 (growth period of 20 days)
# this followed a 4-week acclimation period.
hourly_seq_jevne <- seq(ymd_hms("2014-06-01 00:00:00"),ymd_hms("2014-06-20 00:00:00"), by="hour") # hourly sequence of POSIXct values
times_jevne <- seq(0, 19*24, 1) #20 days stepped hourly


###### Environmental forcing data #####

N_files <- c(D1_N="D1_N", D4_N="D4_N", S1="S1_N", S4_N="S4_N")
N_filepaths <- N_files %>% map_chr(\(x) paste0("./validation_data/jevne2020/", x, ".csv"))

jevne_N <- read_csv(N_filepaths, id = "trt", col_names = c("date", "N"), col_types = "cd") %>% 
  mutate(date = ymd_hms(paste0(date, ":00")))  %>%
  mutate(trt = tools::file_path_sans_ext(basename(trt)),
         trt = as_factor(str_sub(trt,1,2)),
        date = round_date(date, unit="day"),
        N = N/14/10^6) #N is imported in µg/L, we want mol/L

PAR_files <- c(D1_PAR="D1_PAR", D4_PAR="D4_PAR", S1_PAR="S1_PAR", S4_PAR="S4_PAR")
PAR_filepaths <- PAR_files %>% map_chr(\(x) paste0("./validation_data/jevne2020/", x, ".csv"))

jevne_PAR <- read_csv(PAR_filepaths, id = "trt", col_names = c("date", "PAR"), col_types = "cd") %>% 
  mutate(date = ymd_hms(paste0(date, ":00")))  %>%
  mutate(trt = tools::file_path_sans_ext(basename(trt)),
         trt = as_factor(str_sub(trt,1,2)),
         date = round_date(date, unit="day"))

temp_files <- c(D1_temp="D1_temp", D4_temp="D4_temp", S1_temp="S1_temp", S4_temp="S4_temp")
temp_filepaths <- temp_files %>% map_chr(\(x) paste0("./validation_data/jevne2020/", x, ".csv"))

jevne_temp <- read_csv(temp_filepaths, id = "trt", col_names = c("date", "temp"), col_types = "cd") %>% 
  mutate(date = ymd_hms(paste0(date, ":00")))   %>%
  mutate(trt = tools::file_path_sans_ext(basename(trt)),
         trt = as_factor(str_sub(trt,1,2)),
         date = round_date(date, unit="day"))

jevne_N<- jevne_N %>% arrange(date) %>% group_by(trt, date) %>% summarise(N = mean(N))
jevne_temp<- jevne_temp %>% arrange(date) %>% group_by(trt, date) %>% summarise(temp = mean(temp))
jevne_PAR <- jevne_PAR %>% arrange(date) %>% group_by(trt, date) %>% summarise(PAR = mean(PAR))

env_data_jevne <- jevne_temp %>% ungroup() %>% 
  full_join(jevne_N %>% ungroup()) %>% 
  full_join(jevne_PAR %>% ungroup()) %>% 
  full_join(data.frame(date=hourly_seq_jevne)) %>%
  arrange(date) %>% 
  complete(date, trt) %>% 
  filter(date < as_datetime("2014-06-20: 01:00:00") & date > as_datetime("2014-05-21: 00:00:00")) %>%
  arrange(date) %>% 
  filter(!is.na(trt))
  
env_data_jevne <- env_data_jevne %>%
  group_by(trt) %>% 
  mutate(across(c(temp, N, PAR), ~na.approx(.x, rule=2))) %>% 
  mutate(temp_K=temp+273.15) %>% 
  filter(date > as_datetime("2014-05-31: 23:00:00"))

###### Initial conditions #####
# The kelp started the acclimation period with an avg length of 45.8 ± 2.4 cm
# Assuming they grew to be about as large as 1 SD above the mean at the beginning of the acclimation period would give a mean around 48.5
View(output_matsson %>% filter(L_allometric < 49 & L_allometric > 48) %>% select(m_EC, m_EN, M_V, W, L_allometric, level))
#when length is around 48.5 cm, m_EC around 0.313, m_EN = 0.000918, M_V=0.0295, W=1.16
state_jevne <- c(m_EC = 0.313, #Reserve density of C reserve (initial mass of C reserve per initial mass of structure)
              m_EN = 0.000918,  #Reserve density of N reserve (initial mass of N reserve per initial mass of structure)
              M_V = 1.16/(w_V+0.000918*w_EN+0.313*w_EC)) #molM_V #initial mass of structure
W <- 1.16 #initial biomass for conversions

####### ~ S1 #####
S1 <- env_data_jevne %>% filter(trt=="S1")
# Irradiance forcing function
I_field <- approxfun(x = seq(from = 0, to = 19*24, by = 1), y = S1$PAR*3600*1e-6, method = "linear", rule = 2) 
# Temperature forcing function
T_field <- approxfun(x = seq(from = 0, to = 19*24, by = 1), y = S1$temp_K, method = "linear", rule = 2)
# Nitrate forcing function
N_field <- approxfun(x = seq(from = 0, to = 19*24, by = 1), y = S1$N, method = "linear", rule = 2) 

S1_output <- params_nested %>% mutate(std_L = future_map(data, function(df) {
              temp_params <- params_Lo
              temp_params[c("T_A", "T_H", "T_AH")] <- c(df$T_A, df$T_H, df$T_AH)
              ode_output <- ode(y = state_jevne, t = times_jevne, func = rates_Lo, parms = temp_params)
              ode_output <- as.data.frame(ode_output) #convert deSolve output into data frame
        ode_output })) %>% 
     select(-data) %>%
     unnest(cols=std_L) %>%
     ungroup() %>%
     group_by(level) %>%
     mutate(date=hourly_seq_jevne, trt="S1", #add date and trt columns
            .before=1) #putting new columns at the beginning for readability

####### ~ D1 #####
D1 <- env_data_jevne %>% filter(trt=="D1")
I_field <- approxfun(x = seq(from = 0, to = 19*24, by = 1), y = D1$PAR*3600*1e-6, method = "linear", rule = 2) 
T_field <- approxfun(x = seq(from = 0, to = 19*24, by = 1), y = D1$temp_K, method = "linear", rule = 2)
N_field <- approxfun(x = seq(from = 0, to = 19*24, by = 1), y = D1$N, method = "linear", rule = 2) 

D1_output <- params_nested %>% mutate(std_L = future_map(data, function(df) {
  temp_params <- params_Lo
  temp_params[c("T_A", "T_H", "T_AH")] <- c(df$T_A, df$T_H, df$T_AH)
  ode_output <- ode(y = state_jevne, t = times_jevne, func = rates_Lo, parms = temp_params)
  ode_output <- as.data.frame(ode_output) #convert deSolve output into data frame
  ode_output }))  %>% 
  select(-data) %>%
  unnest(cols=std_L) %>%
  ungroup() %>%
  group_by(level) %>%
  mutate(date=hourly_seq_jevne, trt="D1", #add date and trt columns
         .before=1) #putting new columns at the beginning for readability

######## ~ S4 #####
S4 <- env_data_jevne %>% filter(trt=="S4")
I_field <- approxfun(x = seq(from = 0, to = 19*24, by = 1), y = S4$PAR*3600*1e-6, method = "linear", rule = 2) 
T_field <- approxfun(x = seq(from = 0, to = 19*24, by = 1), y = S4$temp_K, method = "linear", rule = 2)
N_field <- approxfun(x = seq(from = 0, to = 19*24, by = 1), y = S4$N, method = "linear", rule = 2) 

S4_output <- params_nested %>% mutate(std_L = future_map(data, function(df) {
  temp_params <- params_Lo
  temp_params[c("T_A", "T_H", "T_AH")] <- c(df$T_A, df$T_H, df$T_AH)
  ode_output <- ode(y = state_jevne, t = times_jevne, func = rates_Lo, parms = temp_params)
  ode_output <- as.data.frame(ode_output) #convert deSolve output into data frame
  ode_output }))  %>% 
  select(-data) %>%
  unnest(cols=std_L) %>%
  ungroup() %>%
  group_by(level) %>%
  mutate(date=hourly_seq_jevne, trt="S4", #add date and trt columns
         .before=1) #putting new columns at the beginning for readability

####### ~ D4 #####
D4 <- env_data_jevne %>% filter(trt=="D4")
I_field <- approxfun(x = seq(from = 0, to = 19*24, by = 1), y = D4$PAR*3600*1e-6, method = "linear", rule = 2) 
T_field <- approxfun(x = seq(from = 0, to = 19*24, by = 1), y = D4$temp_K, method = "linear", rule = 2)
N_field <- approxfun(x = seq(from = 0, to = 19*24, by = 1), y = D4$N, method = "linear", rule = 2) 

D4_output <- params_nested %>% mutate(std_L = future_map(data, function(df) {
  temp_params <- params_Lo
  temp_params[c("T_A", "T_H", "T_AH")] <- c(df$T_A, df$T_H, df$T_AH)
  ode_output <- ode(y = state_jevne, t = times_jevne, func = rates_Lo, parms = temp_params)
  ode_output <- as.data.frame(ode_output) #convert deSolve output into data frame
  ode_output }))  %>% 
  select(-data) %>%
  unnest(cols=std_L) %>%
  ungroup() %>%
  group_by(level) %>%
  mutate(date=hourly_seq_jevne, trt="D4", #add date and trt columns
         .before=1) #putting new columns at the beginning for readability

jevne_output <- bind_rows(S1_output, D1_output, S4_output, D4_output)


###### Import observed data and combine with model output #####
jevne_growth <- read_delim("./validation_data/jevne2020/jevne2020growthN.CSV", delim=";", show_col_types = FALSE) %>% select(Sampleday, Treatment, Replicate, Growth_mean)

jevne_growth <- jevne_growth %>% mutate(date = dmy(Sampleday),
                                        trt = as.factor(Treatment),
                                        rep = as.factor(Replicate),
                                        growth = as.double(str_replace(Growth_mean,",", ".")), .keep="unused")
jevne_growth_means <- jevne_growth %>% 
  group_by(date, trt) %>% 
  summarise(mean_growth = mean(growth),
            sd = sd(growth))

jevne_model_growth <- jevne_output %>% 
  filter(date %in% as_datetime(c("2014-06-01","2014-06-05", "2014-06-10", "2014-06-15", "2014-06-20"))) %>% 
  select(date, trt, level, L_allometric)

jevne_model_growth <- jevne_model_growth %>% 
  arrange(date) %>% group_by(level, trt) %>% 
  mutate(time_elapsed=as.numeric(as.duration(date-lag(date)), "days"),
                              growth = L_allometric-lag(L_allometric),
                              model_growth_rate=growth/time_elapsed,
                              date=date(date))

###### RMSE ####
jevne_rmse <- jevne_growth_means %>% right_join(jevne_model_growth) %>% 
  group_by(level, trt) %>% 
  na.omit() %>% 
  summarise(rmse=rmse(mean_growth, model_growth_rate))

jevne_rmse %>% 
  ungroup() %>% 
  friedman_test(rmse ~ level|trt)

jevne_rmse_matrix <- jevne_rmse %>% 
  ungroup() %>% 
  pivot_wider(names_from = level, values_from = rmse) %>%
  column_to_rownames(var="trt") %>% 
  as.matrix() 

friedman.test(jevne_rmse_matrix)
frdAllPairsNemenyiTest(jevne_rmse_matrix, p.adjust = "bonferroni")

###### Figures ####

# Growth figures
ggplot(data=jevne_output %>% filter(level!="lit"))+
  geom_line(aes(x=date, y=L_allometric, color=trt), linewidth=1)+
  labs(x="Date", y="Kelp frond length (cm)", color=NULL)+
  facet_wrap(~level)+
  theme_classic()+
  theme(text = element_text(size=18),
        axis.title.y = element_text(margin = margin(t = 0, r = 9, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 9, r = 0, b = 0, l = 0)))

ggplot(jevne_growth_means, aes(x=date, y=mean_growth, group=trt, color=trt)) + 
  geom_line()+
  geom_pointrange(aes(ymin=mean_growth-sd, ymax=mean_growth+sd))

ggplot(jevne_growth_means, aes(x=date, y=mean_growth, group=trt, color=trt)) + 
  geom_line()+
  geom_pointrange(aes(ymin=mean_growth-sd, ymax=mean_growth+sd))+
  geom_point(data=jevne_model_growth %>% na.omit(), aes(x=date, y=model_growth_rate))

ggplot(jevne_growth_means, aes(x=date, y=mean_growth)) +
  geom_pointrange(aes(ymin=mean_growth-sd, ymax=mean_growth+sd))+
  geom_point(data=jevne_model_growth %>% na.omit(), aes(x=date, y=model_growth_rate, color=level))+
  facet_wrap(~trt)

### Environmental figures
PAR_plot_jevne<- ggplot(data=env_data_jevne)+
  geom_line(aes(x=date, y=PAR, color=trt), linewidth=1)+
  labs(x=NULL, y=expression(paste("PAR (μmol photons ",  m^-2, " ",s^-1, ")")), color=NULL, linetype=NULL)+
  theme_classic()+
  theme(text = element_text(size=16),
        axis.title.y = element_text(margin = margin(t = 0, r = 9, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 9, r = 0, b = 0, l = 0)))

N_plot_jevne<-ggplot(data=env_data_jevne)+
  geom_line(aes(x=date, y=N*10^6, color=trt), linewidth=1)+
  labs(x=NULL, y=bquote("NO"[3]~" (µM)"), color=NULL, linetype=NULL)+
  theme_classic()+
  theme(text = element_text(size=16),
        axis.title.y = element_text(margin = margin(t = 0, r = 9, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 9, r = 0, b = 0, l = 0)))

temp_plot_jevne<-ggplot(data=env_data_jevne)+
  geom_line(aes(x=date, y=temp, color=trt), linewidth=1)+
  labs(x=NULL, y="Temperature (°C)", color=NULL, linetype=NULL)+
  theme_classic()+
  theme(text = element_text(size=16),
        axis.title.y = element_text(margin = margin(t = 0, r = 9, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 9, r = 0, b = 0, l = 0)))

PAR_plot_jevne+N_plot_jevne+temp_plot_jevne+
  plot_layout(guides = 'collect')