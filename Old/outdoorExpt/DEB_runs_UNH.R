### DEB model runs for outdoor adult sporophyte growth expt ### 

# Created by Ruby Krasnow in November-December 2023

#Import libraries
library(deSolve)
library(tidyverse)
library(lubridate)
library(gridExtra)
library(Metrics)
library(patchwork)

#Required for model runs
source("~/Downloads/SugarKelpDEB_Krasnow/SolveR_R.R") 
source("~/Downloads/SugarKelpDEB_Krasnow/KelpDEB_model.R")
source("./outdoorExpt/outdoorHOBO/outdoor_HOBO.R")
source("./outdoorExpt/outdoor_expt.R")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### Minerals and Organics Section #####
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
#We aren't using the X_N, X_C, or P columns here

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
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Ruby's Note
#Note that * performs element-wise multiplication, NOT matrix multiplication, which is denoted by %*%. 
# This means the value in n_0[1,1] is multiplied by the value in the [1,1] position in the 
# matrix of molecular weights: matrix(c(12, 1, 16, 14), nrow=6, ncol=4, byrow= TRUE).
# I found it more intuitive to achieve the same result by doing: w_O <- t(n_O) %*% c(12,1,16,14)
# Then you can define the molecular weights as w_V <- w_O[3], w_EN <- w_O[4], etc.

# Not entirely sure where the 0.04 N/C comes from in the structure if V is assumed to have the C-mol structure of alginate

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##### Parameters compiled #####
params_Lo_unh <- c(#maximum volume-specific assimilation rate of N before temperature correction
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
  kE_C = 0.02,
  kE_N = 0.04,
  #fraction of rejection flux from growth SU incorporated back into i-reserve
  kappa_Ei = 0.8, #dimensionless
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

####### Initial conditions ############
#Initial conditions of state variables

growth_data %>% filter(date==as_date("2023-06-23")) %>% ggplot()+geom_histogram(aes(x=blade_len))+geom_vline(aes(xintercept=mean(blade_len)), color="red")+geom_vline(aes(xintercept=median(blade_len)), color="blue")

growth_data %>% filter(date==as_date("2023-06-23"))  %>% group_by(trt) %>% summarise(mean_len = mean(blade_len))
growth_data %>% filter(date==as_date("2023-07-13")) %>% group_by(trt) %>% summarise(mean_len = mean(blade_len))

PJN2_meandat %>% filter(Date==as_date("2018-02-14")) %>% pull(mean_length)
sol_all %>% filter(source=="Point Judith Pond N 1", Date == as_datetime("2018-02-14 24:00:00"), params=="orig")  %>% select(L_allometric, M_V, m_EC, m_EN, W)

#Initial conditions of state variables - gives starting length of 35cm, using #L=(DW/0.0155)^(1/1.3587)
state_Lo_UNH <- c(m_EC = 0.01, m_EN = 0.001, M_V = 2/(w_V+0.001*w_EN+0.01*w_EC))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####### Time steps & forcing set-up #######

#going to run from 2023-06-23 01:00:00 to 2023-07-13 23:00:00, the last day of measuring
times_Lo_UNH <- seq(0, 526, 1)

# Irradiance forcing
I_field <- approxfun(x = c(0:526), y = Z$PAR, method = "linear", rule = 2)
# Nitrate forcing
N_field <- approxfun(x = c(0:526), y = Z$nitrate, method = "linear", rule = 2)
# Temperature forcing - control tank
T_field <- approxfun(x = c(0:526), y = Z$highN_temp+273.15, method = "linear", rule = 2) ## HIGH TEMP

W <- 2
#sol_control_tank <- ode(y = state_Lo_UNH, t = times_Lo_UNH, func = rates_Lo, parms = low_params_for_model_cross)
#sol_control_tank <- clean_ode_sol(sol_control_tank, Z$date, times_Lo_UNH, source="tank") %>% mutate(L_allometric)
# 
# ggplot()+
#   geom_smooth(data=sol_control_tank, aes(x=Date, y=L_allometric))


run_model_fun_small <- function(p) {
  ode_out = ode(y = state_Lo_UNH, t = times_Lo_UNH, func = rates_Lo, parms = p)
  ode_out = as.data.frame(ode_out) %>% select(time, W, L_allometric)
  ode_out
}

update_params <- function(T_A, T_H, T_AH) {
  params_temp <- params_Lo_unh
  params_temp[c("T_A", "T_H", "T_AH")] <- c(T_A, T_H, T_AH)
  params_temp
}
 
params_list<- params %>% rowwise() %>%  mutate(
  params_list = list(update_params(T_A, T_H, T_AH)),
  )
  
sols_list<- lapply(params_list$params_list, run_model_fun_small)

params_list_hot <- params_list %>% ungroup() %>% mutate(name=c(1:20)) %>% left_join(enframe(sols_list)) %>% 
  select(-c(name, params_list)) %>% unnest_longer(value) %>% mutate(
    len = value$L_allometric,
    W=value$W,
    date = rep(Z$date, 20),.keep="unused"
  )

final_lens <- params_list_hot %>% group_by(type, level, res) %>% summarise(final_len = max(len))

ggplot()+
  geom_smooth(data=params_list_hot %>% filter(level!="lit", level!="orig", res!="means"), aes(x=date, y=len, color=level))+facet_grid(res~type)+geom_hline(yintercept = 45)



T_field <- approxfun(x = c(0:526), y = Z$control_temp+273.15, method = "linear", rule = 2)
sols_list_cold<- lapply(params_list$params_list, run_model_fun_small)

params_list_cold <- params_list %>% ungroup() %>% mutate(name=c(1:20)) %>% left_join(enframe(sols_list_cold)) %>% 
  select(-c(name, params_list)) %>% unnest_longer(value) %>% mutate(
    len = value$L_allometric,
    W=value$W,
    date = rep(Z$date, 20),.keep="unused"
  )

hot_plot<-ggplot()+
  geom_smooth(data=params_list_hot %>% filter(level!="lit", level!="orig", res=="cross"), aes(x=date, y=len, color=level))+facet_wrap(~type)+labs(x=NULL, y="calibrated under heat")

cold_plot<-ggplot()+
  geom_smooth(data=params_list_cold %>% filter(level!="lit", level!="orig", res=="cross"), aes(x=date, y=len, color=level))+facet_wrap(~type)+labs(x=NULL, y="calibrated under normal")

hot_plot/cold_plot
