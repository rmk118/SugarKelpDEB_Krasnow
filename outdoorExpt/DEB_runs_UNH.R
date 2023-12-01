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
source("SolveR_R.R")
source("KelpDEB_model.R")
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

# New (refit from literature)
new_params_for_model <- params_Lo
new_params_for_model[c("T_A", "T_H", "T_AH")]<-c(new_T_A, new_T_H, new_T_AH)

# Means
high_params_for_model <- params_Lo
med_params_for_model <- params_Lo
low_params_for_model <- params_Lo

high_params_for_model[c("T_A", "T_H", "T_AH")]<-params %>% filter(type=="ctrl", level=="high", res=="means") %>% select(T_A, T_H, T_AH)
med_params_for_model[c("T_A", "T_H", "T_AH")]<-params %>% filter(type=="ctrl", level=="med", res=="means") %>% select(T_A, T_H, T_AH)
low_params_for_model[c("T_A", "T_H", "T_AH")]<-params %>% filter(type=="ctrl", level=="low", res=="means") %>% select(T_A, T_H, T_AH)

# Crosses
high_params_for_model_cross <- params_Lo
med_params_for_model_cross <- params_Lo
low_params_for_model_cross <- params_Lo

high_params_for_model_cross[c("T_A", "T_H", "T_AH")]<-params %>% filter(type=="ctrl", level=="high", res=="all") %>% select(T_A, T_H, T_AH)
med_params_for_model_cross[c("T_A", "T_H", "T_AH")]<-params %>% filter(type=="ctrl", level=="med", res=="all") %>% select(T_A, T_H, T_AH)
low_params_for_model_cross[c("T_A", "T_H", "T_AH")]<-params %>% filter(type=="ctrl", level=="low", res=="all") %>% select(T_A, T_H, T_AH)

# Replicates
high_params_for_model_rep <- params_Lo
med_params_for_model_rep <- params_Lo
low_params_for_model_rep <- params_Lo

high_params_for_model_rep[c("T_A", "T_H", "T_AH")]<-params %>% filter(type=="ctrl", level=="high", res=="all") %>% select(T_A, T_H, T_AH)
med_params_for_model_rep[c("T_A", "T_H", "T_AH")]<-params %>% filter(type=="ctrl", level=="med", res=="all") %>% select(T_A, T_H, T_AH)
low_params_for_model_rep[c("T_A", "T_H", "T_AH")]<-params %>% filter(type=="ctrl", level=="low", res=="all") %>% select(T_A, T_H, T_AH)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

####### Initial conditions ############
#Initial conditions of state variables

growth_data %>% filter(date==as_date("2023-06-23")) %>% ggplot()+geom_histogram(aes(x=blade_len))+geom_vline(aes(xintercept=mean(blade_len)), color="red")+geom_vline(aes(xintercept=median(blade_len)), color="blue")

growth_data %>% filter(date==as_date("2023-06-23")) %>% summarise(median_len = median(blade_len))
PJN2_meandat %>% filter(Date==as_date("2018-02-14")) %>% pull(mean_length)
sol_all %>% filter(source=="Point Judith Pond N 1", Date == as_datetime("2018-02-14 24:00:00"), params=="orig")  %>% select(L_allometric, M_V, m_EC, m_EN, W)

#Initial conditions of state variables - estimated
start_Wd <- 0.00387*33^(1.469)  
state_Lo_UNH <- c(m_EC = 0.08, #mol C/molM_V  #Reserve density of C reserve (initial mass of C reserve per initial mass of structure)
  m_EN = 0.008, #mol N/molM_V #Reserve density of N reserve (initial mass of N reserve per initial mass of structure)
  # M_V = Wd/(w_V+M_EN*w_EN+M_EC*w_EC)
  M_V = start_Wd/(w_V+0.008*w_EN+0.08*w_EC)) #mol M_V #initial mass of structure

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####### Time steps & forcing set-up #######

#going to run from 2023-06-23 01:00:00 to 2023-07-13 23:00:00, the last day of measuring
times_Lo_UNH <- seq(0, 519, 1)

# Irradiance forcing
I_field <- approxfun(x = seq(from = 0, to = 519, by = 1), y = Z$PAR, method = "linear", rule = 2)
# Nitrate forcing
N_field <- approxfun(x = seq(from = 0, to = 519, by = 1), y = Z$nitrate, method = "linear", rule = 2)