
#Required for Calibration Code
source("SolveR_R.R")
source("KelpDEB_model.R")
source("N_uptake_Calibration.R")
source("Photosynthesis_Calibration.R")


state_Johansson <- c(m_EC = 0.3, #mol C/molM_V  #Reserve density of C reserve (initial mass of C reserve per initial mass of structure)
                     m_EN = 0.01, #mol N/molM_V #Reserve density of N reserve (initial mass of N reserve per initial mass of structure)
                     M_V = 0.005/(w_V+0.01*w_EN+0.3*w_EC)) #molM_V #initial mass of structure

#### Literature data for comparison/Calibration ####
###### Nitrate uptake #######
#Espinoza and Chapman (1983) and Ahn et al. (1998)
EC1983_9C_Nuptake_StM <- read.csv("EspinozaChapman1983_Nuptake_9C_StMargaretsBay.csv", header = TRUE, fileEncoding="UTF-8-BOM")
EC1983_18C_Nuptake_StM <- read.csv("EspinozaChapman1983_Nuptake_18C_StMargaretsBay.csv", header = TRUE, fileEncoding="UTF-8-BOM")

#conversions 9C
EC1983_9C_Nuptake_StM$N <- EC1983_9C_Nuptake_StM$ResidualNitrateConcentration
EC1983_9C_Nuptake_StM$N <- round(EC1983_9C_Nuptake_StM$N, digits = 2)
EC1983_9C_Nuptake_StM$N <- EC1983_9C_Nuptake_StM$N/1000000 #microM to M
EC1983_9C_Nuptake_StM$NuptakeRate <- EC1983_9C_Nuptake_StM$NuptakeRate/1000000/w_EN #convert micro g N gDW–1 h–1 to mol N gDW–1 h–1

#conversions 18C
EC1983_18C_Nuptake_StM$N <- EC1983_18C_Nuptake_StM$ResidualNitrateConcentration
EC1983_18C_Nuptake_StM$N <- round(EC1983_18C_Nuptake_StM$N, digits = 2)
EC1983_18C_Nuptake_StM$N <- EC1983_18C_Nuptake_StM$N/1000000 #microM to M
EC1983_18C_Nuptake_StM$NuptakeRate <- EC1983_18C_Nuptake_StM$NuptakeRate/1000000/w_EN

#testing rounding
sol_EspinozaChapman1983_N_9$N <- round(sol_EspinozaChapman1983_N_9$N*1000000, digits = 3)/1000000
sol_EspinozaChapman1983_N_18$N <- round(sol_EspinozaChapman1983_N_18$N*1000000, digits = 3)/1000000

N_calibration <- ggplot() +
  geom_line(data = sol_EspinozaChapman1983_N_9, mapping = aes(x = N*1000000, y = J_EN_A*1000000, color = "Model of Espinoza and Chapman (1983) at 9°C")) +
  geom_line(data = sol_EspinozaChapman1983_N_18, mapping = aes(x = N*1000000, y = J_EN_A*1000000, color = "Model of Espinoza and Chapman (1983) at 18°C")) +
  geom_point(data = EC1983_9C_Nuptake_StM, mapping = aes(x = N*1000000, y = NuptakeRate*1000000, color="Espinoza and Chapman (1983), St. Margaret's Bay, 9°C"), size=3) +
  geom_point(data = EC1983_18C_Nuptake_StM, mapping = aes(x = N*1000000, y = NuptakeRate*1000000, color="Espinoza and Chapman (1983), St. Margaret's Bay, 18°C"), shape = 23, fill = 'grey', size=3) +
  xlim(0, 80) +
  scale_color_manual(values = c("gray60", "gray0", "gray60", "gray0")) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  labs(x= bquote('N Concentration (μmol' ~NO[3]^{"-"}~ 'L'^"-1"*')'), y = bquote('N uptake (μmol' ~NO[3]^{"-"}~ 'g DW'^"-1"*' h'^"-1"*')')) +
  ggtitle('a)')

#Error calculations
er9 <- merge(EC1983_9C_Nuptake_StM, sol_EspinozaChapman1983_N_9, all.x = TRUE)
rmse(er9$NuptakeRate, er9$J_EN_A) #3.683799e-07

er18 <- merge(EC1983_18C_Nuptake_StM, sol_EspinozaChapman1983_N_18, all.x = TRUE)
rmse(er18$NuptakeRate, er18$J_EN_A) #2.606024e-07

######## Photosynthesis related ####
#Johansson2002
Johansson2002 <- read.csv("Johansson2002.csv", header = TRUE, fileEncoding="UTF-8-BOM")
#conversions
Johansson2002$Irradiance <- Johansson2002$Irradiance*3600*1e-6 #micromol photons m-2 s-1 to mol photons m-2 h-1
Johansson2002$O2production <- Johansson2002$O2production/1e+6*32/1000*3600 #micromol O2 kg DW-1 s-1 to g O2/g/h
Johansson2002$O2productionSHIFT <- Johansson2002$O2production + 0.001720976 #from net to gross

Photosynthesis_calibration <- ggplot(data = Johansson2002) +
  geom_line(data = sol_Johansson2002, mapping = aes(x = I, y = J_O*1000)) +
  geom_point(mapping = aes(x = Irradiance, y = O2productionSHIFT*1000), size = 3) +
  scale_color_grey() +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  labs(x= bquote('Irradiance (mol γ m'^"-2"*' h'^"-1)"), y = bquote('Oxygen production (mg' ~O[2]~ 'g DW'^"-1"*' h'^"-1"*')')) +
  ggtitle('b)')

#error calculations
Johansson2002$I <- round(Johansson2002$Irradiance, digits = 6)
sol_Johansson2002$I <- round(sol_Johansson2002$I, digits = 6)
erPhoto <- merge(Johansson2002, sol_Johansson2002, all.x = TRUE)
rmse(erPhoto$O2productionSHIFT, erPhoto$J_O)

######## Combine calibration plot (Figure 5) #######
#Figure 5
grid.arrange(N_calibration, Photosynthesis_calibration, ncol=2)