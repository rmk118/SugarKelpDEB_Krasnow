#Recalibration of Arrhenius temperature response curve for DEB kelp model developed by Venolia et al., 2020
#A complete description of the Venolia et al., (2020) model is 
# available at https://doi.org/10.1016/j.ecolmodel.2020.109151

# This code was written by Ruby Krasnow between November-December 2023
# Last updated: Dec 30, 2023

#Packages needed
library(tidyverse)
library(patchwork)
library(lubridate)
library(tseries)
library(minpack.lm)
library(gt)

#### Original literature data (cited in Venolia et al., 2020) ###############################################################

lit_data <- read.csv("orig_lit_data.csv") %>% 
  select(c(paper, paper_full, temp, rate, std_rate, extra_info)) %>% 
  filter(temp!=23) %>% 
  mutate(extra_info = if_else(paper!="bl182" & extra_info=="Germany", NA, extra_info),
         paper_full=as_factor(paper_full)) %>% 
  mutate(paper_full = fct_relevel(paper_full, c("bl1982_France", "bl1982_Norway", "bl1982_Germany", "bl1982_UK", "fl1980","dd1987",  "d1987_0", "d1987_5", "d1987_10", "d1987_15", "d1987_20")),
         temp_K = temp+273.15,
         level="lit") 

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

lit_data_fit <- lit_data %>% 
  mutate(y=exp((T_A/T_0)-(T_A/temp_K)) * (1+exp((T_AL/T_0)-(T_AL/T_L)) + exp((T_AH/T_H)-(T_AH/T_0))) * ((1+exp((T_AL/temp_K)-(T_AL/T_L))+exp((T_AH/T_H)-(T_AH/temp_K)))^-1))

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

# Import new data, from studies not used in the original paper
new_lit_df <- read_csv("new_lit.csv", col_types = "cdddccc") %>% #specify column types
  mutate(temp_K=temp+273.15, #add column with temperature in Kelvin
         level = case_when(level=="high"~ "warm", 
                           level=="low" ~"cold",
                           .default=level)) %>% 
  na.omit()

all_lit_data_new <- bind_rows(lit_data, new_lit_df)
ggplot()+geom_point(data=all_lit_data_new, aes(x=temp, y=std_rate, color=level))

# NLS -------------------------------------------------------------------------

nls_fun_simple <- function(df, x) {
  #level can be "warm" or "cold"
  #df should be a data frame that includes the original literature values (in order to fit the colder portion of the curve), along with the new points you want to include in the nls fitting
  
  data <- df %>% ungroup() %>% filter(level=="lit" | level==x)

  nls_res <- nlsLM(formula=std_rate ~ exp((T_a/T_0)-(T_a/temp_K)) *
                     (1+exp((T_AL/T_0)-(T_AL/T_L)) + exp((T_ah/T_h)-(T_ah/T_0))) *
                     ((1+exp((T_AL/temp_K)-(T_AL/T_L))+exp((T_ah/T_h)-(T_ah/temp_K)))^-1),
                   start=list(T_a=T_a, T_h=T_h, T_ah=T_ah),
                   lower=c(6000,282.15,12000),
                   upper=c(12500,292.15,25000),
                   data = data)

  temp_T_A <- nls_res$m$getPars()[[1]]
  temp_T_H <- nls_res$m$getPars()[[2]]
  temp_T_AH <- nls_res$m$getPars()[[3]]

  vector_smooth <- exp((temp_T_A/T_0)-(temp_T_A/temp_smooth))*(1+exp((T_AL/T_0)-(T_AL/T_L))+exp((temp_T_AH/temp_T_H)-(temp_T_AH/T_0))) * ((1+exp((T_AL/temp_smooth)-(T_AL/T_L))+exp((temp_T_AH/temp_T_H)-(temp_T_AH/temp_smooth)))^-1)

  df_out <- data.frame(T_A=temp_T_A, T_H=temp_T_H, T_AH=temp_T_AH,level=x, temp_smooth=temp_smooth, std_rate=vector_smooth)
  df_out
}

new_lit <- c(warm="warm", cold="cold") %>% 
  map(~nls_fun_simple(all_lit_data_new, .x)) %>% 
  bind_rows(.id="level2")

params_new_lit <- new_lit %>% select(T_A, T_H, T_AH,level) %>% distinct() %>% 
  add_row(T_A=T_A, T_H=T_H,T_AH=T_AH,level="orig")%>% 
  add_row(T_A=new_T_A, T_H=new_T_H,T_AH=new_T_AH, level="lit")

params_nested <- params_new_lit %>% nest(data = c(T_A, T_H, T_AH))

ggplot()+
  geom_point(data=all_lit_data_new, aes(x=temp, y=std_rate, color=level), size=1)+
  geom_line(data=new_lit, aes(x=temp_smooth-273.15, y=std_rate, color=level), linewidth=1.2)+
  geom_line(aes(x=temp_smooth-273.15, y=y_smooth_venolia, color="original"), linewidth=1.2 )+
  theme_classic()+
  labs(x="Temperature (°C)", y="Standardized rate", color=NULL)+
  scale_color_manual(values=c("cold"='#0f85a0',"warm"="#dd4124",  "lit"="black","original"="black"),
                     breaks=c("cold","warm","original"),
                     labels=c("warm"="Warm", "cold"="Cold" ,"original"="Original"))+
  theme(text = element_text(size=15),
        axis.title.y = element_text(margin = margin(t = 0, r = 9, b = 0, l = 0), face="bold"),
        axis.title.x = element_text(margin = margin(t = 9, r = 0, b = 0, l = 0), face="bold"))


gt(params_new_lit %>% mutate(T_H=T_H-273.15)) %>% 
  fmt_number(columns = "T_H", decimals = 2) %>% 
  fmt_number(columns = c("T_AH", "T_A"), decimals = 0) %>% 
  fmt(columns = level, fns=str_to_title) %>% 
  cols_move_to_start(level) %>% 
  cols_label(T_A = html("T<sub>A</sub> (K)"),
             T_H = html("T<sub>H</sub> (°C)"),
             T_AH = html("T<sub>AH</sub> (K)"),
             level="") %>% 
  cols_width(everything() ~ px(80)) 
#%>%
  #gtsave(filename="params.docx")

