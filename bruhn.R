# Bruhn et al., 2016 model comparison

# BRUHN 2016 --------------------------------------------------------------
# Stations overview and timeline:
# 1. Deployment date: Dec 6, 2011, Sampling 2 date: Jun 12, 2012
# - Odby Bay = monitoring station VIB3702
# - Riisgarde Broad = monitoring station VIB3727

# 2. Deployment date: Oct 28, 2011, Sampling 2 date: May 25, 2012
# - Færker Vig = monitoring station VIB3708

# Sampling 1 date was Apr 11, 2012 for all 3 stations
# Lysen Broad and Fur Sund are not included because of their distance from the monitoring stations

###### Time steps ######
# Experimental period for Odby Bay (OB) and Riisgarde Broad (RB) was Dec 6, 2011 to Jun 12, 2012 
# (growth period of 190 days)
hourly_seq_bruhn1 <- seq(ymd_hms("2011-12-06 00:00:00"),ymd_hms("2012-06-12 00:00:00"), by="hour") # hourly sequence of POSIXct values
times_bruhn1 <- seq(0, 189*24, 1) #190 days stepped hourly

# Experimental period for Færker Vig (FV) was Oct 28, 2011 to May 25, 2012 
# (growth period of 211 days)
hourly_seq_bruhn2 <- seq(ymd_hms("2011-10-28 00:00:00"),ymd_hms("2012-05-25 00:00:00"), by="hour") # hourly sequence of POSIXct values
times_bruhn2 <- seq(0, 210*24, 1) #211 days stepped hourly

###### Environmental forcing data #####
stations_list <- c(VIB3702="VIB3702",VIB3708="VIB3708", VIB3727="VIB3727")

# Temperature
temp_filepaths <- stations_list %>% map_chr(\(x) paste0("./validation_data/bruhn2016/", x, "_temp.csv"))

bruhn_temp <- read_csv(temp_filepaths, id = "station", col_names = c("date", "temp"), col_types = "cd") %>% 
  mutate(date = ymd_hms(paste0(date, " 00:00:00")))  %>%
  mutate(station = tools::file_path_sans_ext(basename(station)),
         station = as_factor(str_sub(station,1,7)))

# Nitrate
N_filepaths <- stations_list %>% map_chr(\(x) paste0("./validation_data/bruhn2016/", x, "_N.csv"))

bruhn_N <- read_csv(N_filepaths, id = "station", col_names = c("date", "N"), col_types = "cd") %>% 
  mutate(date = ymd_hms(paste0(date, " 00:00:00")))  %>%
  mutate(station = tools::file_path_sans_ext(basename(station)),
         station = as_factor(str_sub(station,1,7)),
         N = N/10^6) #N is imported in µmol/L, the model needs it in mol/L

tempN_data_bruhn1 <- bruhn_temp %>% 
  full_join(bruhn_N) %>% 
  filter(station %in% c("VIB3702", "VIB3727")) %>% 
  mutate(station=fct_drop(station)) %>% 
  full_join(data.frame(date=hourly_seq_bruhn1)) %>% 
  arrange(date) %>% 
  complete(date, station) %>% 
  filter(!is.na(station)) %>%
  group_by(station) %>% 
  mutate(across(c(temp, N), ~na.approx(.x, rule=2))) %>% 
  filter(date %in% hourly_seq_bruhn1) %>% mutate(temp_K=temp+273.15)

tempN_data_bruhn2 <- bruhn_temp %>% 
  full_join(bruhn_N) %>% 
  filter(station=="VIB3708") %>% 
  full_join(data.frame(date=hourly_seq_bruhn2)) %>% 
  arrange(date) %>%
  mutate(station=as_factor("VIB3708")) %>% 
  mutate(across(c(temp, N), ~na.approx(.x, rule=2))) %>% 
  filter(date %in% hourly_seq_bruhn2) %>% mutate(temp_K=temp+273.15)

# PAR
PAR_filepaths_1.5 <- stations_list %>% map_chr(\(x) paste0("./validation_data/bruhn2016/", x, "_1.5.csv"))

bruhn_PAR_1.5 <- read_csv(PAR_filepaths_1.5, id = "station", col_names = c("date", "PAR"), col_types = "cd") %>% 
  mutate(date = ymd_hms(paste0(date, " 00:00:00")))  %>%
  mutate(station = as.factor(tools::file_path_sans_ext(basename(station))))

PAR_filepaths_2.5 <- stations_list %>% map_chr(\(x) paste0("./validation_data/bruhn2016/", x, "_2.5.csv"))

bruhn_PAR_2.5 <- read_csv(PAR_filepaths_2.5, id = "station", col_names = c("date", "PAR"), col_types = "cd") %>% 
  mutate(date = ymd_hms(paste0(date, " 00:00:00")))  %>%
  mutate(station = as_factor(tools::file_path_sans_ext(basename(station))))

bruhn_PAR <- bruhn_PAR_1.5 %>% bind_rows(bruhn_PAR_2.5) %>% 
  mutate(depth = as_factor(str_split_i(station, "_", 2)),
         station = as_factor(str_split_i(station, "_", 1)))

PAR_data_bruhn1 <- bruhn_PAR %>% 
  filter(station %in% c("VIB3702", "VIB3727")) %>% 
  mutate(station=fct_drop(station)) %>% 
  full_join(data.frame(date=hourly_seq_bruhn1)) %>% 
  arrange(date) %>% 
  complete(date, station, depth)  %>% filter(!is.na(station) & !is.na(depth)) %>%
  group_by(depth, station) %>% 
  mutate(across(c(PAR), ~na.approx(.x, rule=2))) %>% 
  filter(date %in% hourly_seq_bruhn1) %>% 
  mutate(PAR_mol = PAR*3600*1e-6)

PAR_data_bruhn2 <- bruhn_PAR %>% 
  filter(station=="VIB3708") %>% 
  full_join(data.frame(date=hourly_seq_bruhn2)) %>% 
  arrange(date) %>% 
  complete(date, depth)  %>% mutate(station=as_factor("VIB3708")) %>% 
  filter(!is.na(depth)) %>%
  group_by(depth, station) %>% 
  mutate(across(c(PAR), ~na.approx(.x, rule=2))) %>% 
  filter(date %in% hourly_seq_bruhn2) %>% 
  mutate(PAR_mol = PAR*3600*1e-6)

###### Initial conditions #####
# The kelp started the acclimation period with an avg length of 1mm, using the original length-weight relationship we can find W=0.00014
W <- 0.00014 #initial biomass for conversions
state_bruhn <- c(m_EC = 0.002, m_EN = 0.001,
                 M_V = W/(w_V+0.001*w_EN+0.002*w_EC))

###### Model runs #####
####### ~ Odby Bay #####
T_field <- approxfun(x = times_bruhn1, y = tempN_data_bruhn1 %>% filter(station=="VIB3702") %>% pull(temp_K), method = "linear", rule = 2)

N_field <- approxfun(x = times_bruhn1, y = (tempN_data_bruhn1 %>% filter(station=="VIB3702") %>% pull(N)), method = "linear", rule = 2) 

####### ~ 1.5 m 
I_field <- approxfun(x = times_bruhn1, y = PAR_data_bruhn1 %>% filter(station=="VIB3702" & depth=="1.5") %>% pull(PAR_mol), method = "linear", rule = 2) 

VIB3702_1.5_output <- params_nested %>% mutate(std_L = future_map(data, function(df) {
  temp_params <- params_Lo
  temp_params[c("T_A", "T_H", "T_AH")] <- c(df$T_A, df$T_H, df$T_AH)
  ode_output <- ode(y = state_bruhn, t = times_bruhn1, func = rates_Lo, parms = temp_params)
  ode_output <- as.data.frame(ode_output) #convert deSolve output into data frame
  ode_output }))  %>% 
  select(-data) %>%
  unnest(cols=std_L) %>%
  ungroup() %>%
  group_by(level) %>%
  mutate(date=hourly_seq_bruhn1, station="VIB3702", depth=1.5, #add date, station, and depth columns
         .before=1) #putting new columns at the beginning for readability

####### ~ 2.5 m 
I_field <- approxfun(x = times_bruhn1, y = PAR_data_bruhn1 %>% filter(station=="VIB3702" & depth=="2.5") %>% pull(PAR_mol), method = "linear", rule = 2) 

VIB3702_2.5_output <- params_nested %>% mutate(std_L = future_map(data, function(df) {
  temp_params <- params_Lo
  temp_params[c("T_A", "T_H", "T_AH")] <- c(df$T_A, df$T_H, df$T_AH)
  ode_output <- ode(y = state_bruhn, t = times_bruhn1, func = rates_Lo, parms = temp_params)
  ode_output <- as.data.frame(ode_output) #convert deSolve output into data frame
  ode_output }))  %>% 
  select(-data) %>%
  unnest(cols=std_L) %>%
  ungroup() %>%
  group_by(level) %>%
  mutate(date=hourly_seq_bruhn1, station="VIB3702", depth=2.5, #add date, station, and depth columns
         .before=1) #putting new columns at the beginning for readability

######## ~ Riisgarde Broad #####
T_field <- approxfun(x = times_bruhn1, y = tempN_data_bruhn1 %>% filter(station=="VIB3727") %>% pull(temp_K), method = "linear", rule = 2)

N_field <- approxfun(x = times_bruhn1, y = (tempN_data_bruhn1 %>% filter(station=="VIB3727") %>% pull(N)), method = "linear", rule = 2) 

####### ~ 1.5 m 
I_field <- approxfun(x = times_bruhn1, y = PAR_data_bruhn1 %>% filter(station=="VIB3727" & depth=="1.5") %>% pull(PAR_mol), method = "linear", rule = 2) 

VIB3727_1.5_output <- params_nested %>% mutate(std_L = future_map(data, function(df) {
  temp_params <- params_Lo
  temp_params[c("T_A", "T_H", "T_AH")] <- c(df$T_A, df$T_H, df$T_AH)
  ode_output <- ode(y = state_bruhn, t = times_bruhn1, func = rates_Lo, parms = temp_params)
  ode_output <- as.data.frame(ode_output) #convert deSolve output into data frame
  ode_output }))  %>% 
  select(-data) %>%
  unnest(cols=std_L) %>%
  ungroup() %>%
  group_by(level) %>%
  mutate(date=hourly_seq_bruhn1, station="VIB3727", depth=1.5, #add date, station, and depth columns
         .before=1) #putting new columns at the beginning for readability

####### ~ 2.5 m 
I_field <- approxfun(x = times_bruhn1, y = PAR_data_bruhn1 %>% filter(station=="VIB3727" & depth=="2.5") %>% pull(PAR_mol), method = "linear", rule = 2) 

VIB3727_2.5_output <- params_nested %>% mutate(std_L = future_map(data, function(df) {
  temp_params <- params_Lo
  temp_params[c("T_A", "T_H", "T_AH")] <- c(df$T_A, df$T_H, df$T_AH)
  ode_output <- ode(y = state_bruhn, t = times_bruhn1, func = rates_Lo, parms = temp_params)
  ode_output <- as.data.frame(ode_output) #convert deSolve output into data frame
  ode_output }))  %>% 
  select(-data) %>%
  unnest(cols=std_L) %>%
  ungroup() %>%
  group_by(level) %>%
  mutate(date=hourly_seq_bruhn1, station="VIB3727", depth=2.5, #add date, station, and depth columns
         .before=1) #putting new columns at the beginning for readability

######## ~ Færker Vig #####
T_field <- approxfun(x = times_bruhn2, y = tempN_data_bruhn2 %>% filter(station=="VIB3708") %>% pull(temp_K), method = "linear", rule = 2)

N_field <- approxfun(x = times_bruhn2, y = (tempN_data_bruhn2 %>% filter(station=="VIB3708") %>% pull(N)), method = "linear", rule = 2) 

####### ~ 1.5 m 
I_field <- approxfun(x = times_bruhn2, y = PAR_data_bruhn2 %>% filter(station=="VIB3708" & depth=="1.5") %>% pull(PAR_mol), method = "linear", rule = 2) 

VIB3708_1.5_output <- params_nested %>% mutate(std_L = future_map(data, function(df) {
  temp_params <- params_Lo
  temp_params[c("T_A", "T_H", "T_AH")] <- c(df$T_A, df$T_H, df$T_AH)
  ode_output <- ode(y = state_bruhn, t = times_bruhn2, func = rates_Lo, parms = temp_params)
  ode_output <- as.data.frame(ode_output) #convert deSolve output into data frame
  ode_output }))  %>% 
  select(-data) %>%
  unnest(cols=std_L) %>%
  ungroup() %>%
  group_by(level) %>%
  mutate(date=hourly_seq_bruhn2, station="VIB3708", depth=1.5, #add date, station, and depth columns
         .before=1) #putting new columns at the beginning for readability

####### ~ 2.5 m 
I_field <- approxfun(x = times_bruhn2, y = PAR_data_bruhn2 %>% filter(station=="VIB3708" & depth=="2.5") %>% pull(PAR_mol), method = "linear", rule = 2) 

VIB3708_2.5_output <- params_nested %>% mutate(std_L = future_map(data, function(df) {
  temp_params <- params_Lo
  temp_params[c("T_A", "T_H", "T_AH")] <- c(df$T_A, df$T_H, df$T_AH)
  ode_output <- ode(y = state_bruhn, t = times_bruhn2, func = rates_Lo, parms = temp_params)
  ode_output <- as.data.frame(ode_output) #convert deSolve output into data frame
  ode_output }))  %>% 
  select(-data) %>%
  unnest(cols=std_L) %>%
  ungroup() %>%
  group_by(level) %>%
  mutate(date=hourly_seq_bruhn2, station="VIB3708", depth=2.5, #add date, station, and depth columns
         .before=1) #putting new columns at the beginning for readability

###### Import observed data and combine with model output #####
bruhn_growth<- read_csv("./validation_data/bruhn2016/bruhn_length.csv",col_names = c("label", "len"), col_types = "cd") %>% 
  mutate(station=case_when(str_starts(label,"ob")~"VIB3702",
                           str_starts(label,"fv")~"VIB3708",
                           str_starts(label,"rb")~"VIB3727"),
         depth=case_when(str_detect(label,"1.5")~1.5,
                         str_detect(label,"2.5")~2.5),
         date = case_when(str_ends(label,"april")~ymd_hms("2012-04-11 00:00:00"),
                          str_ends(label,"june")~ymd_hms("2012-06-12 00:00:00"),
                          str_ends(label,"may")~ymd_hms("2012-05-25 00:00:00")))


bruhn_output <- VIB3702_2.5_output %>% bind_rows(VIB3702_1.5_output, VIB3727_1.5_output, VIB3727_2.5_output, VIB3708_1.5_output, VIB3708_2.5_output) %>% 
  filter(level!="lit")

bruhn_samples <- bruhn_output %>% filter(date %in% as_datetime(c("2012-04-11", "2012-05-25", "2012-06-12"))) %>% 
  select(date, station, depth, level, W, L_allometric) %>% 
  rename(model=L_allometric) 

bruhn_samples <- bruhn_growth %>% left_join(bruhn_samples) %>% 
  pivot_longer(cols=c(len, model), names_to = "type", values_to = "length") %>% 
  mutate(sample_time = if_else(date==ymd_hms("2012-05-25 00:00:00") | date==ymd_hms("2012-06-12 00:00:00"),"june", "april"))

###### SGR #####

bruhn_sgr<- read_csv("./validation_data/bruhn2016/bruhnSGR.csv", col_types = "cdc") %>% 
  mutate(station=case_when(str_starts(label,"ob")~"VIB3702",
                           str_starts(label,"fv")~"VIB3708",
                           str_starts(label,"rb")~"VIB3727"),
         depth=case_when(str_detect(label,"1.5")~1.5,
                         str_detect(label,"2.5")~2.5),
         date = case_when(str_ends(label,"april")~ymd_hms("2012-04-11 00:00:00"),
                          str_ends(label,"june")~ymd_hms("2012-06-12 00:00:00"),
                          str_ends(label,"may")~ymd_hms("2012-05-25 00:00:00")))

bruhn_model_sgr_1<- bruhn_samples %>% filter(date==ymd_hms("2012-04-11 00:00:00"), type=="model") %>% 
  mutate(model_SGR=case_when(
    station=="VIB3708" & date==ymd_hms("2012-04-11 00:00:00") ~ 100*log(W/0.00014)/166,
    station!="VIB3708" & date==ymd_hms("2012-04-11 00:00:00") ~100*log(W/0.00014)/127))

bruhn_model_sgr_2<- bruhn_samples %>% filter(type=="model") %>% 
  arrange(date) %>% 
  group_by(station, depth, level) %>% 
  mutate(numerator = log(W/lag(W)),
         time_elapsed=as.numeric(as.duration(date-lag(date)), "days"),
         model_SGR = 100*numerator/time_elapsed) %>% na.omit()

bruhn_sgr_wide <- bruhn_sgr %>% pivot_wider(names_from = "var", values_from = "SGR")  %>% 
  mutate(sample_time = if_else(date==ymd_hms("2012-05-25 00:00:00") | date==ymd_hms("2012-06-12 00:00:00"),"june", "april"))

###### Figures #####

#ggplot(data=bruhn_sgr_wide %>% filter(date==ymd_hms("2012-04-11 00:00:00")))+ 
ggplot(data=bruhn_sgr_wide)+ 
  geom_bar(aes(x=station, y=mean, fill=station), stat="identity")+
  geom_errorbar(aes(x=station, ymin=low, ymax=high), width=0.4)+
  labs(x="Date", y="SGR")+
  theme_classic()+
  facet_grid(sample_time~depth)+
  geom_point(data=bruhn_model_sgr_1, aes(x=station, y=model_SGR, shape=level))+
  geom_point(data=bruhn_model_sgr_2, aes(x=station, y=model_SGR, shape=level))+
  theme(text = element_text(size=18),
        axis.title.y = element_text(margin = margin(t = 0, r = 9, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 9, r = 0, b = 0, l = 0)))

### Growth figure
ggplot()+ 
  geom_point(data=bruhn_growth %>% filter(date<=ymd_hms("2012-04-11 00:00:00")), aes(x=date, y=len))+
  geom_line(data= bruhn_output %>% filter(date<=ymd_hms("2012-04-11 00:00:00")), aes(x=date, y=L_allometric, color=level), linewidth=1)+
  geom_line(data= bruhn_output %>% filter(date<=ymd_hms("2012-04-11 00:00:00"), level=="warm"), aes(x=date, y=(W*1000/0.729)^(1/1.911)), linewidth=0.5, linetype="dotted")+
  labs(x="Date", y="Kelp frond length (cm)", color=NULL)+
  theme_classic()+
  facet_grid(depth~station, scales="free")+
  theme(text = element_text(size=18),
        axis.title.y = element_text(margin = margin(t = 0, r = 9, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 9, r = 0, b = 0, l = 0)))+
  ylim(0,NA)

ggplot()+ 
  geom_point(data=bruhn_growth, aes(x=date, y=len))+
  geom_line(data= bruhn_output%>% filter(level=="warm"), aes(x=date, y=L_allometric, color="orig"), linewidth=1)+
  geom_line(data= bruhn_output %>% filter(level=="warm"), aes(x=date, y=(W*1000/0.729)^(1/1.911), color="scrosati"), linewidth=1)+
  geom_line(data= bruhn_output %>% filter(level=="warm"), aes(x=date, y=(W/0.0155)^(1/1.358), color="stagnol"), linewidth=1)+
  labs(x="Date", y="Kelp frond length (cm)", color=NULL)+
  theme_classic()+
  facet_grid(depth~station, scales="free")+
  theme(text = element_text(size=18),
        axis.title.y = element_text(margin = margin(t = 0, r = 9, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 9, r = 0, b = 0, l = 0)))+
  ylim(0,NA)


ggplot()+ 
  geom_point(data=bruhn_growth, aes(x=date, y=len))+
  geom_line(data= bruhn_output %>% filter(level=="warm"), aes(x=date, y=L_allometric, color="orig"), linewidth=1)+
  geom_line(data= bruhn_output %>% filter(level=="warm"), aes(x=date, y=(W*1000/0.729)^(1/1.911), color="scrosati"), linewidth=1)+
  geom_line(data= bruhn_output %>% filter(level=="warm"), aes(x=date, y=(W/0.0155)^(1/1.358), color="stagnol"), linewidth=1)+
  labs(x="Date", y="Kelp frond length (cm)", color=NULL)+
  theme_classic()+
  facet_grid(depth~station, scales="free")+
  theme(text = element_text(size=18),
        axis.title.y = element_text(margin = margin(t = 0, r = 9, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 9, r = 0, b = 0, l = 0)))+
  ylim(0,NA)


### Environmental figures
PAR_plot_bruhn<- ggplot(data=PAR_data_bruhn1 %>% bind_rows(PAR_data_bruhn2))+
  geom_line(aes(x=date, y=PAR, color=station, linetype=depth), linewidth=1)+
  labs(x=NULL, y=expression(paste("PAR (μmol photons ",  m^-2, " ",s^-1, ")")), color="Station", linetype="Depth")+
  theme_classic()+
  theme(text = element_text(size=16),
        axis.title.y = element_text(margin = margin(t = 0, r = 9, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 9, r = 0, b = 0, l = 0)))

N_plot_bruhn<-ggplot(data=tempN_data_bruhn1 %>% bind_rows(tempN_data_bruhn2))+
  geom_line(aes(x=date, y=N*10^6, group=station, color=station), linewidth=1)+
  labs(x=NULL, y=bquote("NO"[3]~" (µM)"), color=NULL)+
  theme_classic()+
  theme(text = element_text(size=16),
        axis.title.y = element_text(margin = margin(t = 0, r = 9, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 9, r = 0, b = 0, l = 0)))+
  guides(color="none")

temp_plot_bruhn<-ggplot(data=tempN_data_bruhn1 %>% bind_rows(tempN_data_bruhn2))+
  geom_line(aes(x=date, y=temp, group=station, color=station), linewidth=1)+
  labs(x=NULL, y="Temperature (°C)", color=NULL, linetype=NULL)+
  theme_classic()+
  theme(text = element_text(size=16),
        axis.title.y = element_text(margin = margin(t = 0, r = 9, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 9, r = 0, b = 0, l = 0)))+
  guides(color="none")

PAR_plot_bruhn+N_plot_bruhn+temp_plot_bruhn+
  plot_layout(guides = 'collect')

