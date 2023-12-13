# testing new literature

new_lit_df <- read_csv("new_lit.csv", col_types = "cdddccc") %>% 
  mutate(temp_K=temp+273.15, level = case_when(level=="high"~ "warm", 
                                               level=="low" ~"cold",
                                               .default=level)) %>% 
  na.omit()


all_lit_data_new <- bind_rows(lit_data %>% select(-c(ctrl_group, stress_group)) %>% mutate(level="lit"),
                              new_lit_df)

ggplot()+geom_point(data=all_lit_data_new, aes(x=temp, y=std_rate, color=level))

# -------------------------------------------------------------------------



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

# -------------------------------------------------------------------------


# -------------------------------------------------------------------------


  temp_T_A <- nls_res$m$getPars()[[1]]
  temp_T_H <- nls_res$m$getPars()[[2]]
  temp_T_AH <- nls_res$m$getPars()[[3]]

  vector_smooth <- exp((temp_T_A/T_0)-(temp_T_A/temp_smooth))*(1+exp((T_AL/T_0)-(T_AL/T_L))+exp((temp_T_AH/temp_T_H)-(temp_T_AH/T_0))) * ((1+exp((T_AL/temp_smooth)-(T_AL/T_L))+exp((temp_T_AH/temp_T_H)-(temp_T_AH/temp_smooth)))^-1)

  df_out <- data.frame(T_A=temp_T_A, T_H=temp_T_H, T_AH=temp_T_AH,level=x, temp_smooth=temp_smooth, std_rate=vector_smooth)
  df_out
  #data
}

new_lit <- c(warm="warm", cold="cold") %>% map(~nls_fun_simple(all_lit_data_new, .x)) %>% 
 bind_rows(.id="level2")

params_new_lit<-new_lit %>% select(T_A, T_H, T_AH,level) %>% distinct() %>% 
  add_row(T_A=T_A, T_H=T_H,T_AH=T_AH,level="orig")%>% 
  add_row(T_A=new_T_A, T_H=new_T_H,T_AH=new_T_AH, level="lit")
params_new_lit

ggplot()+
  # geom_line(data=new_lit, aes(x=temp_smooth-273.15, y=std_rate, color=level), linewidth=2)+
  geom_point(data=all_lit_data_new, aes(x=temp, y=std_rate, color=level))+
  #geom_smooth(data=all_lit_data_new, aes(x=temp, y=std_rate, color=level))+
  theme_bw()+
  #ylim(0,5)+
  labs(x="Temperature (°C)", y="Standardized rate", color=NULL)

ggplot()+
  geom_line(data=new_lit, aes(x=temp_smooth-273.15, y=std_rate, color=level), linewidth=1.5)+
  geom_point(data=all_lit_data_new, aes(x=temp, y=std_rate, color=level), size=1)+
  geom_line(aes(x=temp_smooth-273.15, y=y_smooth_venolia, color="original"), linewidth=1 )+
  theme_classic()+
  labs(x="Temperature (°C)", y="Standardized rate", color=NULL)+
  scale_color_manual(values=c("cold"='#0f85a0',"warm"="#dd4124",  "lit"="#edd746","original"="black"),
                     breaks=c("cold","warm",  "lit","original"),
                     labels=c("warm"="Warm", "cold"="Cold" ,"lit"="Literature","original"="Original"))+
  theme(text = element_text(size=18),
        axis.title.y = element_text(margin = margin(t = 0, r = 9, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 9, r = 0, b = 0, l = 0)))

