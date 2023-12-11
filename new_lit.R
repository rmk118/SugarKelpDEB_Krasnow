# testing new literature

new_lit_df <- read_csv("new_lit.csv", col_types = "cdddccc") %>% mutate(temp_K=temp+273.15)

gd1988high <- tibble(temp=c(8,15,18,20),
                    paper="gd1988", level="high",
                    rate=c(4.70, 3.97, 3.45, 0.91)) %>% mutate(std_rate=rate/0.91, temp_K=temp+273.15)

p2015 <- tibble(temp=c(5,10,15,17.5,20,5,10,15,17.5,20,22.5 ), paper="p2015", level="high", std_rate=c(2.735,2.265,2.13,1.515,1, 3.302,3.597,3.268,2.151,1,0.008), temp_K=temp+273.15)
#nd2019 <- tibble(temp=c(15,18,21), paper="nd2019", level="high", std_rate=c(11.65,5.53,1), temp_K=temp+273.15)

# low_lit_df <- tibble(temp=c(8,15,18,20),
#                      paper="gd1988",
#                      level="low",
#                      rate=c(3.84, 3.66, 2.99, 0.91)) %>% mutate(std_rate=rate/0.91, temp_K=temp+273.15)

low_lit_df <- new_lit_df %>% filter(str_detect(paper_full, "Bergen")) %>% filter(str_detect(paper_full, "_20_", TRUE))

all_lit_data_new <- #bind_rows(lit_data %>% select(-c(ctrl_group, stress_group)) %>% mutate(level="lit"),
                             # new_lit_df, gd1988high)
  bind_rows(lit_data %>% select(-c(ctrl_group, stress_group)) %>% mutate(level="lit"),
                              low_lit_df, p2015)

ggplot()+geom_point(data=all_lit_data_new, aes(x=temp, y=std_rate, color=level))


nls_fun_simple <- function(df, x) {
  #level can be "high" or "low"
  #df should be a data frame that includes the original literature values (in order to fit the lower portion of the curve), along with the new points you want to include in the nls fitting
  
  data <- df %>% ungroup() %>% filter(level=="lit" | level==x)


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

  df_out <- data.frame(T_A=temp_T_A, T_H=temp_T_H, T_AH=temp_T_AH,level=x, temp_smooth=temp_smooth, std_rate=vector_smooth)
  df_out
  #data
}

new_lit <- c(high="high", low="low") %>% map(~nls_fun_simple(all_lit_data_new, .x)) %>% 
 bind_rows(.id="level2")




params_new_lit<-new_lit %>% select(T_A, T_H, T_AH,level) %>% distinct() %>% 
  add_row(T_A=T_A, T_H=T_H,T_AH=T_AH,level="orig")%>% 
  add_row(T_A=new_T_A, T_H=new_T_H,T_AH=new_T_AH, level="lit")
params_new_lit


high_nls<- nlsLM(formula=std_rate ~ exp((T_a/T_0)-(T_a/temp_K)) *
        (1+exp((T_AL/T_0)-(T_AL/T_L)) + exp((T_ah/T_h)-(T_ah/T_0))) *
        ((1+exp((T_AL/temp_K)-(T_AL/T_L))+exp((T_ah/T_h)-(T_ah/temp_K)))^-1),
      start=list(T_a=T_a, T_h=T_h, T_ah=T_ah),
      lower=c(6000,285.15,12000),
      upper=c(12500,294.15,25000),
      data = all_lit_data_new %>% filter(level!="low"))
summary(high_nls)

low_nls<-nlsLM(formula=std_rate ~ exp((T_a/T_0)-(T_a/temp_K)) *
        (1+exp((T_AL/T_0)-(T_AL/T_L)) + exp((T_ah/T_h)-(T_ah/T_0))) *
        ((1+exp((T_AL/temp_K)-(T_AL/T_L))+exp((T_ah/T_h)-(T_ah/temp_K)))^-1),
      start=list(T_a=T_a, T_h=T_h, T_ah=T_ah),
      lower=c(6000,282.15,12000),
      upper=c(12500,290.15,25000),
      data = all_lit_data_new %>% filter(level!="high"))
summary(low_nls)

ggplot()+
  # geom_line(data=new_lit, aes(x=temp_smooth-273.15, y=std_rate, color=level), linewidth=2)+
  geom_point(data=all_lit_data_new, aes(x=temp, y=std_rate, color=level))+
  #geom_smooth(data=all_lit_data_new, aes(x=temp, y=std_rate, color=level))+
  theme_bw()+
  #ylim(0,5)+
  labs(x="Temperature (°C)", y="Standardized rate", color=NULL)

ggplot()+
  geom_line(data=new_lit, aes(x=temp_smooth-273.15, y=std_rate, color=level), linewidth=2)+
  geom_point(data=all_lit_data_new, aes(x=temp, y=std_rate, color=level))+
  annotate(geom='line', x=temp_smooth-273.15,y=y_smooth_venolia)+
  theme_bw()+
  #ylim(0,5)+
  labs(x="Temperature (°C)", y="Standardized rate", color=NULL)

