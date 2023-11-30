
cross144 <- c(growth_rates_no_flags %>% filter(cross==144, trt!="control") %>% mutate(rgr= if_else(rgr<0 | is.na(rgr), 0, rgr)) %>% summarise(mean(rgr)))

ref_cross <- growth_rates_no_flags %>% group_by(cross, date, trt) %>%
  filter(growth_hp > 0) %>% 
  summarise(mean_growth = mean(growth_hp, na.rm=TRUE),
            mean_rgr = mean(growth_rate_hp, na.rm=TRUE),
            mean_rgr = mean(rgr, na.rm=TRUE),
            max_rgr = max(rgr, na.rm=TRUE)) %>% add_week_fun() %>% 
  left_join(weekly_means_degC) %>% 
  filter(!(week %in% c(1,5)))  %>% 
  group_by(cross) %>% 
  filter(round(mean_temp,0)==20) %>% 
  summarize(ref_rgr = mean(max_rgr)) %>% 
  select(cross, ref_rgr)  %>%
  add_row(cross=as.factor(144), ref_rgr=cross144[[1]])


unh_cross<-growth_rates_no_flags %>% 
  group_by(cross, date, trt, mean_temp) %>%
  summarize(rgr = mean(rgr, na.rm=TRUE)) %>% 
  left_join(ref_cross) %>% 
  mutate(std_rate = rgr/ref_rgr) %>% 
  ungroup() %>%
  mutate(ctrl_group = case_when(
    as.character(cross) %in% levels(high_ctrl$cross) ~ "high",
    as.character(cross) %in% levels(low_ctrl$cross) ~ "low",
    .default = "med"),
    stress_group = case_when(
      as.character(cross) %in% levels(high_stress$cross) ~ "high",
      as.character(cross) %in% levels(low_stress$cross) ~ "low",
      .default = "med") ) %>%
  mutate(ctrl_group = as_factor(ctrl_group),
         stress_group = as_factor(stress_group)) %>% 
  mutate(ctrl_group = fct_relevel(ctrl_group, c("high", "med", "low")),stress_group = fct_relevel(stress_group, c("high", "med", "low"))) %>% 
  ungroup() %>% 
  mutate(paper="unh_cross2023", paper_full = paste(paper, stress_group, sep="_"), temp=mean_temp, rate=rgr, extra_info=paste(trt, cross,sep="_"), temp_K= temp + 273.15) %>%  #adding columns with identifying info
  select(paper, paper_full, temp, rate,std_rate, extra_info, temp_K, ctrl_group, stress_group)

stress_cross_df <- map(.x=c(high="high", med="med", low="low"), .f=nls_fun, df=unh_cross, type="stress") %>% bind_rows(.id="level") %>% mutate(res="cross")

ctrl_cross_df <- map(.x=c(high="high", med="med", low="low"), .f=nls_fun, df=unh_cross, type="ctrl") %>% bind_rows(.id="level") %>% mutate(res="cross")

all_calibrations_updated <- rbind(stress_means_df, stress_all_df, ctrl_means_df, ctrl_all_df, stress_cross_df, ctrl_cross_df) %>% 
  mutate(level = as_factor(level),
         level = fct_expand(level, "lit"),
         level = fct_relevel(level, c("high", "med", "low", "lit")))


all_calibrations_updated %>% select(T_A, T_H, T_AH,type,level,res) %>% distinct() %>% 
  add_row(T_A=T_A, T_H=T_H,T_AH=T_AH,type="orig", level="orig", res="orig")%>% 
  add_row(T_A=new_T_A, T_H=new_T_H,T_AH=new_T_AH,type="lit", level="lit", res="lit")

ggplot()+
  geom_line(data=all_calibrations_updated %>% group_by(type, level, res), aes(x=temp_smooth-273.15, y=std_rate, color=level), linewidth=2)+
  geom_point(data=all_lit_data, aes(x=temp, y=std_rate, color=level))+
  theme_bw()+
  facet_grid(type~res)+
  #ylim(0,5)+
  labs(x="Temperature (°C)", y="Standardized rate", color=NULL)
