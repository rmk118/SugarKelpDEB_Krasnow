#testing halves instead of thirds

growth_no_zero <- growth_data %>% 
  group_by(id) %>% 
  mutate(time_elapsed = as.integer(date-lag(date)), #days between sampling events
         rgr=((log(hp)-log(lag(hp)))/time_elapsed)*100) %>% #relative growth rate, in % per day
  mutate(rgr=if_else(rgr<0,0,rgr)) %>% na.omit() %>% 
  add_trt_fun() #%>% filter(trt!="lowN")

crosses<-growth_no_zero %>% group_by(cross, date, trt) %>%
  summarise(mean_rgr = mean(rgr, na.rm=TRUE),
            max_rgr = max(rgr, na.rm=TRUE)) %>% add_week_fun() %>% 
  left_join(weekly_means_degC %>% mutate(mean_temp=round(mean_temp,0))) %>% 
  group_by(cross) %>% 
  filter(mean_temp==20) %>% 
  summarize(ref_rgr = max(max_rgr)) %>% 
  select(cross, ref_rgr) %>% 
  filter(ref_rgr>0.7)

mean_growth_by_cross <- growth_no_zero %>% 
  filter(cross %in% c(fct_drop(crosses$cross))) %>% 
  group_by(cross) %>% 
  summarise(gr_hp=mean(rgr)) %>% 
  arrange(-gr_hp) %>% 
  mutate(cross = fct_drop(cross))

top_half <- mean_growth_by_cross %>% slice_max(gr_hp, n=11) %>% mutate(cross = fct_drop(cross))
bottom_half <- mean_growth_by_cross %>% slice_min(gr_hp, n=11) %>% mutate(cross = fct_drop(cross))


growth_rates_halves <- growth_no_zero %>%
  filter(cross %in% c(fct_drop(crosses$cross))) %>% 
  mutate(group = case_when(
    as.character(cross) %in% levels(top_half$cross) ~ "top",
    as.character(cross) %in% levels(bottom_half$cross) ~ "bottom")) %>% 
  ungroup()

# calculate weekly mean and maximum relative growth rates for each half
weekly_halves <- growth_rates_halves %>% 
  group_by(group, trt, date) %>%
  summarise(mean_rgr = mean(rgr, na.rm=TRUE),
            max_rgr = max(rgr, na.rm=TRUE)) %>% 
  add_week_fun() %>% #add week labels
  filter(!(week %in% c(1,5))) #remove weeks without growth rates

ref_half <- weekly_halves %>% 
  group_by(group) %>% 
  left_join(weekly_means_degC) %>% 
  filter(round(mean_temp,0)==20) %>% 
  mutate(ref_rgr = max_rgr) %>% 
 filter(trt!="lowN") %>% 
  select(group, ref_rgr)


# Finding standardized RGR at the group level
unh_half <- weekly_halves %>% 
  left_join(ref_half) %>%  #add a column to the weekly means (stress grouping) with reference RGRs
  left_join(weekly_means_degC) %>% 
  mutate(std_rgr = max_rgr/ref_rgr) %>%  # standardized max RGR (max RGR divided by ref temp RGR)
        # std_rgr_mean = mean_rgr/ref_rgr, # standardized mean RGR (mean RGR divided by ref temp RGR)
        # std_rgr = case_when(week==2 & trt!="control"~ std_rgr_max, #use max for week 2 non-control
                             #.default= std_rgr_mean)) %>% 
  ungroup() %>%
  mutate(paper="unh_half2023", paper_full=paste(paper, group, sep="_"), temp=mean_temp, rate=mean_rgr, extra_info=NA, temp_K=temp+273.15) %>% #adding columns with identifying info
  select(paper, paper_full, temp, rate, std_rgr, extra_info, temp_K, group)


##### Cross level ####################################################################

half_test_cross<-growth_no_zero %>%
  group_by(cross, trt, date) %>%
  summarise(mean_rgr = mean(rgr, na.rm=TRUE)) %>% 
  add_week_fun() %>% #add week labels
  filter(!(week %in% c(1,5))) %>% 
  left_join(crosses) %>% 
  left_join(weekly_means_degC) %>% 
  na.omit() %>% 
  mutate(std_rate=mean_rgr/ref_rgr) %>%
  mutate(group = case_when(
    as.character(cross) %in% levels(top_half$cross) ~ "top",
    as.character(cross) %in% levels(bottom_half$cross) ~ "bottom")) %>%
  ungroup() %>% 
  mutate(paper="unh_cross2023", paper_full = paste(paper, group, sep="_"), temp=mean_temp, rate=mean_rgr, extra_info=paste(trt, cross,sep="_"), temp_K= temp + 273.15) %>%  #adding columns with identifying info
  select( paper, paper_full, temp, rate,std_rate, extra_info, temp_K, group)

#lit_data <- lit_data %>% select(-c(stress_group, ctrl_group)) %>% mutate(group="lit")

#### Combine new and lit data #### 
lit_data_plus_half <- bind_rows(lit_data, unh_half %>% mutate(std_rate=std_rgr, .keep="unused")) %>% mutate(res="means")
lit_data_plus_cross_half <- bind_rows(lit_data, half_test_cross) %>%  mutate(res="cross")


all_lit_data_half <- rbind(lit_data_plus_half, lit_data_plus_cross_half) %>% 
  mutate(level = case_when(str_ends(paper_full, "top") ~ "top",
                           str_ends(paper_full, "bottom") ~ "bottom",
                           .default = "lit"),
         level = as_factor(level),
         level = fct_relevel(level, c("top", "bottom", "lit")))

ggplot()+geom_point(data=all_lit_data_half, aes(x=temp, y=std_rate, color=level))+facet_wrap(~res)

nls_fun_simple <- function(df, level) {
  #level can be "high", "med", or "low"
  #df should be a data frame that includes the original literature values (in order to fit the lower portion of the curve), along with the new points you want to include in the nls fitting
  
  relevant_grouping <- paste("group")
  data <- df %>% filter(.data[[relevant_grouping]]=="lit" | .data[[relevant_grouping]]==level)
  
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
  
  df_out <- data.frame(T_A=temp_T_A, T_H=temp_T_H, T_AH=temp_T_AH,level=level, temp_smooth=temp_smooth, std_rate=vector_smooth)
  df_out
}

#### means
means_half_df <- map(.x=c(top="top", bottom="bottom"), .f=nls_fun_simple, df=lit_data_plus_half) %>% bind_rows(.id="level") %>% mutate(res="means")

#### crosses
cross_half_df <- map(.x=c(top="top", bottom="bottom"), .f=nls_fun_simple, df=lit_data_plus_cross_half) %>% bind_rows(.id="level") %>% mutate(res="cross")


halves_calibrations <- rbind(means_half_df, cross_half_df) %>% 
 # halves_calibrations <- means_half_df %>% 
  mutate(level = as_factor(level),
         level = fct_expand(level, "lit"),
         level = fct_relevel(level, c("top", "bottom", "lit")))

params_halves<-halves_calibrations %>% select(T_A, T_H, T_AH,level,res) %>% distinct() %>% 
  add_row(T_A=T_A, T_H=T_H,T_AH=T_AH,level="orig", res="orig")%>% 
  add_row(T_A=new_T_A, T_H=new_T_H,T_AH=new_T_AH,level="lit", res="lit")

ggplot()+
  geom_line(data=halves_calibrations %>% ungroup(), aes(x=temp_smooth-273.15, y=std_rate, color=level), linewidth=2)+
  theme_bw()+
  facet_wrap(~res)+
  annotate(geom='line', x=temp_smooth-273.15,y=y_smooth_venolia)+
  labs(x="Temperature (°C)", y="Standardized rate", color=NULL)+
  scale_color_manual(values=c("top"="#dd4124", "bottom"='#0f85a0'),
                     breaks=c("top", "bottom"),
                     labels=c("top"="Top half", "bottom"="Bottom half"))

