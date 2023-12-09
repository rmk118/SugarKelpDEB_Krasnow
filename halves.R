#testing halves instead of thirds

top_half_ctrl <- mean_growth_by_cross_ctrl %>% slice_max(gr_hp, n=15) %>% mutate(cross = fct_drop(cross))
bottom_half_ctrl <- mean_growth_by_cross_ctrl %>% slice_min(gr_hp, n=15) %>% mutate(cross = fct_drop(cross))

top_half_stress <- mean_growth_by_cross_stress %>% slice_max(gr_hp, n=15) %>% mutate(cross = fct_drop(cross))
bottom_half_stress <- mean_growth_by_cross_stress %>% slice_min(gr_hp, n=15) %>% mutate(cross = fct_drop(cross))

growth_rates_halves <- growth_rates_no_flags %>%
  mutate(ctrl_group = case_when(
    as.character(cross) %in% levels(top_half_ctrl$cross) ~ "top",
    as.character(cross) %in% levels(bottom_half_ctrl$cross) ~ "bottom"), 
    ctrl_group = as.factor(ctrl_group),
    stress_group = case_when(
      as.character(cross) %in% levels(top_half_stress$cross) ~ "top",
      as.character(cross) %in% levels(bottom_half_stress$cross) ~ "bottom"), 
    stress_group = as.factor(stress_group)) %>%
  ungroup()


# calculate weekly mean and maximum relative growth rates for each half, where grouping was done based on growth under non-control temperatures
weekly_halves_stress <- growth_rates_halves %>% 
  group_by(stress_group, trt, date, mean_temp) %>%
  summarise(mean_growth = mean(growth_hp, na.rm=TRUE),
            mean_rgr = mean(growth_rate_hp, na.rm=TRUE),
            mean_rgr = mean(rgr, na.rm=TRUE),
            max_rgr = max(rgr, na.rm=TRUE)) %>% 
  add_week_fun() %>% #add week labels
  filter(!(week %in% c(1,5))) #remove weeks without growth rates

# calculate weekly mean and maximum relative growth rates for each of the 3 groups, where grouping was done based on growth under control temperatures
weekly_halves_ctrl <- growth_rates_halves %>% 
  group_by(ctrl_group, trt, date, mean_temp) %>% 
  summarise(mean_growth = mean(growth_hp, na.rm=TRUE),
            mean_rgr = mean(growth_rate_hp, na.rm=TRUE),
            mean_rgr = mean(rgr, na.rm=TRUE),
            max_rgr = max(rgr, na.rm=TRUE)) %>% add_week_fun() %>%
  filter(!(week %in% c(1,5)))


ref_stress_half <- weekly_halves_stress %>% 
  group_by(stress_group) %>% 
  filter(round(mean_temp,1)==20) %>% 
  mutate(ref_rgr = max_rgr) %>% 
  select(stress_group, ref_rgr)

ref_ctrl_half <- weekly_halves_ctrl %>% 
  group_by(ctrl_group) %>% 
  filter(round(mean_temp,1)==20) %>% 
  mutate(ref_rgr = max_rgr) %>% 
  select(ctrl_group, ref_rgr)

# Finding standardized RGR at the group level
unh_stress_half <- weekly_halves_stress %>% 
  left_join(ref_stress_half) %>%  #add a column to the weekly means (stress grouping) with reference RGRs
  mutate(across(c(mean_rgr, max_rgr),
                ~if_else(.x<0 | is.na(.x), 0, .x))) %>% # Replace negative or missing RGRs with 0
  mutate(std_rgr_max = max_rgr/ref_rgr, # standardized max RGR (max RGR divided by ref temp RGR)
         std_rgr_mean = mean_rgr/ref_rgr, # standardized mean RGR (mean RGR divided by ref temp RGR)
         std_rgr = case_when(week==2 & trt!="control"~ std_rgr_max, #use max for week 2 non-control
                             .default= std_rgr_mean)) %>% 
  ungroup() %>%
  mutate(paper="unh_stress_all2023", paper_full=paste(paper, stress_group, sep="_"), temp=mean_temp, rate=mean_rgr, extra_info=NA, temp_K=temp+273.15) %>% #adding columns with identifying info
  select(paper, paper_full, temp, rate, std_rgr, extra_info, temp_K, stress_group)

unh_ctrl_half <- weekly_halves_ctrl %>% 
  left_join(ref_ctrl_half) %>%  #add a column to the weekly means (ctrl grouping) with reference RGRs
  mutate(across(c(mean_rgr, max_rgr),
                ~if_else(.x<0 | is.na(.x), 0, .x))) %>% # Replace negative or missing RGRs with 0
  mutate(std_rgr_max = max_rgr/ref_rgr, # standardized max RGR (max RGR divided by ref temp RGR)
         std_rgr_mean = mean_rgr/ref_rgr, # standardized mean RGR (mean RGR divided by ref temp RGR)
         std_rgr = case_when(week==2 & trt!="control"~ std_rgr_max, #use max for week 2 non-control
                             .default= std_rgr_mean)) %>% 
  ungroup() %>%
  mutate(paper="unh_ctrl_all2023", paper_full=paste(paper, ctrl_group, sep="_"), temp=mean_temp, rate=mean_rgr, extra_info=NA, temp_K=temp+273.15) %>% #adding columns with identifying info
  select(paper, paper_full, temp, rate, std_rgr, extra_info, temp_K, ctrl_group)


##### Cross level ####################################################################

###### Stress #####
unh_cross_stress_half <-growth_rates_halves %>% 
  group_by(cross, date, trt, mean_temp) %>% #group by cross, week, and treatment
  summarize(rgr = mean(rgr, na.rm=TRUE)) %>%  #find mean RGR for each cross under each trt in each week
  mutate(stress_group = case_when(
    as.character(cross) %in% levels(top_half_stress$cross) ~ "top",
    as.character(cross) %in% levels(bottom_half_stress$cross) ~ "bottom"), 
    stress_group = as.factor(stress_group)) %>% # add stress grouping labels
  left_join(ref_stress_half) %>% # add reference RGRs for stress groupings
  mutate(std_rate = rgr/ref_rgr) %>% # find standardized rate by dividing RGR at non-20°C temps by the reference rate, which is the stress group RGR at 20°C
  mutate(std_rate = if_else(std_rate<0 | is.na(std_rate), 0, std_rate)) %>% #replace negative rates with 0
  mutate(ctrl_group = case_when(
    as.character(cross) %in% levels(top_half_ctrl$cross) ~ "top",
    as.character(cross) %in% levels(bottom_half_ctrl$cross) ~ "bottom"), 
    ctrl_group = as.factor(ctrl_group)) %>%  #add control grouping labels
  ungroup() %>% 
  mutate(paper="unh_cross_stress2023", paper_full = paste(paper, stress_group, sep="_"), temp=mean_temp, rate=rgr, extra_info=paste(trt, cross,sep="_"), temp_K= temp + 273.15) %>%  #adding columns with identifying info
  select(paper, paper_full, temp, rate,std_rate, extra_info, temp_K, ctrl_group, stress_group)

###### Control #####
unh_cross_ctrl_half <-growth_rates_halves %>% 
  group_by(cross, date, trt, mean_temp) %>% #group by cross, week, and treatment
  summarize(rgr = mean(rgr, na.rm=TRUE)) %>%  #find mean RGR for each cross under each trt in each week
  mutate(ctrl_group = case_when(
    as.character(cross) %in% levels(top_half_ctrl$cross) ~ "top",
    as.character(cross) %in% levels(bottom_half_ctrl$cross) ~ "bottom"), 
    ctrl_group = as.factor(ctrl_group)) %>% # add ctrl grouping labels
  left_join(ref_ctrl_half) %>% # add reference RGRs for ctrl groupings
  mutate(std_rate = rgr/ref_rgr) %>% # find standardized rate by dividing RGR at non-20°C temps by the reference rate, which is the ctrl group RGR at 20°C
  mutate(std_rate = if_else(std_rate<0 | is.na(std_rate), 0, std_rate)) %>% #replace negative rates with 0
  mutate(stress_group = case_when(
    as.character(cross) %in% levels(top_half_stress$cross) ~ "top",
    as.character(cross) %in% levels(bottom_half_stress$cross) ~ "bottom"), 
    stress_group = as.factor(stress_group)) %>%  #add stress grouping labels
  ungroup() %>% 
  mutate(paper="unh_cross_ctrl2023", paper_full = paste(paper, ctrl_group, sep="_"), temp=mean_temp, rate=rgr, extra_info=paste(trt, cross,sep="_"), temp_K= temp + 273.15) %>%  #adding columns with identifying info
  select(paper, paper_full, temp, rate,std_rate, extra_info, temp_K, ctrl_group, stress_group)


#### Combine new and lit data #### 
lit_data_plus_half <- bind_rows(lit_data, 
                           unh_stress_half %>% mutate(std_rate=std_rgr,ctrl_group = NA,.keep="unused")) %>% 
  mutate(type="stress", res="means")

lit_data_plus_cross_half <- bind_rows(lit_data,
                               unh_cross_stress_half) %>%  mutate(type="stress", res="cross")

lit_data_plus_ctrl_half <- bind_rows(lit_data,
                                unh_ctrl_half %>% mutate(std_rate=std_rgr, stress_group = NA,.keep="unused")) %>%
  mutate(type="ctrl", res="means")

lit_data_plus_cross_ctrl_half <- bind_rows(lit_data,
                                    unh_cross_ctrl_half) %>% mutate(type="ctrl", res="cross")

all_lit_data_half <- rbind(lit_data_plus_half, lit_data_plus_cross_half, lit_data_plus_ctrl_half, lit_data_plus_cross_ctrl_half) %>% distinct() %>% # removes duplicate rows (i.e., the original literature data)
  mutate(level = case_when(str_ends(paper_full, "top") ~ "top",
                           str_ends(paper_full, "bottom") ~ "bottom",
                           .default = "lit"),
         level = as_factor(level),
         level = fct_relevel(level, c("top", "bottom", "lit")))

ggplot()+geom_point(data=all_lit_data_half, aes(x=temp, y=std_rate, color=level))+facet_grid(type~res)


#### Stress means
stress_means_half_df <- map(.x=c(top="top", bottom="bottom"), .f=nls_fun, df=lit_data_plus_half, type="stress") %>% bind_rows(.id="level") %>% mutate(res="means")

#### Stress crosses
stress_cross_half_df <- map(.x=c(top="top", bottom="bottom"), .f=nls_fun, df=lit_data_plus_cross_half, type="stress") %>% bind_rows(.id="level") %>% mutate(res="cross")

##### Control means
ctrl_means_half_df <- map(.x=c(top="top", bottom="bottom"), .f=nls_fun, df=lit_data_plus_ctrl_half, type="ctrl") %>% bind_rows(.id="level") %>% mutate(res="means")

#### Control crosses
ctrl_cross_half_df <- map(.x=c(top="top", bottom="bottom"), .f=nls_fun, df=lit_data_plus_cross_ctrl_half, type="ctrl") %>% bind_rows(.id="level") %>% mutate(res="cross")

halves_calibrations <- rbind(stress_means_half_df, ctrl_means_half_df, stress_cross_half_df, ctrl_cross_half_df) %>% 
  mutate(level = as_factor(level),
         level = fct_expand(level, "lit"),
         level = fct_relevel(level, c("top", "bottom", "lit")))

params_halves<-all_calibrations %>% select(T_A, T_H, T_AH,type,level,res) %>% distinct() %>% 
  add_row(T_A=T_A, T_H=T_H,T_AH=T_AH,type="orig", level="orig", res="orig")%>% 
  add_row(T_A=new_T_A, T_H=new_T_H,T_AH=new_T_AH,type="lit", level="lit", res="lit")

ggplot()+
  geom_line(data=halves_calibrations %>% ungroup(), aes(x=temp_smooth-273.15, y=std_rate, color=level), linewidth=2)+
  theme_bw()+
  facet_grid(type~res)+
  annotate(geom='line', x=temp_smooth-273.15,y=y_smooth_venolia)+
  labs(x="Temperature (°C)", y="Standardized rate", color=NULL)+
  scale_color_manual(values=c("top"="#dd4124", "bottom"='#0f85a0'),
                     breaks=c("top", "bottom"),
                     labels=c("top"="Top half", "bottom"="Bottom half"))

