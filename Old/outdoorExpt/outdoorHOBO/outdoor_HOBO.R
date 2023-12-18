#Organizing outdoor tank sporophyte heat stress HOBO data
#Ruby Krasnow
#Last updated: Nov 30, 2023

library(tidyverse)
library(patchwork)
library(lubridate)
library(clipr)
library(tseries)
library(zoo)

### Import data ####
control_files <- list.files("~/Downloads/SugarKelpDEB_Krasnow/Old/outdoorExpt/outdoorHOBO/outdoor_HOBO_data/highN_control", full.names = TRUE)
control <- suppressWarnings(read_csv(control_files,
              col_select = c("num","dateTime", "temp","lux"),
             col_types = c("d", "c", "d", "d"))) %>% 
  na.omit() %>%
  mutate(datetime = parse_date_time(dateTime, c("mdy IMS p", "mdy HM")),
         trt="control",.keep="unused")

lowN_files <- list.files("~/Downloads/SugarKelpDEB_Krasnow/Old/outdoorExpt/outdoorHOBO/outdoor_HOBO_data/lowN_temp_stress", full.names = TRUE)
lowN<-suppressWarnings(read_csv(lowN_files,
         col_select = c("num","dateTime", "temp","lux"),
         col_types = c("d", "c", "d", "d"))) %>%
  na.omit() %>%
  mutate(datetime = parse_date_time(dateTime, "mdy HM"),
         trt="lowN", .keep="unused")

highN_files <- list.files("~/Downloads/SugarKelpDEB_Krasnow/Old/outdoorExpt/outdoorHOBO/outdoor_HOBO_data/highN_temp_stress", full.names = TRUE)
highN<-suppressWarnings(read_csv(highN_files,
                                col_select = c("num","dateTime", "temp","lux"),
                                col_types = c("d", "c", "d", "d")))%>%
  na.omit() %>%
  mutate(datetime = parse_date_time(dateTime, "mdy HM"),
         trt="highN", .keep="unused")

dat_full <- bind_rows(control, lowN, highN) %>% #combine dataframes
  select(c(trt, datetime,temp, lux)) %>% #select relevant columns
  mutate(trt = as.factor(trt)) #convert treatment to factor

ggplot(data=dat_full %>% filter(temp<30))+
  geom_line(aes(x=datetime, y=temp, color=trt))+
  scale_x_datetime(breaks="4 days", date_labels = "%b %d")+
  labs(x=NULL, y="Temperature (°C)", color="Treatment")+
  theme_light()+
  theme(text = element_text(size=14),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title = element_text(size=13), axis.text = element_text(size=14))+
  scale_color_hue(labels = c("Control + N", "Heat + N", "Heat"))+
  geom_vline(xintercept = c(as_datetime("2023-06-19"), as_datetime("2023-06-30")))

all_hourly <- dat_full %>% 
  mutate(date = round_date(datetime, "hour")) %>% 
  group_by(trt, date) %>% 
  summarize(temp_hourly = mean(temp, na.rm=TRUE),
            #PAR_hourly = mean(lux*0.0185, na.rm=TRUE)) %>% 
            PAR_hourly = mean(lux, na.rm=TRUE)) %>% 
  filter(date > ymd_hm("2023-06-11 23:59") & date < ymd_hm("2023-07-18 23:59"))

week1 <- interval(ymd_hm("2023-06-11 23:59"), ymd("2023-06-20"))
week2 <- interval(ymd("2023-06-20"), ymd("2023-06-27")) #6/28 growth reflects Week 2
week3 <- interval(ymd("2023-06-27"), ymd("2023-07-04")) #7/05 growth reflects Week 3
week4 <- interval(ymd("2023-07-04"), ymd("2023-07-11")) #7/13 growth reflects Week 4
week5 <- interval(ymd("2023-07-11"), ymd_hm("2023-07-18 23:59"))

all_hourly %>% #week 2
  filter(date %within% week2) %>% 
  filter(trt!="control")%>% summarise(mean(temp_hourly)) #highN = 19.9

all_hourly <- all_hourly %>% mutate(week = case_when(
  date %within% week1 ~ 1,
  date %within% week2 ~ 2,
  date %within% week3 ~ 3,
  date %within% week4 ~ 4,
  date %within% week5 ~ 5
))

weekly_means_degC <- all_hourly %>% group_by(week, trt) %>% 
  summarise(mean_temp = mean(temp_hourly)) %>% 
  mutate(mean_temp = round(mean_temp,2)) %>% arrange(trt)

# Temperature
ggplot(data=all_hourly %>% filter(date < ymd("2023-07-18")))+
  geom_vline(xintercept = c(as_datetime(c("2023-06-12","2023-06-20","2023-06-27","2023-07-04","2023-07-11","2023-07-18"))), linetype="dashed")+
  geom_line(aes(x=date, y=temp_hourly, color=trt))+
  scale_x_datetime(breaks="7 days", date_labels = "%b %d")+
  labs(x=NULL, y="Temperature (°C)", color="Treatment")+
  theme_light()+
  theme(text = element_text(size=14),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title = element_text(size=13), axis.text = element_text(size=14))+
  scale_color_hue(labels = c("Control + N", "Heat + N", "Heat"))

# PAR
ggplot(data=all_hourly %>% filter(date < ymd("2023-07-14") & date > ymd("2023-06-22"), PAR_hourly<750) %>% group_by(date) %>% summarise(PAR=mean(PAR_hourly)))+
  geom_vline(xintercept = c(as_datetime(c("2023-06-12","2023-06-20","2023-06-27","2023-07-04","2023-07-11"))), linetype="dashed")+
  geom_line(aes(x=date, y=PAR*60*60/10^5.6))+
  #geom_line(aes(x=date, y=PAR))+
  scale_x_datetime(breaks="7 days", date_labels = "%b %d")+
  labs(x=NULL, y="PAR", color="Treatment")+
  theme_light()+
  theme(text = element_text(size=14),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title = element_text(size=13), axis.text = element_text(size=14))

PAR_forcing <- all_hourly %>% filter(date < ymd("2023-07-14") & date > ymd("2023-06-22"), PAR_hourly<750) %>% group_by(date) %>% summarise(PAR=mean(PAR_hourly)*60*60/10^5.6*0.5)

temp_forcing <- all_hourly %>% filter(date < ymd("2023-07-14") & date > ymd("2023-06-22")) %>% select(date, trt, temp_hourly) %>% pivot_wider(names_from = trt, values_from = temp_hourly) %>% rename(control_temp=control, highN_temp=highN, lowN_temp=lowN)
