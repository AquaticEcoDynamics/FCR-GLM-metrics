library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(stringr)

setwd(".../FCR-GLM-metrics/Uncertainty_analysis/Run1")

obs_3<- read.csv("glm3_reweight_ies.3.obs.csv")

#MOM
obs_mom <- read.csv("model_files/field_data/mom_observed.csv") %>% 
  mutate(DateTime = as.Date(DateTime, format="%Y-%m-%d"))%>%
  filter(DateTime > "2016-12-01") 


#MOM posterior 
obs_mom_post <- obs_3 %>% 
  dplyr:: select(grep("mom", names(obs_3))) %>%
  mutate(realizations = obs_3$real_name) %>%
  select(realizations, everything())

obs_mom_post_melt <- reshape2::melt(obs_mom_post) 

obs_mom_post_melt$realizations <- as.numeric(obs_mom_post_melt$realizations)
arrange_post <- arrange(obs_mom_post_melt, realizations) %>%
  mutate(DateTime = rep(obs_mom$DateTime, times=233)) %>%
  mutate(DateTime = as.Date(DateTime, format="%Y-%m-%d"))


means_by_variable <- arrange_post %>%
  group_by(variable) %>%
  summarize(mean_value = mean(value)) %>%
  mutate(DateTime = obs_mom$DateTime, observations = obs_mom$deviation) %>%
  mutate(Difference = abs(mean_value - observations))

  noise_mom <- max(means_by_variable$Difference)/2

#TD

obs_td <- read.csv("model_files/field_data/obs_td.csv") %>%
  mutate(DateTime = as.Date(DateTime, format="%Y-%m-%d")) %>%
  filter(DateTime > "2016-12-01" & DateTime < "2020-01-01") %>%
  mutate(month = lubridate::month(DateTime)) %>%
  filter(between(month, 4, 9)) %>%
  mutate(year = lubridate::year(DateTime)) %>%
  group_by(year) %>%
  filter(year == 2017 | year== 2018 | year==2019)

#TD posterior 
obs_td_post <- obs_3 %>% 
  dplyr:: select(grep("td", names(obs_3))) %>%
  mutate(realizations = obs_3$real_name) %>%
  select(realizations, everything())

obs_td_post_melt <- reshape2::melt(obs_td_post) 

obs_td_post_melt$realizations <- as.numeric(obs_td_post_melt$realizations)
arrange_td_post <- arrange(obs_td_post_melt, realizations) %>%
  mutate(DateTime = rep(obs_td$DateTime, times=233)) %>%
  mutate(DateTime = as.Date(DateTime, format="%Y-%m-%d")) %>%
  mutate(month = lubridate::month(DateTime)) %>%
  filter(between(month, 4, 9)) %>%
  mutate(year = lubridate::year(DateTime)) %>%
  group_by(year) %>%
  mutate(col= 1)

means_by_variable_TD <- arrange_td_post %>%
  group_by(variable) %>%
  summarize(mean_value = mean(value)) %>%
  mutate(DateTime = obs_td$DateTime, observations = obs_td$thermo.depth) %>%
  mutate(Difference = abs(mean_value - observations))

  noise_td <- max(means_by_variable_TD$Difference)/2

#Temperature
obs_temp <- read.csv("model_files/field_data/CleanedObsTemp.csv") %>%
  mutate(DateTime = as.Date(DateTime, format="%Y-%m-%d")) %>%
  filter(DateTime > "2016-12-01" & DateTime < "2020-01-01") %>%
  arrange(Depth) %>%
  filter(Depth != 9.2) 

#Temperature posterior
obs_temp_post <- obs_3 %>% 
  dplyr:: select(grep("_t_", names(obs_3))) %>%
  mutate(realizations = obs_3$real_name) %>%
  select(realizations, everything())

obs_temp_post_melt <- reshape2::melt(obs_temp_post) 

depth0.1_post <- obs_temp_post_melt[grepl('^wq0.1', obs_temp_post_melt$variable), ] %>%
  mutate(Depth = str_sub(variable, start=3, end=5))

other_depths_post<- obs_temp_post_melt[!grepl('^wq0.1', obs_temp_post_melt$variable), ]%>%
  mutate(Depth = str_sub(variable, start=3, end=3))

temp_df_post <- rbind(depth0.1_post, other_depths_post)
temp_df_post$realizations <- as.numeric(temp_df_post$realizations)

means_by_variable_Temp <- temp_df_post %>%
  group_by(variable) %>%
  summarize(mean_value = mean(value)) %>%
  mutate(DateTime = obs_temp$DateTime, observations = obs_temp$temp) %>%
  mutate(Difference = abs(mean_value - observations))

  noise_temp <- max(means_by_variable_Temp$Difference)/2

#Oxygen
obs_oxy <- read.csv("model_files/field_data/CleanedObsOxy.csv") %>%
  mutate(DateTime = as.Date(DateTime, format="%Y-%m-%d")) %>%
  filter(DateTime > "2016-12-01" & DateTime < "2020-01-01") %>%
  arrange(Depth) %>%
  filter(Depth != 9.2) 

#Oxy posterior
obs_oxy_post <- obs_3 %>% 
  dplyr:: select(grep("_ox_", names(obs_3))) %>%
  mutate(realizations = obs_3$real_name) %>%
  select(realizations, everything())

obs_oxy_post_melt <- reshape2::melt(obs_oxy_post) 

depth0.1_oxy_post <- obs_oxy_post_melt[grepl('^wq0.1', obs_oxy_post_melt$variable), ] %>%
  mutate(Depth = str_sub(variable, start=3, end=5))

other_depths_oxy_post<- obs_oxy_post_melt[!grepl('^wq0.1', obs_oxy_post_melt$variable), ]%>%
  mutate(Depth = str_sub(variable, start=3, end=3))

oxy_df_post <- rbind(depth0.1_oxy_post, other_depths_oxy_post)
oxy_df_post$realizations <- as.numeric(oxy_df_post$realizations)

means_by_variable_Oxy <- oxy_df_post %>%
  group_by(variable) %>%
  summarize(mean_value = mean(value)) %>%
  mutate(DateTime = obs_oxy$DateTime, observations = obs_oxy$OXY_oxy) %>%
  mutate(Difference = abs(mean_value - observations))

  noise_oxy <- max(means_by_variable_Oxy$Difference)/2


#Schmidt Stability
#SS observations
obs_SS <- read.csv("model_files/field_data/Obs_SS.csv") %>% 
  mutate(datetime = as.Date(datetime, format="%Y-%m-%d"))%>%
  filter(datetime > "2016-12-01" & datetime < "2020-01-01") 


#SS posterior
obs_ss_post <- obs_3 %>% 
  dplyr:: select(grep("ss", names(obs_3))) %>%
  mutate(realizations = obs_3$real_name) %>%
  select(realizations, everything())

obs_ss_post_melt <- reshape2::melt(obs_ss_post)

obs_ss_post_melt$realizations <- as.numeric(obs_ss_post_melt$realizations)
arrange_post_ss <- arrange(obs_ss_post_melt, realizations) %>%
  mutate(DateTime = rep(obs_SS$datetime, times=233)) %>%
  mutate(DateTime = as.Date(DateTime, format="%Y-%m-%d"))


means_by_variable_ss <- arrange_post_ss %>%
  group_by(variable) %>%
  summarize(mean_value = mean(value)) %>%
  mutate(DateTime = obs_SS$datetime, observations = obs_SS$schmidt.stability) %>%
  mutate(Difference = abs(mean_value - observations))

  noise_ss <- max(means_by_variable_ss$Difference)/2

#Number of anoxic layers
#Anoxic layers observations
obs_al <- read.csv("model_files/field_data/anoxia_observed.csv") %>% 
  mutate(DateTime = as.Date(DateTime, format="%Y-%m-%d"))%>%
  filter(DateTime > "2016-12-01" & DateTime < "2020-01-01") 

#Anoxia posterior
obs_al_post <- obs_3 %>% 
  dplyr:: select(grep("af", names(obs_3))) %>%
  mutate(realizations = obs_3$real_name) %>%
  select(realizations, everything())

obs_al_post_melt <- reshape2::melt(obs_al_post)

obs_al_post_melt$realizations <- as.numeric(obs_al_post_melt$realizations)
arrange_post_al <- arrange(obs_al_post_melt, realizations) %>%
  mutate(DateTime = rep(obs_al$DateTime, times=233)) %>%
  mutate(DateTime = as.Date(DateTime, format="%Y-%m-%d"))


means_by_variable_al <- arrange_post_al %>%
  group_by(variable) %>%
  summarize(mean_value = mean(value)) %>%
  mutate(DateTime = obs_al$DateTime, observations = obs_al$Count) %>%
  mutate(Difference = abs(mean_value - observations))

noise_al <- max(means_by_variable_al$Difference)/2