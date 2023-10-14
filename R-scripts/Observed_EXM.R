#Packages
library(ncdf4)
library(tidyr)
library(dplyr)
library(glmtools)
library(GLMr)
library(rLakeAnalyzer)
library(lubridate)
library(ggplot2)
library(reshape2)
library(stats)
library(utils)
library(ggpubr)
library(patchwork)

#Working directory
setwd(".../FCR-GLM-metrics/observations")

#Observed temperature data
obs_temp<-read.csv('CleanedObsTemp.csv')
obs_temp$DateTime <- as.Date(obs_temp$DateTime, format="%Y-%m-%d")

#Filtering out days with missing data and selecting modeling period
obs_temp <- subset(obs_temp, DateTime != "2019-04-29" & DateTime != "2019-05-30") %>%
  filter(DateTime > "2015-07-07" & DateTime < "2020-01-01")

#Interpolation of observed temperature data over depth
estimate_temp_by_date <- function(target_date, target_depth) {
  data_for_date <- obs_temp %>% 
    filter(DateTime == target_date) %>%
    arrange(Depth)
  approx(data_for_date$Depth, data_for_date$temp, xout = target_depth)$y
}

temp_interp_depth <- crossing(
  tibble(DateTime = unique(obs_temp$DateTime)),
  tibble(Depth = c(0.1, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 9.2))
) %>%
  group_by(DateTime) %>%
  mutate(temp = estimate_temp_by_date(DateTime[1], Depth))  

#Creating format required by rLakeAnalyzer::ts.thermo.depth and rLakeAnalyzer::ts.schmidt.stability
temp_wtr <- pivot_wider(temp_interp_depth, names_from=Depth, names_prefix='wtr_', values_from=temp)
#write.csv(temp_wtr, obs_temp_wtr.csv)

#Calculating thermocline depth (from 1 April to 30 September)
thermo_depth_df <- temp_wtr %>%
  mutate(month = lubridate::month(DateTime)) %>%
  filter(between(month, 4, 9)) %>%
  mutate(year = lubridate::year(DateTime)) %>%
  group_by(year)

#Deleting last two columns to have the right format for rLakeAnalyzer::ts.thermo.depth 
thermo_depth_df <- thermo_depth_df[, -c(13, 14)]
thermo_depth<- ts.thermo.depth(thermo_depth_df, Smin = 0.1, na.rm=TRUE, seasonal=FALSE)
#write.csv(thermo_depth, "obs_td.csv")

#Calculating Schmidt stability
bathy <- read.csv('bathymetry.csv')
schmidt_stability <- ts.schmidt.stability(temp_wtr, bathy, na.rm=TRUE)
#write.csv(schmidt_stability, "obs_SS.csv")

#Interpolation of observed temp data over depth by 0.1 m using estimate_temp_by_date function
temp_interp_depth_obs <- crossing(
  tibble(DateTime = unique(obs_temp$DateTime)),
  tibble(Depth = seq(0.1, 9.2, by=0.1))
) %>%
  group_by(DateTime) %>%
  mutate(temp = estimate_temp_by_date(DateTime[1], Depth))

temp_interp_depth_obs$DateTime <- as.Date(temp_interp_depth_obs$DateTime, format="%Y-%m-%d")

#interpolation over time
estimate_temp_by_time <- function(target_depth, target_date) {
  data_for_depth <- temp_interp_depth_obs %>% 
    filter(Depth == target_depth) %>%
    arrange(DateTime)
  
  approx(data_for_depth$DateTime, data_for_depth$temp, xout = target_date)$y
}

temp_obs_interp <- crossing(
  tibble(DateTime = seq(as.Date("2016-12-02", format="%Y-%m-%d"), as.Date("2019-12-31", format="%Y-%m-%d"), by = 1)),
  tibble(Depth = unique(temp_interp_depth_obs$Depth))
) %>%
  group_by(Depth) %>%
  mutate(temp = estimate_temp_by_time(Depth[1], DateTime))

#Save interpolated observed temperature dataset: 'obs_temp_interpolated.csv'
#write.csv(temp_obs_interp, 'obs_temp_interpolated.csv', row.names=FALSE)

#Observed oxygen data
obs_oxy<-read.csv('CleanedObsOxy.csv') %>%
  filter(DateTime > "2015-07-07" & DateTime < "2020-01-01") %>%
  mutate(DateTime = as.Date(DateTime, format="%Y-%m-%d"))

#MOM
epi_oxy_obs <- filter(obs_oxy, Depth==1)
hypo_oxy_obs <- filter(obs_oxy, Depth==8)
met_oxy_obs <- filter(obs_oxy, Depth==4)
obs_mom <- merge(epi_oxy_obs, hypo_oxy_obs, by="DateTime") %>%
  dplyr::rename(epi_oxy = OXY_oxy.x, hypo_oxy = OXY_oxy.y)

obs_mom$exp_oxy_obs <- (obs_mom$epi_oxy+obs_mom$hypo_oxy)/2

obs_mom <- merge(obs_mom, met_oxy_obs, by="DateTime")%>%
  dplyr::rename(met_oxy = OXY_oxy)
obs_mom$deviation <- obs_mom$met_oxy - obs_mom$exp_oxy_obs 

#write.csv(obs_mom[, c("DateTime", "deviation")], "mom_observed.csv", row.names=FALSE)

#Interpolation of observed oxygen data over depth by 0.1 m
estimate_oxy_by_date <- function(target_date, target_depth) {
  data_for_date <- obs_oxy %>% 
    filter(DateTime == target_date) %>%
    arrange(Depth)
  approx(data_for_date$Depth, data_for_date$OXY_oxy, xout = target_depth)$y
}

oxy_interp_depth_obs <- crossing(
  tibble(DateTime = unique(obs_oxy$DateTime)),
  tibble(Depth = seq(0.1, 9.2, by=0.1))
) %>%
  group_by(DateTime) %>%
  mutate(OXY_oxy = estimate_oxy_by_date(DateTime[1], Depth))

# Gettind rid of NAs
oxy_interp_depth_obs <-na.omit(oxy_interp_depth_obs)

#Interpolation of observed oxygen over time
estimate_oxy_by_depth <- function(target_depth, target_date) {
  data_for_depth <- oxy_interp_depth_obs %>% 
    filter(Depth == target_depth) %>%
    arrange(DateTime)
  approx(data_for_depth$DateTime, data_for_depth$OXY_oxy, xout = target_date)$y
}

oxy_obs_interp <- crossing(
  #interpolating in the calibration period
  tibble(DateTime = seq(as.Date("2015-07-09", format="%Y-%m-%d"), as.Date("2019-12-06", format="%Y-%m-%d"), by = 1)),
  tibble(Depth = unique(oxy_interp_depth_obs$Depth))
) %>%
  group_by(Depth) %>%
  mutate(OXY_oxy = estimate_oxy_by_depth(Depth[1], DateTime))

#Save interpolated observed oxy dataset: 'obs_oxy_interpolated.csv'
#write.csv(oxy_obs_interp, 'obs_oxy_interpolated.csv', row.names=FALSE)

#creating anoxia dataset: unit conversion from mmol/m3 to mg/L, selecting oxygen concentrations < 1mg/L, filer between May and November
obs_anoxia <- oxy_obs_interp %>%
  mutate(OXY_oxy=OXY_oxy*32/1000) %>%
  mutate(OXY_oxy=ifelse(OXY_oxy<=1, 1, 0)) %>%
  na.omit() %>%
  mutate(month = lubridate::month(DateTime)) %>%
  filter(between(month, 5, 11))

#Creating empty dataframe for loop
uniqueDates  <- unique(obs_anoxia$DateTime)

newData_obs  <- data.frame(
  DateTime = unique(obs_anoxia$DateTime),
  Count = length(uniqueDates)
)

#Calculating number of anoxic layers each day
for(i in 1:length(uniqueDates)){
  
  filteredData  <- filter(obs_anoxia, DateTime==uniqueDates[i] & OXY_oxy==1)
  newData_obs$Count[i] <- nrow(filteredData)
}

#write.csv(newData_obs, 'anoxia_observed.csv', row.names=FALSE)

