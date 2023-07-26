library(zoo)
library(dplyr)


setwd('.../FCR-GLM-metrics')
obs_td <- read.csv('Uncertainty_analysis/model_files/field_data/obs_td.csv') %>%
  filter(DateTime > "2016-12-01", DateTime < "2020-01-01")
obs_mom <- read.csv('Uncertainty_analysis/model_files/field_data/mom_observed.csv') %>%
  filter(DateTime > "2016-12-01", DateTime < "2020-01-01")
obs_ss <- read.csv('Uncertainty_analysis/model_files/field_data/Obs_SS.csv') %>%
  rename(DateTime = datetime) %>%
  filter(DateTime > "2016-12-01", DateTime < "2020-01-01")
obs_oxy <- read.csv('Uncertainty_analysis/model_files/field_data/CleanedObsOxy.csv') %>%
  filter(DateTime > "2016-12-01", DateTime < "2020-01-01", Depth != 9.2)
obs_temp <- read.csv('Uncertainty_analysis/model_files/field_data/CleanedObsTemp.csv') %>%
  filter(DateTime > "2016-12-01", DateTime < "2020-01-01", Depth != 9.2)
obs_af <- read.csv('Uncertainty_analysis/model_files/field_data/anoxia_observed.csv') %>%
  filter(DateTime > "2016-12-01", DateTime < "2020-01-01")

#Linear approximation
obs_mom$DateTime <- as.Date(obs_mom$DateTime, format="%Y-%m-%d")
known_dates <- obs_mom$DateTime
known_values <- obs_mom$deviation

# Generate sequence of dates between the sampling dates
interpolation_dates <- seq(from = min(known_dates), to = max(known_dates), by = "day")

# Perform linear interpolation for each date
interpolated_values <- approx(known_dates, known_values, xout = interpolation_dates)$y

# Combine the dates and interpolated values into a data frame
interpolated_data_mom <- data.frame(DateTime = interpolation_dates, Value = interpolated_values)


#MOM
interpolated_data_mom <- interpolated_data_mom %>%
  mutate(rolling_avg = rollmean(Value, k=7, fill=NA, align='right')) %>%
  mutate(Difference = Value - rolling_avg)
new_mom_sd <- sd(interpolated_data_mom$Difference, na.rm=TRUE)

#TD
#Linear approximation
obs_td$DateTime <- as.Date(obs_td$DateTime, format="%Y-%m-%d")
known_dates <- obs_td$DateTime
known_values <- obs_td$thermo.depth

# Generate sequence of dates between the sampling dates
interpolation_dates <- seq(from = min(known_dates), to = max(known_dates), by = "day")

# Perform linear interpolation for each date
interpolated_values <- approx(known_dates, known_values, xout = interpolation_dates)$y

# Combine the dates and interpolated values into a data frame
interpolated_data_td <- data.frame(DateTime = interpolation_dates, Value = interpolated_values)

#TD
interpolated_data_td <- interpolated_data_td %>%
  mutate(rolling_avg = rollmean(Value, k=7, fill=NA, align='right')) %>%
  mutate(Difference = Value - rolling_avg)
new_td_sd <- sd(interpolated_data_td$Difference, na.rm=TRUE)

#SS
#Linear approximation
obs_ss$DateTime <- as.Date(obs_ss$DateTime, format="%Y-%m-%d")
known_dates <- obs_ss$DateTime
known_values <- obs_ss$schmidt.stability

# Generate sequence of dates between the sampling dates
interpolation_dates <- seq(from = min(known_dates), to = max(known_dates), by = "day")

# Perform linear interpolation for each date
interpolated_values <- approx(known_dates, known_values, xout = interpolation_dates)$y

# Combine the dates and interpolated values into a data frame
interpolated_data_ss <- data.frame(DateTime = interpolation_dates, Value = interpolated_values)

interpolated_data_ss <- interpolated_data_ss %>%
  mutate(rolling_avg = rollmean(Value, k=7, fill=NA, align='right')) %>%
  mutate(Difference = Value - rolling_avg)
new_ss_sd <- sd(interpolated_data_ss$Difference, na.rm=TRUE)

#AF
#Linear approximation
obs_af$DateTime <- as.Date(obs_af$DateTime, format="%Y-%m-%d")
known_dates <- obs_af$DateTime
known_values <- obs_af$Count

# Generate sequence of dates between the sampling dates
interpolation_dates <- seq(from = min(known_dates), to = max(known_dates), by = "day")

# Perform linear interpolation for each date
interpolated_values <- approx(known_dates, known_values, xout = interpolation_dates)$y

# Combine the dates and interpolated values into a data frame
interpolated_data_af <- data.frame(DateTime = interpolation_dates, Value = interpolated_values)

interpolated_data_af <- interpolated_data_af %>%
  mutate(rolling_avg = rollmean(Value, k=7, fill=NA, align='right')) %>%
  mutate(Difference = Value - rolling_avg)
new_af_sd <- sd(interpolated_data_af$Difference, na.rm=TRUE)


#OXY
depths <- c(0.1, 1, 2, 3, 4, 5, 6, 7, 8, 9)

match_noise_oxy <- data.frame()
#match <- data.frame()

for(i in 1:length(depths)){
  
  
    match <- obs_oxy %>%
    dplyr::filter(Depth == depths[i])
    known_dates <- as.Date(match$DateTime, format="%Y-%m-%d")
    known_values <- match$OXY_oxy
      
    # Generate sequence of dates between the sampling dates
    interpolation_dates <- seq(from = min(known_dates), to = max(known_dates), by = "day")
      
    # Perform linear interpolation for each date
    interpolated_values <- approx(known_dates, known_values, xout = interpolation_dates)$y
      
    # Combine the dates and interpolated values into a data frame
    interpolated_data <- data.frame(DateTime = interpolation_dates, Value = interpolated_values) %>%
      mutate(rolling_avg = rollmean(Value, k=7, fill=NA, align='right')) %>%
      mutate(Difference = Value - rolling_avg) %>%
      mutate(Depth = round(depths[i], 1))
      
    match_noise_oxy <- rbind(match_noise_oxy, interpolated_data)
  
}
match_noise_oxy$DateTime <- as.Date(match_noise_oxy$DateTime, format="%Y-%m-%d")
obs_oxy$DateTime <- as.Date(obs_oxy$DateTime, format="%Y-%m-%d")
new_oxy_sd <- sd(match_noise_oxy$Difference, na.rm = TRUE)


#TEMP
depths <- c(0.1, 1, 2, 3, 4, 5, 6, 7, 8, 9)

match_noise_temp <- data.frame()

for(i in 1:length(depths)){
  
  
  match <- obs_temp %>%
    dplyr::filter(Depth == depths[i])
  known_dates <- as.Date(match$DateTime, format="%Y-%m-%d")
  known_values <- match$temp
  
  # Generate sequence of dates between the sampling dates
  interpolation_dates <- seq(from = min(known_dates), to = max(known_dates), by = "day")
  
  # Perform linear interpolation for each date
  interpolated_values <- approx(known_dates, known_values, xout = interpolation_dates)$y
  
  # Combine the dates and interpolated values into a data frame
  interpolated_data <- data.frame(DateTime = interpolation_dates, Value = interpolated_values) %>%
    mutate(rolling_avg = rollmean(Value, k=7, fill=NA, align='right')) %>%
    mutate(Difference = Value - rolling_avg) %>%
    mutate(Depth = round(depths[i], 1))
  
  match_noise_temp <- rbind(match_noise_temp, interpolated_data)
  
}
match_noise_temp$DateTime <- as.Date(match_noise_temp$DateTime, format="%Y-%m-%d")
obs_temp$DateTime <- as.Date(obs_temp$DateTime, format="%Y-%m-%d")
new_temp_sd <- sd(match_noise_temp$Difference, na.rm = TRUE)

