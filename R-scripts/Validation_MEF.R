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

#Deep mixing 2
setwd(".../FCR-GLM-metrics")
sim_folder <- getwd()
error <- read.csv('Observations/error_stats.csv')

#oxygen
var="OXY_oxy"
obs_oxy<-read.csv('Observations/CleanedObsOxy.csv') %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) 
obs_oxy$DateTime <- as.Date(obs_oxy$DateTime, format="%Y-%m-%d")

depths<- c(0.1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9.2) 

#Deep mixing 2 naive
new<- file.path(sim_folder, 'Validation/valid_naive_deepm2/output/output.nc')

new_oxy <- get_var(new, "OXY_oxy", reference="surface", z_out=depths) %>%
  pivot_longer(cols=starts_with("OXY_oxy_"), names_to="Depth", names_prefix="OXY_oxy_", values_to = "OXY_oxy") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) 
new_oxy$DateTime <- as.Date(new_oxy$DateTime, format="%Y-%m-%d")

oxygen <- merge(obs_oxy, new_oxy, by=c("DateTime","Depth")) %>%
  dplyr::rename(obsoxy = OXY_oxy.x, mod_oxy = OXY_oxy.y) 
oxygen$DateTime <- as.Date(oxygen$DateTime, format="%Y-%m-%d")

#MEF
for (i in 1:nrow(oxygen)) {
  oxygen$MEF_1[i]<- ((oxygen$mod_oxy[i]- oxygen$obsoxy[i])^2)
  oxygen$MEF_2[i]<- ((oxygen$obsoxy[i]-mean(oxygen$obsoxy))^2)
  MEF_oxy<- 1-(sum(oxygen$MEF_1)/sum(oxygen$MEF_2))
}

#deep mixing 2 naive
error[error$metric=="oxy" & error$calibration=="PEST_N", "Validation.deepm2"] <- MEF_oxy

#Deep mixing 2 exm w1
new_w1 <- file.path(sim_folder, 'Validation/valid_exm_deepm2_w1/output/output.nc')

new_oxy_w1 <- get_var(new_w1, "OXY_oxy", reference="surface", z_out=depths) %>%
  pivot_longer(cols=starts_with("OXY_oxy_"), names_to="Depth", names_prefix="OXY_oxy_", values_to = "OXY_oxy") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) 
new_oxy_w1$DateTime <- as.Date(new_oxy_w1$DateTime, format="%Y-%m-%d")

oxygen_w1 <- merge(obs_oxy, new_oxy_w1, by=c("DateTime","Depth")) %>%
  dplyr::rename(obsoxy = OXY_oxy.x, mod_oxy = OXY_oxy.y) 
oxygen_w1$DateTime <- as.Date(oxygen_w1$DateTime, format="%Y-%m-%d")

#MEF
for (i in 1:nrow(oxygen_w1)) {
  oxygen_w1$MEF_1[i]<- ((oxygen_w1$mod_oxy[i]- oxygen_w1$obsoxy[i])^2)
  oxygen_w1$MEF_2[i]<- ((oxygen_w1$obsoxy[i]-mean(oxygen_w1$obsoxy))^2)
  MEF_oxy_w1<- 1-(sum(oxygen_w1$MEF_1)/sum(oxygen_w1$MEF_2))
}

error[error$metric=="oxy" & error$calibration=="PEST_exm_w1", "Validation.deepm2"] <- MEF_oxy_w1

#Deep mixing 2 exm w2
new_w2 <- file.path(sim_folder, 'Validation/valid_exm_deepm2_w2/output/output.nc')

new_oxy_w2 <- get_var(new_w2, "OXY_oxy", reference="surface", z_out=depths) %>%
  pivot_longer(cols=starts_with("OXY_oxy_"), names_to="Depth", names_prefix="OXY_oxy_", values_to = "OXY_oxy") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) 
new_oxy_w2$DateTime <- as.Date(new_oxy_w2$DateTime, format="%Y-%m-%d")

oxygen_w2 <- merge(obs_oxy, new_oxy_w2, by=c("DateTime","Depth")) %>%
  dplyr::rename(obsoxy = OXY_oxy.x, mod_oxy = OXY_oxy.y) 
oxygen_w2$DateTime <- as.Date(oxygen_w2$DateTime, format="%Y-%m-%d")

#MEF
for (i in 1:nrow(oxygen_w2)) {
  oxygen_w2$MEF_1[i]<- ((oxygen_w2$mod_oxy[i]- oxygen_w2$obsoxy[i])^2)
  oxygen_w2$MEF_2[i]<- ((oxygen_w2$obsoxy[i]-mean(oxygen_w2$obsoxy))^2)
  MEF_oxy_w2<- 1-(sum(oxygen_w2$MEF_1)/sum(oxygen_w2$MEF_2))
}

error[error$metric=="oxy" & error$calibration=="PEST_exm_w2", "Validation.deepm2"] <- MEF_oxy_w2

#Deep mixing 2 exm w3
new_w3 <- file.path(sim_folder, 'Validation/valid_exm_deepm2_w3/output/output.nc')

new_oxy_w3 <- get_var(new_w3, "OXY_oxy", reference="surface", z_out=depths) %>%
  pivot_longer(cols=starts_with("OXY_oxy_"), names_to="Depth", names_prefix="OXY_oxy_", values_to = "OXY_oxy") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) 
new_oxy_w3$DateTime <- as.Date(new_oxy_w3$DateTime, format="%Y-%m-%d")

oxygen_w3 <- merge(obs_oxy, new_oxy_w3, by=c("DateTime","Depth")) %>%
  dplyr::rename(obsoxy = OXY_oxy.x, mod_oxy = OXY_oxy.y) 
oxygen_w3$DateTime <- as.Date(oxygen_w3$DateTime, format="%Y-%m-%d")

#MEF
for (i in 1:nrow(oxygen_w3)) {
  oxygen_w3$MEF_1[i]<- ((oxygen_w3$mod_oxy[i]- oxygen_w3$obsoxy[i])^2)
  oxygen_w3$MEF_2[i]<- ((oxygen_w3$obsoxy[i]-mean(oxygen_w3$obsoxy))^2)
  MEF_oxy_w3<- 1-(sum(oxygen_w3$MEF_1)/sum(oxygen_w3$MEF_2))
}

error[error$metric=="oxy" & error$calibration=="PEST_exm_w3", "Validation.deepm2"] <- MEF_oxy_w3

#Temperature
var="temp"
obs_temp<-read.csv('Observations/CleanedObsTemp.csv') %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) 
obs_temp$DateTime <- as.Date(obs_temp$DateTime, format="%Y-%m-%d")

depths<- c(0.1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9.2) 

#Deep mixing 2 naive
new<- file.path(sim_folder, 'Validation/valid_naive_deepm2/output/output.nc')

new_temp <- get_var(new, "temp", reference="surface", z_out=depths) %>%
  pivot_longer(cols=starts_with("temp_"), names_to="Depth", names_prefix="temp_", values_to = "temp") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) 
new_temp$DateTime <- as.Date(new_temp$DateTime, format="%Y-%m-%d")

temp <- merge(obs_temp, new_temp, by=c("DateTime","Depth")) %>%
  dplyr::rename(obtemp = temp.x, modtemp = temp.y) 
temp$DateTime <- as.Date(temp$DateTime, format="%Y-%m-%d")

for (i in 1:nrow(temp)) {
  temp$MEF_1[i]<- ((temp$modtemp[i]- temp$obtemp[i])^2)
  temp$MEF_2[i]<- ((temp$obtemp[i]-mean(temp$obtemp))^2)
  MEF_temp<- 1-(sum(temp$MEF_1)/sum(temp$MEF_2))
}

error[error$metric=="temp" & error$calibration=="PEST_N", "Validation.deepm2"] <- MEF_temp

#Deep mixing 2 exm w1
new_w1 <- file.path(sim_folder, 'Validation/valid_exm_deepm2_w1/output/output.nc')

new_temp_w1 <- get_var(new_w1, "temp", reference="surface", z_out=depths) %>%
  pivot_longer(cols=starts_with("temp_"), names_to="Depth", names_prefix="temp_", values_to = "temp") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) 
new_temp_w1$DateTime <- as.Date(new_temp_w1$DateTime, format="%Y-%m-%d")

temp_w1 <- merge(obs_temp, new_temp_w1, by=c("DateTime","Depth")) %>%
  dplyr::rename(obtemp = temp.x, modtemp = temp.y) 
temp_w1$DateTime <- as.Date(temp_w1$DateTime, format="%Y-%m-%d")

for (i in 1:nrow(temp_w1)) {
  temp_w1$MEF_1[i]<- ((temp_w1$modtemp[i]- temp_w1$obtemp[i])^2)
  temp_w1$MEF_2[i]<- ((temp_w1$obtemp[i]-mean(temp_w1$obtemp))^2)
  MEF_temp_w1<- 1-(sum(temp_w1$MEF_1)/sum(temp_w1$MEF_2))
}

error[error$metric=="temp" & error$calibration=="PEST_exm_w1", "Validation.deepm2"] <- MEF_temp_w1

#Deep mixing 2 exm w2
new_w2 <- file.path(sim_folder, 'Validation/valid_exm_deepm2_w2/output/output.nc')

new_temp_w2 <- get_var(new_w2, "temp", reference="surface", z_out=depths) %>%
  pivot_longer(cols=starts_with("temp_"), names_to="Depth", names_prefix="temp_", values_to = "temp") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) 
new_temp_w2$DateTime <- as.Date(new_temp_w2$DateTime, format="%Y-%m-%d")

temp_w2 <- merge(obs_temp, new_temp_w2, by=c("DateTime","Depth")) %>%
  dplyr::rename(obtemp = temp.x, modtemp = temp.y) 
temp_w2$DateTime <- as.Date(temp_w2$DateTime, format="%Y-%m-%d")

for (i in 1:nrow(temp_w2)) {
  temp_w2$MEF_1[i]<- ((temp_w2$modtemp[i]- temp_w2$obtemp[i])^2)
  temp_w2$MEF_2[i]<- ((temp_w2$obtemp[i]-mean(temp_w2$obtemp))^2)
  MEF_temp_w2<- 1-(sum(temp_w2$MEF_1)/sum(temp_w2$MEF_2))
}

error[error$metric=="temp" & error$calibration=="PEST_exm_w2", "Validation.deepm2"] <- MEF_temp_w2

#Deep mixing 2 exm w3
new_w3 <- file.path(sim_folder, 'Validation/valid_exm_deepm2_w3/output/output.nc')

new_temp_w3 <- get_var(new_w3, "temp", reference="surface", z_out=depths) %>%
  pivot_longer(cols=starts_with("temp_"), names_to="Depth", names_prefix="temp_", values_to = "temp") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) 
new_temp_w3$DateTime <- as.Date(new_temp_w3$DateTime, format="%Y-%m-%d")

temp_w3 <- merge(obs_temp, new_temp_w3, by=c("DateTime","Depth")) %>%
  dplyr::rename(obtemp = temp.x, modtemp = temp.y) 
temp_w3$DateTime <- as.Date(temp_w3$DateTime, format="%Y-%m-%d")

for (i in 1:nrow(temp_w3)) {
  temp_w3$MEF_1[i]<- ((temp_w3$modtemp[i]- temp_w3$obtemp[i])^2)
  temp_w3$MEF_2[i]<- ((temp_w3$obtemp[i]-mean(temp_w3$obtemp))^2)
  MEF_temp_w3<- 1-(sum(temp_w3$MEF_1)/sum(temp_w3$MEF_2))
}

error[error$metric=="temp" & error$calibration=="PEST_exm_w3", "Validation.deepm2"] <- MEF_temp_w3

#TD
#Observed thermocline deph
obs_TD<- read.csv('Observations/obs_td.csv') %>%
  mutate(DateTime = as.Date(DateTime, format="%Y-%m-%d")) %>%
  mutate(month = lubridate::month(DateTime)) %>%
  filter(between(month, 4, 9)) %>%
  mutate(year = lubridate::year(DateTime)) %>%
  group_by(year) %>%
  filter(DateTime  > "2015-07-19" & DateTime < "2016-12-02")

obs_TD$DateTime <- as.Date(obs_TD$DateTime, format="%Y-%m-%d")

#Schmidt stability
bathy <- read.csv('observations/bathymetry.csv')

#Schmidt stability observed
schmidt_stability_obs<- read.csv("Observations/Obs_SS.csv") %>%
  mutate(datetime = as.Date(datetime, format = "%Y-%m-%d")) %>%
  filter(datetime > "2015-07-19" & datetime < "2016-12-02")
schmidt_stability_obs$datetime<- as.Date(schmidt_stability_obs$datetime, format='%Y-%m-%d')

#MOM Observed
obs_mom<-read.csv('observations/mom_observed.csv') %>%
  filter(DateTime > "2015-07-19" & DateTime < "2016-12-02")
obs_mom$DateTime <- as.Date(obs_mom$DateTime, format = "%Y-%m-%d")
depths<- c(1, 4, 8) 

###############################################################################################################
#Deepm2 naive
PEST_calib <- file.path(sim_folder, 'Validation/valid_naive_deepm2/output/output.nc')

#Thermocline depth
#Modelled thermocline depth Deep2 naive
temp<- get_var(PEST_calib, var_name="temp", reference="surface") 
colClean <- function(x){ colnames(x) <- gsub("temp", "wtr", colnames(x)); x } 
colClean(temp)
temp$DateTime <- as.Date(temp$DateTime, format="%Y-%m-$d")

thermo_depth_model <- ts.thermo.depth(temp, Smin = 0.1, na.rm=TRUE, seasonal=FALSE)  %>% 
  dplyr::rename(td_model = thermo.depth, DateTime = datetime) %>%
  mutate(month = lubridate::month(DateTime)) %>%
  filter(between(month, 4, 9)) %>%
  mutate(year = lubridate::year(DateTime)) %>%
  group_by(year) 

thermo_depth_model$DateTime <- as.Date(thermo_depth_model$DateTime, format="%Y-%m-%d")

#Calculating error (MEF) of thermocline depth for Figure 7
td_merge <- merge(thermo_depth_model, obs_TD, by="DateTime")

for (i in 1:nrow(td_merge)) {
  td_merge$MEFF_1[i]<- ((td_merge$td_model[i]- td_merge$thermo.depth[i])^2)
  td_merge$MEFF_2[i]<- ((td_merge$thermo.depth[i]-mean(td_merge$thermo.depth))^2)
  MEFF_TD<- 1-(sum(td_merge$MEFF_1)/sum(td_merge$MEFF_2))
}

#Adding calculated Deep mixing 2 naive model TD MEF to error table
error[error$metric=="TD" & error$calibration=="PEST_N", "Validation.deepm2"] <- MEFF_TD

#Schmidt stability
#Modelled Schmidt stability (deepm2 naive model)
schmidt_stability <- ts.schmidt.stability(temp, bathy, na.rm=TRUE) %>% 
  dplyr::rename(ss_PEST = schmidt.stability)

schmidt_stability$datetime <- as.Date(schmidt_stability$datetime, format="%Y-%m-%d")

#Calculating error (MEF) of schmidt stability for Figure 7
SS_merge <- merge(schmidt_stability_obs, schmidt_stability)

for (i in 1:nrow(SS_merge)) {
  SS_merge$MEFF_1[i]<- ((SS_merge$ss_PEST[i]- SS_merge$schmidt.stability[i])^2)
  SS_merge$MEFF_2[i]<- ((SS_merge$schmidt.stability[i]-mean(SS_merge$schmidt.stability))^2)
  MEFF_SS<- 1-(sum(SS_merge$MEFF_1)/sum(SS_merge$MEFF_2))
}

#Adding calculated Deep mixing 2 naive model SS MEF to error table
error[error$metric=="SS" & error$calibration=="PEST_N", "Validation.deepm2"] <- MEFF_SS

#MOM
#Oxygen deepm2 naive model
oxy <- get_var(PEST_calib, "OXY_oxy", reference="surface", z_out=depths) %>%
  pivot_longer(cols=starts_with("OXY_oxy_"), names_to="Depth", names_prefix="OXY_oxy_", values_to = "OXY_oxy") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) 
oxy$DateTime <- as.Date(oxy$DateTime, format="%Y-%m-%d")

#Merge into one dataset
epi_oxy <- filter(oxy, Depth==1)
hypo_oxy <- filter(oxy, Depth==8)
met_oxy <- filter(oxy, Depth==4)
merge_mod<- merge(epi_oxy, hypo_oxy, by="DateTime") %>%
  dplyr::rename(epi_oxy = OXY_oxy.x, hypo_oxy = OXY_oxy.y)
merge_mod$exp_oxy <- (merge_mod$epi_oxy+merge_mod$hypo_oxy)/2

#Calculate deviation between modelled met oxy and expected met oxy
merge_mod <- merge(merge_mod, met_oxy[, c("DateTime", "OXY_oxy")], by="DateTime")%>%
  dplyr::rename(met_oxy = OXY_oxy)
merge_mod$deviation <- merge_mod$met_oxy - merge_mod$exp_oxy

obs_mom$DateTime <- as.Date(obs_mom$DateTime, format="%Y-%m-%d")
merge_mod$DateTime <- as.Date(merge_mod$DateTime, format="%Y-%m-%d")

#Calculating error (MEF) of metalimnetic oxygen minima for Figure 7
merge_mom <- merge(obs_mom, merge_mod, by="DateTime")

for (i in 1:nrow(merge_mom)) {
  merge_mom$MEFF_1[i]<- ((merge_mom$deviation.y[i]- merge_mom$deviation.x[i])^2)
  merge_mom$MEFF_2[i]<- ((merge_mom$deviation.x[i]-mean(merge_mom$deviation.x))^2)
  MEFF_MOM<- 1-(sum(merge_mom$MEFF_1)/sum(merge_mom$MEFF_2))
}

#Adding calculated Deep mixing 2 naive model MOM MEF to error table
error[error$metric=="MOM" & error$calibration=="PEST_N", "Validation.deepm2"] <- MEFF_MOM

#############################################################################################################

#Deepm2 exm w1
PEST_calib <- file.path(sim_folder, 'Validation/valid_exm_deepm2_w1/output/output.nc')

#Thermocline depth
#Modelled thermocline depth Deepm2 exm w1
temp<- get_var(PEST_calib, var_name="temp", reference="surface") 
colClean <- function(x){ colnames(x) <- gsub("temp", "wtr", colnames(x)); x } 
colClean(temp)
temp$DateTime <- as.Date(temp$DateTime, format="%Y-%m-$d")

thermo_depth_model <- ts.thermo.depth(temp, Smin = 0.1, na.rm=TRUE, seasonal=FALSE)  %>% 
  dplyr::rename(td_model = thermo.depth, DateTime = datetime) %>%
  mutate(month = lubridate::month(DateTime)) %>%
  filter(between(month, 4, 9)) %>%
  mutate(year = lubridate::year(DateTime)) %>%
  group_by(year) 

thermo_depth_model$DateTime <- as.Date(thermo_depth_model$DateTime, format="%Y-%m-%d")

#Calculating error (MEF) of thermocline depth for Figure 7
td_merge <- merge(thermo_depth_model, obs_TD, by="DateTime")

for (i in 1:nrow(td_merge)) {
  td_merge$MEFF_1[i]<- ((td_merge$td_model[i]- td_merge$thermo.depth[i])^2)
  td_merge$MEFF_2[i]<- ((td_merge$thermo.depth[i]-mean(td_merge$thermo.depth))^2)
  MEFF_TD<- 1-(sum(td_merge$MEFF_1)/sum(td_merge$MEFF_2))
}

#Adding calculated MEF of deepm2 w1 to error table
error[error$metric=="TD" & error$calibration=="PEST_exm_w1", "Validation.deepm2"] <- MEFF_TD

#Schmidt stability
#Modelled Schmidt stability (deepm2 exm w1 model)
schmidt_stability <- ts.schmidt.stability(temp, bathy, na.rm=TRUE) %>% 
  dplyr::rename(ss_PEST = schmidt.stability)

schmidt_stability$datetime <- as.Date(schmidt_stability$datetime, format="%Y-%m-%d")

#Calculating error (MEF) of schmidt stability for Figure 7
SS_merge <- merge(schmidt_stability_obs, schmidt_stability)

for (i in 1:nrow(SS_merge)) {
  SS_merge$MEFF_1[i]<- ((SS_merge$ss_PEST[i]- SS_merge$schmidt.stability[i])^2)
  SS_merge$MEFF_2[i]<- ((SS_merge$schmidt.stability[i]-mean(SS_merge$schmidt.stability))^2)
  MEFF_SS<- 1-(sum(SS_merge$MEFF_1)/sum(SS_merge$MEFF_2))
}

#Adding calculated Deep mixing 2 exm w1 model SS MEF to error table
error[error$metric=="SS" & error$calibration=="PEST_exm_w1", "Validation.deepm2"] <- MEFF_SS

#MOM
#Oxygen Deep mixing 2 exm w1 model 
oxy <- get_var(PEST_calib, "OXY_oxy", reference="surface", z_out=depths) %>%
  pivot_longer(cols=starts_with("OXY_oxy_"), names_to="Depth", names_prefix="OXY_oxy_", values_to = "OXY_oxy") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) 
oxy$DateTime <- as.Date(oxy$DateTime, format="%Y-%m-%d")

#Merge into one dataset
epi_oxy <- filter(oxy, Depth==1)
hypo_oxy <- filter(oxy, Depth==8)
met_oxy <- filter(oxy, Depth==4)
merge_mod<- merge(epi_oxy, hypo_oxy, by="DateTime") %>%
  dplyr::rename(epi_oxy = OXY_oxy.x, hypo_oxy = OXY_oxy.y)
merge_mod$exp_oxy <- (merge_mod$epi_oxy+merge_mod$hypo_oxy)/2

#Calculate deviation between modelled met oxy and expected met oxy
merge_mod <- merge(merge_mod, met_oxy[, c("DateTime", "OXY_oxy")], by="DateTime")%>%
  dplyr::rename(met_oxy = OXY_oxy)
merge_mod$deviation <- merge_mod$met_oxy - merge_mod$exp_oxy

obs_mom$DateTime <- as.Date(obs_mom$DateTime, format="%Y-%m-%d")
merge_mod$DateTime <- as.Date(merge_mod$DateTime, format="%Y-%m-%d")

#Calculating error (MEF) of metalimnetic oxygen minima for Figure 7
merge_mom <- merge(obs_mom, merge_mod, by="DateTime")

for (i in 1:nrow(merge_mom)) {
  merge_mom$MEFF_1[i]<- ((merge_mom$deviation.y[i]- merge_mom$deviation.x[i])^2)
  merge_mom$MEFF_2[i]<- ((merge_mom$deviation.x[i]-mean(merge_mom$deviation.x))^2)
  MEFF_MOM<- 1-(sum(merge_mom$MEFF_1)/sum(merge_mom$MEFF_2))
}

#Adding calculated Deep mixing 2 exm w1 model MOM MEF to error table
error[error$metric=="MOM" & error$calibration=="PEST_exm_w1", "Validation.deepm2"] <- MEFF_MOM
######################################################################################################################

#Deepm2 exm w2
PEST_calib <- file.path(sim_folder, 'Validation/valid_exm_deepm2_w2/output/output.nc')

#Thermocline depth
#Modelled thermocline depths
temp<- get_var(PEST_calib, var_name="temp", reference="surface") 
colClean <- function(x){ colnames(x) <- gsub("temp", "wtr", colnames(x)); x } 
colClean(temp)
temp$DateTime <- as.Date(temp$DateTime, format="%Y-%m-$d")

thermo_depth_model <- ts.thermo.depth(temp, Smin = 0.1, na.rm=TRUE, seasonal=FALSE)  %>% 
  dplyr::rename(td_model = thermo.depth, DateTime = datetime) %>%
  mutate(month = lubridate::month(DateTime)) %>%
  filter(between(month, 4, 9)) %>%
  mutate(year = lubridate::year(DateTime)) %>%
  group_by(year) 

thermo_depth_model$DateTime <- as.Date(thermo_depth_model$DateTime, format="%Y-%m-%d")

#Calculating error (MEF) of thermocline depth for Figure 7
td_merge <- merge(thermo_depth_model, obs_TD, by="DateTime")

for (i in 1:nrow(td_merge)) {
  td_merge$MEFF_1[i]<- ((td_merge$td_model[i]- td_merge$thermo.depth[i])^2)
  td_merge$MEFF_2[i]<- ((td_merge$thermo.depth[i]-mean(td_merge$thermo.depth))^2)
  MEFF_TD<- 1-(sum(td_merge$MEFF_1)/sum(td_merge$MEFF_2))
}

#Adding calculated MEF of deepm2 w2 to error table
error[error$metric=="TD" & error$calibration=="PEST_exm_w2", "Validation.deepm2"] <- MEFF_TD

#Schmidt stability
#Modelled Schmidt stability (deepm2 exm w2 model)
schmidt_stability <- ts.schmidt.stability(temp, bathy, na.rm=TRUE) %>% 
  dplyr::rename(ss_PEST = schmidt.stability)

schmidt_stability$datetime <- as.Date(schmidt_stability$datetime, format="%Y-%m-%d")

#Calculating error (MEF) of schmidt stability for Figure 7
SS_merge <- merge(schmidt_stability_obs, schmidt_stability)

for (i in 1:nrow(SS_merge)) {
  SS_merge$MEFF_1[i]<- ((SS_merge$ss_PEST[i]- SS_merge$schmidt.stability[i])^2)
  SS_merge$MEFF_2[i]<- ((SS_merge$schmidt.stability[i]-mean(SS_merge$schmidt.stability))^2)
  MEFF_SS<- 1-(sum(SS_merge$MEFF_1)/sum(SS_merge$MEFF_2))
}

#Adding calculated Deep mixing 2 exm w2 model SS MEF to error table
error[error$metric=="SS" & error$calibration=="PEST_exm_w2", "Validation.deepm2"] <- MEFF_SS

#MOM
#Oxygen Deep mixing 2 exm w2 model 
oxy <- get_var(PEST_calib, "OXY_oxy", reference="surface", z_out=depths) %>%
  pivot_longer(cols=starts_with("OXY_oxy_"), names_to="Depth", names_prefix="OXY_oxy_", values_to = "OXY_oxy") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) 
oxy$DateTime <- as.Date(oxy$DateTime, format="%Y-%m-%d")

#Merge into one dataset
epi_oxy <- filter(oxy, Depth==1)
hypo_oxy <- filter(oxy, Depth==8)
met_oxy <- filter(oxy, Depth==4)
merge_mod<- merge(epi_oxy, hypo_oxy, by="DateTime") %>%
  dplyr::rename(epi_oxy = OXY_oxy.x, hypo_oxy = OXY_oxy.y)
merge_mod$exp_oxy <- (merge_mod$epi_oxy+merge_mod$hypo_oxy)/2

#Calculate deviation between modelled met oxy and expected met oxy
merge_mod <- merge(merge_mod, met_oxy[, c("DateTime", "OXY_oxy")], by="DateTime")%>%
  dplyr::rename(met_oxy = OXY_oxy)
merge_mod$deviation <- merge_mod$met_oxy - merge_mod$exp_oxy

obs_mom$DateTime <- as.Date(obs_mom$DateTime, format="%Y-%m-%d")
merge_mod$DateTime <- as.Date(merge_mod$DateTime, format="%Y-%m-%d")

#Calculating error (MEF) of metalimnetic oxygen minima for Figure 7
merge_mom <- merge(obs_mom, merge_mod, by="DateTime")

for (i in 1:nrow(merge_mom)) {
  merge_mom$MEFF_1[i]<- ((merge_mom$deviation.y[i]- merge_mom$deviation.x[i])^2)
  merge_mom$MEFF_2[i]<- ((merge_mom$deviation.x[i]-mean(merge_mom$deviation.x))^2)
  MEFF_MOM<- 1-(sum(merge_mom$MEFF_1)/sum(merge_mom$MEFF_2))
}

#Adding calculated Deep mixing 2 exm w2 model MOM MEF to error table
error[error$metric=="MOM" & error$calibration=="PEST_exm_w2", "Validation.deepm2"] <- MEFF_MOM
##########################################################################################################################

#Deepm2 exm w3
PEST_calib <- file.path(sim_folder, 'Validation/valid_exm_deepm2_w3/output/output.nc')

#Thermocline depth
#Modelled thermocline depths
temp<- get_var(PEST_calib, var_name="temp", reference="surface") 
colClean <- function(x){ colnames(x) <- gsub("temp", "wtr", colnames(x)); x } 
colClean(temp)
temp$DateTime <- as.Date(temp$DateTime, format="%Y-%m-$d")

thermo_depth_model <- ts.thermo.depth(temp, Smin = 0.1, na.rm=TRUE, seasonal=FALSE)  %>% 
  dplyr::rename(td_model = thermo.depth, DateTime = datetime) %>%
  mutate(month = lubridate::month(DateTime)) %>%
  filter(between(month, 4, 9)) %>%
  mutate(year = lubridate::year(DateTime)) %>%
  group_by(year) 

thermo_depth_model$DateTime <- as.Date(thermo_depth_model$DateTime, format="%Y-%m-%d")

#Calculating error (MEF) of thermocline depth for Figure 7
td_merge <- merge(thermo_depth_model, obs_TD, by="DateTime")

for (i in 1:nrow(td_merge)) {
  td_merge$MEFF_1[i]<- ((td_merge$td_model[i]- td_merge$thermo.depth[i])^2)
  td_merge$MEFF_2[i]<- ((td_merge$thermo.depth[i]-mean(td_merge$thermo.depth))^2)
  MEFF_TD<- 1-(sum(td_merge$MEFF_1)/sum(td_merge$MEFF_2))
}

#Adding calculated MEF of deepm2 w3 to error table
error[error$metric=="TD" & error$calibration=="PEST_exm_w3", "Validation.deepm2"] <- MEFF_TD

#Schmidt stability
#Modelled Schmidt stability (deepm2 exm w3 model)
schmidt_stability <- ts.schmidt.stability(temp, bathy, na.rm=TRUE) %>% 
  dplyr::rename(ss_PEST = schmidt.stability)

schmidt_stability$datetime <- as.Date(schmidt_stability$datetime, format="%Y-%m-%d")

#Calculating error (MEF) of schmidt stability for Figure 7
SS_merge <- merge(schmidt_stability_obs, schmidt_stability)

for (i in 1:nrow(SS_merge)) {
  SS_merge$MEFF_1[i]<- ((SS_merge$ss_PEST[i]- SS_merge$schmidt.stability[i])^2)
  SS_merge$MEFF_2[i]<- ((SS_merge$schmidt.stability[i]-mean(SS_merge$schmidt.stability))^2)
  MEFF_SS<- 1-(sum(SS_merge$MEFF_1)/sum(SS_merge$MEFF_2))
}

#Adding calculated Deep mixing 2 exm w3 model SS MEF to error table
error[error$metric=="SS" & error$calibration=="PEST_exm_w3", "Validation.deepm2"] <- MEFF_SS

#MOM
#Oxygen Deep mixing 2 exm w3 model 
oxy <- get_var(PEST_calib, "OXY_oxy", reference="surface", z_out=depths) %>%
  pivot_longer(cols=starts_with("OXY_oxy_"), names_to="Depth", names_prefix="OXY_oxy_", values_to = "OXY_oxy") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) 
oxy$DateTime <- as.Date(oxy$DateTime, format="%Y-%m-%d")

#Merge into one dataset
epi_oxy <- filter(oxy, Depth==1)
hypo_oxy <- filter(oxy, Depth==8)
met_oxy <- filter(oxy, Depth==4)
merge_mod<- merge(epi_oxy, hypo_oxy, by="DateTime") %>%
  dplyr::rename(epi_oxy = OXY_oxy.x, hypo_oxy = OXY_oxy.y)
merge_mod$exp_oxy <- (merge_mod$epi_oxy+merge_mod$hypo_oxy)/2

#Calculate deviation between modelled met oxy and expected met oxy
merge_mod <- merge(merge_mod, met_oxy[, c("DateTime", "OXY_oxy")], by="DateTime")%>%
  dplyr::rename(met_oxy = OXY_oxy)
merge_mod$deviation <- merge_mod$met_oxy - merge_mod$exp_oxy

obs_mom$DateTime <- as.Date(obs_mom$DateTime, format="%Y-%m-%d")
merge_mod$DateTime <- as.Date(merge_mod$DateTime, format="%Y-%m-%d")

#Calculating error (MEF) of metalimnetic oxygen minima for Figure 7
merge_mom <- merge(obs_mom, merge_mod, by="DateTime")

for (i in 1:nrow(merge_mom)) {
  merge_mom$MEFF_1[i]<- ((merge_mom$deviation.y[i]- merge_mom$deviation.x[i])^2)
  merge_mom$MEFF_2[i]<- ((merge_mom$deviation.x[i]-mean(merge_mom$deviation.x))^2)
  MEFF_MOM<- 1-(sum(merge_mom$MEFF_1)/sum(merge_mom$MEFF_2))
}

#Adding calculated Deep mixing 2 exm w3 model MOM MEF to error table
error[error$metric=="MOM" & error$calibration=="PEST_exm_w3", "Validation.deepm2"] <- MEFF_MOM
##############################################################################################################################

#Anoxia
#Deepm2 naive
output <- nc_open('Validation/valid_naive_deepm2/output/output.nc')

oxy<- ncvar_get(output, "OXY_oxy")
depth<- ncvar_get(output, "z")
depth[depth >= 100] <- NA
tallest_layer <- ncvar_get(output, "NS")

# Set every column of depth dataset ascending- result: 'out'
out <- apply(depth, 2, sort, decreasing=TRUE, na.last=TRUE)

#Loop for depth reference surface- result: 'new'
new <- out
for (i in 1:nrow(new)) {
  for (j in 1:ncol(new)) {
    new[i-1, j] <- out[1, j] - out[i, j]
  }
}

#Maximum depth
for (i in 1:dim(tallest_layer)) {
  new[tallest_layer[i], i] <- out[1, i]
}

#Minimum depth (surface)
new <- rbind(seq(0, 0, length.out = ncol(new)), new)

# Reverse every column of oxy dataset- result: 'oxy_out'
oxy_out <- apply(oxy, 2, rev)

# Putting NAs at the end of each column 
na_fun <- function(x) {
  num.na <- sum(is.na(x))
  x <- x[!is.na(x)]
  x <- c(x, rep(NA, num.na))
  return(x)
}

# apply na_fun over each column in the oxy_out dataframe
oxy_out <- apply(oxy_out, 2, na_fun)
oxy_out <-rbind(oxy_out[1,], oxy_out)

#Melting data into one column, creating dataframe depth, oxy
co<- melt(oxy_out, na.rm=TRUE)
co <- na.omit(co)
ct <- melt(new, na.rm=TRUE)
ct<- na.omit(ct)
df <- cbind(co$Var2, co$value, ct$value)

#Creating dataframe for time 
time <- data.frame(seq(as.Date("2015-07-20"), as.Date("2016-11-30"), by="day"))
ID <- seq.int(1:500)
time <- cbind(ID, time)
colnames(time) <- c("ID", "DateTime")
colnames(df) <- c("ID", "Oxy", "Depth")

#Merge time, depth, oxy according to ID
merge <- merge(time, df, all=TRUE)
merge$DateTime <- as.Date(merge$DateTime, format="%Y-%m-%d")

#Oxygen estimator function at any depth on certain date
estimate_oxy_by_date <- function(target_date, target_depth) {
  data_for_date <- merge %>% 
    filter(DateTime == target_date) %>%
    arrange(Depth)
  
  approx(data_for_date$Depth, data_for_date$Oxy, xout = target_depth)$y
}

#Interpolate DO in 0.02m increments in the water column
oxy_interp_depth <- crossing(
  tibble(DateTime = unique(merge$DateTime)),
  tibble(Depth = seq(0, 9.2, by = 0.1))
) %>%
  group_by(DateTime) %>%
  mutate(Oxy = estimate_oxy_by_date(DateTime[1], Depth))

#Model anoxia 
anoxia <- oxy_interp_depth
anoxia$Oxy <- anoxia$Oxy * 32/1000
anoxia$Oxy<- ifelse(anoxia$Oxy<=1, 1, 0)
anoxia <- na.omit(anoxia)
anoxia$DateTime <- as.Date(anoxia$DateTime, format="%Y-%m-%d")

#Calculating error (MEF) of number of anoxic layers per day for Figure 7
#Observations
obs_anoxic_layers <- read.csv('Observations/anoxia_observed.csv') %>%
  mutate(DateTime=as.Date(DateTime, format="%Y-%m-%d"))

#Modelled
anoxia <- mutate(anoxia, month = lubridate::month(DateTime)) %>%
  filter(between(month, 5, 11))

#Creating empty dataframe for loop
uniqueDates  <- unique(anoxia$DateTime)

newData_mod  <- data.frame(
  DateTime = unique(anoxia$DateTime),
  Count = length(uniqueDates)
)

#Calculating number of anoxic layers each day
for(i in 1:length(uniqueDates)){
  
  filteredData  <- filter(anoxia, DateTime==uniqueDates[i] & Oxy==1)
  newData_mod$Count[i] <- nrow(filteredData)
}

newData_mod$DateTime <- as.Date(newData_mod$DateTime, format="%Y-%m-%d")
merge_anoxia <- merge(obs_anoxic_layers, newData_mod, by="DateTime")

for (i in 1:nrow(merge_anoxia)) {
  merge_anoxia$MEFF_1[i]<- ((merge_anoxia$Count.y[i]- merge_anoxia$Count.x[i])^2)
  merge_anoxia$MEFF_2[i]<- ((merge_anoxia$Count.x[i]-mean(merge_anoxia$Count.x))^2)
  MEFF_anoxia<- 1-(sum(merge_anoxia$MEFF_1)/sum(merge_anoxia$MEFF_2))
}

#Adding calculated MEF to error table
#Deepm2 naive
error[error$metric=="A" & error$calibration=="PEST_N", "Validation.deepm2"] <- MEFF_anoxia
##########################################################################################################################

#Deepm2 exm w1
output <- nc_open('Validation/valid_exm_deepm2_w1/output/output.nc')

oxy<- ncvar_get(output, "OXY_oxy")
depth<- ncvar_get(output, "z")
depth[depth >= 100] <- NA
tallest_layer <- ncvar_get(output, "NS")

# Set every column of depth dataset ascending- result: 'out'
out <- apply(depth, 2, sort, decreasing=TRUE, na.last=TRUE)

#Loop for depth reference surface- result: 'new'
new <- out
for (i in 1:nrow(new)) {
  for (j in 1:ncol(new)) {
    new[i-1, j] <- out[1, j] - out[i, j]
  }
}

#Maximum depth
for (i in 1:dim(tallest_layer)) {
  new[tallest_layer[i], i] <- out[1, i]
}

#Minimum depth (surface)
new <- rbind(seq(0, 0, length.out = ncol(new)), new)

# Reverse every column of oxy dataset- result: 'oxy_out'
oxy_out <- apply(oxy, 2, rev)

# Putting NAs at the end of each column 
na_fun <- function(x) {
  num.na <- sum(is.na(x))
  x <- x[!is.na(x)]
  x <- c(x, rep(NA, num.na))
  return(x)
}

# apply na_fun over each column in the oxy_out dataframe
oxy_out <- apply(oxy_out, 2, na_fun)
oxy_out <-rbind(oxy_out[1,], oxy_out)

#Melting data into one column, creating dataframe depth, oxy
co<- melt(oxy_out, na.rm=TRUE)
co <- na.omit(co)
ct <- melt(new, na.rm=TRUE)
ct<- na.omit(ct)
df <- cbind(co$Var2, co$value, ct$value)

#Creating dataframe for time 
time <- data.frame(seq(as.Date("2015-07-20"), as.Date("2016-11-30"), by="day"))
ID <- seq.int(1:502)
time <- cbind(ID, time)
colnames(time) <- c("ID", "DateTime")
colnames(df) <- c("ID", "Oxy", "Depth")

#Merge time, depth, oxy according to ID
merge <- merge(time, df, all=TRUE)
merge$DateTime <- as.Date(merge$DateTime, format="%Y-%m-%d")

#Oxygen estimator function at any depth on certain date
estimate_oxy_by_date <- function(target_date, target_depth) {
  data_for_date <- merge %>% 
    filter(DateTime == target_date) %>%
    arrange(Depth)
  
  approx(data_for_date$Depth, data_for_date$Oxy, xout = target_depth)$y
}

#Interpolate DO in 0.02m increments in the water column
oxy_interp_depth <- crossing(
  tibble(DateTime = unique(merge$DateTime)),
  tibble(Depth = seq(0, 9.2, by = 0.1))
) %>%
  group_by(DateTime) %>%
  mutate(Oxy = estimate_oxy_by_date(DateTime[1], Depth))

#Model anoxia 
anoxia <- oxy_interp_depth
anoxia$Oxy <- anoxia$Oxy * 32/1000
anoxia$Oxy<- ifelse(anoxia$Oxy<=1, 1, 0)
anoxia <- na.omit(anoxia)
anoxia$DateTime <- as.Date(anoxia$DateTime, format="%Y-%m-%d")

#Calculating error (MEF) of number of anoxic layers per day for Figure 7
#Observations
obs_anoxic_layers <- read.csv('Observations/anoxia_observed.csv') %>%
  mutate(DateTime=as.Date(DateTime, format="%Y-%m-%d"))

#Modelled
anoxia <- mutate(anoxia, month = lubridate::month(DateTime)) %>%
  filter(between(month, 5, 11))

#Creating empty dataframe for loop
uniqueDates  <- unique(anoxia$DateTime)

newData_mod  <- data.frame(
  DateTime = unique(anoxia$DateTime),
  Count = length(uniqueDates)
)

#Calculating number of anoxic layers each day
for(i in 1:length(uniqueDates)){
  
  filteredData  <- filter(anoxia, DateTime==uniqueDates[i] & Oxy==1)
  newData_mod$Count[i] <- nrow(filteredData)
}

newData_mod$DateTime <- as.Date(newData_mod$DateTime, format="%Y-%m-%d")
merge_anoxia <- merge(obs_anoxic_layers, newData_mod, by="DateTime")

for (i in 1:nrow(merge_anoxia)) {
  merge_anoxia$MEFF_1[i]<- ((merge_anoxia$Count.y[i]- merge_anoxia$Count.x[i])^2)
  merge_anoxia$MEFF_2[i]<- ((merge_anoxia$Count.x[i]-mean(merge_anoxia$Count.x))^2)
  MEFF_anoxia<- 1-(sum(merge_anoxia$MEFF_1)/sum(merge_anoxia$MEFF_2))
}

#Deepm2 exm w1
error[error$metric=="A" & error$calibration=="PEST_exm_w1", "Validation.deepm2"] <- MEFF_anoxia
###############################################################################################################################

#Deepm2 exm w2
output <- nc_open('Validation/valid_exm_deepm2_w2/output/output.nc')

oxy<- ncvar_get(output, "OXY_oxy")
depth<- ncvar_get(output, "z")
depth[depth >= 100] <- NA
tallest_layer <- ncvar_get(output, "NS")

# Set every column of depth dataset ascending- result: 'out'
out <- apply(depth, 2, sort, decreasing=TRUE, na.last=TRUE)

#Loop for depth reference surface- result: 'new'
new <- out
for (i in 1:nrow(new)) {
  for (j in 1:ncol(new)) {
    new[i-1, j] <- out[1, j] - out[i, j]
  }
}

#Maximum depth
for (i in 1:dim(tallest_layer)) {
  new[tallest_layer[i], i] <- out[1, i]
}

#Minimum depth (surface)
new <- rbind(seq(0, 0, length.out = ncol(new)), new)

# Reverse every column of oxy dataset- result: 'oxy_out'
oxy_out <- apply(oxy, 2, rev)

# Putting NAs at the end of each column 
na_fun <- function(x) {
  num.na <- sum(is.na(x))
  x <- x[!is.na(x)]
  x <- c(x, rep(NA, num.na))
  return(x)
}

# apply na_fun over each column in the oxy_out dataframe
oxy_out <- apply(oxy_out, 2, na_fun)
oxy_out <-rbind(oxy_out[1,], oxy_out)

#Melting data into one column, creating dataframe depth, oxy
co<- melt(oxy_out, na.rm=TRUE)
co <- na.omit(co)
ct <- melt(new, na.rm=TRUE)
ct<- na.omit(ct)
df <- cbind(co$Var2, co$value, ct$value)

#Creating dataframe for time 
time <- data.frame(seq(as.Date("2015-07-20"), as.Date("2016-11-30"), by="day"))
ID <- seq.int(1:502)
time <- cbind(ID, time)
colnames(time) <- c("ID", "DateTime")
colnames(df) <- c("ID", "Oxy", "Depth")

#Merge time, depth, oxy according to ID
merge <- merge(time, df, all=TRUE)
merge$DateTime <- as.Date(merge$DateTime, format="%Y-%m-%d")

#Oxygen estimator function at any depth on certain date
estimate_oxy_by_date <- function(target_date, target_depth) {
  data_for_date <- merge %>% 
    filter(DateTime == target_date) %>%
    arrange(Depth)
  
  approx(data_for_date$Depth, data_for_date$Oxy, xout = target_depth)$y
}

#Interpolate DO in 0.02m increments in the water column
oxy_interp_depth <- crossing(
  tibble(DateTime = unique(merge$DateTime)),
  tibble(Depth = seq(0, 9.2, by = 0.1))
) %>%
  group_by(DateTime) %>%
  mutate(Oxy = estimate_oxy_by_date(DateTime[1], Depth))

#Model anoxia 
anoxia <- oxy_interp_depth
anoxia$Oxy <- anoxia$Oxy * 32/1000
anoxia$Oxy<- ifelse(anoxia$Oxy<=1, 1, 0)
anoxia <- na.omit(anoxia)
anoxia$DateTime <- as.Date(anoxia$DateTime, format="%Y-%m-%d")

#Calculating error (MEF) of number of anoxic layers per day for Figure 7
#Observations
obs_anoxic_layers <- read.csv('Observations/anoxia_observed.csv') %>%
  mutate(DateTime=as.Date(DateTime, format="%Y-%m-%d"))

#Modelled
anoxia <- mutate(anoxia, month = lubridate::month(DateTime)) %>%
  filter(between(month, 5, 11))

#Creating empty dataframe for loop
uniqueDates  <- unique(anoxia$DateTime)

newData_mod  <- data.frame(
  DateTime = unique(anoxia$DateTime),
  Count = length(uniqueDates)
)

#Calculating number of anoxic layers each day
for(i in 1:length(uniqueDates)){
  
  filteredData  <- filter(anoxia, DateTime==uniqueDates[i] & Oxy==1)
  newData_mod$Count[i] <- nrow(filteredData)
}

newData_mod$DateTime <- as.Date(newData_mod$DateTime, format="%Y-%m-%d")
merge_anoxia <- merge(obs_anoxic_layers, newData_mod, by="DateTime")

for (i in 1:nrow(merge_anoxia)) {
  merge_anoxia$MEFF_1[i]<- ((merge_anoxia$Count.y[i]- merge_anoxia$Count.x[i])^2)
  merge_anoxia$MEFF_2[i]<- ((merge_anoxia$Count.x[i]-mean(merge_anoxia$Count.x))^2)
  MEFF_anoxia<- 1-(sum(merge_anoxia$MEFF_1)/sum(merge_anoxia$MEFF_2))
}

#Deepm2 exm w2
error[error$metric=="A" & error$calibration=="PEST_exm_w2", "Validation.deepm2"] <- MEFF_anoxia
####################################################################################################

#Deepm2 exm w3
output <- nc_open('Validation/valid_exm_deepm2_w3/output/output.nc')

oxy<- ncvar_get(output, "OXY_oxy")
depth<- ncvar_get(output, "z")
depth[depth >= 100] <- NA
tallest_layer <- ncvar_get(output, "NS")

# Set every column of depth dataset ascending- result: 'out'
out <- apply(depth, 2, sort, decreasing=TRUE, na.last=TRUE)

#Loop for depth reference surface- result: 'new'
new <- out
for (i in 1:nrow(new)) {
  for (j in 1:ncol(new)) {
    new[i-1, j] <- out[1, j] - out[i, j]
  }
}

#Maximum depth
for (i in 1:dim(tallest_layer)) {
  new[tallest_layer[i], i] <- out[1, i]
}

#Minimum depth (surface)
new <- rbind(seq(0, 0, length.out = ncol(new)), new)

# Reverse every column of oxy dataset- result: 'oxy_out'
oxy_out <- apply(oxy, 2, rev)

# Putting NAs at the end of each column 
na_fun <- function(x) {
  num.na <- sum(is.na(x))
  x <- x[!is.na(x)]
  x <- c(x, rep(NA, num.na))
  return(x)
}

# apply na_fun over each column in the oxy_out dataframe
oxy_out <- apply(oxy_out, 2, na_fun)
oxy_out <-rbind(oxy_out[1,], oxy_out)

#Melting data into one column, creating dataframe depth, oxy
co<- melt(oxy_out, na.rm=TRUE)
co <- na.omit(co)
ct <- melt(new, na.rm=TRUE)
ct<- na.omit(ct)
df <- cbind(co$Var2, co$value, ct$value)

#Creating dataframe for time 
time <- data.frame(seq(as.Date("2015-07-20"), as.Date("2016-11-30"), by="day"))
ID <- seq.int(1:502)
time <- cbind(ID, time)
colnames(time) <- c("ID", "DateTime")
colnames(df) <- c("ID", "Oxy", "Depth")

#Merge time, depth, oxy according to ID
merge <- merge(time, df, all=TRUE)
merge$DateTime <- as.Date(merge$DateTime, format="%Y-%m-%d")

#Oxygen estimator function at any depth on certain date
estimate_oxy_by_date <- function(target_date, target_depth) {
  data_for_date <- merge %>% 
    filter(DateTime == target_date) %>%
    arrange(Depth)
  
  approx(data_for_date$Depth, data_for_date$Oxy, xout = target_depth)$y
}

#Interpolate DO in 0.02m increments in the water column
oxy_interp_depth <- crossing(
  tibble(DateTime = unique(merge$DateTime)),
  tibble(Depth = seq(0, 9.2, by = 0.1))
) %>%
  group_by(DateTime) %>%
  mutate(Oxy = estimate_oxy_by_date(DateTime[1], Depth))

#Model anoxia 
anoxia <- oxy_interp_depth
anoxia$Oxy <- anoxia$Oxy * 32/1000
anoxia$Oxy<- ifelse(anoxia$Oxy<=1, 1, 0)
anoxia <- na.omit(anoxia)
anoxia$DateTime <- as.Date(anoxia$DateTime, format="%Y-%m-%d")

#Calculating error (MEF) of number of anoxic layers per day for Figure 7
#Observations
obs_anoxic_layers <- read.csv('Observations/anoxia_observed.csv') %>%
  mutate(DateTime=as.Date(DateTime, format="%Y-%m-%d"))

#Modelled
anoxia <- mutate(anoxia, month = lubridate::month(DateTime)) %>%
  filter(between(month, 5, 11))

#Creating empty dataframe for loop
uniqueDates  <- unique(anoxia$DateTime)

newData_mod  <- data.frame(
  DateTime = unique(anoxia$DateTime),
  Count = length(uniqueDates)
)

#Calculating number of anoxic layers each day
for(i in 1:length(uniqueDates)){
  
  filteredData  <- filter(anoxia, DateTime==uniqueDates[i] & Oxy==1)
  newData_mod$Count[i] <- nrow(filteredData)
}

newData_mod$DateTime <- as.Date(newData_mod$DateTime, format="%Y-%m-%d")
merge_anoxia <- merge(obs_anoxic_layers, newData_mod, by="DateTime")

for (i in 1:nrow(merge_anoxia)) {
  merge_anoxia$MEFF_1[i]<- ((merge_anoxia$Count.y[i]- merge_anoxia$Count.x[i])^2)
  merge_anoxia$MEFF_2[i]<- ((merge_anoxia$Count.x[i]-mean(merge_anoxia$Count.x))^2)
  MEFF_anoxia<- 1-(sum(merge_anoxia$MEFF_1)/sum(merge_anoxia$MEFF_2))
}

#Deepm2 exm w3
error[error$metric=="A" & error$calibration=="PEST_exm_w3", "Validation.deepm2"] <- MEFF_anoxia

write.csv(error, 'Observations/error_stats.csv', row.names=FALSE)

#############################################################################################################################
#Deep mixing 1
setwd(".../FCR-GLM-metrics")
sim_folder <- getwd()
error <- read.csv('Observations/error_stats.csv')

#oxygen
var="OXY_oxy"
obs_oxy<-read.csv('Observations/CleanedObsOxy.csv') %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) 
obs_oxy$DateTime <- as.Date(obs_oxy$DateTime, format="%Y-%m-%d")

depths<- c(0.1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9.2) 

#Deep mixing 1 naive
new<- file.path(sim_folder, 'Validation/valid_naive_deepm1/output/output.nc')

new_oxy <- get_var(new, "OXY_oxy", reference="surface", z_out=depths) %>%
  pivot_longer(cols=starts_with("OXY_oxy_"), names_to="Depth", names_prefix="OXY_oxy_", values_to = "OXY_oxy") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) 
new_oxy$DateTime <- as.Date(new_oxy$DateTime, format="%Y-%m-%d")

oxygen <- merge(obs_oxy, new_oxy, by=c("DateTime","Depth")) %>%
  dplyr::rename(obsoxy = OXY_oxy.x, mod_oxy = OXY_oxy.y) 
oxygen$DateTime <- as.Date(oxygen$DateTime, format="%Y-%m-%d")

#MEF
for (i in 1:nrow(oxygen)) {
  oxygen$MEF_1[i]<- ((oxygen$mod_oxy[i]- oxygen$obsoxy[i])^2)
  oxygen$MEF_2[i]<- ((oxygen$obsoxy[i]-mean(oxygen$obsoxy))^2)
  MEF_oxy<- 1-(sum(oxygen$MEF_1)/sum(oxygen$MEF_2))
}

#deep mixing 1 naive
error[error$metric=="oxy" & error$calibration=="PEST_N", "Validation.deepm1"] <- MEF_oxy

#Deep mixing 1 exm w1
new_w1 <- file.path(sim_folder, 'Validation/valid_exm_deepm1_w1/output/output.nc')

new_oxy_w1 <- get_var(new_w1, "OXY_oxy", reference="surface", z_out=depths) %>%
  pivot_longer(cols=starts_with("OXY_oxy_"), names_to="Depth", names_prefix="OXY_oxy_", values_to = "OXY_oxy") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) 
new_oxy_w1$DateTime <- as.Date(new_oxy_w1$DateTime, format="%Y-%m-%d")

oxygen_w1 <- merge(obs_oxy, new_oxy_w1, by=c("DateTime","Depth")) %>%
  dplyr::rename(obsoxy = OXY_oxy.x, mod_oxy = OXY_oxy.y) 
oxygen_w1$DateTime <- as.Date(oxygen_w1$DateTime, format="%Y-%m-%d")

#MEF
for (i in 1:nrow(oxygen_w1)) {
  oxygen_w1$MEF_1[i]<- ((oxygen_w1$mod_oxy[i]- oxygen_w1$obsoxy[i])^2)
  oxygen_w1$MEF_2[i]<- ((oxygen_w1$obsoxy[i]-mean(oxygen_w1$obsoxy))^2)
  MEF_oxy_w1<- 1-(sum(oxygen_w1$MEF_1)/sum(oxygen_w1$MEF_2))
}

error[error$metric=="oxy" & error$calibration=="PEST_exm_w1", "Validation.deepm1"] <- MEF_oxy_w1

#Deep mixing 1 exm w2
new_w2 <- file.path(sim_folder, 'Validation/valid_exm_deepm1_w2/output/output.nc')

new_oxy_w2 <- get_var(new_w2, "OXY_oxy", reference="surface", z_out=depths) %>%
  pivot_longer(cols=starts_with("OXY_oxy_"), names_to="Depth", names_prefix="OXY_oxy_", values_to = "OXY_oxy") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) 
new_oxy_w2$DateTime <- as.Date(new_oxy_w2$DateTime, format="%Y-%m-%d")

oxygen_w2 <- merge(obs_oxy, new_oxy_w2, by=c("DateTime","Depth")) %>%
  dplyr::rename(obsoxy = OXY_oxy.x, mod_oxy = OXY_oxy.y) 
oxygen_w2$DateTime <- as.Date(oxygen_w2$DateTime, format="%Y-%m-%d")

#MEF
for (i in 1:nrow(oxygen_w2)) {
  oxygen_w2$MEF_1[i]<- ((oxygen_w2$mod_oxy[i]- oxygen_w2$obsoxy[i])^2)
  oxygen_w2$MEF_2[i]<- ((oxygen_w2$obsoxy[i]-mean(oxygen_w2$obsoxy))^2)
  MEF_oxy_w2<- 1-(sum(oxygen_w2$MEF_1)/sum(oxygen_w2$MEF_2))
}

error[error$metric=="oxy" & error$calibration=="PEST_exm_w2", "Validation.deepm1"] <- MEF_oxy_w2

#Deep mixing 1 exm w3
new_w3 <- file.path(sim_folder, 'Validation/valid_exm_deepm1_w3/output/output.nc')

new_oxy_w3 <- get_var(new_w3, "OXY_oxy", reference="surface", z_out=depths) %>%
  pivot_longer(cols=starts_with("OXY_oxy_"), names_to="Depth", names_prefix="OXY_oxy_", values_to = "OXY_oxy") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) 
new_oxy_w3$DateTime <- as.Date(new_oxy_w3$DateTime, format="%Y-%m-%d")

oxygen_w3 <- merge(obs_oxy, new_oxy_w3, by=c("DateTime","Depth")) %>%
  dplyr::rename(obsoxy = OXY_oxy.x, mod_oxy = OXY_oxy.y) 
oxygen_w3$DateTime <- as.Date(oxygen_w3$DateTime, format="%Y-%m-%d")

#MEF
for (i in 1:nrow(oxygen_w3)) {
  oxygen_w3$MEF_1[i]<- ((oxygen_w3$mod_oxy[i]- oxygen_w3$obsoxy[i])^2)
  oxygen_w3$MEF_2[i]<- ((oxygen_w3$obsoxy[i]-mean(oxygen_w3$obsoxy))^2)
  MEF_oxy_w3<- 1-(sum(oxygen_w3$MEF_1)/sum(oxygen_w3$MEF_2))
}

error[error$metric=="oxy" & error$calibration=="PEST_exm_w3", "Validation.deepm1"] <- MEF_oxy_w3

#Temperature
var="temp"
obs_temp<-read.csv('Observations/CleanedObsTemp.csv') %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) 
obs_temp$DateTime <- as.Date(obs_temp$DateTime, format="%Y-%m-%d")

depths<- c(0.1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9.2) 

#Deep mixing 1 naive
new<- file.path(sim_folder, 'Validation/valid_naive_deepm1/output/output.nc')

new_temp <- get_var(new, "temp", reference="surface", z_out=depths) %>%
  pivot_longer(cols=starts_with("temp_"), names_to="Depth", names_prefix="temp_", values_to = "temp") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) 
new_temp$DateTime <- as.Date(new_temp$DateTime, format="%Y-%m-%d")

temp <- merge(obs_temp, new_temp, by=c("DateTime","Depth")) %>%
  dplyr::rename(obtemp = temp.x, modtemp = temp.y) 
temp$DateTime <- as.Date(temp$DateTime, format="%Y-%m-%d")

for (i in 1:nrow(temp)) {
  temp$MEF_1[i]<- ((temp$modtemp[i]- temp$obtemp[i])^2)
  temp$MEF_2[i]<- ((temp$obtemp[i]-mean(temp$obtemp))^2)
  MEF_temp<- 1-(sum(temp$MEF_1)/sum(temp$MEF_2))
}

error[error$metric=="temp" & error$calibration=="PEST_N", "Validation.deepm1"] <- MEF_temp

#Deep mixing 1 exm w1
new_w1 <- file.path(sim_folder, 'Validation/valid_exm_deepm1_w1/output/output.nc')

new_temp_w1 <- get_var(new_w1, "temp", reference="surface", z_out=depths) %>%
  pivot_longer(cols=starts_with("temp_"), names_to="Depth", names_prefix="temp_", values_to = "temp") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) 
new_temp_w1$DateTime <- as.Date(new_temp_w1$DateTime, format="%Y-%m-%d")

temp_w1 <- merge(obs_temp, new_temp_w1, by=c("DateTime","Depth")) %>%
  dplyr::rename(obtemp = temp.x, modtemp = temp.y) 
temp_w1$DateTime <- as.Date(temp_w1$DateTime, format="%Y-%m-%d")

for (i in 1:nrow(temp_w1)) {
  temp_w1$MEF_1[i]<- ((temp_w1$modtemp[i]- temp_w1$obtemp[i])^2)
  temp_w1$MEF_2[i]<- ((temp_w1$obtemp[i]-mean(temp_w1$obtemp))^2)
  MEF_temp_w1<- 1-(sum(temp_w1$MEF_1)/sum(temp_w1$MEF_2))
}

error[error$metric=="temp" & error$calibration=="PEST_exm_w1", "Validation.deepm1"] <- MEF_temp_w1

#Deep mixing 1 exm w2
new_w2 <- file.path(sim_folder, 'Validation/valid_exm_deepm1_w2/output/output.nc')

new_temp_w2 <- get_var(new_w2, "temp", reference="surface", z_out=depths) %>%
  pivot_longer(cols=starts_with("temp_"), names_to="Depth", names_prefix="temp_", values_to = "temp") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) 
new_temp_w2$DateTime <- as.Date(new_temp_w2$DateTime, format="%Y-%m-%d")

temp_w2 <- merge(obs_temp, new_temp_w2, by=c("DateTime","Depth")) %>%
  dplyr::rename(obtemp = temp.x, modtemp = temp.y) 
temp_w2$DateTime <- as.Date(temp_w2$DateTime, format="%Y-%m-%d")

for (i in 1:nrow(temp_w2)) {
  temp_w2$MEF_1[i]<- ((temp_w2$modtemp[i]- temp_w2$obtemp[i])^2)
  temp_w2$MEF_2[i]<- ((temp_w2$obtemp[i]-mean(temp_w2$obtemp))^2)
  MEF_temp_w2<- 1-(sum(temp_w2$MEF_1)/sum(temp_w2$MEF_2))
}

error[error$metric=="temp" & error$calibration=="PEST_exm_w2", "Validation.deepm1"] <- MEF_temp_w2

#Deep mixing 1 exm w3
new_w3 <- file.path(sim_folder, 'Validation/valid_exm_deepm1_w3/output/output.nc')

new_temp_w3 <- get_var(new_w3, "temp", reference="surface", z_out=depths) %>%
  pivot_longer(cols=starts_with("temp_"), names_to="Depth", names_prefix="temp_", values_to = "temp") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) 
new_temp_w3$DateTime <- as.Date(new_temp_w3$DateTime, format="%Y-%m-%d")

temp_w3 <- merge(obs_temp, new_temp_w3, by=c("DateTime","Depth")) %>%
  dplyr::rename(obtemp = temp.x, modtemp = temp.y) 
temp_w3$DateTime <- as.Date(temp_w3$DateTime, format="%Y-%m-%d")

for (i in 1:nrow(temp_w3)) {
  temp_w3$MEF_1[i]<- ((temp_w3$modtemp[i]- temp_w3$obtemp[i])^2)
  temp_w3$MEF_2[i]<- ((temp_w3$obtemp[i]-mean(temp_w3$obtemp))^2)
  MEF_temp_w3<- 1-(sum(temp_w3$MEF_1)/sum(temp_w3$MEF_2))
}

error[error$metric=="temp" & error$calibration=="PEST_exm_w3", "Validation.deepm1"] <- MEF_temp_w3

#TD
#Observed thermocline deph
obs_TD<- read.csv('Observations/obs_td.csv') %>%
  mutate(DateTime = as.Date(DateTime, format = "%Y-%m-%d"))%>%
  mutate(month = lubridate::month(DateTime)) %>%
  filter(between(month, 4, 9)) %>%
  mutate(year = lubridate::year(DateTime)) %>%
  group_by(year) %>%
  filter(DateTime  > "2015-07-19" & DateTime < "2016-12-02")

obs_TD$DateTime <- as.Date(obs_TD$DateTime, format="%Y-%m-%d")

#Schmidt stability
bathy <- read.csv('Observations/bathymetry.csv')

#Schmidt stability observed
schmidt_stability_obs<- read.csv("Observations/Obs_SS.csv") %>%
  mutate(datetime = as.Date(datetime, format = '%Y-%m-%d')) %>%
  filter(datetime > "2015-07-19" & datetime < "2016-12-02")
schmidt_stability_obs$datetime<- as.Date(schmidt_stability_obs$datetime, format='%Y-%m-%d')

#MOM Observed
obs_mom<-read.csv('Observations/mom_observed.csv') %>%
  mutate(DateTime = as.Date(DateTime, format = "%Y-%m-%d"))%>%
  filter(DateTime > "2015-07-19" & DateTime < "2016-12-02")
obs_mom$DateTime <- as.Date(obs_mom$DateTime, format = "%Y-%m-%d")
depths<- c(1, 4, 8) 

###############################################################################################################
#Deepm1 naive
PEST_calib <- file.path(sim_folder, 'Validation/valid_naive_deepm1/output/output.nc')

#Thermocline depth
#Modelled thermocline depth Deep2 naive
temp<- get_var(PEST_calib, var_name="temp", reference="surface") 
colClean <- function(x){ colnames(x) <- gsub("temp", "wtr", colnames(x)); x } 
colClean(temp)
temp$DateTime <- as.Date(temp$DateTime, format="%Y-%m-$d")

thermo_depth_model <- ts.thermo.depth(temp, Smin = 0.1, na.rm=TRUE, seasonal=FALSE)  %>% 
  dplyr::rename(td_model = thermo.depth, DateTime = datetime) %>%
  mutate(month = lubridate::month(DateTime)) %>%
  filter(between(month, 4, 9)) %>%
  mutate(year = lubridate::year(DateTime)) %>%
  group_by(year) 

thermo_depth_model$DateTime <- as.Date(thermo_depth_model$DateTime, format="%Y-%m-%d")

#Calculating error (MEF) of thermocline depth for Figure 7
td_merge <- merge(thermo_depth_model, obs_TD, by="DateTime")

for (i in 1:nrow(td_merge)) {
  td_merge$MEFF_1[i]<- ((td_merge$td_model[i]- td_merge$thermo.depth[i])^2)
  td_merge$MEFF_2[i]<- ((td_merge$thermo.depth[i]-mean(td_merge$thermo.depth))^2)
  MEFF_TD<- 1-(sum(td_merge$MEFF_1)/sum(td_merge$MEFF_2))
}

#Adding calculated Deep mixing 1 naive model TD MEF to error table
error[error$metric=="TD" & error$calibration=="PEST_N", "Validation.deepm1"] <- MEFF_TD

#Schmidt stability
#Modelled Schmidt stability (deepm2 naive model)
schmidt_stability <- ts.schmidt.stability(temp, bathy, na.rm=TRUE) %>% 
  dplyr::rename(ss_PEST = schmidt.stability)

schmidt_stability$datetime <- as.Date(schmidt_stability$datetime, format="%Y-%m-%d")

#Calculating error (MEF) of schmidt stability for Figure 7
SS_merge <- merge(schmidt_stability_obs, schmidt_stability)

for (i in 1:nrow(SS_merge)) {
  SS_merge$MEFF_1[i]<- ((SS_merge$ss_PEST[i]- SS_merge$schmidt.stability[i])^2)
  SS_merge$MEFF_2[i]<- ((SS_merge$schmidt.stability[i]-mean(SS_merge$schmidt.stability))^2)
  MEFF_SS<- 1-(sum(SS_merge$MEFF_1)/sum(SS_merge$MEFF_2))
}

#Adding calculated Deep mixing 1 naive model SS MEF to error table
error[error$metric=="SS" & error$calibration=="PEST_N", "Validation.deepm1"] <- MEFF_SS

#MOM
#Oxygen deepm1 naive model
oxy <- get_var(PEST_calib, "OXY_oxy", reference="surface", z_out=depths) %>%
  pivot_longer(cols=starts_with("OXY_oxy_"), names_to="Depth", names_prefix="OXY_oxy_", values_to = "OXY_oxy") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) 
oxy$DateTime <- as.Date(oxy$DateTime, format="%Y-%m-%d")

#Merge into one dataset
epi_oxy <- filter(oxy, Depth==1)
hypo_oxy <- filter(oxy, Depth==8)
met_oxy <- filter(oxy, Depth==4)
merge_mod<- merge(epi_oxy, hypo_oxy, by="DateTime") %>%
  dplyr::rename(epi_oxy = OXY_oxy.x, hypo_oxy = OXY_oxy.y)
merge_mod$exp_oxy <- (merge_mod$epi_oxy+merge_mod$hypo_oxy)/2

#Calculate deviation between modelled met oxy and expected met oxy
merge_mod <- merge(merge_mod, met_oxy[, c("DateTime", "OXY_oxy")], by="DateTime")%>%
  dplyr::rename(met_oxy = OXY_oxy)
merge_mod$deviation <- merge_mod$met_oxy - merge_mod$exp_oxy

obs_mom$DateTime <- as.Date(obs_mom$DateTime, format="%Y-%m-%d")
merge_mod$DateTime <- as.Date(merge_mod$DateTime, format="%Y-%m-%d")

#Calculating error (MEF) of metalimnetic oxygen minima for Figure 7
merge_mom <- merge(obs_mom, merge_mod, by="DateTime")

for (i in 1:nrow(merge_mom)) {
  merge_mom$MEFF_1[i]<- ((merge_mom$deviation.y[i]- merge_mom$deviation.x[i])^2)
  merge_mom$MEFF_2[i]<- ((merge_mom$deviation.x[i]-mean(merge_mom$deviation.x))^2)
  MEFF_MOM<- 1-(sum(merge_mom$MEFF_1)/sum(merge_mom$MEFF_2))
}

#Adding calculated Deep mixing 1 naive model MOM MEF to error table
error[error$metric=="MOM" & error$calibration=="PEST_N", "Validation.deepm1"] <- MEFF_MOM

#############################################################################################################

#Deepm1 exm w1
PEST_calib <- file.path(sim_folder, 'Validation/valid_exm_deepm1_w1/output/output.nc')

#Thermocline depth
#Modelled thermocline depth Deepm1 exm w1
temp<- get_var(PEST_calib, var_name="temp", reference="surface") 
colClean <- function(x){ colnames(x) <- gsub("temp", "wtr", colnames(x)); x } 
colClean(temp)
temp$DateTime <- as.Date(temp$DateTime, format="%Y-%m-$d")

thermo_depth_model <- ts.thermo.depth(temp, Smin = 0.1, na.rm=TRUE, seasonal=FALSE)  %>% 
  dplyr::rename(td_model = thermo.depth, DateTime = datetime) %>%
  mutate(month = lubridate::month(DateTime)) %>%
  filter(between(month, 4, 9)) %>%
  mutate(year = lubridate::year(DateTime)) %>%
  group_by(year) 

thermo_depth_model$DateTime <- as.Date(thermo_depth_model$DateTime, format="%Y-%m-%d")

#Calculating error (MEF) of thermocline depth for Figure 7
td_merge <- merge(thermo_depth_model, obs_TD, by="DateTime")

for (i in 1:nrow(td_merge)) {
  td_merge$MEFF_1[i]<- ((td_merge$td_model[i]- td_merge$thermo.depth[i])^2)
  td_merge$MEFF_2[i]<- ((td_merge$thermo.depth[i]-mean(td_merge$thermo.depth))^2)
  MEFF_TD<- 1-(sum(td_merge$MEFF_1)/sum(td_merge$MEFF_2))
}

#Adding calculated MEF of deepm1 w1 to error table
error[error$metric=="TD" & error$calibration=="PEST_exm_w1", "Validation.deepm1"] <- MEFF_TD

#Schmidt stability
#Modelled Schmidt stability (deepm2 exm w1 model)
schmidt_stability <- ts.schmidt.stability(temp, bathy, na.rm=TRUE) %>% 
  dplyr::rename(ss_PEST = schmidt.stability)

schmidt_stability$datetime <- as.Date(schmidt_stability$datetime, format="%Y-%m-%d")

#Calculating error (MEF) of schmidt stability for Figure 7
SS_merge <- merge(schmidt_stability_obs, schmidt_stability)

for (i in 1:nrow(SS_merge)) {
  SS_merge$MEFF_1[i]<- ((SS_merge$ss_PEST[i]- SS_merge$schmidt.stability[i])^2)
  SS_merge$MEFF_2[i]<- ((SS_merge$schmidt.stability[i]-mean(SS_merge$schmidt.stability))^2)
  MEFF_SS<- 1-(sum(SS_merge$MEFF_1)/sum(SS_merge$MEFF_2))
}

#Adding calculated Deep mixing 2 exm w1 model SS MEF to error table
error[error$metric=="SS" & error$calibration=="PEST_exm_w1", "Validation.deepm1"] <- MEFF_SS

#MOM
#Oxygen Deep mixing 1 exm w1 model 
oxy <- get_var(PEST_calib, "OXY_oxy", reference="surface", z_out=depths) %>%
  pivot_longer(cols=starts_with("OXY_oxy_"), names_to="Depth", names_prefix="OXY_oxy_", values_to = "OXY_oxy") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) 
oxy$DateTime <- as.Date(oxy$DateTime, format="%Y-%m-%d")

#Merge into one dataset
epi_oxy <- filter(oxy, Depth==1)
hypo_oxy <- filter(oxy, Depth==8)
met_oxy <- filter(oxy, Depth==4)
merge_mod<- merge(epi_oxy, hypo_oxy, by="DateTime") %>%
  dplyr::rename(epi_oxy = OXY_oxy.x, hypo_oxy = OXY_oxy.y)
merge_mod$exp_oxy <- (merge_mod$epi_oxy+merge_mod$hypo_oxy)/2

#Calculate deviation between modelled met oxy and expected met oxy
merge_mod <- merge(merge_mod, met_oxy[, c("DateTime", "OXY_oxy")], by="DateTime")%>%
  dplyr::rename(met_oxy = OXY_oxy)
merge_mod$deviation <- merge_mod$met_oxy - merge_mod$exp_oxy

obs_mom$DateTime <- as.Date(obs_mom$DateTime, format="%Y-%m-%d")
merge_mod$DateTime <- as.Date(merge_mod$DateTime, format="%Y-%m-%d")

#Calculating error (MEF) of metalimnetic oxygen minima for Figure 7
merge_mom <- merge(obs_mom, merge_mod, by="DateTime")

for (i in 1:nrow(merge_mom)) {
  merge_mom$MEFF_1[i]<- ((merge_mom$deviation.y[i]- merge_mom$deviation.x[i])^2)
  merge_mom$MEFF_2[i]<- ((merge_mom$deviation.x[i]-mean(merge_mom$deviation.x))^2)
  MEFF_MOM<- 1-(sum(merge_mom$MEFF_1)/sum(merge_mom$MEFF_2))
}

#Adding calculated Deep mixing 2 exm w1 model MOM MEF to error table
error[error$metric=="MOM" & error$calibration=="PEST_exm_w1", "Validation.deepm1"] <- MEFF_MOM
######################################################################################################################

#Deepm1 exm w2
PEST_calib <- file.path(sim_folder, 'Validation/valid_exm_deepm1_w2/output/output.nc')

#Thermocline depth
#Modelled thermocline depths
temp<- get_var(PEST_calib, var_name="temp", reference="surface") 
colClean <- function(x){ colnames(x) <- gsub("temp", "wtr", colnames(x)); x } 
colClean(temp)
temp$DateTime <- as.Date(temp$DateTime, format="%Y-%m-$d")

thermo_depth_model <- ts.thermo.depth(temp, Smin = 0.1, na.rm=TRUE, seasonal=FALSE)  %>% 
  dplyr::rename(td_model = thermo.depth, DateTime = datetime) %>%
  mutate(month = lubridate::month(DateTime)) %>%
  filter(between(month, 4, 9)) %>%
  mutate(year = lubridate::year(DateTime)) %>%
  group_by(year) 

thermo_depth_model$DateTime <- as.Date(thermo_depth_model$DateTime, format="%Y-%m-%d")

#Calculating error (MEF) of thermocline depth for Figure 7
td_merge <- merge(thermo_depth_model, obs_TD, by="DateTime")

for (i in 1:nrow(td_merge)) {
  td_merge$MEFF_1[i]<- ((td_merge$td_model[i]- td_merge$thermo.depth[i])^2)
  td_merge$MEFF_2[i]<- ((td_merge$thermo.depth[i]-mean(td_merge$thermo.depth))^2)
  MEFF_TD<- 1-(sum(td_merge$MEFF_1)/sum(td_merge$MEFF_2))
}

#Adding calculated MEF of deepm1 w2 to error table
error[error$metric=="TD" & error$calibration=="PEST_exm_w2", "Validation.deepm1"] <- MEFF_TD

#Schmidt stability
#Modelled Schmidt stability (deepm1 exm w2 model)
schmidt_stability <- ts.schmidt.stability(temp, bathy, na.rm=TRUE) %>% 
  dplyr::rename(ss_PEST = schmidt.stability)

schmidt_stability$datetime <- as.Date(schmidt_stability$datetime, format="%Y-%m-%d")

#Calculating error (MEF) of schmidt stability for Figure 7
SS_merge <- merge(schmidt_stability_obs, schmidt_stability)

for (i in 1:nrow(SS_merge)) {
  SS_merge$MEFF_1[i]<- ((SS_merge$ss_PEST[i]- SS_merge$schmidt.stability[i])^2)
  SS_merge$MEFF_2[i]<- ((SS_merge$schmidt.stability[i]-mean(SS_merge$schmidt.stability))^2)
  MEFF_SS<- 1-(sum(SS_merge$MEFF_1)/sum(SS_merge$MEFF_2))
}

#Adding calculated Deep mixing 1 exm w2 model SS MEF to error table
error[error$metric=="SS" & error$calibration=="PEST_exm_w2", "Validation.deepm1"] <- MEFF_SS

#MOM
#Oxygen Deep mixing 1 exm w2 model 
oxy <- get_var(PEST_calib, "OXY_oxy", reference="surface", z_out=depths) %>%
  pivot_longer(cols=starts_with("OXY_oxy_"), names_to="Depth", names_prefix="OXY_oxy_", values_to = "OXY_oxy") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) 
oxy$DateTime <- as.Date(oxy$DateTime, format="%Y-%m-%d")

#Merge into one dataset
epi_oxy <- filter(oxy, Depth==1)
hypo_oxy <- filter(oxy, Depth==8)
met_oxy <- filter(oxy, Depth==4)
merge_mod<- merge(epi_oxy, hypo_oxy, by="DateTime") %>%
  dplyr::rename(epi_oxy = OXY_oxy.x, hypo_oxy = OXY_oxy.y)
merge_mod$exp_oxy <- (merge_mod$epi_oxy+merge_mod$hypo_oxy)/2

#Calculate deviation between modelled met oxy and expected met oxy
merge_mod <- merge(merge_mod, met_oxy[, c("DateTime", "OXY_oxy")], by="DateTime")%>%
  dplyr::rename(met_oxy = OXY_oxy)
merge_mod$deviation <- merge_mod$met_oxy - merge_mod$exp_oxy

obs_mom$DateTime <- as.Date(obs_mom$DateTime, format="%Y-%m-%d")
merge_mod$DateTime <- as.Date(merge_mod$DateTime, format="%Y-%m-%d")

#Calculating error (MEF) of metalimnetic oxygen minima for Figure 7
merge_mom <- merge(obs_mom, merge_mod, by="DateTime")

for (i in 1:nrow(merge_mom)) {
  merge_mom$MEFF_1[i]<- ((merge_mom$deviation.y[i]- merge_mom$deviation.x[i])^2)
  merge_mom$MEFF_2[i]<- ((merge_mom$deviation.x[i]-mean(merge_mom$deviation.x))^2)
  MEFF_MOM<- 1-(sum(merge_mom$MEFF_1)/sum(merge_mom$MEFF_2))
}

#Adding calculated Deep mixing 1 exm w2 model MOM MEF to error table
error[error$metric=="MOM" & error$calibration=="PEST_exm_w2", "Validation.deepm1"] <- MEFF_MOM
##########################################################################################################################

#Deepm1 exm w3
PEST_calib <- file.path(sim_folder, 'Validation/valid_exm_deepm1_w3/output/output.nc')

#Thermocline depth
#Modelled thermocline depths
temp<- get_var(PEST_calib, var_name="temp", reference="surface") 
colClean <- function(x){ colnames(x) <- gsub("temp", "wtr", colnames(x)); x } 
colClean(temp)
temp$DateTime <- as.Date(temp$DateTime, format="%Y-%m-$d")

thermo_depth_model <- ts.thermo.depth(temp, Smin = 0.1, na.rm=TRUE, seasonal=FALSE)  %>% 
  dplyr::rename(td_model = thermo.depth, DateTime = datetime) %>%
  mutate(month = lubridate::month(DateTime)) %>%
  filter(between(month, 4, 9)) %>%
  mutate(year = lubridate::year(DateTime)) %>%
  group_by(year) 

thermo_depth_model$DateTime <- as.Date(thermo_depth_model$DateTime, format="%Y-%m-%d")

#Calculating error (MEF) of thermocline depth for Figure 7
td_merge <- merge(thermo_depth_model, obs_TD, by="DateTime")

for (i in 1:nrow(td_merge)) {
  td_merge$MEFF_1[i]<- ((td_merge$td_model[i]- td_merge$thermo.depth[i])^2)
  td_merge$MEFF_2[i]<- ((td_merge$thermo.depth[i]-mean(td_merge$thermo.depth))^2)
  MEFF_TD<- 1-(sum(td_merge$MEFF_1)/sum(td_merge$MEFF_2))
}

#Adding calculated MEF of deepm1 w3 to error table
error[error$metric=="TD" & error$calibration=="PEST_exm_w3", "Validation.deepm1"] <- MEFF_TD

#Schmidt stability
#Modelled Schmidt stability (deepm1 exm w3 model)
schmidt_stability <- ts.schmidt.stability(temp, bathy, na.rm=TRUE) %>% 
  dplyr::rename(ss_PEST = schmidt.stability)

schmidt_stability$datetime <- as.Date(schmidt_stability$datetime, format="%Y-%m-%d")

#Calculating error (MEF) of schmidt stability for Figure 7
SS_merge <- merge(schmidt_stability_obs, schmidt_stability)

for (i in 1:nrow(SS_merge)) {
  SS_merge$MEFF_1[i]<- ((SS_merge$ss_PEST[i]- SS_merge$schmidt.stability[i])^2)
  SS_merge$MEFF_2[i]<- ((SS_merge$schmidt.stability[i]-mean(SS_merge$schmidt.stability))^2)
  MEFF_SS<- 1-(sum(SS_merge$MEFF_1)/sum(SS_merge$MEFF_2))
}

#Adding calculated Deep mixing 1 exm w3 model SS MEF to error table
error[error$metric=="SS" & error$calibration=="PEST_exm_w3", "Validation.deepm1"] <- MEFF_SS

#MOM
#Oxygen Deep mixing 1 exm w3 model 
oxy <- get_var(PEST_calib, "OXY_oxy", reference="surface", z_out=depths) %>%
  pivot_longer(cols=starts_with("OXY_oxy_"), names_to="Depth", names_prefix="OXY_oxy_", values_to = "OXY_oxy") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) 
oxy$DateTime <- as.Date(oxy$DateTime, format="%Y-%m-%d")

#Merge into one dataset
epi_oxy <- filter(oxy, Depth==1)
hypo_oxy <- filter(oxy, Depth==8)
met_oxy <- filter(oxy, Depth==4)
merge_mod<- merge(epi_oxy, hypo_oxy, by="DateTime") %>%
  dplyr::rename(epi_oxy = OXY_oxy.x, hypo_oxy = OXY_oxy.y)
merge_mod$exp_oxy <- (merge_mod$epi_oxy+merge_mod$hypo_oxy)/2

#Calculate deviation between modelled met oxy and expected met oxy
merge_mod <- merge(merge_mod, met_oxy[, c("DateTime", "OXY_oxy")], by="DateTime")%>%
  dplyr::rename(met_oxy = OXY_oxy)
merge_mod$deviation <- merge_mod$met_oxy - merge_mod$exp_oxy

obs_mom$DateTime <- as.Date(obs_mom$DateTime, format="%Y-%m-%d")
merge_mod$DateTime <- as.Date(merge_mod$DateTime, format="%Y-%m-%d")

#Calculating error (MEF) of metalimnetic oxygen minima for Figure 7
merge_mom <- merge(obs_mom, merge_mod, by="DateTime")

for (i in 1:nrow(merge_mom)) {
  merge_mom$MEFF_1[i]<- ((merge_mom$deviation.y[i]- merge_mom$deviation.x[i])^2)
  merge_mom$MEFF_2[i]<- ((merge_mom$deviation.x[i]-mean(merge_mom$deviation.x))^2)
  MEFF_MOM<- 1-(sum(merge_mom$MEFF_1)/sum(merge_mom$MEFF_2))
}

#Adding calculated Deep mixing 1 exm w3 model MOM MEF to error table
error[error$metric=="MOM" & error$calibration=="PEST_exm_w3", "Validation.deepm1"] <- MEFF_MOM
##############################################################################################################################

#Anoxia
#Deepm1 naive
output <- nc_open('Validation/valid_naive_deepm1/output/output.nc')

oxy<- ncvar_get(output, "OXY_oxy")
depth<- ncvar_get(output, "z")
depth[depth >= 100] <- NA
tallest_layer <- ncvar_get(output, "NS")

# Set every column of depth dataset ascending- result: 'out'
out <- apply(depth, 2, sort, decreasing=TRUE, na.last=TRUE)

#Loop for depth reference surface- result: 'new'
new <- out
for (i in 1:nrow(new)) {
  for (j in 1:ncol(new)) {
    new[i-1, j] <- out[1, j] - out[i, j]
  }
}

#Maximum depth
for (i in 1:dim(tallest_layer)) {
  new[tallest_layer[i], i] <- out[1, i]
}

#Minimum depth (surface)
new <- rbind(seq(0, 0, length.out = ncol(new)), new)

# Reverse every column of oxy dataset- result: 'oxy_out'
oxy_out <- apply(oxy, 2, rev)

# Putting NAs at the end of each column 
na_fun <- function(x) {
  num.na <- sum(is.na(x))
  x <- x[!is.na(x)]
  x <- c(x, rep(NA, num.na))
  return(x)
}

# apply na_fun over each column in the oxy_out dataframe
oxy_out <- apply(oxy_out, 2, na_fun)
oxy_out <-rbind(oxy_out[1,], oxy_out)

#Melting data into one column, creating dataframe depth, oxy
co<- melt(oxy_out, na.rm=TRUE)
co <- na.omit(co)
ct <- melt(new, na.rm=TRUE)
ct<- na.omit(ct)
df <- cbind(co$Var2, co$value, ct$value)

#Creating dataframe for time 
time <- data.frame(seq(as.Date("2015-07-20"), as.Date("2016-11-30"), by="day"))
ID <- seq.int(1:500)
time <- cbind(ID, time)
colnames(time) <- c("ID", "DateTime")
colnames(df) <- c("ID", "Oxy", "Depth")

#Merge time, depth, oxy according to ID
merge <- merge(time, df, all=TRUE)
merge$DateTime <- as.Date(merge$DateTime, format="%Y-%m-%d")

#Oxygen estimator function at any depth on certain date
estimate_oxy_by_date <- function(target_date, target_depth) {
  data_for_date <- merge %>% 
    filter(DateTime == target_date) %>%
    arrange(Depth)
  
  approx(data_for_date$Depth, data_for_date$Oxy, xout = target_depth)$y
}

#Interpolate DO in 0.02m increments in the water column
oxy_interp_depth <- crossing(
  tibble(DateTime = unique(merge$DateTime)),
  tibble(Depth = seq(0, 9.2, by = 0.1))
) %>%
  group_by(DateTime) %>%
  mutate(Oxy = estimate_oxy_by_date(DateTime[1], Depth))

#Model anoxia 
anoxia <- oxy_interp_depth
anoxia$Oxy <- anoxia$Oxy * 32/1000
anoxia$Oxy<- ifelse(anoxia$Oxy<=1, 1, 0)
anoxia <- na.omit(anoxia)
anoxia$DateTime <- as.Date(anoxia$DateTime, format="%Y-%m-%d")

#Calculating error (MEF) of number of anoxic layers per day for Figure 7
#Observations
obs_anoxic_layers <- read.csv('Observations/anoxia_observed.csv') %>%
  mutate(DateTime=as.Date(DateTime, format="%Y-%m-%d"))

#Modelled
anoxia <- mutate(anoxia, month = lubridate::month(DateTime)) %>%
  filter(between(month, 5, 11))

#Creating empty dataframe for loop
uniqueDates  <- unique(anoxia$DateTime)

newData_mod  <- data.frame(
  DateTime = unique(anoxia$DateTime),
  Count = length(uniqueDates)
)

#Calculating number of anoxic layers each day
for(i in 1:length(uniqueDates)){
  
  filteredData  <- filter(anoxia, DateTime==uniqueDates[i] & Oxy==1)
  newData_mod$Count[i] <- nrow(filteredData)
}

newData_mod$DateTime <- as.Date(newData_mod$DateTime, format="%Y-%m-%d")
merge_anoxia <- merge(obs_anoxic_layers, newData_mod, by="DateTime")

for (i in 1:nrow(merge_anoxia)) {
  merge_anoxia$MEFF_1[i]<- ((merge_anoxia$Count.y[i]- merge_anoxia$Count.x[i])^2)
  merge_anoxia$MEFF_2[i]<- ((merge_anoxia$Count.x[i]-mean(merge_anoxia$Count.x))^2)
  MEFF_anoxia<- 1-(sum(merge_anoxia$MEFF_1)/sum(merge_anoxia$MEFF_2))
}

#Adding calculated MEF to error table
#Deepm1 naive
error[error$metric=="A" & error$calibration=="PEST_N", "Validation.deepm1"] <- MEFF_anoxia
##########################################################################################################################

#Deepm1 exm w1
output <- nc_open('Validation/valid_exm_deepm1_w1/output/output.nc')

oxy<- ncvar_get(output, "OXY_oxy")
depth<- ncvar_get(output, "z")
depth[depth >= 100] <- NA
tallest_layer <- ncvar_get(output, "NS")

# Set every column of depth dataset ascending- result: 'out'
out <- apply(depth, 2, sort, decreasing=TRUE, na.last=TRUE)

#Loop for depth reference surface- result: 'new'
new <- out
for (i in 1:nrow(new)) {
  for (j in 1:ncol(new)) {
    new[i-1, j] <- out[1, j] - out[i, j]
  }
}

#Maximum depth
for (i in 1:dim(tallest_layer)) {
  new[tallest_layer[i], i] <- out[1, i]
}

#Minimum depth (surface)
new <- rbind(seq(0, 0, length.out = ncol(new)), new)

# Reverse every column of oxy dataset- result: 'oxy_out'
oxy_out <- apply(oxy, 2, rev)

# Putting NAs at the end of each column 
na_fun <- function(x) {
  num.na <- sum(is.na(x))
  x <- x[!is.na(x)]
  x <- c(x, rep(NA, num.na))
  return(x)
}

# apply na_fun over each column in the oxy_out dataframe
oxy_out <- apply(oxy_out, 2, na_fun)
oxy_out <-rbind(oxy_out[1,], oxy_out)

#Melting data into one column, creating dataframe depth, oxy
co<- melt(oxy_out, na.rm=TRUE)
co <- na.omit(co)
ct <- melt(new, na.rm=TRUE)
ct<- na.omit(ct)
df <- cbind(co$Var2, co$value, ct$value)

#Creating dataframe for time 
time <- data.frame(seq(as.Date("2015-07-20"), as.Date("2016-11-30"), by="day"))
ID <- seq.int(1:500)
time <- cbind(ID, time)
colnames(time) <- c("ID", "DateTime")
colnames(df) <- c("ID", "Oxy", "Depth")

#Merge time, depth, oxy according to ID
merge <- merge(time, df, all=TRUE)
merge$DateTime <- as.Date(merge$DateTime, format="%Y-%m-%d")

#Oxygen estimator function at any depth on certain date
estimate_oxy_by_date <- function(target_date, target_depth) {
  data_for_date <- merge %>% 
    filter(DateTime == target_date) %>%
    arrange(Depth)
  
  approx(data_for_date$Depth, data_for_date$Oxy, xout = target_depth)$y
}

#Interpolate DO in 0.02m increments in the water column
oxy_interp_depth <- crossing(
  tibble(DateTime = unique(merge$DateTime)),
  tibble(Depth = seq(0, 9.2, by = 0.1))
) %>%
  group_by(DateTime) %>%
  mutate(Oxy = estimate_oxy_by_date(DateTime[1], Depth))

#Model anoxia 
anoxia <- oxy_interp_depth
anoxia$Oxy <- anoxia$Oxy * 32/1000
anoxia$Oxy<- ifelse(anoxia$Oxy<=1, 1, 0)
anoxia <- na.omit(anoxia)
anoxia$DateTime <- as.Date(anoxia$DateTime, format="%Y-%m-%d")

#Calculating error (MEF) of number of anoxic layers per day for Figure 7
#Observations
obs_anoxic_layers <- read.csv('Observations/anoxia_observed.csv') %>%
  mutate(DateTime=as.Date(DateTime, format="%Y-%m-%d"))

#Modelled
anoxia <- mutate(anoxia, month = lubridate::month(DateTime)) %>%
  filter(between(month, 5, 11))

#Creating empty dataframe for loop
uniqueDates  <- unique(anoxia$DateTime)

newData_mod  <- data.frame(
  DateTime = unique(anoxia$DateTime),
  Count = length(uniqueDates)
)

#Calculating number of anoxic layers each day
for(i in 1:length(uniqueDates)){
  
  filteredData  <- filter(anoxia, DateTime==uniqueDates[i] & Oxy==1)
  newData_mod$Count[i] <- nrow(filteredData)
}

newData_mod$DateTime <- as.Date(newData_mod$DateTime, format="%Y-%m-%d")
merge_anoxia <- merge(obs_anoxic_layers, newData_mod, by="DateTime")

for (i in 1:nrow(merge_anoxia)) {
  merge_anoxia$MEFF_1[i]<- ((merge_anoxia$Count.y[i]- merge_anoxia$Count.x[i])^2)
  merge_anoxia$MEFF_2[i]<- ((merge_anoxia$Count.x[i]-mean(merge_anoxia$Count.x))^2)
  MEFF_anoxia<- 1-(sum(merge_anoxia$MEFF_1)/sum(merge_anoxia$MEFF_2))
}

#Deepm1 exm w1
error[error$metric=="A" & error$calibration=="PEST_exm_w1", "Validation.deepm1"] <- MEFF_anoxia
###############################################################################################################################

#Deepm1 exm w2
output <- nc_open('Validation/valid_exm_deepm1_w2/output/output.nc')

oxy<- ncvar_get(output, "OXY_oxy")
depth<- ncvar_get(output, "z")
depth[depth >= 100] <- NA
tallest_layer <- ncvar_get(output, "NS")

# Set every column of depth dataset ascending- result: 'out'
out <- apply(depth, 2, sort, decreasing=TRUE, na.last=TRUE)

#Loop for depth reference surface- result: 'new'
new <- out
for (i in 1:nrow(new)) {
  for (j in 1:ncol(new)) {
    new[i-1, j] <- out[1, j] - out[i, j]
  }
}

#Maximum depth
for (i in 1:dim(tallest_layer)) {
  new[tallest_layer[i], i] <- out[1, i]
}

#Minimum depth (surface)
new <- rbind(seq(0, 0, length.out = ncol(new)), new)

# Reverse every column of oxy dataset- result: 'oxy_out'
oxy_out <- apply(oxy, 2, rev)

# Putting NAs at the end of each column 
na_fun <- function(x) {
  num.na <- sum(is.na(x))
  x <- x[!is.na(x)]
  x <- c(x, rep(NA, num.na))
  return(x)
}

# apply na_fun over each column in the oxy_out dataframe
oxy_out <- apply(oxy_out, 2, na_fun)
oxy_out <-rbind(oxy_out[1,], oxy_out)

#Melting data into one column, creating dataframe depth, oxy
co<- melt(oxy_out, na.rm=TRUE)
co <- na.omit(co)
ct <- melt(new, na.rm=TRUE)
ct<- na.omit(ct)
df <- cbind(co$Var2, co$value, ct$value)

#Creating dataframe for time 
time <- data.frame(seq(as.Date("2015-07-20"), as.Date("2016-11-30"), by="day"))
ID <- seq.int(1:500)
time <- cbind(ID, time)
colnames(time) <- c("ID", "DateTime")
colnames(df) <- c("ID", "Oxy", "Depth")

#Merge time, depth, oxy according to ID
merge <- merge(time, df, all=TRUE)
merge$DateTime <- as.Date(merge$DateTime, format="%Y-%m-%d")

#Oxygen estimator function at any depth on certain date
estimate_oxy_by_date <- function(target_date, target_depth) {
  data_for_date <- merge %>% 
    filter(DateTime == target_date) %>%
    arrange(Depth)
  
  approx(data_for_date$Depth, data_for_date$Oxy, xout = target_depth)$y
}

#Interpolate DO in 0.02m increments in the water column
oxy_interp_depth <- crossing(
  tibble(DateTime = unique(merge$DateTime)),
  tibble(Depth = seq(0, 9.2, by = 0.1))
) %>%
  group_by(DateTime) %>%
  mutate(Oxy = estimate_oxy_by_date(DateTime[1], Depth))

#Model anoxia 
anoxia <- oxy_interp_depth
anoxia$Oxy <- anoxia$Oxy * 32/1000
anoxia$Oxy<- ifelse(anoxia$Oxy<=1, 1, 0)
anoxia <- na.omit(anoxia)
anoxia$DateTime <- as.Date(anoxia$DateTime, format="%Y-%m-%d")

#Calculating error (MEF) of number of anoxic layers per day for Figure 7
#Observations
obs_anoxic_layers <- read.csv('Observations/anoxia_observed.csv') %>%
  mutate(DateTime=as.Date(DateTime, format="%Y-%m-%d"))

#Modelled
anoxia <- mutate(anoxia, month = lubridate::month(DateTime)) %>%
  filter(between(month, 5, 11))

#Creating empty dataframe for loop
uniqueDates  <- unique(anoxia$DateTime)

newData_mod  <- data.frame(
  DateTime = unique(anoxia$DateTime),
  Count = length(uniqueDates)
)

#Calculating number of anoxic layers each day
for(i in 1:length(uniqueDates)){
  
  filteredData  <- filter(anoxia, DateTime==uniqueDates[i] & Oxy==1)
  newData_mod$Count[i] <- nrow(filteredData)
}

newData_mod$DateTime <- as.Date(newData_mod$DateTime, format="%Y-%m-%d")
merge_anoxia <- merge(obs_anoxic_layers, newData_mod, by="DateTime")

for (i in 1:nrow(merge_anoxia)) {
  merge_anoxia$MEFF_1[i]<- ((merge_anoxia$Count.y[i]- merge_anoxia$Count.x[i])^2)
  merge_anoxia$MEFF_2[i]<- ((merge_anoxia$Count.x[i]-mean(merge_anoxia$Count.x))^2)
  MEFF_anoxia<- 1-(sum(merge_anoxia$MEFF_1)/sum(merge_anoxia$MEFF_2))
}

#Deepm1 exm w2
error[error$metric=="A" & error$calibration=="PEST_exm_w2", "Validation.deepm1"] <- MEFF_anoxia
####################################################################################################

#Deepm1 exm w3
output <- nc_open('Validation/valid_exm_deepm1_w3/output/output.nc')

oxy<- ncvar_get(output, "OXY_oxy")
depth<- ncvar_get(output, "z")
depth[depth >= 100] <- NA
tallest_layer <- ncvar_get(output, "NS")

# Set every column of depth dataset ascending- result: 'out'
out <- apply(depth, 2, sort, decreasing=TRUE, na.last=TRUE)

#Loop for depth reference surface- result: 'new'
new <- out
for (i in 1:nrow(new)) {
  for (j in 1:ncol(new)) {
    new[i-1, j] <- out[1, j] - out[i, j]
  }
}

#Maximum depth
for (i in 1:dim(tallest_layer)) {
  new[tallest_layer[i], i] <- out[1, i]
}

#Minimum depth (surface)
new <- rbind(seq(0, 0, length.out = ncol(new)), new)

# Reverse every column of oxy dataset- result: 'oxy_out'
oxy_out <- apply(oxy, 2, rev)

# Putting NAs at the end of each column 
na_fun <- function(x) {
  num.na <- sum(is.na(x))
  x <- x[!is.na(x)]
  x <- c(x, rep(NA, num.na))
  return(x)
}

# apply na_fun over each column in the oxy_out dataframe
oxy_out <- apply(oxy_out, 2, na_fun)
oxy_out <-rbind(oxy_out[1,], oxy_out)

#Melting data into one column, creating dataframe depth, oxy
co<- melt(oxy_out, na.rm=TRUE)
co <- na.omit(co)
ct <- melt(new, na.rm=TRUE)
ct<- na.omit(ct)
df <- cbind(co$Var2, co$value, ct$value)

#Creating dataframe for time 
time <- data.frame(seq(as.Date("2015-07-20"), as.Date("2016-11-30"), by="day"))
ID <- seq.int(1:500)
time <- cbind(ID, time)
colnames(time) <- c("ID", "DateTime")
colnames(df) <- c("ID", "Oxy", "Depth")

#Merge time, depth, oxy according to ID
merge <- merge(time, df, all=TRUE)
merge$DateTime <- as.Date(merge$DateTime, format="%Y-%m-%d")

#Oxygen estimator function at any depth on certain date
estimate_oxy_by_date <- function(target_date, target_depth) {
  data_for_date <- merge %>% 
    filter(DateTime == target_date) %>%
    arrange(Depth)
  
  approx(data_for_date$Depth, data_for_date$Oxy, xout = target_depth)$y
}

#Interpolate DO in 0.02m increments in the water column
oxy_interp_depth <- crossing(
  tibble(DateTime = unique(merge$DateTime)),
  tibble(Depth = seq(0, 9.2, by = 0.1))
) %>%
  group_by(DateTime) %>%
  mutate(Oxy = estimate_oxy_by_date(DateTime[1], Depth))

#Model anoxia 
anoxia <- oxy_interp_depth
anoxia$Oxy <- anoxia$Oxy * 32/1000
anoxia$Oxy<- ifelse(anoxia$Oxy<=1, 1, 0)
anoxia <- na.omit(anoxia)
anoxia$DateTime <- as.Date(anoxia$DateTime, format="%Y-%m-%d")

#Calculating error (MEF) of number of anoxic layers per day for Figure 7
#Observations
obs_anoxic_layers <- read.csv('Observations/anoxia_observed.csv') %>%
  mutate(DateTime=as.Date(DateTime, format="%Y-%m-%d"))

#Modelled
anoxia <- mutate(anoxia, month = lubridate::month(DateTime)) %>%
  filter(between(month, 5, 11))

#Creating empty dataframe for loop
uniqueDates  <- unique(anoxia$DateTime)

newData_mod  <- data.frame(
  DateTime = unique(anoxia$DateTime),
  Count = length(uniqueDates)
)

#Calculating number of anoxic layers each day
for(i in 1:length(uniqueDates)){
  
  filteredData  <- filter(anoxia, DateTime==uniqueDates[i] & Oxy==1)
  newData_mod$Count[i] <- nrow(filteredData)
}

newData_mod$DateTime <- as.Date(newData_mod$DateTime, format="%Y-%m-%d")
merge_anoxia <- merge(obs_anoxic_layers, newData_mod, by="DateTime")

for (i in 1:nrow(merge_anoxia)) {
  merge_anoxia$MEFF_1[i]<- ((merge_anoxia$Count.y[i]- merge_anoxia$Count.x[i])^2)
  merge_anoxia$MEFF_2[i]<- ((merge_anoxia$Count.x[i]-mean(merge_anoxia$Count.x))^2)
  MEFF_anoxia<- 1-(sum(merge_anoxia$MEFF_1)/sum(merge_anoxia$MEFF_2))
}

#Deepm1 exm w3
error[error$metric=="A" & error$calibration=="PEST_exm_w3", "Validation.deepm1"] <- MEFF_anoxia

write.csv(error, 'Observations/error_stats.csv', row.names=FALSE)



