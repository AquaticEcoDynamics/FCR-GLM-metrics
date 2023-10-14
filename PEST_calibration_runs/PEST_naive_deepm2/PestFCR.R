library(dplyr)
library(remotes)
library(GLMr)
library(glmtools)
library(tidyr)
library(base)
library(utils)

#Set working directory
setwd(".../FCR-GLM-metrics/Calibrated_models/PEST_naive_deepm2")
sim_folder <- getwd()
nc_file <- file.path(sim_folder, 'output/output.nc')
depths<- c(0.1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9.2) 

obstemp<-read.csv('CleanedObsTemp.csv') %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST")))

obs_oxy<-read.csv('CleanedObsOxy.csv') %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST")))

modtemp <- get_temp(nc_file, reference="surface", z_out=depths) %>%
  pivot_longer(cols=starts_with("temp_"), names_to="Depth", names_prefix="temp_", values_to = "temp") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) 

watertemp<-merge(modtemp, obstemp, by=c("DateTime","Depth")) %>%
  dplyr::rename(modtemp = temp.x, obstemp = temp.y)

match <- vector('list')

for(i in 1:length(unique(watertemp$Depth))){
  
  
  match[[i]] <- watertemp %>%
    dplyr::filter(Depth == depths[i])
  write.csv(match[[i]]$modtemp, paste0(".../FCR-GLM-metrics/Calibrated_models/PEST_naive_deepm2/",  #file path for saving csv files, should be set to the working directory
                                       "match", depths[i],
                                       ".csv"), row.names=FALSE)
}



mod_oxy <- get_var(nc_file, "OXY_oxy", reference="surface", z_out=depths) %>%
  pivot_longer(cols=starts_with("OXY_oxy_"), names_to="Depth", names_prefix="OXY_oxy_", values_to = "OXY_oxy") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) 

oxy_compare <-merge(mod_oxy, obs_oxy, by=c("DateTime","Depth")) %>%
  dplyr::rename(mod_oxy = OXY_oxy.x, obs_oxy = OXY_oxy.y)

#Loop for writing csv files for PEST
match_oxy <- vector('list')

for(i in 1:length(unique(oxy_compare$Depth))){
  
  
  match_oxy[[i]] <- oxy_compare %>%
    filter(Depth == depths[i])
  write.csv(match_oxy[[i]]$mod_oxy, paste0(".../FCR-GLM-metrics/Calibrated_models/PEST_naive_deepm2/", #file path for saving csv files, should be set to the working directory
                                                    "matchoxy", depths[i],
                                                    ".csv"), row.names=FALSE)
  
}



