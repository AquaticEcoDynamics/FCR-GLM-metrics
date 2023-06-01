#Loading packages
library(dplyr)
library(remotes)
library(GLMr)
library(glmtools)
library(tidyr)
library(rLakeAnalyzer)
library(base)
library(utils)
library(ncdf4)
library(reshape2)

#Set working directory
setwd(".../FCR-GLM-metrics/Calibrated_models/PEST_EXM_weight1_deepm1")
sim_folder <- getwd()
nc_file <- file.path(sim_folder, 'output/output.nc')
depths<- c(0.1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9.2) 

obstemp<-read.csv('field_data/CleanedObsTemp.csv') %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST")))

obs_oxy<-read.csv('field_data/CleanedObsOxy.csv') %>%
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
  write.csv(match[[i]]$modtemp, paste0(".../FCR-GLM-metrics/Calibrated_models/PEST_EXM_weight1_deepm1/", #file path for saving csv files, should be set to the working directory
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
  write.csv(match_oxy[[i]]$mod_oxy, paste0(".../FCR-GLM-metrics/Calibrated_models/PEST_EXM_weight1_deepm1/", #file path for saving csv files, should be set to the working directory
                                                    "matchoxy", depths[i],
                                                    ".csv"), row.names=FALSE)
  
}

#Temperature 
temp<- get_var(nc_file, var_name="temp", reference="surface") 
colClean <- function(x){ colnames(x) <- gsub("temp", "wtr", colnames(x)); x } 
colClean(temp)
temp$DateTime <- as.Date(temp$DateTime, format="%Y-%m-$d")

obs_TD<- read.csv('field_data/obs_td.csv') %>%
  mutate(month = lubridate::month(DateTime)) %>%
  filter(between(month, 4, 9))

obs_TD$DateTime <- as.Date(obs_TD$DateTime, format="%Y-%m-%d")

#Modelled thermocline depths
thermo_depth_model <- ts.thermo.depth(temp, Smin = 0.1, na.rm=TRUE, seasonal=FALSE)  %>% 
  dplyr::rename(td_model = thermo.depth, DateTime = datetime) %>%
  mutate(month = lubridate::month(DateTime)) %>%
  filter(between(month, 4, 9)) 

thermo_depth_model$DateTime <- as.Date(thermo_depth_model$DateTime, format="%Y-%m-%d")

#Merge
td_merge <- merge(thermo_depth_model, obs_TD, by="DateTime", sort=FALSE)

#Need to figure out what to do with NaN
write.csv(td_merge$td_model, "match_td.csv", row.names=FALSE)

#Observed schmidt stability plus bathy data 
obs_SS<- read.csv('field_data/Obs_SS.csv') 
obs_SS$datetime <- as.Date(obs_SS$datetime, format="%Y-%m-%d") 

depths<- c(0,	0.4,	0.7,	1,	1.3,	1.6,	1.9,	2.3,	2.6,	2.9,	3.2,	3.5,	3.8,	4.1,	4.4,	4.7,	5,	5.3,	5.6,	5.9,	6.2,	6.5,	6.8,	7.1,	7.4,	7.7,	8,	8.3,	8.7,	9,	9.3)
areas<- c(119880.9164,	111302.1604,	103030.4489,	95068.47603,	85430.25429,	76424.14457,	68146.39437,	59666.85885,	51179.89109,	42851.13714,	36269.3312,	31155.95008,	28442.99667,	25790.70893,	22583.1399,	19631.05422,	16834.20941,	14484.22802,	12399.67051,	10811.30792,	9469.324081,	8228.697419,	6929.077352,	5637.911458,	4358.358439,	3239.620513,	2179.597283,	1201.23579,	494.615572,	61.408883,	0)
bathy_new<- cbind(depths, areas) 
bathy_new<- as.data.frame(bathy_new)

schmidt_stability_model <- rLakeAnalyzer::ts.schmidt.stability(temp, bathy_new, na.rm=TRUE)
schmidt_stability_model$datetime <- as.Date(schmidt_stability_model$datetime, format="%Y-%m-%d")

#Merge
ss_merge <- merge(obs_SS, schmidt_stability_model, by='datetime')
write.csv(ss_merge$schmidt.stability.y, "match_ss.csv", row.names=FALSE)

#MOM
mom_obs<- read.csv("field_data/mom_observed.csv") %>%
  mutate(DateTime = as.Date(DateTime, format = "%Y-%m-%d")) 

epi_oxy <- filter(mod_oxy, Depth==1)
hypo_oxy <- filter(mod_oxy, Depth==8)
met_oxy <- filter(mod_oxy, Depth==4)

mom_mod<- merge(epi_oxy, hypo_oxy, by="DateTime") %>%
  dplyr::rename(epi_oxy = OXY_oxy.x, hypo_oxy = OXY_oxy.y) %>%
  mutate(exp_oxy = (epi_oxy + hypo_oxy)/2) %>%
  merge(met_oxy, by="DateTime") %>%
  dplyr::rename(met_oxy = OXY_oxy) %>%
  mutate(deviation = met_oxy - exp_oxy, DateTime = as.Date(DateTime, format = "%Y-%m-%d"))
  
merge_mom <- merge(mom_mod, mom_obs, by="DateTime") %>%
  dplyr::rename(mom_mod = deviation.x, mom_obs = deviation.y)

write.csv(merge_mom$mom_mod, "match_mom.csv", row.names=FALSE)

#Anoxic factor

output <- nc_open("output/output.nc")
oxy<- ncvar_get(output, "OXY_oxy")
depth<- ncvar_get(output, "z") 
depth[depth >= 100] <- NA
tallest_layer <- ncvar_get(output, "NS")

# Set every column of depth dataset ascending- result: 'out'
depth <- apply(depth, 2, sort, decreasing=TRUE, na.last=TRUE)

#Loop for depth reference surface- result: 'new'
new <- depth
for (i in 1:nrow(new)-1) {
  for (j in 1:ncol(new)) {
    new[i, j] <- depth[1, j] - depth[i+1, j]
  }
}

#Maximum depth
for (i in 1:dim(tallest_layer)) {
  new[tallest_layer[i], i] <- depth[1, i]
}

#Minimum depth (surface)
new <- rbind(seq(0, 0, length.out = ncol(new)), new)

# Reverse every column of oxy dataset- result: 'oxy_out'
oxy <- apply(oxy, 2, rev)

# Putting NAs at the end of each column 
beetroot <- function(x) {
  # count NA
  num.na <- sum(is.na(x))
  # remove NA
  x <- x[!is.na(x)]
  # glue the number of NAs at the end
  x <- c(x, rep(NA, num.na))
  return(x)
}

# apply beetroot over each column in the oxy_out dataframe
oxy <- apply(oxy, 2, beetroot)
oxy <-rbind(oxy[1,], oxy)

#Melting data into one column, creating dataframe depth, oxy
oxy <- melt(oxy, na.rm=TRUE) %>%
  na.omit()

new <- melt(new, na.rm=TRUE) %>%
  na.omit()

combined <- cbind(oxy$Var2, oxy$value, new$value)
colnames(combined) <- c("ID", "Oxy", "Depth")

#Creating dataframe for time 
time <- data.frame(DateTime=seq(as.Date("2016-12-02"), as.Date("2019-12-31"), by="day")) %>%
  mutate(ID = seq.int(1:1125))

combined <- merge(time, combined, all=TRUE)
combined$DateTime <- as.Date(combined$DateTime, format="%Y-%m-%d")

#Oxygen estimator function at any depth on certain date
estimate_oxy_by_date <- function(target_date, target_depth) {
  data_for_date <- combined %>% 
    filter(DateTime == target_date) %>%
    arrange(Depth)
  
  approx(data_for_date$Depth, data_for_date$Oxy, xout = target_depth)$y
}

oxy_interp_depth <- crossing(
  tibble(DateTime = unique(combined$DateTime)),
  tibble(Depth = seq(0.1, 9.2, by = 0.1))
) %>%
  group_by(DateTime) %>%
  mutate(Oxy = estimate_oxy_by_date(DateTime[1], Depth))


#Model anoxia 
anoxia <- oxy_interp_depth %>%
  mutate(Oxy=Oxy*32/1000) %>%
  mutate(Oxy=ifelse(Oxy<=1, 1, 0)) %>%
  #na.omit() %>%
  mutate(month = lubridate::month(DateTime)) %>%
  filter(between(month, 5, 11))

anoxia$DateTime <- as.Date(anoxia$DateTime, format="%Y=%m-%d")
#anoxia<- na.omit(anoxia)

uniqueDates  <- unique(anoxia$DateTime)

newData  <- data.frame(
  DateTime = unique(anoxia$DateTime),
  Count = length(uniqueDates)
)

for(i in 1:length(uniqueDates)){
  
  filteredData  <- filter(anoxia, DateTime==uniqueDates[i] & Oxy==1)
  newData$Count[i] <- nrow(filteredData)
}

write.csv(newData$Count, "match_af.csv", row.names=FALSE)


