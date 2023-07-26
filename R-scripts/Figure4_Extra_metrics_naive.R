#Packages
library(tidyr)
library(dplyr)
library(glmtools)
library(GLMr)
library(rLakeAnalyzer)
library(lubridate)
library(reshape2)
library(ggplot2)
library(ncdf4)
library(reshape2)
library(stats)
library(utils)
library(ggplot2)
library(ggpubr)

#Set working directory
setwd(".../FCR-GLM-metrics")
sim_folder <- getwd()
PEST_calib <- file.path(sim_folder, 'Calibrated_models/Deepm2_naive/output/output.nc')
error <- read.csv("observations/error_stats.csv")

#Observed thermocline deph
obs_TD<- read.csv('observations/obs_td.csv') %>%
  mutate(month = lubridate::month(DateTime)) %>%
  filter(between(month, 4, 9)) %>%
  mutate(year = lubridate::year(DateTime)) %>%
  group_by(year) %>%
  filter(DateTime  > "2016-12-01" & DateTime < "2020-01-01")

obs_TD$DateTime <- as.Date(obs_TD$DateTime, format="%Y-%m-%d")

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

#Adding calculated MEF to error table
#error[error$metric=="TD" & error$calibration=="PEST_N", "Calibration.deepm2"] <- MEFF_TD

#Ice on and off

#Observations
obs_ice<- read.csv("observations/Ice_Data_2013_2022.csv") %>%
  dplyr::rename(DateTime = Date) %>%
  #filter(DateTime > "2015-07-12" & DateTime < "2016-12-02")
  filter(DateTime > "2016-12-01" & DateTime < "2020-01-01")
obs_ice$DateTime <- as.Date(obs_ice$DateTime, format="%Y-%m-%d")

#Ice model
ice<-glmtools::get_var(PEST_calib, var_name="white_ice_thickness")
iceblue<-glmtools::get_var(PEST_calib, var_name="blue_ice_thickness")
icesnow <- glmtools::get_var(PEST_calib, var_name="snow_thickness")
mod_ice_PEST<- cbind(ice, iceblue$blue_ice_thickness, icesnow$snow_thickness) %>%
  dplyr::rename(blue_ice_thickness = `iceblue$blue_ice_thickness`, snow_thickness = `icesnow$snow_thickness`) %>%
  mutate(IceOn_mod = ifelse(white_ice_thickness>0 | blue_ice_thickness >0 | snow_thickness >0, 1, 0)) %>%
  mutate(DateTime = as.Date(DateTime, format="%Y-%m-%d")) %>%
  mutate(dayofyear = lubridate::yday(DateTime)) %>%
  mutate(month = lubridate::month(DateTime)) %>%
  filter(month==12 | month==1 | month==2) %>%
  mutate(year = lubridate::year(DateTime)) %>%
  group_by(year) 

mod_ice_PEST$DateTime <- as.Date(mod_ice_PEST$DateTime, format="%Y-%m-%d")

#Merge modelled and observed ice datasets
merge_ice<- merge(obs_ice, mod_ice_PEST, all.y=TRUE) %>%
  mutate(dayofyear = lubridate::yday(DateTime)) %>%
  replace_na(list(IceOn = 0)) %>%
  mutate(month = lubridate::month(DateTime)) %>%
  filter(month==12 | month==1 |month==2)

merge_ice_new<- merge_ice

#Mapping 0 to 4 and 1 to 1 for plotting reasons
merge_ice_new$IceOn[merge_ice_new$IceOn == 0] <- 4
merge_ice_new$IceOn[merge_ice_new$IceOn == 1] <- 1

merge_ice_new$IceOn_mod[merge_ice_new$IceOn_mod == 0] <- 4
merge_ice_new$IceOn_mod[merge_ice_new$IceOn_mod == 1] <- 1


#Combined time series of thermocline depth and ice cover
plot_2 <- thermo_depth_model %>%
  ggplot2::ggplot(ggplot2::aes(x = DateTime, y = td_model, colour="PEST_N TD", group=year)) +
  ggplot2::geom_line()+
  geom_point(data=merge_ice_new, aes(x=DateTime, y=IceOn_mod, colour="PEST_N ice"), pch=0, size=1.8, stroke=0.5)+
  geom_point(data=merge_ice_new, aes(x=DateTime, y=IceOn, colour="Obs. ice"), pch=20, size=1)+
  ggplot2::geom_point(data=obs_TD[-1,], aes(x=DateTime, y=thermo.depth, colour="Obs. TD"), pch=10)+
  ggplot2::labs(x = "Date", y = "Depth (m)")+
  geom_rect(ymin = 0, ymax = -6, 
            xmin = as.Date("2016-11-30"), xmax = as.Date("2017-03-02"), colour='blue', fill='NA') +
  geom_rect(ymin = 0, ymax = -6, 
            xmin = as.Date("2017-11-30"), xmax = as.Date("2018-03-02"), colour='blue', fill='NA') +
  geom_rect(ymin = 0, ymax = -6, 
            xmin = as.Date("2018-11-30"), xmax = as.Date("2019-03-02"), colour='blue', fill='NA') +
  geom_rect(ymin = 0, ymax = -6, 
            xmin = as.Date("2019-11-30"), xmax = as.Date("2020-01-02"), colour='blue', fill='NA') +
  ggplot2::scale_y_reverse(limits = c(6, -0.5), breaks= c(0, 1, 2, 3, 4, 5, 6), labels= c("0", "1", "2", "3", "4", "5", "6"), sec.axis = sec_axis(~./1, name = "Ice cover", breaks=c(0, 1, 2, 3, 4, 5, 6), labels=c(" ", "On ", " ", "", "Off", " ", " ")))+
  ggplot2::scale_x_date(expand=c(0.01,0.01), date_breaks = "6 months")+
  ggplot2::scale_colour_manual(name="", values=c("Obs. TD"="black", "PEST_N TD"="#FF61CC", "PEST_N ice"="skyblue", "Obs. ice"="black"), guide=guide_legend(override.aes=list(linetype=c(NA, 1, NA, NA), size=c(1.5, 0.8, 2, 2), shape=c(10, NA, 0, 20))))+
  ggplot2::theme_light() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face= "bold", size = 12),
    axis.title.y = ggplot2::element_text(face="bold", size= 10),
    axis.title.x = ggplot2::element_text(face="bold", size= 10),
    axis.text.x = ggplot2::element_text(size=9),
    axis.text.y = ggplot2::element_text(size=9),
    legend.text = ggplot2::element_text(size= 7),
    legend.title = ggplot2::element_blank(),
    legend.position=c(0.5, 0.96),
    legend.direction="horizontal",
    legend.key.height = unit(2, "mm")
   )

plot_2


#Schmidt stability
bathy <- read.csv('observations/bathymetry.csv')

#Schmidt stability observed
schmidt_stability_obs<- read.csv("observations/Obs_SS.csv") %>%
  filter(datetime > "2016-12-01" & datetime < "2020-01-01")
#filter(datetime > "2015-07-12" & datetime < "2016-12-01")
schmidt_stability_obs$datetime<- as.Date(schmidt_stability_obs$datetime, format='%Y-%m-%d')

#Schmidt stability dataset (using same temp dataset as for thermoclinde depth)
#temp<- get_var(PEST_calib, var_name="temp", reference="surface")
#colClean <- function(x){ colnames(x) <- gsub("temp", "wtr", colnames(x)); x } 
#temp <- colClean(temp)
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


#Adding calculated MEF to error table
#error[error$metric=="SS" & error$calibration=="PEST_N", "Calibration.deepm2"] <- MEFF_SS

SS_plot <- ggplot(data=schmidt_stability_obs, aes(x=datetime, y=schmidt.stability, colour="Obs. SS")) +
  geom_point(pch=10)+
  ylab("Schmidt stability")+
  xlab("Date")+
  ylim(c(0, 70))+
  geom_line(data=schmidt_stability, aes(x=datetime, y=ss_PEST, colour="PEST_N SS"))+
  scale_x_date(expand=c(0.01, 20))+
  scale_colour_manual(values=c("Obs. SS"="black", "PEST_N SS" ="#00BFC4"), guide=guide_legend(override.aes=list(linetype=c(NA, 1), shape=c(10, NA))))+
  ggplot2::theme_light() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face= "bold", size = 12),
    axis.title.y = ggplot2::element_text(face="bold", size= 10),
    axis.title.x = ggplot2::element_text(face="bold", size= 10),
    legend.text = ggplot2::element_text(size= 7),
    axis.text.x = ggplot2::element_text(size=9),
    axis.text.y = ggplot2::element_text(size=9),
    legend.title = ggplot2::element_blank(),
    legend.position = c(0.5, 0.92),
    legend.direction="horizontal",
    legend.key.height = unit(2, "mm")
  )

SS_plot

#MOM
#Observed
obs_mom<-read.csv('observations/mom_observed.csv') %>%
  #filter(DateTime > "2015-07-12" & DateTime < "2017-01-01")
  filter(DateTime > "2016-12-02" & DateTime < "2020-01-01")
obs_mom$DateTime <- as.Date(obs_mom$DateTime, format = "%Y-%m-%d")

#Modelled
#Get modelled oxy at 1, 4, 8 m depths
depths<- c(1, 4, 8) 
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


#Adding calculated MEF to error table
#error[error$metric=="MOM" & error$calibration=="PEST_N", "Calibration.deepm2"] <- MEFF_MOM

#plot MOM
plot_MOM <-ggplot(data=merge_mod, aes(x=DateTime, y=deviation, colour="PEST_N MOM"))+
  geom_line()+
  geom_point(data=obs_mom, aes(x=DateTime, y=deviation, colour="Obs. MOM"), pch=10)+
  geom_hline(yintercept=0, lty=3)+
  xlab("Date")+
  ylab(expression(bold(MOM~(mmol/m^{3}))))+
  scale_x_date(expand=c(0.01,0.01))+
  ggplot2::scale_colour_manual(name="Legend", values=c("Obs. MOM"="black", "PEST_N MOM" ="#CD9600"), guide=guide_legend(override.aes=list(linetype=c(NA, 1), shape=c(10, NA))))+
  ggplot2::theme_light() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face= "bold", size = 12),
    axis.title.y = ggplot2::element_text(face="bold", size= 10),
    axis.title.x = ggplot2::element_text(face="bold", size= 10),
    legend.text = ggplot2::element_text(size= 7), 
    axis.text.x = ggplot2::element_text(size=9),
    axis.text.y = ggplot2::element_text(size=9),
    legend.title = ggplot2::element_blank(),
    legend.position = c(0.35, 0.92),
    legend.direction="horizontal",
    plot.margin = unit(c(5.5,10,5.5,5.5), "pt"),
    legend.key.height = unit(2, "mm")
    
  )
plot_MOM

#Anoxia
output <- nc_open('Calibrated_models/Deepm2_naive/output/output.nc')
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
time <- data.frame(seq(as.Date("2016-12-02"), as.Date("2019-12-31"), by="day"))
ID <- seq.int(1:1125)
#time <- data.frame(seq(as.Date("2015-07-13"), as.Date("2016-12-01"), by="day"))
#ID <- seq.int(1:508)
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

#Observed
obs_oxy<-read.csv('observations/CleanedObsOxy.csv') %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) %>%
  #filter(DateTime > "2015-07-12" & DateTime < "2016-12-31")
  filter(DateTime > "2016-12-01" & DateTime < "2020-01-01")

obs_anoxia <- obs_oxy %>%
  mutate(OXY_oxy=OXY_oxy*32/1000) %>%
  mutate(OXY_oxy=ifelse(OXY_oxy<=1, 1, 0)) %>%
  na.omit() 
obs_anoxia$DateTime <- as.Date(obs_anoxia$DateTime, format="%Y-%m-%d")

#Anoxia plot 1 model
PEST_anoxia<- ggplot() +
  geom_raster(anoxia, mapping =aes(DateTime, Depth, fill = as.character(Oxy))) +
  ylab("Depth (m)")+
  xlab("Date")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_y_reverse(expand=c(0,0)) +
  scale_x_date(expand=c(0,0))+
  geom_point(filter(obs_anoxia, obs_anoxia$OXY_oxy==1), mapping = aes(x = DateTime, y = Depth, colour = "Observed anoxia"), pch=4)+
  scale_fill_manual(
    labels = c("Anoxic water (DO <=1mg/L)", "'Oxic' water (DO>1mg/L)"),
    values = c("1"="#F8766D","0"="blue")
  )+
  scale_colour_manual(
    values = c("Observed anoxia"="black")
  ) +
  guides(fill = guide_legend(order = 1)) +
  theme(
    legend.background = element_rect(fill="white"),
    legend.spacing = unit(-6, "pt"),
    plot.title = ggplot2::element_text(face= "bold", size = 12),
    axis.title.y = ggplot2::element_text(face="bold",size= 10),
    axis.title.x = ggplot2::element_text(face="bold", size= 10),
    legend.text = ggplot2::element_text(size= 7), 
    axis.text.x = ggplot2::element_text(size=9),
    axis.text.y = ggplot2::element_text(size=9),
    legend.title = ggplot2::element_blank(),
    legend.position = c(0.5, 0.92),
    legend.direction="horizontal",
    legend.box = "horizontal",
    legend.key.size = unit(3, 'mm'),
    legend.key.height = unit(2, "mm")#change legend key size
  )

PEST_anoxia

#Calculating error (MEF) of number of anoxic layers per day for Figure 7
#Observations
obs_anoxic_layers <- read.csv('observations/anoxia_observed.csv') %>%
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

merge_anoxia <- merge(obs_anoxic_layers, newData_mod, by="DateTime")

for (i in 1:nrow(merge_anoxia)) {
  merge_anoxia$MEFF_1[i]<- ((merge_anoxia$Count.y[i]- merge_anoxia$Count.x[i])^2)
  merge_anoxia$MEFF_2[i]<- ((merge_anoxia$Count.x[i]-mean(merge_anoxia$Count.x))^2)
  MEFF_anoxia<- 1-(sum(merge_anoxia$MEFF_1)/sum(merge_anoxia$MEFF_2))
}

#Adding calculated MEF to error table
#error[error$metric=="A" & error$calibration=="PEST_N", "Calibration.deepm2"] <- MEFF_anoxia
#write.csv(error, 'observations/error_stats.csv')

#Sediment temperature obs
obs_temp<-read.csv('observations/CleanedObsTemp.csv') %>%
  filter(DateTime >= "2016-12-01" & DateTime <= "2020-01-01")
obs_temp$DateTime <- as.Date(obs_temp$DateTime, format="%Y-%m-%d")

#Sediment temperature modelled
#Zone 2

#water temp
temp_z2 <- get_var(PEST_calib, var_name="temp", reference="surface", z_out=5) %>%
  mutate(doy = lubridate::yday(DateTime)) %>%
  rename(water_temp_z2 =temp_5)

#read nml file
nml_file <- read_nml(nml_file = file.path(sim_folder, 'Calibrated_models/Deepm2_naive/glm3.nml'))
sed_temp_mean <- get_nml_value(nml_file, 'sed_temp_mean')
sed_temp_amp <- get_nml_value(nml_file, 'sed_temp_amplitude')
sed_temp_doy <- get_nml_value(nml_file, 'sed_temp_peak_doy')

#get parameter values form nml file for zone 2
for(i in 1:nrow(temp_z2)){
  temp_z2$sed_temp_z2[i] <- sed_temp_mean[2] + sed_temp_amp[2]*cos((2*pi/365)*(temp_z2$doy[i]-sed_temp_doy[2]))
}
temp_z2$DateTime <- as.Date(temp_z2$DateTime, format="%Y-%m-%d")

zone2<- ggplot(data=temp_z2, aes(x=DateTime, y=water_temp_z2, colour="PEST_N water temp 5 m"))+
  geom_line()+
  geom_line(aes(x=DateTime, y=sed_temp_z2, colour="PEST_N sed temp z2"))+
  geom_point(data=filter(obs_temp, Depth==5), aes(x=DateTime, y=temp, colour="Obs. water temp 5 m"), pch=10)+
  labs(x = "Date", y = "Temperature (°C)")+
  scale_x_date(expand=c(0.01,0.01))+
  ylim(c(0, 30))+
  scale_colour_manual(name="Legend", values=c("Obs. water temp 5 m"="black", "PEST_N water temp 5 m" ="black", "PEST_N sed temp z2"="green"), guide=guide_legend(override.aes=list(linetype=c(NA, 1, 1), shape=c(10, NA, NA))))+
  theme_light() +
    theme(
      plot.title = ggplot2::element_text(face= "bold", size = 12),
      axis.title.y = ggplot2::element_text(face="bold", size= 10, margin = margin(t = 0, r = 10, b = 0, l = 0)),
      axis.title.x = ggplot2::element_text(face="bold", size= 10),
      legend.text = ggplot2::element_text(size= 7),
      axis.text.x = ggplot2::element_text(size=9),
      axis.text.y = ggplot2::element_text(size=9),
      legend.title = ggplot2::element_blank(),
      legend.position= c(0.49, 0.92),
      legend.direction="horizontal",
      plot.margin = unit(c(5.5,10,5.5,10), "pt"),
      legend.key.width = unit(5, "mm"),
      legend.margin = margin(2, 0, 2, 0)

  )
zone2
mix <- ggarrange(plot_2,                                                 
                 ggarrange(SS_plot, zone2, ncol=2, labels=c("b", "c")),
                 ggarrange(PEST_anoxia, plot_MOM, ncol = 2, labels = c("d", "e")), 
                 nrow = 3, 
                 labels = "a"
) 
mix
ggsave("Results/Figure4.png",
       plot = mix,
       dpi=500,
       width = 246.2, 
       height = 210, 
       units = "mm")