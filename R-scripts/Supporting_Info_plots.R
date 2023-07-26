library(ncdf4)
library(glmtools)
library(reshape2)
library(dplyr)
library(stats)
library(utils)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(GLMr)
library(rLakeAnalyzer)
library(lubridate)
library(patchwork)

#Set working directory
setwd(".../FCR-GLM-metrics")

#Figure B1

#Modelled oxy
sim_folder <- getwd()
output <- nc_open("Calibrated_models/Deepm1_naive/output/output.nc")
oxy<- ncvar_get(output, "OXY_oxy")
depth<- ncvar_get(output, "z")
depth[depth >= 100] <- NA
tallest_layer <- ncvar_get(output, "NS")
oxyunit <- expression(bold(~mmol/m^{3}))

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

#Add surface row
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

# apply function over each column in the oxy_out dataframe
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

#Interpolate DO in 0.1m increments in the water column
mod_oxy_int <- crossing(
  tibble(DateTime = unique(merge$DateTime)),
  tibble(Depth = seq(0.1, 9.2, by = 0.1))
) %>%
  group_by(DateTime) %>%
  mutate(Oxy = estimate_oxy_by_date(DateTime[1], Depth))

#Contour plot for modelled oxygen
modelled <-ggplot() + 
  geom_raster(mod_oxy_int, mapping =aes(DateTime, Depth, fill = Oxy))+
  scale_y_reverse(expand=c(0.01,0.01), limits=c(9.3, 0))+
  scale_x_date(expand=c(0.025,0.025), limits=c(as.Date("2016-12-01"), as.Date("2019-12-06")))+
  guides(fill = guide_colourbar(ticks.colour= "black"))+
  scale_fill_distiller(name = oxyunit, palette = 'Spectral', direction = -1, na.value = "grey90", guide="colourbar", limits = c(0, 500))+
  ylab("Depth (m)")+
  ggtitle("Dissolved Oxygen", subtitle = "Modelled DO")+
  theme(
    legend.background = element_rect(fill="white"),
    legend.spacing.y = unit(3, "mm"),
    plot.title = ggplot2::element_text(face= "bold", size = 12),
    plot.subtitle = ggplot2::element_text(size = 10),
    axis.title.y = ggplot2::element_text(face="bold",size= 10),
    axis.title.x = ggplot2::element_text(face="bold", size= 10),
    legend.text = ggplot2::element_text(size= 8), 
    axis.text.x = ggplot2::element_text(size=9),
    axis.text.y = ggplot2::element_text(size=9),
    legend.title = ggplot2::element_text(face="bold", size=10),
    legend.key.size = unit(7, 'mm')
  )
modelled

#Observations
obs_oxy <- read.csv('observations/CleanedObsOxy.csv') %>%
  mutate(DateTime = as.Date(DateTime, format="%Y-%m-%d"))
obs_oxy_int <- read.csv('observations/obs_oxy_interpolated.csv')%>%
  mutate(DateTime = as.Date(DateTime, format="%Y-%m-%d"))

#Plot observed
observed<- ggplot() + 
  geom_raster(obs_oxy_int, mapping =aes(DateTime, Depth, fill = OXY_oxy))+
  scale_y_reverse(expand=c(0.01,0.01), limits=c(9.3, 0))+
  scale_x_date(expand=c(0.025,0.025), limits=c(as.Date("2016-12-01"), as.Date("2019-12-06")))+
  guides(fill = guide_colourbar(ticks.colour= "black"))+
  ggtitle(" ", subtitle = "Observed DO")+
  scale_fill_distiller(name = oxyunit, palette = 'Spectral', direction = -1, na.value = "grey90", guide="colourbar", limits=c(0, 500))+
  geom_point(obs_oxy, mapping = aes(x = DateTime, y = Depth), colour="black", pch=4, size=0.5)+
  ylab("Depth (m)")+
  labs(fill = "Legend", colour = NULL)+
  theme(
    legend.background = element_rect(fill="white"),
    legend.spacing.y = unit(3, "mm"),
    plot.title = ggplot2::element_text(face= "bold", size = 12),
    plot.subtitle = ggplot2::element_text(size = 10),
    axis.title.y = ggplot2::element_text(face="bold",size= 10),
    axis.title.x = ggplot2::element_text(face="bold", size= 10),
    legend.text = ggplot2::element_text(size= 8),
    axis.text.x = ggplot2::element_text(size=9),
    axis.text.y = ggplot2::element_text(size=9),
    legend.title = ggplot2::element_text(face="bold", size=10),
    legend.key.size = unit(7, 'mm')
  )
observed

#Difference
merge_int <- merge(obs_oxy_int, mod_oxy_int, all=TRUE)
#modelled - observed
merge_int$difference <- merge_int$Oxy-merge_int$OXY_oxy
merge_int$DateTime <- as.Date(merge_int$DateTime, format = "%Y-%m-%d")

difference <- ggplot() + 
  geom_raster(merge_int, mapping =aes(DateTime, Depth, fill = difference))+
  scale_y_reverse(expand=c(0.01,0.01), limits=c(9.3, 0))+
  scale_x_date(expand=c(0.025,0.025), limits=c(as.Date("2016-12-01"), as.Date("2019-12-06")))+
  ggtitle(" ", subtitle = "Difference DO (Modelled - Observed)")+
  #guides(fill=guide_legend(reverse=TRUE, byrow=TRUE))+
  guides(fill = guide_colourbar(ticks.colour= "black"))+
  ylab("Depth (m)")+
  scale_fill_distiller(name = oxyunit, palette = 'Spectral', direction = -1, na.value = "grey90", guide='colourbar')+
  theme(
    legend.background = element_rect(fill="white"),
    legend.spacing.y = unit(3, "mm"),
    plot.title = ggplot2::element_text(face= "bold", size = 12),
    plot.subtitle = ggplot2::element_text(size = 10),
    axis.title.y = ggplot2::element_text(face="bold",size= 10),
    axis.title.x = ggplot2::element_text(face="bold", size= 10),
    legend.text = ggplot2::element_text(size= 8), 
    legend.title = ggplot2::element_text(face="bold", size=10),
    axis.text.x = ggplot2::element_text(size=9),
    axis.text.y = ggplot2::element_text(size=9),
    #legend.position = c(0.98, 0.85),
    #legend.direction="vertical",
    #legend.box = "horizontal",
    legend.key.size = unit(7, 'mm'),
    #legend.key.height = unit(2, "mm")#change legend key size
  )
difference

#CONTOUR PLOTS TEMPERATURE

#Modelled
sim_folder <- getwd()
output <- nc_open("Calibrated_models/Deepm1_naive/output/output.nc")
temp<- ncvar_get(output, "temp")
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

# Reverse every column of temp dataset- result: 'temp_out'
temp_out <- apply(temp, 2, rev)

# apply function over each column in the temp_out dataframe
temp_out <- apply(temp_out, 2, na_fun)
temp_out <-rbind(temp_out[1,], temp_out)

#Melting data into one column, creating dataframe depth, temp
melt_t<- melt(temp_out, na.rm=TRUE)
melt_t <- na.omit(melt_t)
melt_d <- melt(new, na.rm=TRUE)
melt_d<- na.omit(melt_d)
df <- cbind(melt_t$Var2, melt_t$value, melt_d$value)

#Creating dataframe for time 
time <- data.frame(seq(as.Date("2016-12-02"), as.Date("2019-12-31"), by="day"))
ID <- seq.int(1:1125)
#time <- data.frame(seq(as.Date("2015-07-13"), as.Date("2016-12-01"), by="day"))
#ID <- seq.int(1:508)
time <- cbind(ID, time)
colnames(time) <- c("ID", "DateTime")
colnames(df) <- c("ID", "Temp", "Depth")

#Merge time, depth, temp according to ID
merge <- merge(time, df, all=TRUE)
merge$DateTime <- as.Date(merge$DateTime, format="%Y-%m-%d")

#Temperature estimator function at any depth on certain date
estimate_temp_by_date <- function(target_date, target_depth) {
  data_for_date <- merge %>% 
    filter(DateTime == target_date) %>%
    arrange(Depth)
  
  approx(data_for_date$Depth, data_for_date$Temp, xout = target_depth)$y
}
#Interpolate temp in 0.1m increments in the water column
temp_interp_depth <- crossing(
  tibble(DateTime = unique(merge$DateTime)),
  tibble(Depth = seq(0.1, 9.2, by = 0.1))
) %>%
  group_by(DateTime) %>%
  mutate(Temp = estimate_temp_by_date(DateTime[1], Depth))

modelledt <-ggplot() + 
  geom_raster(temp_interp_depth, mapping =aes(DateTime, Depth, fill = Temp))+
  scale_y_reverse(expand=c(0.01,0.01), limits=c(9.3, 0))+
  scale_x_date(expand=c(0.025,0.025), limits=c(as.Date("2016-12-01"), as.Date("2019-12-06")))+
  #guides(fill=guide_legend(title="Temp", reverse=TRUE, byrow=TRUE))+
  guides(fill = guide_colourbar(ticks.colour= "black"))+
  scale_fill_distiller(name = "°C", palette = 'Spectral', direction = -1, na.value = "grey90", guide="colourbar", limits = c(-1, 30))+
  ylab("Depth (m)")+
  ggtitle("Temperature", subtitle = "Modelled temp")+
  theme(
    legend.background = element_rect(fill="white"),
    legend.spacing.y = unit(3, "mm"),
    plot.title = ggplot2::element_text(face= "bold", size = 12),
    plot.subtitle = ggplot2::element_text(size = 10),
    axis.title.y = ggplot2::element_text(face="bold",size= 10),
    axis.title.x = ggplot2::element_text(face="bold", size= 10),
    legend.text = ggplot2::element_text(size= 8), 
    legend.title = ggplot2::element_text(face="bold", size=10),
    axis.text.x = ggplot2::element_text(size=9),
    axis.text.y = ggplot2::element_text(size=9),
    legend.key.size = unit(7, 'mm'),
  )
modelledt

#Observed temp
obs_temp <- read.csv('observations/CleanedObsTemp.csv') %>%
  mutate(DateTime = as.Date(DateTime, format="%Y-%m-%d"))
obs_temp_int <- read.csv('observations/obs_temp_interpolated.csv')%>%
  mutate(DateTime = as.Date(DateTime, format="%Y-%m-%d"))

#Plot observed temp

observedt<- ggplot() + 
  geom_raster(obs_temp_int, mapping =aes(DateTime, Depth, fill = temp))+
  scale_y_reverse(expand=c(0.01,0.01), limits=c(9.3, 0))+
  scale_x_date(expand=c(0.025,0.025), limits=c(as.Date("2016-12-01"), as.Date("2019-12-06")))+
  guides(fill = guide_colourbar(ticks.colour= "black"))+
  ggtitle(" ", subtitle = "Observed temp")+
  scale_fill_distiller(name = "°C", palette = 'Spectral', direction = -1, na.value = "grey90", guide="colourbar", limits=c(0, 30))+
  geom_point(obs_temp, mapping = aes(x = DateTime, y = Depth), colour="black", pch=4, size=0.5)+
  ylab("Depth (m)")+
  labs(fill = "Legend", colour = NULL)+
  theme(
    legend.background = element_rect(fill="white"),
    legend.spacing.y = unit(3, "mm"),
    plot.title = ggplot2::element_text(face= "bold", size = 12),
    plot.subtitle = ggplot2::element_text(size = 10),
    axis.title.y = ggplot2::element_text(face="bold",size= 10),
    axis.title.x = ggplot2::element_text(face="bold", size= 10),
    legend.text = ggplot2::element_text(size= 8), 
    legend.title = ggplot2::element_text(face="bold", size=10),
    axis.text.x = ggplot2::element_text(size=9),
    axis.text.y = ggplot2::element_text(size=9),
    legend.key.size = unit(7, 'mm'),
  )
observedt

#Difference
merge_temp <- merge(obs_temp_int, temp_interp_depth, all=TRUE)
#modelled - observed
merge_temp$difference <- merge_temp$Temp-merge_temp$temp

#difference plot
differencet <- ggplot() + 
  geom_raster(merge_temp, mapping =aes(DateTime, Depth, fill = difference))+
  scale_y_reverse(expand=c(0.01,0.01), limits=c(9.3, 0))+
  scale_x_date(expand=c(0.025,0.025), limits=c(as.Date("2016-12-01"), as.Date("2019-12-06")))+
  ggtitle(" ", subtitle = "Difference temp (Modelled - Observed)")+
  guides(fill = guide_colourbar(ticks.colour= "black"))+
  ylab("Depth (m)")+
  scale_fill_distiller(name = "°C", palette = 'Spectral', direction = -1, na.value = "grey90", guide="colourbar") +
  theme(
    legend.background = element_rect(fill="white"),
    legend.spacing.y = unit(3, "mm"),
    plot.title = ggplot2::element_text(face= "bold", size = 12),
    plot.subtitle = ggplot2::element_text(size = 10),
    axis.title.y = ggplot2::element_text(face="bold",size= 10),
    axis.title.x = ggplot2::element_text(face="bold", size= 10),
    legend.text = ggplot2::element_text(size= 8), 
    legend.title = ggplot2::element_text(face="bold", size=10),
    axis.text.x = ggplot2::element_text(size=9),
    axis.text.y = ggplot2::element_text(size=9),
    #panel.background = ggplot2::element_rect(fill = 'white', colour="black"),
    #panel.border = ggplot2::element_rect(colour = "black"),
    #legend.position = c(0.98, 0.85),
    #legend.direction="vertical",
    #legend.box = "horizontal",
    legend.key.size = unit(7, 'mm'),
    #legend.key.height = unit(2, "mm")#change legend key size
  )
differencet


mix <- ggarrange(modelledt, modelled, observedt, observed, differencet, difference, ncol=2, nrow=3, labels=c("a", "b", "c", "d", "e", "f"))
mix

ggsave("Results/Figure_S1.png",
       plot = mix,
       width = 246.2, 
       height = 210, 
       units = "mm")

#Figure B2
sim_folder <- getwd()
PEST_calib <- file.path(sim_folder, 'Calibrated_models/Deepm1_naive/output/output.nc')

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
output <- nc_open('Calibrated_models/Deepm1_naive/output/output.nc')
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
nml_file <- read_nml(nml_file = file.path(sim_folder, 'Calibrated_models/Deepm1_naive/glm3.nml'))
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
ggsave("Results/Figure_S2.png",
       plot = mix,
       dpi=500,
       width = 246.2, 
       height = 210, 
       units = "mm")

#Figure B3

sim_folder <- getwd()
#weight 1
w1 <- file.path(sim_folder, 'Calibrated_models/Deepm1_exm_weight1/output/output.nc')
#weight 2
w2 <- file.path(sim_folder, 'Calibrated_models/Deepm1_exm_weight2/output/output.nc')
#weight 3
w3 <- file.path(sim_folder, 'Calibrated_models/Deepm1_exm_weight3/output/output.nc')
error <- read.csv("observations/error_stats.csv")

#oxygen
var="OXY_oxy"
obs_oxy<-read.csv('observations/CleanedObsOxy.csv') %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) %>%
  filter(DateTime > "2016-12-01")
obs_oxy$DateTime <- as.Date(obs_oxy$DateTime, format="%Y-%m-%d")

depths<- c(0.1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9.2) 

oxy_w1 <- get_var(w1, "OXY_oxy", reference="surface", z_out=depths) %>%
  pivot_longer(cols=starts_with("OXY_oxy_"), names_to="Depth", names_prefix="OXY_oxy_", values_to = "OXY_oxy") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) 
oxy_w1$DateTime <- as.Date(oxy_w1$DateTime, format="%Y-%m-%d")

oxy_w2 <- get_var(w2, "OXY_oxy", reference="surface", z_out=depths) %>%
  pivot_longer(cols=starts_with("OXY_oxy_"), names_to="Depth", names_prefix="OXY_oxy_", values_to = "OXY_oxy") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) 
oxy_w2$DateTime <- as.Date(oxy_w2$DateTime, format="%Y-%m-%d")

oxy_w3 <- get_var(w3, "OXY_oxy", reference="surface", z_out=depths) %>%
  pivot_longer(cols=starts_with("OXY_oxy_"), names_to="Depth", names_prefix="OXY_oxy_", values_to = "OXY_oxy") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) 
oxy_w3$DateTime <- as.Date(oxy_w3$DateTime, format="%Y-%m-%d")

oxygen <- merge(obs_oxy, oxy_w1, by=c("DateTime","Depth")) %>%
  dplyr::rename(obsoxy = OXY_oxy.x, mod_oxy_w1 = OXY_oxy.y) %>%
  merge(oxy_w2, by=c("DateTime","Depth")) %>%
  dplyr::rename(mod_oxy_w2 = OXY_oxy) %>%
  merge(oxy_w3, by=c("DateTime","Depth")) %>%
  dplyr::rename(mod_oxy_w3 = OXY_oxy)
oxygen$DateTime <- as.Date(oxygen$DateTime, format="%Y-%m-%d")

#Plots
plot_depths <- c(1, 4, 9)
plot <- vector('list')

for(i in 1:length(plot_depths)){
  
  
  plot[[i]] <- oxygen %>%
    filter(Depth == plot_depths[i]) %>%
    ggplot2::ggplot(ggplot2::aes(x = DateTime, y = obsoxy, colour="Obs. oxy")) +
    ggplot2::geom_point(pch=10)+
    geom_line(data=filter(oxy_w1, Depth==plot_depths[i]), aes(x=DateTime, y=OXY_oxy, colour="PEST_w1 oxy"), lty=1, size=0.5)+
    geom_line(data=filter(oxy_w2, Depth==plot_depths[i]), aes(x=DateTime, y=OXY_oxy, colour="PEST_w2 oxy"), lty=2, size=0.5)+
    geom_line(data=filter(oxy_w3, Depth==plot_depths[i]), aes(x=DateTime, y=OXY_oxy, colour="PEST_w3 oxy"), lty=3, size=0.8)+
    ggplot2::ggtitle(" ", subtitle=" ")+
    xlab("Date")+
    ylab(expression(bold(Oxygen~(mmol/m^{3}))))+
    ylim(c(0, 420))+
    ggplot2::scale_colour_manual(name="Legend", values=c("Obs. oxy"="black", "PEST_w1 oxy"="#00BA38", "PEST_w2 oxy"="#00BA38", "PEST_w3 oxy"="#00BA38"), guide=guide_legend(override.aes=list(linetype=c(NA, 1, 2, 3), shape=c(10, NA, NA, NA))))+
    ggplot2::theme_light() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face= "bold", size = 12),
      axis.title.y = ggplot2::element_text(face = "bold", size= 10),
      axis.title.x = ggplot2::element_text(face = "bold", size= 10),
      legend.text = ggplot2::element_text(size= 8),
      axis.text.x = ggplot2::element_text(size=9),
      axis.text.y = ggplot2::element_text(size=9),
      plot.subtitle=element_text(size=10),
      legend.title = ggplot2::element_blank()
      
    )
  
  
}
plot

combinedPlot <- patchwork::wrap_plots(plot, ncol = 1) +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

combinedPlot

#Temperature
obs_temp<-read.csv('observations/CleanedObsTemp.csv') %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) %>%
  filter(DateTime > "2016-12-01" & DateTime < "2020-01-01")
obs_temp$DateTime <- as.Date(obs_temp$DateTime, format="%Y-%m-%d")

temp_w1 <- get_var(w1, "temp", reference="surface", z_out=depths) %>%
  pivot_longer(cols=starts_with("temp_"), names_to="Depth", names_prefix="temp_", values_to = "temp") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) 
temp_w1$DateTime <- as.Date(temp_w1$DateTime, format="%Y-%m-%d")

temp_w2 <- get_var(w2, "temp", reference="surface", z_out=depths) %>%
  pivot_longer(cols=starts_with("temp_"), names_to="Depth", names_prefix="temp_", values_to = "temp") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) 
temp_w2$DateTime <- as.Date(temp_w2$DateTime, format="%Y-%m-%d")

temp_w3 <- get_var(w3, "temp", reference="surface", z_out=depths) %>%
  pivot_longer(cols=starts_with("temp_"), names_to="Depth", names_prefix="temp_", values_to = "temp") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) 
temp_w3$DateTime <- as.Date(temp_w3$DateTime, format="%Y-%m-%d")

temp <- merge(obs_temp, temp_w1, by=c("DateTime","Depth")) %>%
  dplyr::rename(obtemp = temp.x, modtemp = temp.y) %>%
  merge(temp_w2, by=c("DateTime","Depth")) %>%
  dplyr::rename(mod_temp_w2 = temp) %>%
  merge(temp_w3, by=c("DateTime","Depth")) %>%
  dplyr::rename(mod_temp_w3 = temp)
temp$DateTime <- as.Date(temp$DateTime, format="%Y-%m-%d")

title<- c("Epilimnion", "Metalimnion", "Hypolimnion")
depth<- c("1 m depth", "4 m depth", "9 m depth")

#Plots
plot1 <- vector('list')

for(i in 1:length(plot_depths)){
  
  
  plot1[[i]] <- temp %>%
    filter(Depth == plot_depths[i]) %>%
    ggplot2::ggplot(ggplot2::aes(x = DateTime, y = obtemp, colour="Obs. temp")) +
    ggplot2::geom_point(pch=10)+
    geom_line(data=filter(temp_w1, Depth==plot_depths[i]), aes(x=DateTime, y=temp, colour="PEST_w1 temp"), lty=1, size=0.5)+
    geom_line(data=filter(temp_w2, Depth==plot_depths[i]), aes(x=DateTime, y=temp, colour="PEST_w2 temp"), lty=2, size=0.5)+
    geom_line(data=filter(temp_w3, Depth==plot_depths[i]), aes(x=DateTime, y=temp, colour="PEST_w3 temp"), lty=3, size=0.8)+
    xlab("Date")+
    ylab("Temperature (°C)")+
    ylim(c(0, 30))+
    ggplot2::ggtitle(title[i], subtitle=depth[i])+
    ggplot2::scale_colour_manual(name="Legend", values=c("Obs. temp"="black", "PEST_w1 temp"="#619CFF", "PEST_w2 temp"="#619CFF", "PEST_w3 temp"="#619CFF"), guide=guide_legend(override.aes=list(linetype=c(NA, 1, 2, 3), shape=c(10, NA, NA, NA))))+
    ggplot2::theme_light() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face= "bold", size = 12),
      axis.title.y = ggplot2::element_text(face = "bold", size= 10),
      axis.title.x = ggplot2::element_text(face = "bold", size= 10),
      legend.text = ggplot2::element_text(size= 8),
      axis.text.x = ggplot2::element_text(size=9),
      axis.text.y = ggplot2::element_text(size=9),
      plot.subtitle=element_text(size=10),
      legend.title = ggplot2::element_blank()
    )
  
  
}
plot1

combinedPlot1 <- patchwork::wrap_plots(plot1, ncol = 1) +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
combinedPlot1

final <- ggarrange(combinedPlot1, combinedPlot, ncol=2, nrow=1)
final
ggsave("Results/Figure_S3.png",
       plot = final,
       width = 246.2, 
       height = 210, 
       units = "mm")

#Figure B4
sim_folder <- getwd()
#weight 1
w1 <- file.path(sim_folder, 'Calibrated_models/Deepm1_exm_weight1/output/output.nc')
#weight 2
w2 <- file.path(sim_folder, 'Calibrated_models/Deepm1_exm_weight2/output/output.nc')
#weight 3
w3 <- file.path(sim_folder, 'Calibrated_models/Deepm1_exm_weight3/output/output.nc')

#Observed thermocline depth
obs_TD<- read.csv('observations/obs_td.csv') %>%
  mutate(month = lubridate::month(DateTime)) %>%
  filter(between(month, 4, 9)) %>%
  mutate(year = lubridate::year(DateTime)) %>%
  group_by(year) %>%
  filter(DateTime  > "2016-12-01" & DateTime < "2020-01-01")

obs_TD$DateTime <- as.Date(obs_TD$DateTime, format="%Y-%m-%d")

#Modelled thermocline depths

#model weight 1
temp_w1<- get_var(w1, var_name="temp", reference="surface") 
colClean <- function(x){ colnames(x) <- gsub("temp", "wtr", colnames(x)); x } 
colClean(temp_w1)
temp_w1$DateTime <- as.Date(temp_w1$DateTime, format="%Y-%m-$d")

thermo_depth_model_w1 <- ts.thermo.depth(temp_w1, Smin = 0.1, na.rm=TRUE, seasonal=FALSE)  %>% 
  dplyr::rename(td_model = thermo.depth, DateTime = datetime) %>%
  mutate(month = lubridate::month(DateTime)) %>%
  filter(between(month, 4, 9)) %>%
  mutate(year = lubridate::year(DateTime)) %>%
  group_by(year) 

thermo_depth_model_w1$DateTime <- as.Date(thermo_depth_model_w1$DateTime, format="%Y-%m-%d")

#model weight 2
temp_w2<- get_var(w2, var_name="temp", reference="surface") 
colClean <- function(x){ colnames(x) <- gsub("temp", "wtr", colnames(x)); x } 
colClean(temp_w2)
temp_w2$DateTime <- as.Date(temp_w2$DateTime, format="%Y-%m-$d")

thermo_depth_model_w2 <- ts.thermo.depth(temp_w2, Smin = 0.1, na.rm=TRUE, seasonal=FALSE)  %>% 
  dplyr::rename(td_model = thermo.depth, DateTime = datetime) %>%
  mutate(month = lubridate::month(DateTime)) %>%
  filter(between(month, 4, 9)) %>%
  mutate(year = lubridate::year(DateTime)) %>%
  group_by(year) 

thermo_depth_model_w2$DateTime <- as.Date(thermo_depth_model_w2$DateTime, format="%Y-%m-%d")

#model weight 3
temp_w3<- get_var(w3, var_name="temp", reference="surface") 
colClean <- function(x){ colnames(x) <- gsub("temp", "wtr", colnames(x)); x } 
colClean(temp_w3)
temp_w3$DateTime <- as.Date(temp_w3$DateTime, format="%Y-%m-$d")

thermo_depth_model_w3 <- ts.thermo.depth(temp_w3, Smin = 0.1, na.rm=TRUE, seasonal=FALSE)  %>% 
  dplyr::rename(td_model = thermo.depth, DateTime = datetime) %>%
  mutate(month = lubridate::month(DateTime)) %>%
  filter(between(month, 4, 9)) %>%
  mutate(year = lubridate::year(DateTime)) %>%
  group_by(year) 

thermo_depth_model_w3$DateTime <- as.Date(thermo_depth_model_w3$DateTime, format="%Y-%m-%d")

plot_2 <- thermo_depth_model_w1 %>%
  ggplot2::ggplot(ggplot2::aes(x = DateTime, y = td_model, colour="PEST_w1 TD", group=year), lty=1, size=0.5) +
  ggplot2::geom_line()+
  ggplot2::geom_line(data=thermo_depth_model_w2, aes(x = DateTime, y = td_model, colour="PEST_w2 TD", group=year), lty=2, size=0.5)+
  facet_wrap(~ year, nrow = 1, scales = "free_x")+
  ggplot2::geom_line(data=thermo_depth_model_w3, aes(x = DateTime, y = td_model, colour="PEST_w3 TD", group=year), lty=3, size=0.8)+
  ggplot2::geom_point(data=obs_TD[-1,], aes(x=DateTime, y=thermo.depth, colour="Obs. TD"), pch=10)+
  ggplot2::labs(x = "Date", y = "Depth (m)")+
  ggplot2::scale_y_reverse(limits = c(6, -0.5), breaks= c(0, 1, 2, 3, 4, 5, 6), labels= c("0", "1", "2", "3", "4", "5", "6"))+
  ggplot2::scale_colour_manual(name="", values=c("Obs. TD"="black", "PEST_w1 TD"="#FF61CC", "PEST_w2 TD"="#FF61CC", "PEST_w3 TD"="#FF61CC"), guide=guide_legend(override.aes=list(linetype=c(NA, 1, 2, 3), shape=c(10, NA, NA, NA))))+
  ggplot2::theme_light() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face= "bold", size = 10),
    axis.title.y = ggplot2::element_text(face="bold", size= 10, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.title.x = ggplot2::element_text(face="bold", size= 10),
    axis.text.x = ggplot2::element_text(size=9),
    axis.text.y = ggplot2::element_text(size=9),
    legend.text = ggplot2::element_text(size= 7),
    legend.title = ggplot2::element_blank(),
    strip.text = element_text(size = 12),
    legend.position=c(0.5, 0.92),
    legend.direction="horizontal",
    legend.key.height = unit(2, "mm"),
    legend.spacing.x = unit(1.5, 'mm')
    
  )

plot_2

#Schmidt stability
bathy <- read.csv('observations/bathymetry.csv')

#Schmidt stability observed
schmidt_stability_obs<- read.csv("observations/Obs_SS.csv") %>%
  filter(datetime > "2016-12-01" & datetime < "2020-01-01")
#filter(datetime > "2015-07-12" & datetime < "2016-12-01")
schmidt_stability_obs$datetime<- as.Date(schmidt_stability_obs$datetime, format='%Y-%m-%d')

#Schmidt stability model weight 1

schmidt_stability_w1 <- ts.schmidt.stability(temp_w1, bathy, na.rm=TRUE) %>% 
  dplyr::rename(ss_w1 = schmidt.stability)
schmidt_stability_w1$datetime <- as.Date(schmidt_stability_w1$datetime, format="%Y-%m-%d")

#Schmidt stability model weight 2

schmidt_stability_w2 <- ts.schmidt.stability(temp_w2, bathy, na.rm=TRUE) %>% 
  dplyr::rename(ss_w2 = schmidt.stability)
schmidt_stability_w2$datetime <- as.Date(schmidt_stability_w2$datetime, format="%Y-%m-%d")

#Schmidt stability model weight 3

schmidt_stability_w3 <- ts.schmidt.stability(temp_w3, bathy, na.rm=TRUE) %>% 
  dplyr::rename(ss_w3 = schmidt.stability)
schmidt_stability_w3$datetime <- as.Date(schmidt_stability_w3$datetime, format="%Y-%m-%d")


SS_plot <- ggplot(data=schmidt_stability_obs, aes(x=datetime, y=schmidt.stability, colour="Obs. SS")) +
  geom_point(pch=10)+
  ylab("Schmidt stability")+
  xlab("Date")+
  ylim(c(0, 70))+
  geom_line(data=schmidt_stability_w1, aes(x=datetime, y=ss_w1, colour="PEST_w1 SS"), lty=1, size=0.5)+
  geom_line(data=schmidt_stability_w2, aes(x=datetime, y=ss_w2, colour="PEST_w2 SS"), lty=2, size=0.5)+
  geom_line(data=schmidt_stability_w3, aes(x=datetime, y=ss_w3, colour="PEST_w3 SS"), lty=3, size=0.8)+
  scale_x_date(expand=c(0.01,0.01))+
  scale_colour_manual(values=c("Obs. SS"="black", "PEST_w1 SS" ="#00BFC4", "PEST_w2 SS" ="#00BFC4", "PEST_w3 SS" ="#00BFC4"), guide=guide_legend(override.aes=list(linetype=c(NA, 1, 2, 3), shape=c(10, NA, NA, NA))))+
  ggplot2::theme_light() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face= "bold", size = 12),
    axis.title.y = ggplot2::element_text(face="bold", size= 10, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.title.x = ggplot2::element_text(face="bold", size= 10),
    legend.text = ggplot2::element_text(size= 7),
    axis.text.x = ggplot2::element_text(size=9),
    axis.text.y = ggplot2::element_text(size=9),
    legend.title = ggplot2::element_blank(),
    legend.position = c(0.5, 0.92),
    legend.direction="horizontal",
    plot.margin = unit(c(5.5,10,5.5,10), "pt"),
    legend.key.height = unit(2, "mm"),
    legend.spacing.x = unit(1.5, 'mm')
    
    
  )
SS_plot

#MOM observed
obs_mom<-read.csv('observations/mom_observed.csv') %>%
  #filter(DateTime > "2015-07-12" & DateTime < "2017-01-01")
  filter(DateTime > "2016-12-02" & DateTime < "2020-01-01")
obs_mom$DateTime <- as.Date(obs_mom$DateTime, format = "%Y-%m-%d")

#Model weight 1 mom
depths<- c(1, 4, 8) 
oxy_w1 <- get_var(w1, "OXY_oxy", reference="surface", z_out=depths) %>%
  pivot_longer(cols=starts_with("OXY_oxy_"), names_to="Depth", names_prefix="OXY_oxy_", values_to = "OXY_oxy") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) 
oxy_w1$DateTime <- as.Date(oxy_w1$DateTime, format="%Y-%m-%d")

epi_oxy <- filter(oxy_w1, Depth==1)
hypo_oxy <- filter(oxy_w1, Depth==8)
met_oxy <- filter(oxy_w1, Depth==4)
merge_w1<- merge(epi_oxy, hypo_oxy, by="DateTime") %>%
  dplyr::rename(epi_oxy = OXY_oxy.x, hypo_oxy = OXY_oxy.y)
merge_w1$exp_oxy <- (merge_w1$epi_oxy+merge_w1$hypo_oxy)/2

merge_w1 <- merge(merge_w1, met_oxy[, c("DateTime", "OXY_oxy")], by="DateTime")%>%
  dplyr::rename(met_oxy = OXY_oxy)
merge_w1$deviation <- merge_w1$met_oxy - merge_w1$exp_oxy

#Model weight 2 mom
oxy_w2 <- get_var(w2, "OXY_oxy", reference="surface", z_out=depths) %>%
  pivot_longer(cols=starts_with("OXY_oxy_"), names_to="Depth", names_prefix="OXY_oxy_", values_to = "OXY_oxy") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) 
oxy_w2$DateTime <- as.Date(oxy_w2$DateTime, format="%Y-%m-%d")

epi_oxy <- filter(oxy_w2, Depth==1)
hypo_oxy <- filter(oxy_w2, Depth==8)
met_oxy <- filter(oxy_w2, Depth==4)
merge_w2<- merge(epi_oxy, hypo_oxy, by="DateTime") %>%
  dplyr::rename(epi_oxy = OXY_oxy.x, hypo_oxy = OXY_oxy.y)
merge_w2$exp_oxy <- (merge_w2$epi_oxy+merge_w2$hypo_oxy)/2

merge_w2 <- merge(merge_w2, met_oxy[, c("DateTime", "OXY_oxy")], by="DateTime")%>%
  dplyr::rename(met_oxy = OXY_oxy)
merge_w2$deviation <- merge_w2$met_oxy - merge_w2$exp_oxy

#Model weight 3 mom
oxy_w3 <- get_var(w3, "OXY_oxy", reference="surface", z_out=depths) %>%
  pivot_longer(cols=starts_with("OXY_oxy_"), names_to="Depth", names_prefix="OXY_oxy_", values_to = "OXY_oxy") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) 
oxy_w3$DateTime <- as.Date(oxy_w3$DateTime, format="%Y-%m-%d")

epi_oxy <- filter(oxy_w3, Depth==1)
hypo_oxy <- filter(oxy_w3, Depth==8)
met_oxy <- filter(oxy_w3, Depth==4)
merge_w3<- merge(epi_oxy, hypo_oxy, by="DateTime") %>%
  dplyr::rename(epi_oxy = OXY_oxy.x, hypo_oxy = OXY_oxy.y)
merge_w3$exp_oxy <- (merge_w3$epi_oxy+merge_w3$hypo_oxy)/2

merge_w3 <- merge(merge_w3, met_oxy[, c("DateTime", "OXY_oxy")], by="DateTime")%>%
  dplyr::rename(met_oxy = OXY_oxy)
merge_w3$deviation <- merge_w3$met_oxy - merge_w3$exp_oxy

obs_mom$DateTime <- as.Date(obs_mom$DateTime, format="%Y-%m-%d")
merge_w1$DateTime <- as.Date(merge_w1$DateTime, format="%Y-%m-%d")
merge_w2$DateTime <- as.Date(merge_w2$DateTime, format="%Y-%m-%d")
merge_w3$DateTime <- as.Date(merge_w3$DateTime, format="%Y-%m-%d")

plot_MOM <-ggplot(data=merge_w1, aes(x=DateTime, y=deviation, colour="PEST_w1 MOM"))+
  geom_line(lty=1, size=0.5)+
  geom_line(data=merge_w2, aes(x=DateTime, y=deviation, colour="PEST_w2 MOM"), lty=2, size=0.5)+
  geom_line(data=merge_w3, aes(x=DateTime, y=deviation, colour="PEST_w3 MOM"), lty=3, size=0.8)+
  geom_point(data=obs_mom, aes(x=DateTime, y=deviation, colour="Obs. MOM"), pch=10)+
  geom_hline(yintercept=0, lty=3)+
  xlab("Date")+
  ylim(c(-250, 250))+
  ylab(expression(bold(MOM~(mmol/m^{3}))))+
  scale_x_date(expand=c(0.01,0.01))+
  ggplot2::scale_colour_manual(name="Legend", values=c("Obs. MOM"="black", "PEST_w1 MOM" ="#CD9600", "PEST_w2 MOM" ="#CD9600", "PEST_w3 MOM" ="#CD9600"), guide=guide_legend(override.aes=list(linetype=c(NA, 1, 2, 3), shape=c(10, NA, NA, NA))))+
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
    plot.margin = unit(c(5.5,10,5.5,5.5), "pt"),
    legend.key.height = unit(2, "mm"),
    legend.spacing.x = unit(1.5, 'mm'),
    legend.margin=margin(2, 2, 2, 2)
    
  )
plot_MOM

#Anoxia
#IMPORTANT: This code is to be run in three different runs, make sure you follow the instructions in the comment sections

# START FIRST MODEL RUN HERE
output <- nc_open('Calibrated_models/Deepm1_exm_weight1/output/output.nc')

#START SECOND AND THIRD RUN HERE
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

# apply beetroot over each column in the oxy_out dataframe
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
#################STOP SECOND RUN HERE, RUN LINE 1274TO LINE 1278  ###########
#################STOP THIRD RUN HERE, RUN LINE 1280 ONWARDS ###########

#Model anoxia 
anoxia <- oxy_interp_depth
anoxia$Oxy <- anoxia$Oxy * 32/1000
anoxia$Oxy<- ifelse(anoxia$Oxy<=1, 1, 0)
#anoxia <- na.omit(anoxia)
anoxia$DateTime <- as.Date(anoxia$DateTime, format="%Y-%m-%d")

#model weight 2 
output <- nc_open('Calibrated_models/Deepm1_exm_weight2/output/output.nc')
#####END OF FIRST RUN, START SECOND RUN FROM LINE 1182, RUN IT UNTIL LINE 1257 #########

# Only second run (Model 2)

anoxia$Oxy1 <- oxy_interp_depth$Oxy
anoxia$Oxy1 <- anoxia$Oxy1 * 32/1000
anoxia$Oxy1<- ifelse(anoxia$Oxy1<=1, 1, 0)
output <- nc_open('Calibrated_models/Deepm1_exm_weight3/output/output.nc')
#####END OF SECOND RUN, START THIRD RUN FROM LINE 1182, RUN UNTIL LINE 1257 #########

# Only third run (Model 3)
anoxia$Oxy2 <- oxy_interp_depth$Oxy
anoxia$Oxy2 <- anoxia$Oxy2 * 32/1000
anoxia$Oxy2<- ifelse(anoxia$Oxy2<=1, 1, 0)

anoxia$sum <- anoxia$Oxy + anoxia$Oxy1 + anoxia$Oxy2

#Observed oxygen anoxia
obs_oxy<-read.csv('observations/CleanedObsOxy.csv') %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) %>%
  #filter(DateTime > "2015-07-12" & DateTime < "2016-12-31")
  filter(DateTime > "2016-12-01" & DateTime < "2020-01-01")

obs_anoxia <- obs_oxy %>%
  mutate(OXY_oxy=OXY_oxy*32/1000) %>%
  mutate(OXY_oxy=ifelse(OXY_oxy<=1, 1, 0)) %>%
  na.omit() 
obs_anoxia$DateTime <- as.Date(obs_anoxia$DateTime, format="%Y-%m-%d")

#Anoxia plot more models
PEST_anoxia<- ggplot() +
  geom_raster(anoxia, mapping =aes(DateTime, Depth, fill = as.character(sum))) +
  ylab("Depth (m)")+
  xlab("Date")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_y_reverse(expand=c(0,0)) +
  scale_x_date(expand=c(0,0))+
  geom_point(filter(obs_anoxia, obs_anoxia$OXY_oxy==1), mapping = aes(x = DateTime, y = Depth, colour = "Observed anoxia"), pch=4)+
  scale_fill_manual(
    labels = c("3 models", "2 models", "1 model", "0 models"),
    values = c("3"="red", "2"="#F8766D", "1"= "bisque", "0"="blue")
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
    legend.key.height = unit(2, "mm")
    
  )

PEST_anoxia
mix1 <- ggarrange(plot_2, SS_plot, PEST_anoxia, plot_MOM, ncol=2, nrow=2, labels=c("a", "b", "c", "d"))
mix1
ggsave("Results/Figure_S4.png",
       plot = mix1,
       width = 246.2, 
       height = 160, 
       #scale=1.6,
       dpi=500,
       units = "mm")
