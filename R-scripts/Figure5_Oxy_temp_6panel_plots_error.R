#Packages
library(tidyr)
library(dplyr)
library(glmtools)
library(GLMr)
library(lubridate)
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
library(patchwork)

#Set working directory
setwd(".../FCR-GLM-metrics")
sim_folder <- getwd()

#Models w1, w2, w3
w1 <- file.path(sim_folder, 'Calibrated_models/Deepm2_exm_weight1/output/output.nc')

oxy_w1 <- get_var(w1, "OXY_oxy", reference="surface", z_out=depths) %>%
  pivot_longer(cols=starts_with("OXY_oxy_"), names_to="Depth", names_prefix="OXY_oxy_", values_to = "OXY_oxy") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) 
oxy_w1$DateTime <- as.Date(oxy_w1$DateTime, format="%Y-%m-%d")

w2 <- file.path(sim_folder, 'Calibrated_models/Deepm2_exm_weight2/output/output.nc')

oxy_w2 <- get_var(w2, "OXY_oxy", reference="surface", z_out=depths) %>%
  pivot_longer(cols=starts_with("OXY_oxy_"), names_to="Depth", names_prefix="OXY_oxy_", values_to = "OXY_oxy") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) 
oxy_w2$DateTime <- as.Date(oxy_w2$DateTime, format="%Y-%m-%d")

w3 <- file.path(sim_folder, 'Calibrated_models/Deepm2_exm_weight3/output/output.nc')

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


#oxygen
var="OXY_oxy"
obs_oxy<-read.csv('Observations/CleanedObsOxy.csv') %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) %>%
  filter(DateTime > "2016-12-01")
obs_oxy$DateTime <- as.Date(obs_oxy$DateTime, format="%Y-%m-%d")

depths<- c(0.1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9.2) 

#Calculating model absolute error over time (prediction - obs)
oxygen <- oxygen %>%
  mutate(error_w1 = abs(mod_oxy_w1 - obsoxy),
         error_w2 = abs(mod_oxy_w2 - obsoxy),
         error_w3 = abs(mod_oxy_w3 - obsoxy),
         obs_dates = 0)


#Error calculation w1
for (i in 1:nrow(oxygen)) {
  oxygen$MEFF_1_w1[i]<- ((oxygen$mod_oxy_w1[i]- oxygen$obsoxy[i])^2)
  oxygen$MEFF_2[i]<- ((oxygen$obsoxy[i]-mean(oxygen$obsoxy))^2)
  MEFF_oxy_w1<- 1-(sum(oxygen$MEFF_1_w1)/sum(oxygen$MEFF_2))
}

#Adding calculated MEF to error table
#error[error$metric=="oxy" & error$calibration=="PEST_exm_w1", "Calibration.deepm2"] <- MEFF_oxy_w1

#Error calculation w2
for (i in 1:nrow(oxygen)) {
  oxygen$MEFF_1_w2[i]<- ((oxygen$mod_oxy_w2[i]- oxygen$obsoxy[i])^2)
  MEFF_oxy_w2<- 1-(sum(oxygen$MEFF_1_w2)/sum(oxygen$MEFF_2))
}

#Adding calculated MEF to error table
#error[error$metric=="oxy" & error$calibration=="PEST_exm_w2", "Calibration.deepm2"] <- MEFF_oxy_w2

#Error calculation w3
for (i in 1:nrow(oxygen)) {
  oxygen$MEFF_1_w3[i]<- ((oxygen$mod_oxy_w3[i]- oxygen$obsoxy[i])^2)
  MEFF_oxy_w3<- 1-(sum(oxygen$MEFF_1_w3)/sum(oxygen$MEFF_2))
}

#Adding calculated MEF to error table
#error[error$metric=="oxy" & error$calibration=="PEST_exm_w3", "Calibration.deepm2"] <- MEFF_oxy_w3

#Oxygen error plot
plot_depths <- c(1, 4, 9)
plot_error_DO <- vector('list')

for(i in 1:length(plot_depths)){
  
  
  plot_error_DO[[i]] <- oxygen %>%
    filter(Depth == plot_depths[i]) %>%
    ggplot2::ggplot(ggplot2::aes(x = DateTime, y = obs_dates, colour="Obs. dates")) +
    ggplot2::geom_point(pch=4)+
    geom_line(data=filter(oxygen, Depth==plot_depths[i]), aes(x=DateTime, y=error_w1, colour="Model w1 DO"), lty=2, linewidth=0.5)+
    geom_line(data=filter(oxygen, Depth==plot_depths[i]), aes(x=DateTime, y=error_w2, colour="Model w2 DO"), lty=1, linewidth=0.5)+
    geom_line(data=filter(oxygen, Depth==plot_depths[i]), aes(x=DateTime, y=error_w3, colour="Model w3 DO"), lty=3, linewidth=0.5)+
    ggplot2::ggtitle(" ", subtitle = "Oxygen")+
    xlab("Date")+
    ylab(expression(bold(Absolute~model~error~(mmol/m^{3}))))+
    ylim(c(0, 250))+
    ggplot2::scale_colour_manual(name="Legend", values=c("Obs. dates" = "black", "Model w1 DO"="#00BA38", "Model w2 DO"="#00BA38", "Model w3 DO"="#00BA38"),
    labels=c("Obs. dates", "Model w1 DO", "Model w2 DO", "Model w3 DO"), guide=guide_legend(override.aes=list(linetype=c(NA, 2, 1, 3), shape=c(4, NA, NA, NA), colour=c("black", "#00BA38", "#00BA38", "#00BA38"))))+
    ggplot2::theme_light() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face= "bold", size = 10),
      axis.title.y = ggplot2::element_text(face = "bold", size= 8),
      axis.title.x = ggplot2::element_text(face = "bold", size= 8),
      legend.text = ggplot2::element_text(size= 8),
      axis.text.x = ggplot2::element_text(size=9),
      axis.text.y = ggplot2::element_text(size=9),
      plot.subtitle=element_text(size=10),
      legend.title = ggplot2::element_blank()
      
    )
  
  
}

plot_error_DO

combinedPlot_DO <- patchwork::wrap_plots(plot_error_DO, ncol = 1) +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

combinedPlot_DO

#Temperature
obs_temp<-read.csv('Observations/CleanedObsTemp.csv') %>%
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
  dplyr::rename(obtemp = temp.x, mod_temp_w1 = temp.y) %>%
  merge(temp_w2, by=c("DateTime","Depth")) %>%
  dplyr::rename(mod_temp_w2 = temp) %>%
  merge(temp_w3, by=c("DateTime","Depth")) %>%
  dplyr::rename(mod_temp_w3 = temp)
temp$DateTime <- as.Date(temp$DateTime, format="%Y-%m-%d")

temp <- temp %>%
  mutate(error_w1 = abs(mod_temp_w1 - obtemp),
         error_w2 = abs(mod_temp_w2 - obtemp),
         error_w3 = abs(mod_temp_w3 - obtemp),
         obs_dates = 0)

#Error calculation w1
for (i in 1:nrow(temp)) {
  temp$MEFF_1_w1[i]<- ((temp$modtemp[i]- temp$obtemp[i])^2)
  temp$MEFF_2[i]<- ((temp$obtemp[i]-mean(temp$obtemp))^2)
  MEFF_temp_w1<- 1-(sum(temp$MEFF_1_w1)/sum(temp$MEFF_2))
}

#Adding calculated MEF to error table
#error[error$metric=="temp" & error$calibration=="PEST_exm_w1", "Calibration.deepm2"] <- MEFF_temp_w1

#Error calculation w2
for (i in 1:nrow(temp)) {
  temp$MEFF_1_w2[i]<- ((temp$mod_temp_w2[i]- temp$obtemp[i])^2)
  MEFF_temp_w2<- 1-(sum(temp$MEFF_1_w2)/sum(temp$MEFF_2))
}

#Adding calculated MEF to error table
#error[error$metric=="temp" & error$calibration=="PEST_exm_w2", "Calibration.deepm2"] <- MEFF_temp_w2

#Error calculation w3
for (i in 1:nrow(temp)) {
  temp$MEFF_1_w3[i]<- ((temp$mod_temp_w3[i]- temp$obtemp[i])^2)
  MEFF_temp_w3<- 1-(sum(temp$MEFF_1_w3)/sum(temp$MEFF_2))
}

#Adding calculated MEF to error table
#error[error$metric=="temp" & error$calibration=="PEST_exm_w3", "Calibration.deepm2"] <- MEFF_temp_w3
#write.csv(error, 'observations/error_stats.csv')

#Temperature error plots 
title<- c("Epilimnion", "Metalimnion", "Hypolimnion")
depth<- c("1 m depth", "4 m depth", "9 m depth")
plot_title <- c("Eplimnion (1 m depth)", "Metalimnion (4 m depth)", "Hypolimnion (9 m depth)")

plot_depths <- c(1, 4, 9)
plot_error_temp <- vector('list')

for(i in 1:length(plot_depths)){
  
  
  plot_error_temp[[i]] <- temp %>%
    filter(Depth == plot_depths[i]) %>%
    ggplot2::ggplot(ggplot2::aes(x = DateTime, y = obs_dates, colour="Obs. dates")) +
    ggplot2::geom_point(pch=4)+
    geom_line(data=filter(temp, Depth==plot_depths[i]), aes(x=DateTime, y=error_w1, colour="Model w1 temp"), lty=2, linewidth=0.5)+
    geom_line(data=filter(temp, Depth==plot_depths[i]), aes(x=DateTime, y=error_w2, colour="Model w2 temp"), lty=1, linewidth=0.5)+
    geom_line(data=filter(temp, Depth==plot_depths[i]), aes(x=DateTime, y=error_w3, colour="Model w3 temp"), lty=3, linewidth=0.5)+
    ggplot2::ggtitle(" ", subtitle=" ")+
    xlab("Date")+
    ylab("Absolute model error (\u00B0C)")+
    ggplot2::ggtitle(plot_title[i], subtitle= "Temperature")+
    ylim(c(0, 5))+
    ggplot2::scale_colour_manual(name="Legend", values=c("Obs. dates" = "black", "Model w1 temp"="#619CFF", "Model w2 temp"="#619CFF", "Model w3 temp"="#619CFF"),
    labels=c("Obs. dates", "Model w1 temp", "Model w2 temp", "Model w3 temp"), guide=guide_legend(override.aes=list(linetype=c(NA, 2, 1, 3), shape=c(4, NA, NA, NA), colour=c("black", "#619CFF", "#619CFF", "#619CFF"))))+
    ggplot2::theme_light() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face= "bold", size = 10),
      axis.title.y = ggplot2::element_text(face = "bold", size= 8),
      axis.title.x = ggplot2::element_text(face = "bold", size= 8),
      legend.text = ggplot2::element_text(size= 8),
      axis.text.x = ggplot2::element_text(size=9),
      axis.text.y = ggplot2::element_text(size=9),
      plot.subtitle=element_text(size=10),
      legend.title = ggplot2::element_blank()
      
    )
  
  
}

plot_error_temp

combinedPlot_temp <- patchwork::wrap_plots(plot_error_temp, ncol = 1) +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

combinedPlot_temp

final <- ggarrange(combinedPlot_temp, combinedPlot_DO, ncol=2, nrow=1)
final
ggsave("/Results/Figure5.png",
       plot = final,
       width = 246.2, 
       height = 210, 
       units = "mm")



