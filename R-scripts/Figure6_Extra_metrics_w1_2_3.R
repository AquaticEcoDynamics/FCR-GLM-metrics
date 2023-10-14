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

setwd(".../FCR-GLM-metrics")
sim_folder <- getwd()
#weight 1
w1 <- file.path(sim_folder, 'Calibrated_models/Deepm2_exm_weight1/output/output.nc')
#weight 2
w2 <- file.path(sim_folder, 'Calibrated_models/Deepm2_exm_weight2/output/output.nc')
#weight 3
w3 <- file.path(sim_folder, 'Calibrated_models/Deepm2_exm_weight3/output/output.nc')
error <- read.csv("observations/error_stats.csv")

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

#Calculating error (MEF) of thermocline depth for Figure 7
td_merge_w1 <- merge(thermo_depth_model_w1, obs_TD, by="DateTime")

for (i in 1:nrow(td_merge_w1)) {
  td_merge_w1$MEFF_1[i]<- ((td_merge_w1$td_model[i]- td_merge_w1$thermo.depth[i])^2)
  td_merge_w1$MEFF_2[i]<- ((td_merge_w1$thermo.depth[i]-mean(td_merge_w1$thermo.depth))^2)
  MEFF_TD_w1<- 1-(sum(td_merge_w1$MEFF_1)/sum(td_merge_w1$MEFF_2))
}

#Adding calculated MEF to error table
#error[error$metric=="TD" & error$calibration=="PEST_exm_w1", "Calibration.deepm2"] <- MEFF_TD_w1

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

#Calculating error (MEF) of thermocline depth for Figure 7
td_merge_w2 <- merge(thermo_depth_model_w2, obs_TD, by="DateTime")

for (i in 1:nrow(td_merge_w2)) {
  td_merge_w2$MEFF_1[i]<- ((td_merge_w2$td_model[i]- td_merge_w2$thermo.depth[i])^2)
  td_merge_w2$MEFF_2[i]<- ((td_merge_w2$thermo.depth[i]-mean(td_merge_w2$thermo.depth))^2)
  MEFF_TD_w2<- 1-(sum(td_merge_w2$MEFF_1)/sum(td_merge_w2$MEFF_2))
}

#Adding calculated MEF to error table
#error[error$metric=="TD" & error$calibration=="PEST_exm_w2", "Calibration.deepm2"] <- MEFF_TD_w2

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

#Calculating error (MEF) of thermocline depth for Figure 7
td_merge_w3 <- merge(thermo_depth_model_w3, obs_TD, by="DateTime")

for (i in 1:nrow(td_merge_w3)) {
  td_merge_w3$MEFF_1[i]<- ((td_merge_w3$td_model[i]- td_merge_w3$thermo.depth[i])^2)
  td_merge_w3$MEFF_2[i]<- ((td_merge_w3$thermo.depth[i]-mean(td_merge_w3$thermo.depth))^2)
  MEFF_TD_w3<- 1-(sum(td_merge_w3$MEFF_1)/sum(td_merge_w3$MEFF_2))
}

#Adding calculated MEF to error table
#error[error$metric=="TD" & error$calibration=="PEST_exm_w3", "Calibration.deepm2"] <- MEFF_TD_w3

plot_2 <- thermo_depth_model_w1 %>%
  ggplot2::ggplot(ggplot2::aes(x = DateTime, y = td_model, colour="Model w1 TD", group=year), lty=1, size=0.5) +
  ggplot2::geom_line()+
  ggplot2::geom_line(data=thermo_depth_model_w2, aes(x = DateTime, y = td_model, colour="Model w2 TD", group=year), lty=2, size=0.5)+
  facet_wrap(~ year, nrow = 1, scales = "free_x")+
  ggplot2::geom_line(data=thermo_depth_model_w3, aes(x = DateTime, y = td_model, colour="Model w3 TD", group=year), lty=3, size=0.8)+
  ggplot2::geom_point(data=obs_TD[-1,], aes(x=DateTime, y=thermo.depth, colour="Obs. TD"), pch=10)+
  ggplot2::labs(x = "Date", y = "Depth (m)")+
  ggplot2::scale_y_reverse(limits = c(6, -0.5), breaks= c(0, 1, 2, 3, 4, 5, 6), labels= c("0", "1", "2", "3", "4", "5", "6"))+
  ggplot2::scale_colour_manual(name="", values=c("Obs. TD"="black", "Model w1 TD"="#FF61CC", "Model w2 TD"="#FF61CC", "Model w3 TD"="#FF61CC"), guide=guide_legend(override.aes=list(linetype=c(NA, 1, 2, 3), shape=c(10, NA, NA, NA))))+
  ggplot2::theme_light() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face= "bold", size = 10),
    axis.title.y = ggplot2::element_text(face="bold", size= 10, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.title.x = ggplot2::element_text(face="bold", size= 10),
    axis.text.x = ggplot2::element_text(size=9),
    axis.text.y = ggplot2::element_text(size=9),
    legend.text = ggplot2::element_text(size= 8),
    legend.title = ggplot2::element_blank(),
    strip.text = element_text(size = 12),
    legend.position=c(0.5, 0.92),
    legend.direction="horizontal",
    legend.key.height = unit(2, "mm"),
    legend.key.width = unit(4, "mm"),
    legend.spacing.x = unit(1.5, 'mm'),
    panel.spacing = unit(0.6, "lines")
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

#Calculating error (MEF) of schmidt stability for Figure 7
SS_merge_w1 <- merge(schmidt_stability_obs, schmidt_stability_w1)

for (i in 1:nrow(SS_merge_w1)) {
  SS_merge_w1$MEFF_1[i]<- ((SS_merge_w1$ss_w1[i]- SS_merge_w1$schmidt.stability[i])^2)
  SS_merge_w1$MEFF_2[i]<- ((SS_merge_w1$schmidt.stability[i]-mean(SS_merge_w1$schmidt.stability))^2)
  MEFF_SS_w1<- 1-(sum(SS_merge_w1$MEFF_1)/sum(SS_merge_w1$MEFF_2))
}


#Adding calculated MEF to error table
#error[error$metric=="SS" & error$calibration=="PEST_exm_w1", "Calibration.deepm2"] <- MEFF_SS_w1

#Schmidt stability model weight 2

schmidt_stability_w2 <- ts.schmidt.stability(temp_w2, bathy, na.rm=TRUE) %>% 
  dplyr::rename(ss_w2 = schmidt.stability)
schmidt_stability_w2$datetime <- as.Date(schmidt_stability_w2$datetime, format="%Y-%m-%d")

#Calculating error (MEF) of schmidt stability for Figure 7
SS_merge_w2 <- merge(schmidt_stability_obs, schmidt_stability_w2)

for (i in 1:nrow(SS_merge_w2)) {
  SS_merge_w2$MEFF_1[i]<- ((SS_merge_w2$ss_w2[i]- SS_merge_w2$schmidt.stability[i])^2)
  SS_merge_w2$MEFF_2[i]<- ((SS_merge_w2$schmidt.stability[i]-mean(SS_merge_w2$schmidt.stability))^2)
  MEFF_SS_w2<- 1-(sum(SS_merge_w2$MEFF_1)/sum(SS_merge_w2$MEFF_2))
}

#Adding calculated MEF to error table
#error[error$metric=="SS" & error$calibration=="PEST_exm_w2", "Calibration.deepm2"] <- MEFF_SS_w2

#Schmidt stability model weight 3

schmidt_stability_w3 <- ts.schmidt.stability(temp_w3, bathy, na.rm=TRUE) %>% 
  dplyr::rename(ss_w3 = schmidt.stability)
schmidt_stability_w3$datetime <- as.Date(schmidt_stability_w3$datetime, format="%Y-%m-%d")

#Calculating error (MEF) of schmidt stability for Figure 7
SS_merge_w3 <- merge(schmidt_stability_obs, schmidt_stability_w3)

for (i in 1:nrow(SS_merge_w3)) {
  SS_merge_w3$MEFF_1[i]<- ((SS_merge_w3$ss_w3[i]- SS_merge_w3$schmidt.stability[i])^2)
  SS_merge_w3$MEFF_2[i]<- ((SS_merge_w3$schmidt.stability[i]-mean(SS_merge_w3$schmidt.stability))^2)
  MEFF_SS_w3<- 1-(sum(SS_merge_w3$MEFF_1)/sum(SS_merge_w3$MEFF_2))
}


#Adding calculated MEF to error table
#error[error$metric=="SS" & error$calibration=="PEST_exm_w3", "Calibration.deepm2"] <- MEFF_SS_w3

SS_plot <- ggplot(data=schmidt_stability_obs, aes(x=datetime, y=schmidt.stability, colour="Obs. SS")) +
  geom_point(pch=10)+
  ylab("Schmidt stability")+
  xlab("Date")+
  ylim(c(0, 70))+
  geom_line(data=schmidt_stability_w1, aes(x=datetime, y=ss_w1, colour="Model w1 SS"), lty=1, size=0.5)+
  geom_line(data=schmidt_stability_w2, aes(x=datetime, y=ss_w2, colour="Model w2 SS"), lty=2, size=0.5)+
  geom_line(data=schmidt_stability_w3, aes(x=datetime, y=ss_w3, colour="Model w3 SS"), lty=3, size=0.8)+
  scale_x_date(expand=c(0.01,0.01))+
  scale_colour_manual(values=c("Obs. SS"="black", "Model w1 SS" ="#00BFC4", "Model w2 SS" ="#00BFC4", "Model w3 SS" ="#00BFC4"), guide=guide_legend(override.aes=list(linetype=c(NA, 1, 2, 3), shape=c(10, NA, NA, NA))))+
  ggplot2::theme_light() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face= "bold", size = 12),
    axis.title.y = ggplot2::element_text(face="bold", size= 10, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.title.x = ggplot2::element_text(face="bold", size= 10),
    legend.text = ggplot2::element_text(size= 8),
    axis.text.x = ggplot2::element_text(size=9),
    axis.text.y = ggplot2::element_text(size=9),
    legend.title = ggplot2::element_blank(),
    legend.position = c(0.5, 0.92),
    legend.direction="horizontal",
    plot.margin = unit(c(5.5,10,5.5,10), "pt"),
    legend.key.width = unit(4, "mm"),
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

#Calculating error (MEF) of metalimnetic oxygen minima for Figure 7
merge_mom_w1 <- merge(obs_mom, merge_w1, by="DateTime")

for (i in 1:nrow(merge_mom_w1)) {
  merge_mom_w1$MEFF_1[i]<- ((merge_mom_w1$deviation.y[i]- merge_mom_w1$deviation.x[i])^2)
  merge_mom_w1$MEFF_2[i]<- ((merge_mom_w1$deviation.x[i]-mean(merge_mom_w1$deviation.x))^2)
  MEFF_MOM_w1<- 1-(sum(merge_mom_w1$MEFF_1)/sum(merge_mom_w1$MEFF_2))
}


#Adding calculated MEF to error table
#error[error$metric=="MOM" & error$calibration=="PEST_exm_w1", "Calibration.deepm2"] <- MEFF_MOM_w1

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

#Calculating error (MEF) of metalimnetic oxygen minima for Figure 7
merge_mom_w2 <- merge(obs_mom, merge_w2, by="DateTime")

for (i in 1:nrow(merge_mom_w2)) {
  merge_mom_w2$MEFF_1[i]<- ((merge_mom_w2$deviation.y[i]- merge_mom_w2$deviation.x[i])^2)
  merge_mom_w2$MEFF_2[i]<- ((merge_mom_w2$deviation.x[i]-mean(merge_mom_w2$deviation.x))^2)
  MEFF_MOM_w2<- 1-(sum(merge_mom_w2$MEFF_1)/sum(merge_mom_w2$MEFF_2))
}


#Adding calculated MEF to error table
#error[error$metric=="MOM" & error$calibration=="PEST_exm_w2", "Calibration.deepm2"] <- MEFF_MOM_w2

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

#Calculating error (MEF) of metalimnetic oxygen minima for Figure 7
merge_mom_w3 <- merge(obs_mom, merge_w3, by="DateTime")

for (i in 1:nrow(merge_mom_w3)) {
  merge_mom_w3$MEFF_1[i]<- ((merge_mom_w3$deviation.y[i]- merge_mom_w3$deviation.x[i])^2)
  merge_mom_w3$MEFF_2[i]<- ((merge_mom_w3$deviation.x[i]-mean(merge_mom_w3$deviation.x))^2)
  MEFF_MOM_w3<- 1-(sum(merge_mom_w3$MEFF_1)/sum(merge_mom_w3$MEFF_2))
}

#Adding calculated MEF to error table
#error[error$metric=="MOM" & error$calibration=="PEST_exm_w3", "Calibration.deepm2"] <- MEFF_MOM_w3

obs_mom$DateTime <- as.Date(obs_mom$DateTime, format="%Y-%m-%d")
merge_w1$DateTime <- as.Date(merge_w1$DateTime, format="%Y-%m-%d")
merge_w2$DateTime <- as.Date(merge_w2$DateTime, format="%Y-%m-%d")
merge_w3$DateTime <- as.Date(merge_w3$DateTime, format="%Y-%m-%d")

plot_MOM <-ggplot(data=merge_w1, aes(x=DateTime, y=deviation, colour="Model w1 MOM"))+
  geom_line(lty=1, size=0.5)+
  geom_line(data=merge_w2, aes(x=DateTime, y=deviation, colour="Model w2 MOM"), lty=2, size=0.5)+
  geom_line(data=merge_w3, aes(x=DateTime, y=deviation, colour="Model w3 MOM"), lty=3, size=0.8)+
  geom_point(data=obs_mom, aes(x=DateTime, y=deviation, colour="Obs. MOM"), pch=10)+
  geom_hline(yintercept=0, lty=3)+
  xlab("Date")+
  ylim(c(-250, 250))+
  ylab(expression(bold(MOM~(mmol/m^{3}))))+
  scale_x_date(expand=c(0.01,0.01))+
  ggplot2::scale_colour_manual(name="Legend", values=c("Obs. MOM"="black", "Model w1 MOM" ="#CD9600", "Model w2 MOM" ="#CD9600", "Model w3 MOM" ="#CD9600"), guide=guide_legend(override.aes=list(linetype=c(NA, 1, 2, 3), shape=c(10, NA, NA, NA))))+
  ggplot2::theme_light() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face= "bold", size = 12),
    axis.title.y = ggplot2::element_text(face="bold", size= 10),
    axis.title.x = ggplot2::element_text(face="bold", size= 10),
    legend.text = ggplot2::element_text(size= 8), 
    axis.text.x = ggplot2::element_text(size=9),
    axis.text.y = ggplot2::element_text(size=9),
    legend.title = ggplot2::element_blank(),
    legend.position = c(0.5, 0.92),
    legend.direction="horizontal",
    plot.margin = unit(c(5.5,10,5.5,5.5), "pt"),
    legend.key.height = unit(2, "mm"),
    legend.key.width = unit(4, "mm"),
    legend.spacing.x = unit(1.5, 'mm'),
    legend.margin = margin(2, 2, 2, 2)
  )
plot_MOM

#Anoxia
#IMPORTANT: This code is to be run in three different runs, make sure you follow the instructions in the comment sections
# START FIRST MODEL RUN HERE
output <- nc_open('Calibrated_models/Deepm2_exm_weight1/output/output.nc')

#All runs (STAR SECOND AND THIRD RUN HERE)
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
#################STOP SECOND RUN HERE, RUN LINE 470 - 476  ###########
#################STOP THIRD RUN HERE, RUN LINE 478 ONWARDS ###########

#Model anoxia 
anoxia <- oxy_interp_depth
anoxia$Oxy <- anoxia$Oxy * 32/1000
anoxia$Oxy<- ifelse(anoxia$Oxy<=1, 1, 0)
#anoxia <- na.omit(anoxia)
anoxia$DateTime <- as.Date(anoxia$DateTime, format="%Y-%m-%d")

#model weight 2 
output <- nc_open('Calibrated_models/Deepm2_exm_weight2/output/output.nc')
#####END OF FIRST RUN, START SECOND RUN FROM LINE 380 - LINE 455 #########

# Only second run (Model 2)

anoxia$Oxy1 <- oxy_interp_depth$Oxy
anoxia$Oxy1 <- anoxia$Oxy1 * 32/1000
anoxia$Oxy1<- ifelse(anoxia$Oxy1<=1, 1, 0)
output <- nc_open('Calibrated_models/Deepm2_exm_weight3/output/output.nc')
#####END OF SECOND RUN, START THIRD RUN FROM LINE 380 - LINE 455 #########

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
  geom_point(filter(obs_anoxia, obs_anoxia$OXY_oxy==1), mapping = aes(x = DateTime, y = Depth, colour = "Obs. anoxia"), pch=4)+
  scale_fill_manual(
    labels = c("3 models", "2 models", "1 model", "0 models"),
    values = c("3"="red", "2"="#F8766D", "1"= "bisque", "0"="blue")
  )+
  scale_colour_manual(
    values = c("Obs. anoxia"="black")
  ) +
  guides(fill = guide_legend(order = 1)) +
  theme(
    legend.background = element_rect(fill="white"),
    legend.spacing = unit(-6, "pt"),
    plot.title = ggplot2::element_text(face= "bold", size = 12),
    axis.title.y = ggplot2::element_text(face="bold",size= 10),
    axis.title.x = ggplot2::element_text(face="bold", size= 10),
    legend.text = ggplot2::element_text(size= 8), 
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
mix1 <- ggarrange(plot_2, SS_plot, PEST_anoxia, plot_MOM, ncol=2, nrow=2, labels=c("a)", "b)", "c)", "d)"), font.label = list(size = 12, color = "black", face= "plain"))
mix1
ggsave("Results/Figure6.png",
       plot = mix1,
       width = 246.2, 
       height = 160, 
       #scale=1.6,
       dpi=500,
       units = "mm")

#Calculating error (MEF) of number of anoxic layers per day for Figure 7
#Observations
obs_anoxic_layers <- read.csv('observations/anoxia_observed.csv') %>%
  mutate(DateTime=as.Date(DateTime, format="%Y-%m-%d"))

#Model w1
anoxia <- mutate(anoxia, month = lubridate::month(DateTime)) %>%
  filter(between(month, 5, 11))

#Creating empty dataframe for loop
uniqueDates  <- unique(anoxia$DateTime)

newData_mod_w1  <- data.frame(
  DateTime = unique(anoxia$DateTime),
  Count = length(uniqueDates)
)

#Calculating number of anoxic layers each day
for(i in 1:length(uniqueDates)){
  
  filteredData  <- filter(anoxia, DateTime==uniqueDates[i] & Oxy==1)
  newData_mod_w1$Count[i] <- nrow(filteredData)
}

merge_anoxia_w1 <- merge(obs_anoxic_layers, newData_mod_w1, by="DateTime")

for (i in 1:nrow(merge_anoxia_w1)) {
  merge_anoxia_w1$MEFF_1[i]<- ((merge_anoxia_w1$Count.y[i]- merge_anoxia_w1$Count.x[i])^2)
  merge_anoxia_w1$MEFF_2[i]<- ((merge_anoxia_w1$Count.x[i]-mean(merge_anoxia_w1$Count.x))^2)
  MEFF_anoxia_w1<- 1-(sum(merge_anoxia_w1$MEFF_1)/sum(merge_anoxia_w1$MEFF_2))
}

#Adding calculated MEF to error table
#error[error$metric=="A" & error$calibration=="PEST_exm_w1", "Calibration.deepm2"] <- MEFF_anoxia_w1

#Model w2
newData_mod_w2  <- data.frame(
  DateTime = unique(anoxia$DateTime),
  Count = length(uniqueDates)
)

#Calculating number of anoxic layers each day
for(i in 1:length(uniqueDates)){
  
  filteredData  <- filter(anoxia, DateTime==uniqueDates[i] & Oxy1==1)
  newData_mod_w2$Count[i] <- nrow(filteredData)
}

merge_anoxia_w2 <- merge(obs_anoxic_layers, newData_mod_w2, by="DateTime")

for (i in 1:nrow(merge_anoxia_w2)) {
  merge_anoxia_w2$MEFF_1[i]<- ((merge_anoxia_w2$Count.y[i]- merge_anoxia_w2$Count.x[i])^2)
  merge_anoxia_w2$MEFF_2[i]<- ((merge_anoxia_w2$Count.x[i]-mean(merge_anoxia_w2$Count.x))^2)
  MEFF_anoxia_w2<- 1-(sum(merge_anoxia_w2$MEFF_1)/sum(merge_anoxia_w2$MEFF_2))
}

#Adding calculated MEF to error table
#error[error$metric=="A" & error$calibration=="PEST_exm_w2", "Calibration.deepm2"] <- MEFF_anoxia_w2

#Model w3
newData_mod_w3  <- data.frame(
  DateTime = unique(anoxia$DateTime),
  Count = length(uniqueDates)
)

#Calculating number of anoxic layers each day
for(i in 1:length(uniqueDates)){
  
  filteredData  <- filter(anoxia, DateTime==uniqueDates[i] & Oxy2==1)
  newData_mod_w3$Count[i] <- nrow(filteredData)
}

merge_anoxia_w3 <- merge(obs_anoxic_layers, newData_mod_w3, by="DateTime")

for (i in 1:nrow(merge_anoxia_w3)) {
  merge_anoxia_w3$MEFF_1[i]<- ((merge_anoxia_w3$Count.y[i]- merge_anoxia_w3$Count.x[i])^2)
  merge_anoxia_w3$MEFF_2[i]<- ((merge_anoxia_w3$Count.x[i]-mean(merge_anoxia_w3$Count.x))^2)
  MEFF_anoxia_w3<- 1-(sum(merge_anoxia_w3$MEFF_1)/sum(merge_anoxia_w3$MEFF_2))
}


#Adding calculated MEF to error table
#error[error$metric=="A" & error$calibration=="PEST_exm_w3", "Calibration.deepm2"] <- MEFF_anoxia_w3

#write.csv(error, 'observations/error_stats.csv', row.names=FALSE)

