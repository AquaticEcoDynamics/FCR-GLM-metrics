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
naive <- file.path(sim_folder, 'Calibrated_models/Deepm2_naive/output/output.nc')
#weight 2 (best performing model system-inspired model)
system_insp <- file.path(sim_folder, 'Calibrated_models/Deepm2_exm_weight2/output/output.nc')
error <- read.csv("Observations/error_stats.csv")

#Observed thermocline depth
obs_TD<- read.csv('Observations/obs_td.csv') %>%
  mutate(DateTime = as.Date(DateTime, format="%Y-%m-%d")) %>%
  mutate(month = lubridate::month(DateTime)) %>%
  filter(between(month, 4, 9)) %>%
  mutate(year = lubridate::year(DateTime)) %>%
  group_by(year) %>%
  filter(DateTime  > "2016-12-01" & DateTime < "2020-01-01")

obs_TD$DateTime <- as.Date(obs_TD$DateTime, format="%Y-%m-%d")

#Modelled thermocline depths

#naive model
temp_naive<- get_var(naive, var_name="temp", reference="surface") 
colClean <- function(x){ colnames(x) <- gsub("temp", "wtr", colnames(x)); x } 
colClean(temp_naive)
temp_naive$DateTime <- as.Date(temp_naive$DateTime, format="%Y-%m-$d")

thermo_depth_model_naive <- ts.thermo.depth(temp_naive, Smin = 0.1, na.rm=TRUE, seasonal=FALSE)  %>% 
  dplyr::rename(td_model = thermo.depth, DateTime = datetime) %>%
  mutate(month = lubridate::month(DateTime)) %>%
  filter(between(month, 4, 9)) %>%
  mutate(year = lubridate::year(DateTime)) %>%
  group_by(year) 

thermo_depth_model_naive$DateTime <- as.Date(thermo_depth_model_naive$DateTime, format="%Y-%m-%d")

#MEF calculation for TD naive model
td_merge <- merge(thermo_depth_model_naive, obs_TD, by="DateTime")

for (i in 1:nrow(td_merge)) {
  td_merge$MEFF_1[i]<- ((td_merge$td_model[i]- td_merge$thermo.depth[i])^2)
  td_merge$MEFF_2[i]<- ((td_merge$thermo.depth[i]-mean(td_merge$thermo.depth))^2)
  MEFF_TD<- 1-(sum(td_merge$MEFF_1)/sum(td_merge$MEFF_2))
}

#Adding calculated MEF to error table
#error[error$metric=="TD" & error$calibration=="PEST_N", "Calibration.deepm2"] <- MEFF_TD

#System-inspired (w2)
temp_system_insp<- get_var(system_insp, var_name="temp", reference="surface") 
colClean <- function(x){ colnames(x) <- gsub("temp", "wtr", colnames(x)); x } 
colClean(temp_system_insp)
temp_system_insp$DateTime <- as.Date(temp_system_insp$DateTime, format="%Y-%m-$d")

thermo_depth_system_insp <- ts.thermo.depth(temp_system_insp, Smin = 0.1, na.rm=TRUE, seasonal=FALSE)  %>% 
  dplyr::rename(td_model = thermo.depth, DateTime = datetime) %>%
  mutate(month = lubridate::month(DateTime)) %>%
  filter(between(month, 4, 9)) %>%
  mutate(year = lubridate::year(DateTime)) %>%
  group_by(year) 

thermo_depth_system_insp$DateTime <- as.Date(thermo_depth_system_insp$DateTime, format="%Y-%m-%d")

plot_2 <- thermo_depth_model_naive %>%
  ggplot2::ggplot(ggplot2::aes(x = DateTime, y = td_model, colour="Naive model TD", group=year)) +
  ggplot2::geom_line(lty=2, linewidth=0.5)+
  ggplot2::geom_line(data=thermo_depth_system_insp, aes(x = DateTime, y = td_model, colour="System-inspired model TD", group=year), lty=1, linewidth=0.5)+
  facet_wrap(~ year, nrow = 1, scales = "free_x")+
  ggplot2::geom_point(data=obs_TD[-1,], aes(x=DateTime, y=thermo.depth, colour="Obs. TD"), pch=10)+
  ggplot2::labs(x = "Date", y = "Depth (m)")+
  ggplot2::scale_y_reverse(limits = c(6, -0.5), breaks= c(0, 1, 2, 3, 4, 5, 6), labels= c("0", "1", "2", "3", "4", "5", "6"))+
  ggplot2::scale_colour_manual(name="", values=c("Obs. TD"="black", "Naive model TD"="#FF61CC", "System-inspired model TD"="#FF61CC"), labels=c("Obs. TD", "Naive model TD", "System-inspired model TD"), guide=guide_legend(override.aes=list(linetype=c(NA, 2, 1), shape=c(10, NA, NA), color = c('black', '#FF61CC', '#FF61CC'))))+
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
bathy <- read.csv('Observations/bathymetry.csv') %>%
  rename(depths = X...depths)

#Schmidt stability observed
schmidt_stability_obs<- read.csv("Observations/Obs_SS.csv") %>%
  mutate(datetime = as.Date(datetime, format="%Y-%m-%d")) %>%
  filter(datetime > "2016-12-01" & datetime < "2020-01-01")
#filter(datetime > "2015-07-12" & datetime < "2016-12-01")
schmidt_stability_obs$datetime<- as.Date(schmidt_stability_obs$datetime, format='%Y-%m-%d')

#Schmidt stability naive

schmidt_stability_naive <- ts.schmidt.stability(temp_naive, bathy, na.rm=TRUE) %>% 
  dplyr::rename(ss_naive = schmidt.stability)
schmidt_stability_naive$datetime <- as.Date(schmidt_stability_naive$datetime, format="%Y-%m-%d")

#Calculating error (MEF) of schmidt stability for Figure 7
SS_merge <- merge(schmidt_stability_obs, schmidt_stability_naive)

for (i in 1:nrow(SS_merge)) {
  SS_merge$MEFF_1[i]<- ((SS_merge$ss_PEST[i]- SS_merge$schmidt.stability[i])^2)
  SS_merge$MEFF_2[i]<- ((SS_merge$schmidt.stability[i]-mean(SS_merge$schmidt.stability))^2)
  MEFF_SS<- 1-(sum(SS_merge$MEFF_1)/sum(SS_merge$MEFF_2))
}

#Adding calculated MEF to error table
#error[error$metric=="SS" & error$calibration=="PEST_N", "Calibration.deepm2"] <- MEFF_SS

#Schmidt stability system-inpired model (w2)

schmidt_stability_system_insp <- ts.schmidt.stability(temp_system_insp, bathy, na.rm=TRUE) %>% 
  dplyr::rename(ss_system_insp = schmidt.stability)
schmidt_stability_system_insp$datetime <- as.Date(schmidt_stability_system_insp$datetime, format="%Y-%m-%d")

SS_plot <- ggplot(data=schmidt_stability_obs, aes(x=datetime, y=schmidt.stability, colour="Obs. SS")) +
  geom_point(pch=10)+
  ylab("Schmidt stability")+
  xlab("Date")+
  ylim(c(0, 70))+
  geom_line(data=schmidt_stability_naive, aes(x=datetime, y=ss_naive, colour="Naive model SS"), lty=2, size=0.5)+
  geom_line(data=schmidt_stability_system_insp, aes(x=datetime, y=ss_system_insp, colour="System-inspired model SS"), lty=1, size=0.5)+
  scale_x_date(expand=c(0.01,0.01))+
  scale_colour_manual(values=c("Obs. SS"="black", "Naive model SS" ="#00BFC4", "System-inspired model SS" ="#00BFC4"), labels=c("Obs. SS", "Naive model SS", "System-inspired model SS"),  guide=guide_legend(override.aes=list(linetype=c(NA, 2, 1), shape=c(10, NA, NA), colour = c("black", "#00BFC4", "#00BFC4"))))+
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
obs_mom<-read.csv('Observations/mom_observed.csv') %>%
  #filter(DateTime > "2015-07-12" & DateTime < "2017-01-01")
  mutate(DateTime = as.Date(DateTime, format="%Y-%m-%d")) %>%
  filter(DateTime > "2016-12-02" & DateTime < "2020-01-01")
obs_mom$DateTime <- as.Date(obs_mom$DateTime, format = "%Y-%m-%d")

#Model naive mom
depths<- c(1, 4, 8) 
oxy_naive <- get_var(naive, "OXY_oxy", reference="surface", z_out=depths) %>%
  pivot_longer(cols=starts_with("OXY_oxy_"), names_to="Depth", names_prefix="OXY_oxy_", values_to = "OXY_oxy") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) 
oxy_naive$DateTime <- as.Date(oxy_naive$DateTime, format="%Y-%m-%d")

epi_oxy <- filter(oxy_naive, Depth==1)
hypo_oxy <- filter(oxy_naive, Depth==8)
met_oxy <- filter(oxy_naive, Depth==4)
merge_naive<- merge(epi_oxy, hypo_oxy, by="DateTime") %>%
  dplyr::rename(epi_oxy = OXY_oxy.x, hypo_oxy = OXY_oxy.y)
merge_naive$exp_oxy <- (merge_naive$epi_oxy+merge_naive$hypo_oxy)/2

merge_naive <- merge(merge_naive, met_oxy[, c("DateTime", "OXY_oxy")], by="DateTime")%>%
  dplyr::rename(met_oxy = OXY_oxy)
merge_naive$deviation <- merge_naive$met_oxy - merge_naive$exp_oxy

#Calculating error (MEF) of metalimnetic oxygen minima for Figure 7
merge_mom <- merge(obs_mom, merge_naive, by="DateTime")

for (i in 1:nrow(merge_mom)) {
  merge_mom$MEFF_1[i]<- ((merge_mom$deviation.y[i]- merge_mom$deviation.x[i])^2)
  merge_mom$MEFF_2[i]<- ((merge_mom$deviation.x[i]-mean(merge_mom$deviation.x))^2)
  MEFF_MOM<- 1-(sum(merge_mom$MEFF_1)/sum(merge_mom$MEFF_2))
}

#Adding calculated MEF to error table
#error[error$metric=="MOM" & error$calibration=="PEST_N", "Calibration.deepm2"] <- MEFF_MOM

#Model system-inspired w2 mom
oxy_system_insp <- get_var(system_insp, "OXY_oxy", reference="surface", z_out=depths) %>%
  pivot_longer(cols=starts_with("OXY_oxy_"), names_to="Depth", names_prefix="OXY_oxy_", values_to = "OXY_oxy") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) 
oxy_system_insp$DateTime <- as.Date(oxy_system_insp$DateTime, format="%Y-%m-%d")

epi_oxy <- filter(oxy_system_insp, Depth==1)
hypo_oxy <- filter(oxy_system_insp, Depth==8)
met_oxy <- filter(oxy_system_insp, Depth==4)
merge_system_insp<- merge(epi_oxy, hypo_oxy, by="DateTime") %>%
  dplyr::rename(epi_oxy = OXY_oxy.x, hypo_oxy = OXY_oxy.y)
merge_system_insp$exp_oxy <- (merge_system_insp$epi_oxy+merge_system_insp$hypo_oxy)/2

merge_system_insp <- merge(merge_system_insp, met_oxy[, c("DateTime", "OXY_oxy")], by="DateTime")%>%
  dplyr::rename(met_oxy = OXY_oxy)
merge_system_insp$deviation <- merge_system_insp$met_oxy - merge_system_insp$exp_oxy

obs_mom$DateTime <- as.Date(obs_mom$DateTime, format="%Y-%m-%d")
merge_naive$DateTime <- as.Date(merge_naive$DateTime, format="%Y-%m-%d")
merge_system_insp$DateTime <- as.Date(merge_system_insp$DateTime, format="%Y-%m-%d")


plot_MOM <-ggplot(data=obs_mom, aes(x=DateTime, y=deviation, colour="Obs. MOM"))+
  geom_point(pch = 10)+
  geom_line(data=merge_naive, aes(x=DateTime, y=deviation, colour="Naive model MOM"),lty=2, linewidth=0.5)+
  geom_line(data=merge_system_insp, aes(x=DateTime, y=deviation, colour="System-inspired model MOM"), lty=1, linewidth=0.5)+
  geom_hline(yintercept=0, lty=3)+
  xlab("Date")+
  ylim(c(-250, 250))+
  ylab(expression(bold(MOM~(mmol/m^{3}))))+
  scale_x_date(expand=c(0.01,0.01))+
  ggplot2::scale_colour_manual(name="Legend", values=c("Obs. MOM"="black", "Naive model MOM" ="#CD9600", "System-inspired model MOM" ="#CD9600"), labels=c("Obs. MOM", "Naive model MOM", "System-inspired model MOM"), guide=guide_legend(override.aes=list(linetype=c(NA, 2, 1), shape=c(10, NA, NA), color=c("black", "#CD9600", "#CD9600"))))+
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
output <- nc_open('Calibrated_models/Deepm2_naive/output/output.nc')

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
#################STOP SECOND RUN HERE  ###########

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

#Observed oxygen anoxia
obs_oxy<-read.csv('Observations/CleanedObsOxy.csv') %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) %>%
  #filter(DateTime > "2015-07-12" & DateTime < "2016-12-31")
  filter(DateTime > "2016-12-01" & DateTime < "2020-01-01")

obs_anoxia <- obs_oxy %>%
  mutate(OXY_oxy=OXY_oxy*32/1000) %>%
  mutate(OXY_oxy=ifelse(OXY_oxy<=1, 1, 0)) %>%
  na.omit() 
obs_anoxia$DateTime <- as.Date(obs_anoxia$DateTime, format="%Y-%m-%d")

anoxia <- anoxia %>%
  mutate(Oxy1 = ifelse(Oxy1 == 1, 2, Oxy1))

anoxia$sum <- anoxia$Oxy + anoxia$Oxy1

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
    labels = c("None", "Naive model", "System-inspired model", "Both"),
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
ggsave("Results/Figure4.png",
       plot = mix1,
       width = 246.2, 
       height = 160, 
       #scale=1.6,
       dpi=500,
       units = "mm")


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

merge_anoxia <- merge(obs_anoxic_layers, newData_mod, by="DateTime")

for (i in 1:nrow(merge_anoxia)) {
  merge_anoxia$MEFF_1[i]<- ((merge_anoxia$Count.y[i]- merge_anoxia$Count.x[i])^2)
  merge_anoxia$MEFF_2[i]<- ((merge_anoxia$Count.x[i]-mean(merge_anoxia$Count.x))^2)
  MEFF_anoxia<- 1-(sum(merge_anoxia$MEFF_1)/sum(merge_anoxia$MEFF_2))
}

#Adding calculated MEF to error table
#error[error$metric=="A" & error$calibration=="PEST_N", "Calibration.deepm2"] <- MEFF_anoxia
#write.csv(error, 'observations/error_stats.csv')

