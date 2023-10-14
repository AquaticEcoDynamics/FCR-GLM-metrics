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

#Set working directory
setwd(".../FCR-GLM-metrics")
sim_folder <- getwd()
error <- read.csv("observations/error_stats.csv")

#Observed oxygen
obs_oxy<-read.csv('observations/CleanedObsOxy.csv') %>%
  filter(DateTime > "2016-12-01")
obs_oxy$DateTime <- as.Date(obs_oxy$DateTime, format="%Y-%m-%d")

depths<- c(0.1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9.2) 

#Deepm2 naive
Deepm2_naive <- file.path(sim_folder, 'Calibrated_models/Deepm2_naive/output/output.nc')

deepm2_naive_oxy <- get_var(Deepm2_naive, "OXY_oxy", reference="surface", z_out=depths) %>%
  pivot_longer(cols=starts_with("OXY_oxy_"), names_to="Depth", names_prefix="OXY_oxy_", values_to = "OXY_oxy")
deepm2_naive_oxy$DateTime <- as.Date(deepm2_naive_oxy$DateTime, format="%Y-%m-%d")

oxygen <- merge(obs_oxy, deepm2_naive_oxy, by=c("DateTime","Depth")) %>%
  dplyr::rename(obsoxy = OXY_oxy.x, mod_oxy = OXY_oxy.y) 
oxygen$DateTime <- as.Date(oxygen$DateTime, format="%Y-%m-%d")

for (i in 1:nrow(oxygen)) {
  oxygen$MEFF_1[i]<- ((oxygen$mod_oxy[i]- oxygen$obsoxy[i])^2)
  oxygen$MEFF_2[i]<- ((oxygen$obsoxy[i]-mean(oxygen$obsoxy))^2)
  MEFF_oxy<- 1-(sum(oxygen$MEFF_1)/sum(oxygen$MEFF_2))
}

#Adding calculated MEF to error table
#error[error$metric=="oxy" & error$calibration=="PEST_N", "Calibration.deepm2"] <- MEFF_oxy

#Deepm1 naive
Deepm1_naive <- file.path(sim_folder, 'Calibrated_models/Deepm1_naive/output/output.nc')

deepm1_naive_oxy <- get_var(Deepm1_naive, "OXY_oxy", reference="surface", z_out=depths) %>%
  pivot_longer(cols=starts_with("OXY_oxy_"), names_to="Depth", names_prefix="OXY_oxy_", values_to = "OXY_oxy")
deepm1_naive_oxy$DateTime <- as.Date(deepm1_naive_oxy$DateTime, format="%Y-%m-%d")

oxygen1 <- merge(obs_oxy, deepm1_naive_oxy, by=c("DateTime","Depth")) %>%
  dplyr::rename(obsoxy = OXY_oxy.x, mod_oxy = OXY_oxy.y) 
oxygen1$DateTime <- as.Date(oxygen1$DateTime, format="%Y-%m-%d")

for (i in 1:nrow(oxygen1)) {
  oxygen1$MEFF_1[i]<- ((oxygen1$mod_oxy[i]- oxygen1$obsoxy[i])^2)
  oxygen1$MEFF_2[i]<- ((oxygen1$obsoxy[i]-mean(oxygen1$obsoxy))^2)
  MEFF_oxy1<- 1-(sum(oxygen1$MEFF_1)/sum(oxygen1$MEFF_2))
}

#Adding calculated MEF to error table
#error[error$metric=="oxy" & error$calibration=="PEST_N", "Calibration.deepm1"] <- MEFF_oxy1

#Observed temperature
obs_temp<-read.csv('observations/CleanedObsTemp.csv') %>%
  filter(DateTime > "2016-12-01")
obs_temp$DateTime <- as.Date(obs_temp$DateTime, format="%Y-%m-%d")

#Deepm2 naive
deepm2_naive_temp <- get_var(Deepm2_naive, "temp", reference="surface", z_out=depths) %>%
  pivot_longer(cols=starts_with("temp_"), names_to="Depth", names_prefix="temp_", values_to = "temp") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) 
deepm2_naive_temp$DateTime <- as.Date(deepm2_naive_temp$DateTime, format="%Y-%m-%d")

temp <- merge(obs_temp, deepm2_naive_temp, by=c("DateTime","Depth")) %>%
  dplyr::rename(obtemp = temp.x, modtemp = temp.y) 
temp$DateTime <- as.Date(temp$DateTime, format="%Y-%m-%d")

for (i in 1:nrow(temp)) {
  temp$MEFF_1[i]<- ((temp$modtemp[i]- temp$obtemp[i])^2)
  temp$MEFF_2[i]<- ((temp$obtemp[i]-mean(temp$obtemp))^2)
  MEFF_t<- 1-(sum(temp$MEFF_1)/sum(temp$MEFF_2))
}

#error[error$metric=="temp" & error$calibration=="PEST_N", "Calibration.deepm2"] <- MEFF_t

#Deepm1 naive
deepm1_naive_temp <- get_var(Deepm1_naive, "temp", reference="surface", z_out=depths) %>%
  pivot_longer(cols=starts_with("temp_"), names_to="Depth", names_prefix="temp_", values_to = "temp") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) 
deepm1_naive_temp$DateTime <- as.Date(deepm1_naive_temp$DateTime, format="%Y-%m-%d")

temp1 <- merge(obs_temp, deepm1_naive_temp, by=c("DateTime","Depth")) %>%
  dplyr::rename(obtemp = temp.x, modtemp = temp.y) 
temp1$DateTime <- as.Date(temp1$DateTime, format="%Y-%m-%d")

for (i in 1:nrow(temp1)) {
  temp1$MEFF_1[i]<- ((temp1$modtemp[i]- temp1$obtemp[i])^2)
  temp1$MEFF_2[i]<- ((temp1$obtemp[i]-mean(temp1$obtemp))^2)
  MEFF_t1<- 1-(sum(temp1$MEFF_1)/sum(temp1$MEFF_2))
}

#error[error$metric=="temp" & error$calibration=="PEST_N", "Calibration.deepm1"] <- MEFF_t1
#write.csv(error, 'observations/error_stats.csv', row.names=FALSE)

#ERROR PLOT
error<- read.csv("observations/error_stats.csv")

#Schmid stabiliy error plot
SS<- ggplot(data=filter(error, metric == "SS"), aes(x=calibration, y=Calibration.deepm1, colour = calibration))+
  geom_rect(ymin = 0.5, ymax = 1.5, 
            xmin = -Inf, xmax = Inf, fill = '#C3D7A4', colour='NA', alpha=0.2) +
  geom_rect(ymin = 0, ymax = 0.5, 
            xmin = -Inf, xmax = Inf, fill = '#F4EDCA', colour= 'NA', alpha=0.2) + 
  geom_rect(ymin = -1, ymax = 0, 
            xmin = -Inf, xmax = Inf, fill = 'salmon', colour='NA', alpha=0.1) +
  geom_point(pch=16, size=2)+  
  ggtitle("Schmidt stability")+
  geom_point(data=filter(error, metric == "SS"), aes(x=calibration, y=Calibration.deepm2), colour="blue", pch=15, size=2)+
  geom_point(data=filter(error, metric == "SS"), aes(x=calibration, y=Validation.deepm2), colour="blue", pch=0, size=2)+
  geom_point(data=filter(error, metric == "SS"), aes(x=calibration, y=Validation.deepm1, colour=calibration), pch=1, size=2) +
  scale_x_discrete(limits = c("reference", "PEST_N", "PEST_exm_w1", "PEST_exm_w2", "PEST_exm_w3"), labels = c("Reference", "Naive", "Model w1", "Model w2", "Model w3"))+
  scale_color_manual(name = "Calibration", 
                     labels = c("Reference", "Naive", "Model w1", "Model w2", "Model w3"),
                     values = c("darkorchid2", "darkorchid2", "darkorchid2", "darkorchid2", "black"))+
  ylab("MEF")+
  xlab("Model")+ 
  ylim(c(-0.5, 1))+
  theme(
    #legend.background = element_rect(fill="white"),
    #legend.spacing = unit(-6, "pt"),
    plot.title = ggplot2::element_text(face= "bold", size = 12),
    axis.title.y = ggplot2::element_text(face="bold",size= 10),
    axis.title.x = ggplot2::element_text(face="bold", size= 10, margin = margin(t = 5, r = 0, b = 0, l = 0)),
    legend.text = ggplot2::element_text(size= 8), 
    legend.title = ggplot2::element_blank(),
    axis.text.x = element_text(size=8),
    axis.text.y = element_text(size=9)
    #legend.position = c(0.5, 0.92),
    #legend.direction="horizontal",
    #legend.box = "horizontal",
    #legend.key.size = unit(3, 'mm'),
    #legend.key.height = unit(2, "mm")#change legend key size
  )

SS

TD<- ggplot(data=filter(error, metric == "TD"), aes(x=calibration, y=Calibration.deepm1, colour = calibration))+
  geom_rect(ymin = 0.5, ymax = 1.5, 
            xmin = -Inf, xmax = Inf, fill = '#C3D7A4', colour='NA', alpha=0.2) +
  geom_rect(ymin = 0, ymax = 0.5, 
            xmin = -Inf, xmax = Inf, fill = '#F4EDCA', colour= 'NA', alpha=0.2) + 
  geom_rect(ymin = -1, ymax = 0, 
            xmin = -Inf, xmax = Inf, fill = 'salmon', colour='NA', alpha=0.1) +
  geom_point(pch=16, size=2)+  
  ggtitle("Thermocline depth")+
  geom_point(data=filter(error, metric == "TD"), aes(x=calibration, y=Calibration.deepm2), colour="blue", pch=15, size=2)+
  geom_point(data=filter(error, metric == "TD"), aes(x=calibration, y=Validation.deepm2), colour="blue", pch=0, size=2)+
  geom_point(data=filter(error, metric == "TD"), aes(x=calibration, y=Validation.deepm1, colour=calibration), pch=1, size=2) +
  scale_x_discrete(limits = c("reference", "PEST_N", "PEST_exm_w1", "PEST_exm_w2", "PEST_exm_w3"), labels = c("Reference", "Naive", "Model w1", "Model w2", "Model w3"))+
  scale_color_manual(name = "Calibration", 
                     labels = c("Reference", "Naive", "Model w1", "Model w2", "Model w3"),
                     values = c("darkorchid2", "darkorchid2", "darkorchid2", "darkorchid2", "black"))+
  ylab("MEF")+
  xlab("Model")+ 
  ylim(c(-0.5, 1))+
  theme(
    #legend.background = element_rect(fill="white"),
    #legend.spacing = unit(-6, "pt"),
    plot.title = ggplot2::element_text(face= "bold", size = 12),
    axis.title.y = ggplot2::element_text(face="bold",size= 10),
    axis.title.x = ggplot2::element_text(face="bold", size= 10, margin = margin(t = 5, r = 0, b = 0, l = 0)),
    legend.text = ggplot2::element_text(size= 8), 
    legend.title = ggplot2::element_blank(),
    axis.text.x = element_text(size=8),
    axis.text.y = element_text(size=9)
    #legend.position = c(0.5, 0.92),
    #legend.direction="horizontal",
    #legend.box = "horizontal",
    #legend.key.size = unit(3, 'mm'),
    #legend.key.height = unit(2, "mm")#change legend key size
  )

TD

MOM<- ggplot(data=filter(error, metric == "MOM"), aes(x=calibration, y=Calibration.deepm1, colour = calibration))+
  geom_rect(ymin = 0.5, ymax = 1.5, 
            xmin = -Inf, xmax = Inf, fill = '#C3D7A4', colour='NA', alpha=0.2) +
  geom_rect(ymin = 0, ymax = 0.5, 
            xmin = -Inf, xmax = Inf, fill = '#F4EDCA', colour= 'NA', alpha=0.2) + 
  geom_rect(ymin = -1, ymax = 0, 
            xmin = -Inf, xmax = Inf, fill = 'salmon', colour='NA', alpha=0.1) +
  geom_point(pch=16, size=2)+  
  ggtitle("Metalimnetic oxygen minimum")+
  geom_point(data=filter(error, metric == "MOM"), aes(x=calibration, y=Calibration.deepm2), colour="blue", pch=15, size=2)+
  geom_point(data=filter(error, metric == "MOM"), aes(x=calibration, y=Validation.deepm2), colour="blue", pch=0, size=2)+
  geom_point(data=filter(error, metric == "MOM"), aes(x=calibration, y=Validation.deepm1, colour=calibration), pch=1, size=2) +
  scale_x_discrete(limits = c("reference", "PEST_N", "PEST_exm_w1", "PEST_exm_w2", "PEST_exm_w3"), labels = c("Reference", "Naive", "Model w1", "Model w2", "Model w3"))+
  scale_color_manual(name = "Calibration", 
                     labels = c("Reference", "Naive", "Model w1", "Model w2", "Model w3"),
                     values = c("darkorchid2", "darkorchid2", "darkorchid2", "darkorchid2", "black"))+
  ylab("MEF")+
  xlab("Model")+ 
  ylim(c(-0.5, 1))+
  theme(
    #legend.background = element_rect(fill="white"),
    #legend.spacing = unit(-6, "pt"),
    plot.title = ggplot2::element_text(face= "bold", size = 12),
    axis.title.y = ggplot2::element_text(face="bold",size= 10),
    axis.title.x = ggplot2::element_text(face="bold", size= 10, margin = margin(t = 5, r = 0, b = 0, l = 0)),
    legend.text = ggplot2::element_text(size= 8), 
    legend.title = ggplot2::element_blank(),
    axis.text.x = element_text(size=8),
    axis.text.y = element_text(size=9)
    #legend.position = c(0.5, 0.92),
    #legend.direction="horizontal",
    #legend.box = "horizontal",
    #legend.key.size = unit(3, 'mm'),
    #legend.key.height = unit(2, "mm")#change legend key size
  )

MOM

AF<- ggplot(data=filter(error, metric == "A"), aes(x=calibration, y=Calibration.deepm1, colour = calibration))+
  geom_rect(ymin = 0.5, ymax = 1.5, 
            xmin = -Inf, xmax = Inf, fill = '#C3D7A4', colour='NA', alpha=0.2) +
  geom_rect(ymin = 0, ymax = 0.5, 
            xmin = -Inf, xmax = Inf, fill = '#F4EDCA', colour= 'NA', alpha=0.2) + 
  geom_rect(ymin = -3.5, ymax = 0, 
            xmin = -Inf, xmax = Inf, fill = 'salmon', colour='NA', alpha=0.1) +
  geom_point(pch=16, size=2)+  
  ggtitle("Number of anoxic layers")+
  geom_point(data=filter(error, metric == "A"), aes(x=calibration, y=Calibration.deepm2), colour="blue", pch=15, size=2)+
  geom_point(data=filter(error, metric == "A"), aes(x=calibration, y=Validation.deepm2), colour="blue", pch=0, size=2)+
  geom_point(data=filter(error, metric == "A"), aes(x=calibration, y=Validation.deepm1, colour=calibration), pch=1, size=2) +
  scale_x_discrete(limits = c("reference", "PEST_N", "PEST_exm_w1", "PEST_exm_w2", "PEST_exm_w3"), labels = c("Reference", "Naive", "Model w1", "Model w2", "Model w3"))+
  scale_color_manual(name = "Calibration", 
                     labels = c("Reference", "Naive", "Model w1", "Model w2", "Model w3"),
                     values = c("darkorchid2", "darkorchid2", "darkorchid2", "darkorchid2", "black"))+
  ylab("MEF")+
  xlab("Model")+ 
  ylim(c(-3.5, 1))+
  theme(
    #legend.background = element_rect(fill="white"),
    #legend.spacing = unit(-6, "pt"),
    plot.title = ggplot2::element_text(face= "bold", size = 12),
    axis.title.y = ggplot2::element_text(face="bold",size= 10),
    axis.title.x = ggplot2::element_text(face="bold", size= 10, margin = margin(t = 5, r = 0, b = 0, l = 0)),
    legend.text = ggplot2::element_text(size= 8), 
    legend.title = ggplot2::element_blank(),
    axis.text.x = element_text(size=8),
    axis.text.y = element_text(size=9)
    #legend.position = c(0.5, 0.92),
    #legend.direction="horizontal",
    #legend.box = "horizontal",
    #legend.key.size = unit(3, 'mm'),
    #legend.key.height = unit(2, "mm")#change legend key size
  )

AF

temp <- ggplot(data=filter(error, metric == "temp"), aes(x=calibration, y=Calibration.deepm1, colour = calibration))+
  geom_rect(ymin = 0.5, ymax = 1.5, 
            xmin = -Inf, xmax = Inf, fill = '#C3D7A4', colour='NA', alpha=0.2) +
  geom_rect(ymin = 0, ymax = 0.5, 
            xmin = -Inf, xmax = Inf, fill = '#F4EDCA', colour= 'NA', alpha=0.2) + 
  geom_rect(ymin = -1, ymax = 0, 
            xmin = -Inf, xmax = Inf, fill = 'salmon', colour='NA', alpha=0.1) +
  geom_point(pch=16, size=2)+  
  ggtitle("Temperature")+
  geom_point(data=filter(error, metric == "temp"), aes(x=calibration, y=Calibration.deepm2), colour="blue", pch=15, size=2)+
  geom_point(data=filter(error, metric == "temp"), aes(x=calibration, y=Validation.deepm2), colour="blue", pch=0, size=2)+
  geom_point(data=filter(error, metric == "temp"), aes(x=calibration, y=Validation.deepm1, colour=calibration), pch=1, size=2) +
  scale_x_discrete(limits = c("reference", "PEST_N", "PEST_exm_w1", "PEST_exm_w2", "PEST_exm_w3"), labels = c("Reference", "Naive", "Model w1", "Model w2", "Model w3"))+
  scale_color_manual(name = "Calibration", 
                     labels = c("Reference", "Naive", "Model w1", "Model w2", "Model w3"),
                     values = c("darkorchid2", "darkorchid2", "darkorchid2", "darkorchid2", "black"))+
  ylab("MEF")+
  xlab("Model")+ 
  ylim(c(-0.5, 1))+
  theme(
    #legend.background = element_rect(fill="white"),
    #legend.spacing = unit(-6, "pt"),
    plot.title = ggplot2::element_text(face= "bold", size = 12),
    axis.title.y = ggplot2::element_text(face="bold",size= 10),
    axis.title.x = ggplot2::element_text(face="bold", size= 10, margin = margin(t = 5, r = 0, b = 0, l = 0)),
    legend.text = ggplot2::element_text(size= 8), 
    legend.title = ggplot2::element_blank(),
    axis.text.x = element_text(size=8),
    axis.text.y = element_text(size=9)
    #legend.position = c(0.5, 0.92),
    #legend.direction="horizontal",
    #legend.box = "horizontal",
    #legend.key.size = unit(3, 'mm'),
    #legend.key.height = unit(2, "mm")#change legend key size
  )

temp

oxy <- ggplot(data=filter(error, metric == "oxy"), aes(x=calibration, y=Calibration.deepm1, colour = calibration))+
  geom_rect(ymin = 0.5, ymax = 1.5, 
            xmin = -Inf, xmax = Inf, fill = '#C3D7A4', colour='NA', alpha=0.2) +
  geom_rect(ymin = 0, ymax = 0.5, 
            xmin = -Inf, xmax = Inf, fill = '#F4EDCA', colour= 'NA', alpha=0.2) + 
  geom_rect(ymin = -1, ymax = 0, 
            xmin = -Inf, xmax = Inf, fill = 'salmon', colour='NA', alpha=0.1) +
  geom_point(pch=16, size=2)+  
  ggtitle("Oxygen")+
  geom_point(data=filter(error, metric == "oxy"), aes(x=calibration, y=Calibration.deepm2), colour="blue", pch=15, size=2)+
  geom_point(data=filter(error, metric == "oxy"), aes(x=calibration, y=Validation.deepm2), colour="blue", pch=0, size=2)+
  geom_point(data=filter(error, metric == "oxy"), aes(x=calibration, y=Validation.deepm1, colour=calibration), pch=1, size=2) +
  scale_x_discrete(limits = c("reference", "PEST_N", "PEST_exm_w1", "PEST_exm_w2", "PEST_exm_w3"), labels = c("Reference", "Naive", "Model w1", "Model w2", "Model w3"))+
  scale_color_manual(name = "Calibration", 
                     labels = c("Reference", "Naive", "Model w1", "Model w2", "Model w3"),
                     values = c("darkorchid2", "darkorchid2", "darkorchid2", "darkorchid2", "black"))+
  ylab("MEF")+
  xlab("Model")+ 
  ylim(c(-0.5, 1))+
  theme(
    #legend.background = element_rect(fill="white"),
    #legend.spacing = unit(-6, "pt"),
    plot.title = ggplot2::element_text(face= "bold", size = 12),
    axis.title.y = ggplot2::element_text(face="bold",size= 10),
    axis.title.x = ggplot2::element_text(face="bold", size= 10, margin = margin(t = 5, r = 0, b = 0, l = 0)),
    legend.text = ggplot2::element_text(size= 8), 
    legend.title = ggplot2::element_blank(),
    axis.text.x = element_text(size=8),
    axis.text.y = element_text(size=9)
    #legend.position = c(0.5, 0.92),
    #legend.direction="horizontal",
    #legend.box = "horizontal",
    #legend.key.size = unit(3, 'mm'),
    #legend.key.height = unit(2, "mm")#change legend key size
  )

oxy

#Legend
legend<- ggplot(data=filter(error, metric == "oxy"), aes(x=calibration, y=Calibration.deepm1, colour = "Calibration DM 1"), pch=16, size=4)+
  geom_point()+
  geom_point(data=filter(error, metric == "oxy"), aes(x=calibration, y=Calibration.deepm2, colour="Calibration DM 2"), pch=15, size=4)+
  geom_point(data=filter(error, metric == "oxy"), aes(x=calibration, y=Validation.deepm2, colour="Validation DM 2"), pch=0, size=4)+
  geom_point(data=filter(error, metric == "oxy"), aes(x=calibration, y=Validation.deepm1, colour="Validation DM 1"), pch=1, size=4)+
  geom_point(data=filter(error, metric == "temp"), aes(x=calibration, y=Calibration.deepm2, colour="Calibration"), pch=16, size=4)+
  geom_point(data=filter(error, metric == "temp"), aes(x=calibration, y=Calibration.deepm1, colour="Validation"), pch=1, size=4)+
  scale_colour_manual(name="Legend", 
                      values=c("Calibration"="black", "Validation"="black", "Calibration DM 1"="darkorchid2", "Validation DM 1"="darkorchid2", "Calibration DM 2" = "blue", "Validation DM 2"="blue"),
                      guide=guide_legend(nrow=1, override.aes=list(shape=c(16, 1, 16, 1, 15, 0), size=c(2, 2, 2, 2, 2, 2))))+
  theme(
    legend.text = ggplot2::element_text(size= 9),
    legend.title = ggplot2::element_blank(),
    legend.direction = "horizontal",
    legend.position = "top",
    #legend.key.size = unit(1, 'mm'),
    legend.key.height = unit(2, "mm"), 
    #legend.spacing.x = unit(0.5, 'mm'),
    #legend.spacing.y = unit(0.5, 'mm'),
    legend.key.width=unit(2,"mm")
    #legend.box.background = element_rect(size=1)
    #legend.margin = margin(t = 2, r = 40, b = 2, l = 40, unit = "mm")
  )
legend

legend_plot<- get_legend(legend)

  
Plot_1 <- ggarrange(temp, oxy, SS, AF, TD, MOM, ncol=2, nrow=3, legend.grob=legend_plot, common.legend = TRUE)+
  theme(legend.position='top', 
        legend.direction= "horizontal")
Plot_1
ggsave("Results/Figure7.png",
       plot = Plot_1,
       dpi=500,
       width = 246.2, 
       height = 210, 
       units = "mm")

