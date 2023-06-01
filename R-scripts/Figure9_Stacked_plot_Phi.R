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
setwd('.../FCR-GLM-metrics/Calibrated_models')
sim_folder <- getwd()

#Phi during the calibration recorded in glm3.iobj file which was converted to csv file for analysis
#Routine deep mixing 2 model
error_routine <- read.csv("Deepm2_routine/glm3.csv")
error_r_1 <- subset(error_routine, select = c(1:5))
error_r <- subset(error_routine, select = -c(1:5))
melted <- melt(error_r, na.rm=TRUE)
iteration <- rep(error_r_1$iteration, times=20)
new <- cbind(iteration, melted)
new1 <- new[order(iteration),]  
colnames(new1) <- c("Iteration", "Obs_group" ,"Phi")

New_obs_group <- c(rep("Temp", times=10), rep("Oxy", times=10))
New_obs_group1 <- rep(New_obs_group, max(error_r_1$iteration)+1)

final_r <- cbind(new1, New_obs_group1)

reorder <- c("Oxy", "Temp")
new_data_r <- aggregate(Phi ~  factor(New_obs_group1, levels=reorder) + Iteration, data=final_r, sum)
colnames(new_data_r) <- c("New_obs_group", "Iteration", "Phi")
r<- ggplot(new_data_r, aes(x=Iteration, y=Phi, fill=New_obs_group)) + 
  geom_area(size=0.1, colour="black")+
  ggtitle("Deep mixing 2 routine") +
  ylim(c(0, 4000))+
  scale_fill_manual(
    name="Legend",
    values = c("Temp"="#619CFF", "Oxy"="#00BA38")
  )+
  ggplot2::theme_light() +
  theme(plot.title = element_text(hjust = 0, size=12, face="bold"),
        axis.title.y = element_text(size=10, face="bold") ,
        axis.title.x = element_text(size=10, face="bold"),
        legend.text = ggplot2::element_text(size= 10),
        plot.subtitle=element_text(size=8),
        legend.title = ggplot2::element_blank(),
        legend.position = "top",
        legend.direction="horizontal",
        legend.key.height = unit(3, "mm"),
        legend.key.size = unit(3, 'mm'),
        legend.key.width=unit(3,"mm"))
r

#Weight 1, deep mixing 2 model
error_w1 <- read.csv("Deepm2_exm_weight1/glm3.csv")
error_w1_1 <- subset(error_w1, select = c(1:5))
error_w1 <- subset(error_w1, select = -c(1:5))

melted <- melt(error_w1, na.rm=TRUE)

iteration <- rep(error_w1_1$iteration, times=24)
new <- cbind(iteration, melted)
new1 <- new[order(iteration),]  
colnames(new1) <- c("Iteration", "Obs_group" ,"Phi")

New_obs_group <- c(rep("Temp", times=10), rep("Oxy", times=10), "TD", "SS", "MOM", "AF")
New_obs_group1 <- rep(New_obs_group, max(error_w1_1$iteration)+1)
final_w1 <- cbind(new1, New_obs_group1)

reorder <- c("Oxy", "Temp", "AF", "MOM", "SS", "TD")
new_data_w1 <- aggregate(Phi ~  factor(New_obs_group1, levels=reorder) + Iteration, data=final_w1, sum)
colnames(new_data_w1) <- c("New_obs_group", "Iteration", "Phi")
w1<- ggplot(new_data_w1, aes(x=Iteration, y=Phi, fill=New_obs_group)) + 
  geom_area(size=0.1, colour="black")+
  ggtitle("Deep mixing 2 w1") +
  ylim(c(0, 4000))+
  scale_fill_manual(
    name="Legend",
    values = c("Temp"="#619CFF", "Oxy"="#00BA38", "TD"="#FF61CC", "SS"="#00BFC4", "MOM"="#CD9600", "AF"="#F8766D"),
    guide=guide_legend(nrow=1)
  )+
  ggplot2::theme_light() +
  theme(plot.title = element_text(hjust = 0, size=12, face="bold"),
        axis.title.y = element_text(size=10, face="bold") ,
        axis.title.x = element_text(size=10, face="bold"),
        legend.text = ggplot2::element_text(size= 10),
        plot.subtitle=element_text(size=8),
        legend.title = ggplot2::element_blank(),
        legend.position = "top",
        legend.direction="horizontal",
        legend.key.height = unit(3, "mm"),
        legend.key.size = unit(3, 'mm'),
        legend.key.width=unit(3,"mm"))
w1

w1_legend <- get_legend(w1)

#Weight 2, deep mixing 2 model

error_w2 <- read.csv("Deepm2_exm_weight2/glm3.csv")
error_w2_1 <- subset(error_w2, select = c(1:5))
error_w2 <- subset(error_w2, select = -c(1:5))

melted <- melt(error_w2, na.rm=TRUE)

iteration <- rep(error_w2_1$iteration, times=24)
new <- cbind(iteration, melted)
new1 <- new[order(iteration),]  
colnames(new1) <- c("Iteration", "Obs_group" ,"Phi")

New_obs_group <- c(rep("Temp", times=10), rep("Oxy", times=10), "TD", "SS", "MOM", "AF")
New_obs_group1 <- rep(New_obs_group, max(error_w2_1$iteration)+1)

final_w2 <- cbind(new1, New_obs_group1)

ggplot(final_w2, aes(x=Iteration, y=Phi, fill=Obs_group)) + 
  geom_area(size=0.1, colour="black")+
  ggtitle("Deep mixing 2 w2")

reorder <- c("Oxy", "Temp", "AF", "MOM", "SS", "TD")
new_data_w2 <- aggregate(Phi ~  factor(New_obs_group1, levels=reorder) + Iteration, data=final_w2, sum)
colnames(new_data_w2) <- c("New_obs_group", "Iteration", "Phi")
w2<- ggplot(new_data_w2, aes(x=Iteration, y=Phi, fill=New_obs_group)) + 
  geom_area(size=0.1, colour="black")+
  ggtitle("Deep mixing 2 w2") +
  ylim(c(0, 4000))+
  scale_fill_manual(
    name="Legend",
    values = c("Temp"="#619CFF", "Oxy"="#00BA38", "TD"="#FF61CC", "SS"="#00BFC4", "MOM"="#CD9600", "AF"="#F8766D"),
    guide=guide_legend(nrow=1)
  )+
  ggplot2::theme_light() +
  theme(plot.title = element_text(hjust = 0, size=12, face="bold"),
        axis.title.y = element_text(size=10, face="bold") ,
        axis.title.x = element_text(size=10, face="bold"),
        legend.text = ggplot2::element_text(size= 10),
        plot.subtitle=element_text(size=8),
        legend.title = ggplot2::element_blank(),
        legend.position = "top",
        legend.direction="horizontal",
        legend.key.height = unit(3, "mm"),
        legend.key.size = unit(3, 'mm'),
        legend.key.width=unit(3,"mm"))
w2

#w3
error_w3 <- read.csv("Deepm2_exm_weight3/glm3.csv")

error_w3_1 <- subset(error_w3, select = c(1:5))
error_w3 <- subset(error_w3, select = -c(1:5))

melted <- melt(error_w3, na.rm=TRUE)

iteration <- rep(error_w3_1$iteration, times=24)
new <- cbind(iteration, melted)
new1 <- new[order(iteration),]  
colnames(new1) <- c("Iteration", "Obs_group" ,"Phi")

New_obs_group <- c(rep("Temp", times=10), rep("Oxy", times=10), "TD", "SS", "MOM", "AF")
New_obs_group1 <- rep(New_obs_group, max(error_w3_1$iteration)+1)

final_w3 <- cbind(new1, New_obs_group1)

reorder <- c("Oxy", "Temp", "AF", "MOM", "SS", "TD")
new_data_w3 <- aggregate(Phi ~  factor(New_obs_group1, levels=reorder) + Iteration, data=final_w3, sum)
colnames(new_data_w3) <- c("New_obs_group", "Iteration", "Phi")

w3<- ggplot(new_data_w3, aes(x=Iteration, y=Phi, fill=New_obs_group)) + 
  geom_area(size=0.1, colour="black")+
  ggtitle("Deep mixing 2 w3") +
  #ylim(c(0, 9000))+
  scale_fill_manual(
    name="Legend",
    values = c("Temp"="#619CFF", "Oxy"="#00BA38", "TD"="#FF61CC", "SS"="#00BFC4", "MOM"="#CD9600", "AF"="#F8766D"),
    guide=guide_legend(nrow=1)
  )+
  ggplot2::theme_light() +
  theme(plot.title = element_text(hjust = 0, size=12, face="bold"),
        axis.title.y = element_text(size=10, face="bold") ,
        axis.title.x = element_text(size=10, face="bold"),
        legend.text = ggplot2::element_text(size= 10),
        plot.subtitle=element_text(size=8),
        legend.title = ggplot2::element_blank(),
        legend.position = "top",
        legend.direction="horizontal",
        legend.key.height = unit(3, "mm"),
        legend.key.size = unit(3, 'mm'),
        legend.key.width=unit(3,"mm"))
w3

stacked<- ggarrange(r, w1, w2, w3, ncol=2, nrow=2, common.legend=TRUE, legend.grob=w1_legend)
stacked

setwd('.../FCR-GLM-metrics')

ggsave("Results/Figure9.png",
       plot = stacked,
       width = 246.2, 
       height = 160, 
       #scale=1.6,
       dpi=500,
       units = "mm")
