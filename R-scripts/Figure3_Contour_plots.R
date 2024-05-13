library(ncdf4)
library(glmtools)
library(reshape2)
library(dplyr)
library(tidyr)
library(stats)
library(utils)
library(ggplot2)
library(ggpubr)

#Set working directory
#setwd(".../FCR-GLM-metrics")
#setwd('/Users/Kamilla/Compile_15_11/Bubbles')
setwd('/Users/Kamilla/Compile_old_GLM1602/Bubbles')
sim_folder <- getwd()
output <- nc_open("output/output.nc")

#Modelled oxy
#sim_folder <- getwd()
#output <- nc_open("Calibrated_models/Deepm2_naive/output/output.nc")
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
time <- data.frame(seq(as.Date("2016-12-20"), as.Date("2019-12-30"), by="day"))
ID <- seq.int(1:1106)
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
output <- nc_open("Calibrated_models/Deepm2_naive/output/output.nc")
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
time <- data.frame(seq(as.Date("2016-12-20"), as.Date("2019-12-30"), by="day"))
ID <- seq.int(1:1106)
#time <- data.frame(seq(as.Date("2016-12-02"), as.Date("2019-12-31"), by="day"))
#ID <- seq.int(1:1125)
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
  

mix <- ggarrange(modelledt, modelled, observedt, observed, differencet, difference, ncol=2, nrow=3, labels=c("a)", "b)", "c)", "d)", "e)", "f)"), font.label = list(size = 12, color = "black", face= "plain"))
mix

ggsave("Results/Figure3.png",
       plot = mix,
       width = 246.2, 
       height = 210, 
       units = "mm")
