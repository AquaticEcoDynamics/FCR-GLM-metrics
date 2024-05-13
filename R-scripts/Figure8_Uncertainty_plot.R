library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)

setwd(".../FCR-GLM-metrics/Uncertainty_analysis/Run2")

#Uncertainty_results
obs_0<- read.csv("glm3_reweight_ies.0.obs.csv")
obs_3<- read.csv("glm3_reweight_ies.3.obs.csv")

#MOM
obs_mom <- read.csv("model_files/field_data/mom_observed.csv") %>% 
  mutate(DateTime = as.Date(DateTime, format="%Y-%m-%d"))%>%
  filter(DateTime > "2016-12-01") 

#MOM prior 
obs_mom_prior <- obs_0 %>% 
  dplyr:: select(grep("mom", names(obs_0))) %>%
  mutate(realizations = obs_0$real_name) %>%
  select(realizations, everything())

obs_mom_prior_melt <- reshape2::melt(obs_mom_prior) 

obs_mom_prior_melt$realizations <- as.numeric(obs_mom_prior_melt$realizations)
arrange <- arrange(obs_mom_prior_melt, realizations) %>%
  mutate(DateTime = rep(obs_mom$DateTime, times=276)) %>%
  mutate(DateTime = as.Date(DateTime, format="%Y-%m-%d"))

#MOM posterior 
obs_mom_post <- obs_3 %>% 
  dplyr:: select(grep("mom", names(obs_3))) %>%
  mutate(realizations = obs_3$real_name) %>%
  select(realizations, everything())

obs_mom_post_melt <- reshape2::melt(obs_mom_post) 

obs_mom_post_melt$realizations <- as.numeric(obs_mom_post_melt$realizations)
arrange_post <- arrange(obs_mom_post_melt, realizations) %>%
  mutate(DateTime = rep(obs_mom$DateTime, times=238)) %>%
  mutate(DateTime = as.Date(DateTime, format="%Y-%m-%d"))

#MOM plot
intervals = 1:19/20
plot1<- ggplot2::ggplot(data=arrange, aes(x=DateTime, y=value))+
  ggfan::geom_fan(intervals=intervals)+
  scale_fill_gradient(name="Prior")+
  ggnewscale::new_scale_fill()+
  ggfan::geom_fan(data=arrange_post, mapping=aes(x=DateTime, y=value), intervals=intervals, alpha=0.6)+
  scale_fill_gradient(name="Posterior", low="red", high="pink") +
  geom_point(data=obs_mom, mapping=aes(x=DateTime, y=deviation, colour="Observations"), pch=16) +
  #geom_line(data=means_by_variable, mapping=aes(x=DateTime, y=mean_value, colour="black"))+
  #ggtitle("Prior vs. Posterior")+
  ylab(expression(bold(MOM~(mmol/m^{3}))))+
  xlab("Date")+
  scale_colour_manual(
    values = c("Observations"="black"), name=""
  ) +
  ggplot2::theme_light() +
  theme(
    plot.title = ggplot2::element_text(face= "bold", size = 10),
    axis.title.y = ggplot2::element_text(face="bold",size= 12),
    axis.title.x = ggplot2::element_text(face="bold", size= 12, margin = margin(t = 5, r = 0, b = 0, l = 0)),
    legend.text = ggplot2::element_text(size= 10), 
    legend.title = ggplot2::element_text(face="bold",size= 12),
    axis.text = element_text(size=12)
  )

plot1

#TD

obs_td <- read.csv("model_files/field_data/obs_td.csv") %>%
  mutate(DateTime = as.Date(DateTime, format="%Y-%m-%d")) %>%
  filter(DateTime > "2016-12-01" & DateTime < "2020-01-01") %>%
  mutate(month = lubridate::month(DateTime)) %>%
  filter(between(month, 4, 9)) %>%
  mutate(year = lubridate::year(DateTime)) %>%
  group_by(year) %>%
  filter(year == 2017 | year== 2018 | year==2019)

#TD prior
obs_td_prior <- obs_0 %>% 
  dplyr:: select(grep("td", names(obs_0))) %>%
  mutate(realizations = obs_0$real_name) %>%
  select(realizations, everything())

obs_td_prior_melt <- reshape2::melt(obs_td_prior) 

obs_td_prior_melt$realizations <- as.numeric(obs_td_prior_melt$realizations)
arrange_td_prior <- arrange(obs_td_prior_melt, realizations) %>%
  mutate(DateTime = rep(obs_td$DateTime, times=276)) %>%
  mutate(DateTime = as.Date(DateTime, format="%Y-%m-%d")) %>%
  mutate(month = lubridate::month(DateTime)) %>%
  filter(between(month, 4, 9)) %>%
  mutate(year = lubridate::year(DateTime)) %>%
  group_by(year) %>%
  mutate(col= 1)

#TD posterior 
obs_td_post <- obs_3 %>% 
  dplyr:: select(grep("td", names(obs_3))) %>%
  mutate(realizations = obs_3$real_name) %>%
  select(realizations, everything())

obs_td_post_melt <- reshape2::melt(obs_td_post) 

obs_td_post_melt$realizations <- as.numeric(obs_td_post_melt$realizations)
arrange_td_post <- arrange(obs_td_post_melt, realizations) %>%
  mutate(DateTime = rep(obs_td$DateTime, times=238)) %>%
  mutate(DateTime = as.Date(DateTime, format="%Y-%m-%d")) %>%
  mutate(month = lubridate::month(DateTime)) %>%
  filter(between(month, 4, 9)) %>%
  mutate(year = lubridate::year(DateTime)) %>%
  group_by(year) %>%
  mutate(col= 1)

#TD plot
intervals = 1:19/20
plot_td_1 <- ggplot2::ggplot(data=arrange_td_prior, aes(x=DateTime, y=value, group=year))+
  facet_wrap(~ year, nrow = 1, scales = "free_x")+
  ggfan::geom_fan(intervals=intervals)+
  scale_fill_gradient(name="Prior")+
  ggnewscale::new_scale_fill()+
  ggfan::geom_fan(data=arrange_td_post, mapping=aes(x=DateTime, y=value, group=year), intervals=intervals, alpha=0.6)+
  scale_fill_gradient(name="Posterior", low="red", high="pink") +
  geom_point(data=obs_td, mapping=aes(x=DateTime, y=thermo.depth, group=year, colour="Observations"), pch=16) +
  ylab("Thermocline depth (m)")+
  #ggtitle("Ensemble 3")+
  xlab("Date")+
  scale_y_reverse(limits=c(10, 0))+
  scale_colour_manual(
    values = c("Observations"="black"), name= substitute(paste(bold("Legend")))
  ) +
  guides(colour = guide_legend(order = 1))+
  ggplot2::theme_light() +
  theme(
    plot.title = ggplot2::element_text(face= "bold", size = 10),
    axis.title.y = ggplot2::element_text(face="bold",size= 12),
    axis.title.x = ggplot2::element_text(face="bold", size= 12, margin = margin(t = 5, r = 0, b = 0, l = 0)),
    legend.text = ggplot2::element_text(size= 10), 
    legend.title = ggplot2::element_text(size= 11),
    axis.text = element_text(size=12),
    panel.spacing = unit(1, "lines"),
    legend.spacing = unit(-4, "pt"),
    legend.key.size = unit(6, 'mm'),
    legend.position = "right",
    strip.text = element_text(size = 16)
  )
plot_td_1

#Schmidt Stability
#SS observations
obs_SS <- read.csv("model_files/field_data/Obs_SS.csv") %>% 
  mutate(datetime = as.Date(datetime, format="%Y-%m-%d"))%>%
  filter(datetime > "2016-12-01" & datetime < "2020-01-01") 

#SS prior 
obs_ss_prior <- obs_0 %>% 
  dplyr:: select(grep("ss", names(obs_0))) %>%
  mutate(realizations = obs_0$real_name) %>%
  select(realizations, everything())

obs_ss_prior_melt <- reshape2::melt(obs_ss_prior)

obs_ss_prior_melt$realizations <- as.numeric(obs_ss_prior_melt$realizations)
arrange_ss <- arrange(obs_ss_prior_melt, realizations) %>%
  mutate(DateTime = rep(obs_SS$datetime, times=276)) %>%
  mutate(DateTime = as.Date(DateTime, format="%Y-%m-%d"))

#SS posterior
obs_ss_post <- obs_3 %>% 
  dplyr:: select(grep("ss", names(obs_3))) %>%
  mutate(realizations = obs_3$real_name) %>%
  select(realizations, everything())

obs_ss_post_melt <- reshape2::melt(obs_ss_post)

obs_ss_post_melt$realizations <- as.numeric(obs_ss_post_melt$realizations)
arrange_post_ss <- arrange(obs_ss_post_melt, realizations) %>%
  mutate(DateTime = rep(obs_SS$datetime, times=238)) %>%
  mutate(DateTime = as.Date(DateTime, format="%Y-%m-%d"))

intervals = 1:19/20
plot_ss<- ggplot2::ggplot(data=arrange_ss, aes(x=DateTime, y=value))+
  ggfan::geom_fan(intervals=intervals)+
  scale_fill_gradient(name="Prior")+
  ggnewscale::new_scale_fill()+
  ggfan::geom_fan(data=arrange_post_ss, mapping=aes(x=DateTime, y=value), intervals=intervals, alpha=0.6)+
  scale_fill_gradient(name="Posterior", low="red", high="pink") +
  geom_point(data=obs_SS, mapping=aes(x=datetime, y=schmidt.stability, colour="Observations"), pch=16) +
  #geom_line(data=means_by_variable, mapping=aes(x=DateTime, y=mean_value, colour="black"))+
  #ggtitle("Prior vs. Posterior")+
  ylab(expression(bold(Schmidt~stability)))+
  xlab("Date")+
  scale_colour_manual(
    values = c("Observations"="black"), name=""
  ) +
  ggplot2::theme_light() +
  theme(
    plot.title = ggplot2::element_text(face= "bold", size = 10),
    axis.title.y = ggplot2::element_text(face="bold",size= 12),
    axis.title.x = ggplot2::element_text(face="bold", size= 12, margin = margin(t = 5, r = 0, b = 0, l = 0)),
    legend.text = ggplot2::element_text(size= 10), 
    legend.title = ggplot2::element_text(face="bold",size= 12),
    axis.text = element_text(size=12)
  )

plot_ss

combined_mom_ss_td <- ggpubr::ggarrange(plot_ss, plot_td_1, plot1, ncol=1, nrow=3, common.legend=TRUE, legend="right", labels=c("a)", "b)", "c)"), font.label = list(color = "black", face= "plain"))
combined_mom_ss_td

ggsave("Results/Figure8.png",
       plot = combined_mom_ss_td,
       width = 300, 
       height = 300,
       dpi=500,
       units = "mm")