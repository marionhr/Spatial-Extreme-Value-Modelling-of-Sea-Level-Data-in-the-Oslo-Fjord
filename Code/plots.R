#!/usr/bin/env Rscript


#plot data with and without trend
plot_data <- function(data, trend_data, location){
  #input: 
  #data: data frame without trend
  #trend_data: data frame with trend
  #location: name of location
  #output:
  #two plots
  ###################################
  
  #plot the data with trend
  trend_data_sub <- subset(trend_data, hour == 0)
  Plot_Trend <- ggplot(trend_data_sub,aes(x=date_to_year(trend_data_sub),y=sealevel))+
    geom_line()+
    geom_smooth(aes(color="trend"))+
    xlab("time (year)")+
    ylab(paste0(location,": sea level height (cm)"))+
    guides(color=guide_legend(" "))+
    ggsave(paste0("/../Plots/Data/",location,"_trend.pdf"),
           path=path, width = 5, height = 3, units = c("in"))

  #plot the data without trend
  data_sub <- subset(data, hour == 0)
  Plot <- ggplot(data_sub,aes(x=date_to_year(data_sub),y=sealevel))+
    geom_line()+
    geom_smooth(aes(color="trend"))+
    xlab("time (year)")+
    ylab(paste0(location,": sea level height (cm)"))+
    guides(color=guide_legend(" "))+
    ggsave(paste0("/../Plots/Data/",location,".pdf"),
           path=path, width = 5, height = 3, units = c("in"))

  return(list("trend"=Plot_Trend, "stationary"=Plot))
}


Helgeroa_plot <- plot_data(Helgeroa_Data_full, Helgeroa_Data_trend, "Helgeroa")
Oscarsborg_plot <- plot_data(Oscarsborg_Data_full, Oscarsborg_Data_trend, "Oscarsborg")
Oslo_plot <- plot_data(Oslo_Data_full, Oslo_Data_trend, "Oslo")
Viker_plot <- plot_data(Viker_Data_full, Viker_Data_trend, "Viker")

Helgeroa_plot$trend + Helgeroa_plot$stationary +
  Oscarsborg_plot$trend + Oscarsborg_plot$stationary +
  Oslo_plot$trend + Oslo_plot$stationary +
  Viker_plot$trend + Viker_plot$stationary +
  plot_layout(ncol=2, guides = "collect") &
  theme(legend.position = "bottom") &
  ggsave(paste0("/../Plots/Data/all_locations.pdf"),
         path=path, width = 10, height = 15, units = c("in"))

##################################################################################

#Plot locations on a map
ggplot() + 
  geom_polygon(data = fhidata::norway_map_counties, mapping = aes(x = long, y = lat, group = group, fill = hole), color = "black")+ 
  scale_fill_manual(values = c("white"))+
  theme_void()+ 
  coord_quickmap() + 
  coord_cartesian(xlim=c(9,12),ylim=c(58.8,60)) +
  geom_point(aes(x=locations$long, y=locations$lat), color="red")+
  geom_text(aes(x=locations$long, y=locations$lat, label=locations$location),hjust=-.1, vjust=1.5, color="red") + 
  guides(fill=FALSE, size=FALSE)+
  ggsave(paste0(path, "/../Plots/Extra_Plots/Locations.pdf"), width = 5, height = 4, units = c("in"))

##################################################################################

#Plot yearly maximum data
Helgeroa_yearly_max_plot <- ggplot(Helgeroa_yearly_max_data) + 
  geom_line(aes(year, yearly_max))+
  xlab("time (year)")+
  ylab(paste0("Helgeroa: sea level height (cm)"))

Oscarsborg_yearly_max_plot <- ggplot(Oscarsborg_yearly_max_data) + 
  geom_line(aes(year, yearly_max))+
  xlab("time (year)")+
  ylab(paste0("Oscarsborg: sea level height (cm)"))


Oslo_yearly_max_plot <- ggplot(Oslo_yearly_max_data) + 
  geom_line(aes(year, yearly_max))+
  xlab("time (year)")+
  ylab(paste0("Oslo: sea level height (cm)"))


Viker_yearly_max_plot <- ggplot(Viker_yearly_max_data) + 
  geom_line(aes(year, yearly_max))+
  xlab("time (year)")+
  ylab(paste0("Viker: sea level height (cm)"))


Helgeroa_yearly_max_plot + Oscarsborg_yearly_max_plot +
  Oslo_yearly_max_plot + Viker_yearly_max_plot +
  plot_layout(ncol=2, guides = "collect") &
  theme(legend.position = "bottom") &
  ggsave(paste0("/../Plots/Data/yearly_max_data.pdf"),
         path=path, width = 10, height = 8, units = c("in"))
  
##################################################################################

#plot dependence between the yearly maximum data for different locations
spatial_dependence_plot <- function(yearly_max_data_1,yearly_max_data_2, loc1, loc2){
  #input:
  #yearly_max_data from two differnt locations
  #loc: names of the locations
  ##############################################
  data_1 <- yearly_max_data_1$yearly_max[is.element(yearly_max_data_1$year, yearly_max_data_2$year)]
  data_2 <- yearly_max_data_2$yearly_max[is.element(yearly_max_data_2$year, yearly_max_data_1$year)]
  plt <- ggplot()+
    geom_point(aes(data_1, data_2))+
    xlab(loc1)+
    ylab(loc2)
  return(plt)
}

p1 <- spatial_dependence_plot(Helgeroa_yearly_max_data, Oscarsborg_yearly_max_data, "Helgeroa", "Oscarsborg")
p2 <- spatial_dependence_plot(Helgeroa_yearly_max_data, Oslo_yearly_max_data, "Helgeroa", "Oslo")
p3 <- spatial_dependence_plot(Helgeroa_yearly_max_data, Viker_yearly_max_data, "Helgeroa", "Viker")

p4 <- spatial_dependence_plot(Oscarsborg_yearly_max_data, Oslo_yearly_max_data, "Oscarsborg", "Oslo")
p5 <- spatial_dependence_plot(Oscarsborg_yearly_max_data, Viker_yearly_max_data, "Oscarsborg", "Viker")

p6 <- spatial_dependence_plot(Oslo_yearly_max_data, Viker_yearly_max_data, "Oslo", "Viker")

p1 + p2 + p3 + p4 + p5 + p6 + 
  plot_layout(ncol=2, guides = "collect") &
  theme(legend.position = "bottom") &
  ggsave(paste0("/../Plots/Data/spatial_dependence_plot.pdf"),
         path=path, width = 10, height = 8, units = c("in"))

##################################################################################

#Plot all the permanent stations

result <- xmlParse(paste0(path,"/../Data/permanente_stasjoner.xml"))
xml_data <- xmlToList(result)
xml_data$stationinfo$location[["name"]]

extract_xml_data <- function(index){
  df <- data.frame(name = xml_data$stationinfo[[index]][["name"]],
                   latitude = xml_data$stationinfo[[index]][["latitude"]],
                   longitude = xml_data$stationinfo[[index]][["longitude"]])
  return(df)
}

station_info <- mapply(extract_xml_data, 1:24)
station_info <- t(station_info)
#rename locations with Norwegian letters
station_info[3,1] <- "Bodø"
station_info[6,1] <- "Heimsjø"
station_info[8,1] <- "Honningsvåg"
station_info[9,1] <- "Kabelvåg"
station_info[12,1] <- "Måløy"
station_info[14,] <- NA
station_info[17,1] <- "Rørvik"
station_info[20,1] <- "Tromsø"
station_info[22,1] <- "Vardø"
station_info[24,1] <- "Ålesund"
station_info <- as.data.frame(station_info)
station_info$latitude <- as.numeric(station_info$latitude)
station_info$longitude <- as.numeric(station_info$longitude)
station_info

#Plot all stations
hjust_perm <- rep(0.6,24)
vjust_perm <- rep(-0.5,24)
hjust_perm[23] <- -0.1
hjust_perm[21] <- -0.05

ggplot() + 
  geom_polygon(data = fhidata::norway_map_counties, mapping = aes(x = long, y = lat, group = group, fill = hole), color = "black")+ 
  scale_fill_manual(values = c("white"))+
  theme_void()+ 
  coord_quickmap() + 
  geom_point(aes(x=station_info$longitude, y=station_info$latitude), color="red")+
  geom_text(aes(x=station_info$longitude, y=station_info$latitude, label=station_info$name), 
            hjust=hjust_perm, 
            vjust=vjust_perm, 
            color="red") + 
  guides(fill=FALSE, size=FALSE)+
  ggsave(paste0(path, "/../Plots/Extra_Plots/permanent_stations.pdf"), width = 10, height = 8, units = c("in"))

##################################################################################

#plot the priors
ggplot()+
  geom_line(aes(x=seq(0,10,0.01),y=dgamma(seq(0,10,0.01),5/4,1/4)))+
  xlab(expression(tau))+
  ylab("Probability Density")+
  ggsave(paste0("/../Plots/Extra_Plots/prior_tau.pdf"), 
         path=path, width = 5, height = 3, units = c("in"))

ggplot()+
  geom_line(aes(x=seq(0,200,0.1),y=dgamma(seq(0,200,0.1),5,1/12)))+
  xlab(expression(rho))+
  ylab("Probability Density")+
  ggsave(paste0("/../Plots/Extra_Plots/prior_rho.pdf"), 
         path=path, width = 5, height = 3, units = c("in"))

ggplot()+
  geom_line(aes(x=seq(0,200,0.1),y=dgamma(seq(0,200,0.1),1.5,1/40)))+
  xlab(expression(rho))+
  ylab("Probability Density")+
  ggsave(paste0("/../Plots/Extra_Plots/different_prior_rho.pdf"), 
         path=path, width = 5, height = 3, units = c("in"))

ggplot()+
  geom_line(aes(x=-0.5+exp(seq(-10,10,0.01))/(1+exp(seq(-10,10,0.01))),y=dnorm(seq(-10,10,0.01),0,10)))+
  xlab(expression(xi))+
  ylab("Probability Density")+
  ggsave(paste0("/../Plots/Extra_Plots/prior_xi.pdf"), 
         path=path, width = 5, height = 3, units = c("in"))

##################################################################################

#plot the correlation function
nu<- 1
dist <- seq(0,200)

range <- 100
kappa <- (8*nu)^(1/2)/range
matern_correlation100 <- besselK(dist * kappa, nu) *
  (kappa * dist) ^ nu * 2^(1 - nu) / gamma(nu)

range <- 40
kappa <- (8*nu)^(1/2)/range
matern_correlation40 <- besselK(dist * kappa, nu) *
  (kappa * dist) ^ nu * 2^(1 - nu) / gamma(nu)

range <- 20
kappa <- (8*nu)^(1/2)/range
matern_correlation20 <- besselK(dist * kappa, nu) *
  (kappa * dist) ^ nu * 2^(1 - nu) / gamma(nu)

range <- 60
kappa <- (8*nu)^(1/2)/range
matern_correlation60 <- besselK(dist * kappa, nu) *
  (kappa * dist) ^ nu * 2^(1 - nu) / gamma(nu)

range <- 80
kappa <- (8*nu)^(1/2)/range
matern_correlation80 <- besselK(dist * kappa, nu) *
  (kappa * dist) ^ nu * 2^(1 - nu) / gamma(nu)


ggplot()+
  geom_line(aes(dist,matern_correlation100, colour="100"))+
  geom_line(aes(dist,matern_correlation80, colour="80"))+
  geom_line(aes(dist,matern_correlation60, colour="60"))+
  geom_line(aes(dist,matern_correlation40, colour="40"))+
  geom_line(aes(dist,matern_correlation20, colour="20"))+
  guides(color=guide_legend("Range"))+
  xlab("Distance")+
  ylab("Matern Correlation Function")+
  ggsave(paste0("/../Plots/Extra_Plots/Matern_Correlation_Function.pdf"), 
         path=path, width = 6, height = 4, units = c("in"))
