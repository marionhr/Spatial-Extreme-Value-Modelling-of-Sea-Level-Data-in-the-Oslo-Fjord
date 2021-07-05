#!/usr/bin/env Rscript

#Load in the data from files:

#Helgeroa
Helgeroa_sealevel_1965_1991 <- read.table(paste0(path,"/../Data/Helgeroa/Helgeroa_observertvannstand_timesverdier_1965-1991.txt"), 
                                      skip = 7, 
                                      col.names = c("year","month","day","hour","min","sealevel"))
Helgeroa_sealevel_1992_2020 <- read.table(paste0(path,"/../Data/Helgeroa/Helgeroa_observertvannstand_timesverdier_1992-2020.txt"), 
                                      skip = 7, 
                                      col.names = c("year","month","day","hour","min","sealevel"))
Helgeroa_tide <- read.table(paste0(path,"/../Data/Helgeroa/Helgeroa_prediksjontidevann_timesverdier_1965-2020.txt"), 
                        skip = 7, 
                        col.names = c("year","month","day","hour","min","tide"))
Nevlunghavn_sealevel_1927_1966 <- read.table(paste0(path,"/../Data/Nevlunghavn_observertvannstand_timesverdier_1927-1966.txt"), 
                                          skip = 7, 
                                          col.names = c("year","month","day","hour","min","sealevel"))


#Add data frames together to one data frame
Helgeroa_sealevel_1965_1991 <- as_tibble(Helgeroa_sealevel_1965_1991)
Helgeroa_sealevel_1992_2020 <- as_tibble(Helgeroa_sealevel_1992_2020)
Helgeroa_tide <- as_tibble(Helgeroa_tide)
Nevlunghavn_sealevel_1927_1966 <- as_tibble(Nevlunghavn_sealevel_1927_1966)

Helgeroa_Data_trend <- full_join(Helgeroa_sealevel_1965_1991, Helgeroa_sealevel_1992_2020)
Helgeroa_Data_trend <- left_join(Helgeroa_Data_trend, Helgeroa_tide)

Helgeroa_Data_trend <- full_join(Nevlunghavn_sealevel_1927_1966, Helgeroa_Data_trend, by=c("year", "month", "day", "hour", "min"))
Helgeroa_Data_trend$sealevel.y[which(!is.na(Helgeroa_Data_trend$sealevel.x) & !is.na(Helgeroa_Data_trend$sealevel.y))] <- rep(NA,length(which(!is.na(Helgeroa_Data_trend$sealevel.x) & !is.na(Helgeroa_Data_trend$sealevel.y))))
Helgeroa_Data_trend <- Helgeroa_Data_trend %>% unite("sealevel", sealevel.x:sealevel.y, na.rm = TRUE, remove = TRUE)
Helgeroa_Data_trend$sealevel <- as.numeric(as.character(Helgeroa_Data_trend$sealevel))
Helgeroa_Data_trend <- Helgeroa_Data_trend %>% arrange(year, month, day, hour)

###################################################################################

#Oscarsborg
Oscarsborg_sealevel_1953_1991 <- read.table(paste0(path,"/../Data/Oscarsborg/Oscarsborg_observertvannstand_timesverdier_1953-1991.txt"), 
                                          skip = 7, 
                                          col.names = c("year","month","day","hour","min","sealevel"))
Oscarsborg_sealevel_1992_2020 <- read.table(paste0(path,"/../Data/Oscarsborg/Oscarsborg_observertvannstand_timesverdier_1992-2020.txt"), 
                                          skip = 7, 
                                          col.names = c("year","month","day","hour","min","sealevel"))
Oscarsborg_tide <- read.table(paste0(path,"/../Data/Oscarsborg/Oscarsborg_prediksjontidevann_timesverdier_1953-2020.txt"), 
                            skip = 7, 
                            col.names = c("year","month","day","hour","min","tide"))



#Add data frames together to one data frame
Oscarsborg_sealevel_1953_1991 <- as_tibble(Oscarsborg_sealevel_1953_1991)
Oscarsborg_sealevel_1992_2020 <- as_tibble(Oscarsborg_sealevel_1992_2020)
Oscarsborg_tide <- as_tibble(Oscarsborg_tide)

Oscarsborg_Data_trend <- full_join(Oscarsborg_sealevel_1953_1991, Oscarsborg_sealevel_1992_2020)
Oscarsborg_Data_trend <- left_join(Oscarsborg_Data_trend, Oscarsborg_tide)

###################################################################################

#Oslo
Oslo_sealevel_1914_1991 <- read.table(paste0(path,"/../Data/Oslo/Oslo_observertvannstand_timesverdier_1914-1991.txt"), 
                                      skip = 7, 
                                      col.names = c("year","month","day","hour","min","sealevel"))
Oslo_sealevel_1992_2020 <- read.table(paste0(path,"/../Data/Oslo/Oslo_observertvannstand_timesverdier_1992-2020.txt"), 
                                      skip = 7, 
                                      col.names = c("year","month","day","hour","min","sealevel"))
Oslo_tide <- read.table(paste0(path,"/../Data/Oslo/Oslo_prediksjontidevann_timesverdier_1914-2020.txt"), 
                        skip = 7, 
                        col.names = c("year","month","day","hour","min","tide"))



#Add data frames together to one data frame
Oslo_sealevel_1914_1991 <- as_tibble(Oslo_sealevel_1914_1991)
Oslo_sealevel_1992_2020 <- as_tibble(Oslo_sealevel_1992_2020)
Oslo_tide <- as_tibble(Oslo_tide)

Oslo_Data_trend <- full_join(Oslo_sealevel_1914_1991, Oslo_sealevel_1992_2020)
Oslo_Data_trend <- left_join(Oslo_Data_trend, Oslo_tide)

###################################################################################

#Viker
Viker_sealevel_1990_2020 <- read.table(paste0(path,"/../Data/Viker/Viker_observertvannstand_timesverdier_1990-2020.txt"), 
                                        skip = 7, 
                                        col.names = c("year","month","day","hour","min","sealevel"))
Viker_tide <- read.table(paste0(path,"/../Data/Viker/Viker_prediksjontidevann_timesverdier_1990-2020.txt"), 
                          skip = 7, 
                          col.names = c("year","month","day","hour","min","tide"))



#Add data frames together to one data frame
Viker_sealevel_1990_2020 <- as_tibble(Viker_sealevel_1990_2020)
Viker_tide <- as_tibble(Viker_tide)

Viker_Data_trend <- left_join(Viker_sealevel_1990_2020, Viker_tide)

###################################################################################

#Find GPS coordinates and distance matrix:

#Location:
locations <- read.table(paste0(path,"/../Data/location.txt"),
                        col.names = c("location","code","latitude","longitude"))
locations <- locations[c(1:3,5),]
locations

#Find distances:
distances <- find_distances(locations)
distances
