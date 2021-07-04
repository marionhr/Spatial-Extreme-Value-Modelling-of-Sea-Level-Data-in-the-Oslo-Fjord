#!/usr/bin/env Rscript

#Find the yearly maximum values:
yearly_maximum_data <- function(data){
  year <- unique(data$year)
  
  yearly_maximum_data <- 1:length(year)*NA
  for (i in 1:length(year)){
    yearly_maximum_data[i] <- max(data$sealevel[data$year==year[i]], na.rm=TRUE)
  }
  
  return(tibble(year=year, yearly_max=yearly_maximum_data))
}


Helgeroa_yearly_max_data <- yearly_maximum_data(Helgeroa_Data)
Oscarsborg_yearly_max_data <- yearly_maximum_data(Oscarsborg_Data)
Oslo_yearly_max_data <- yearly_maximum_data(Oslo_Data)
Viker_yearly_max_data <- yearly_maximum_data(Viker_Data)

ggplot()+
  geom_line(aes(x=1:length(Viker_yearly_max_data$yearly_max), y=Viker_yearly_max_data$yearly_max))

length(Helgeroa_yearly_max_data$year)
length(Oscarsborg_yearly_max_data$year)
length(Oslo_yearly_max_data$year)
length(Viker_yearly_max_data$year)
