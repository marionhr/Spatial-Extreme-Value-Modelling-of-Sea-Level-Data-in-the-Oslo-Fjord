#!/usr/bin/env Rscript

#Find the yearly maximum values:

yearly_maximum_data <- function(data){
  #input: data without a trend and with the years that contain missing data removed
  #output: a data frame with the year and the yearly maximum as two columns
  #################################################################################
  
  year <- unique(data$year)
  
  yearly_maximum_data <- 1:length(year)*NA
  for (i in 1:length(year)){
    yearly_maximum_data[i] <- max(data$sealevel[data$year==year[i]], na.rm=TRUE)
  }
  
  return(tibble(year=year, yearly_max=yearly_maximum_data))
}

#Find the yearly maximum values for the 4 locations
Helgeroa_yearly_max_data <- yearly_maximum_data(Helgeroa_Data)
Oscarsborg_yearly_max_data <- yearly_maximum_data(Oscarsborg_Data)
Oslo_yearly_max_data <- yearly_maximum_data(Oslo_Data)
Viker_yearly_max_data <- yearly_maximum_data(Viker_Data)
