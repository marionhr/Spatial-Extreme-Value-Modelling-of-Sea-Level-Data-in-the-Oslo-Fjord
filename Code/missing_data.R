#!/usr/bin/env Rscript

#Remove the years that are missing more than 3 months of data

remove_years_with_missing_data <- function(data){
  #input: data frame with the trend removed
  #output: data frame without the years with missing data removed
  ################################################################
  
  missing_data <- tibble(year=numeric(),
                         month=numeric(),
                         day=numeric(),
                         number_days=numeric(),
                         number_hours=numeric())
  
  year_rle <- rle(data$year)
  for (year in unique(data$year)){
    if (leap_year(year)){
      days_year <- 366
    }
    else{
      days_year <- 365
    }
    month_rle <- rle(data$month[data$year==year])
    
    #Check if there is any missing data in the year:
    if (!(year_rle$lengths[year_rle$values==year] == 24*days_year)){
      
      for (month in 1:12){
        days_month <- days_in_month(as.Date(paste0(year,"-",month,"-",1)))
        day_rle <- rle(data$day[data$year==year & data$month==month])
        
        #Check if there is any missing data in the month:
        if (!any(month_rle$values==month)){
          missing_data <- missing_data %>% add_row(year=year, month=month, number_days=days_month)
        } 
        else if (!(month_rle$lengths[month_rle$values==month] == 24*days_month)){
          for (day in 1:days_month){
            #Check if there is any missing data in the day:
            if (!any(day_rle$values==day)){
              missing_data <- missing_data %>% add_row(year=year, month=month, day=day, number_hours=24)
            }
            else if (!(day_rle$lengths[day_rle$values==day] == 24)){
              missing_data <- missing_data %>% add_row(year=year, month=month, day=day, number_hours=day_rle$lengths[day_rle$values==day])
            }
          }
        }
      }
    }
  }
  
  year_list <- numeric()
  for (year in unique(missing_data$year)){
    if (leap_year(year)){
      days_year <- 366
    }
    else{
      days_year <- 365
    }
    
    sum_days <- sum(missing_data$number_days[missing_data$year==year], na.rm = TRUE)
    sum_hours <- sum(missing_data$number_hours[missing_data$year==year], na.rm = TRUE)
    
    proportion_missing <- (sum_days*24+sum_hours)/(24*days_year)
    
    if (proportion_missing > (3/12)){
      year_list <- c(year_list, year)
    }
  }
  
  print(year_list)
  
  data <- data[!is.element(data$year,year_list),]
  
  print(unique(data$year))
  return(data)
}

#Find data frames without missing data:
Helgeroa_Data <- remove_years_with_missing_data(Helgeroa_Data_full)
Oscarsborg_Data <- remove_years_with_missing_data(Oscarsborg_Data_full)
Oslo_Data <- remove_years_with_missing_data(Oslo_Data_full)
Viker_Data <- remove_years_with_missing_data(Viker_Data_full)