#!/usr/bin/env Rscript

#Oslo:

#Remove the trend by quadratic regression

remove_trend <- function(data_trend){
  #input: data frame with columns: year, month, day, hour, min, sea level and tide
  #output: the same data frame with the trend in the sea level column removed by quadratic regression
  #####################################################################################################
  df <- data.frame(x = date_to_year(data_trend),
                   y = data_trend$sealevel)
  df$x <- df$x-min(df$x)
  
  reg <- glm(y ~ poly(x, degree=2, raw=TRUE), data = df)
  print(summary(reg))
  
  data <- data_trend
  data$sealevel <- (data$sealevel - 
                           reg$coefficients[1] - 
                           reg$coefficients[2]*df$x - 
                           reg$coefficients[3]*(df$x)^2)
  data$surge <- data$sealevel - data$tide
  
  return(data)
}

#Remove trend at the four locations
Helgeroa_Data_full <- remove_trend(Helgeroa_Data_trend)
Oscarsborg_Data_full <- remove_trend(Oscarsborg_Data_trend)
Oslo_Data_full <- remove_trend(Oslo_Data_trend)
Viker_Data_full <- remove_trend(Viker_Data_trend)
