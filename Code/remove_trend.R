#!/usr/bin/env Rscript

#Oslo:

#Remove the trend by quadratic regression

date_to_year <- function(data){
  return(data$year+1/(as.numeric(strftime(paste0(data$year,"-",12,"-",31),format="%j")))*(as.numeric(strftime(paste0(data$year, "-",data$month,"-",data$day),format="%j"))+data$hour*(1/24)))
}

remove_trend <- function(data_trend){
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


Helgeroa_Data_full <- remove_trend(Helgeroa_Data_trend)
Oscarsborg_Data_full <- remove_trend(Oscarsborg_Data_trend)
Oslo_Data_full <- remove_trend(Oslo_Data_trend)
Viker_Data_full <- remove_trend(Viker_Data_trend)
