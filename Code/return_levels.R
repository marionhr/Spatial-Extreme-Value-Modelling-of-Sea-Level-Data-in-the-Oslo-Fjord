#!/usr/bin/env Rscript

set.seed(10)

#functions:
apply_return_level_function_to_parameter_samples <- function(index_samples, return_period, location_index, estimates){
  #input: 
  #index_samples: which sample
  #return_period: which return period
  #location_index: which location
  #estimates: list of data and parameters
  #output:
  #ret: return level value for one sample and one return period
  #############################################################
  ret <- return_levels(return_period, get(paste0("q.",location_index),estimates$data)[index_samples], get(paste0("lns.",location_index),estimates$data)[index_samples], estimates$data$xi[index_samples], estimates$parameters$alpha, estimates$parameters$beta)
  return(ret)
}

apply_return_periods <- function(return_period, estimates, location_index){
  #input: 
  #return_period: which return period
  #estimates: list of data and parameters
  #location_index: which location
  #output:
  #ret: return level values for all the samples and one return period
  ####################################################################
  ret <- lapply(1:nrow(estimates$data), apply_return_level_function_to_parameter_samples, return_period, location_index, estimates)
  ret <- as.data.frame(ret)
  names(ret) <- 1:nrow(estimates$data)
  ret <- t(ret)
  return(ret)
}

plot_return_levels <- function(median_and_confidence_intervals, return_periods, location, directory, plot_lim, estimates, location_index){
  #input: 
  #median_and_confidence_intervals: matrix containing medians and confidence bounds for the return periods
  #return_periods: list of return periods
  #location: name of location to plot
  #directory: where to save results
  #plot_lim: limits for the plots
  #estimates: list of parameter estimates and parameter values
  #location_index: which location to plot
  #output:
  #plt: plot of return levels for one location
  #######################################################
  plot_lim <- c(min(plot_lim[1],median_and_confidence_intervals[1,3]),max(plot_lim[2],median_and_confidence_intervals[3,]))
  if(location_index!=4){
    empirical_data <- estimates$parameters$x[estimates$parameters$location_indices==location_index]
  }else{
    empirical_data <- plot_lim[1]-10
  }
  
  plt <- ggplot()+
    geom_point(aes(x=length(empirical_data)/1:length(empirical_data), 
                   y=sort(empirical_data,decreasing = TRUE), 
                   colour="empirical \nreturn levels"))+
    geom_line(aes(x=return_periods, y=median_and_confidence_intervals[2,], colour="median\n"))+
    geom_line(aes(x=return_periods, y=median_and_confidence_intervals[1,], colour="95% credible \ninterval\n"),linetype = "dashed")+
    geom_line(aes(x=return_periods, y=median_and_confidence_intervals[3,], colour="95% credible \ninterval\n"),linetype = "dashed")+
    xlab("Return period (years)")+
    ylab(paste0(location,": Return level (cm)"))+
    ylim(plot_lim)

  return(plt)
}


make_return_level_plot <- function(return_periods, estimates, location_index, locations, directory, plot_lim){
  #input:
  #return_periods: list of return periods
  #estimates: list of parameter estimates and parameter values
  #location_index: which location
  #directory: where to save the results
  #locations: list of names for the locations
  #plot_lim: limits for the plots
  #output:
  #plot of return levels for one location
  ##############################################
  
  #find return levels for all samples and all return periods
  ret_lev <- mapply(apply_return_periods, return_periods, MoreArgs=list("estimates"=estimates, "location_index"=location_index))
  
  #find median and confidence bounds
  median_conf_int <- mapply(median_and_confidence_interval, 1:length(return_periods), MoreArgs = list("all_return_levels"=ret_lev))
  
  #plot results
  plt <- plot_return_levels(median_conf_int, return_periods, locations[location_index], directory, plot_lim, estimates, location_index)
  return(plt)
}

plot_all_return_levels <- function(return_periods, ncol_plot, nrow_plot, directory, locations, plot_lim=c(0,1000)){
  #input:
  #return_periods: list of return periods
  #ncol_plot: number of columns in the plot
  #nrow_plot: number of rows in the plot
  #directory: where to save the results
  #locations: list of names for the locations
  #plot_lim: limits for the plots
  #output:
  #plot of return levels
  ##############################################
  
  #find parameter estimates and parameter values
  estimates <- find_data(directory)
  
  #make plots
  nloc <- estimates$parameters$dim
  plots <- map(1:nloc, ~make_return_level_plot(return_periods, estimates, .x, locations, directory, plot_lim))
  ggarrange(plotlist=plots, ncol=ncol_plot, nrow = nrow_plot, common.legend = T, legend="bottom")+
    ggsave(paste0(path, "/../Plots/return_level/",directory,".pdf"),
           width = 10*ncol_plot/4, height = 12*nrow_plot/4+1/2, units = c("in"))
}

####################################################################################################
      
#return level plots
return_periods <- c(1,2,5,10,20,30,40,50,100,150,200)

#Simulated Data 4 locations
directory <- "Simulated_Data_4_locations"
ncol_plot <- 4
nrow_plot <- 1
plot_lim <- c(100,300)
plot_all_return_levels(return_periods, ncol_plot, nrow_plot, directory, 1:4, plot_lim = plot_lim)

#Simulated Data 4 locations common xi
directory <- "Simulated_Data_common_xi"
ncol_plot <- 4
nrow_plot <- 1
plot_lim <- c(100,300)
plot_all_return_levels(return_periods, ncol_plot, nrow_plot, directory, 1:4, plot_lim = plot_lim)

#Simulated Data 8 locations
directory <- "Simulated_Data_8_locations"
ncol_plot <- 4
nrow_plot <- 2
plot_lim <- c(100,300)
plot_all_return_levels(return_periods, ncol_plot, nrow_plot, directory, 1:8, plot_lim = plot_lim)

#Simulated Data 15 locations
directory <- "Simulated_Data_15_locations"
ncol_plot <- 4
nrow_plot <- 4
plot_lim <- c(100,300)
plot_all_return_levels(return_periods, ncol_plot, nrow_plot, directory, 1:15, plot_lim = plot_lim)

#Simulated Data 15 locations different prior
directory <- "Simulated_Data_15_locations_different_prior"
ncol_plot <- 4
nrow_plot <- 4
plot_lim <- c(100,300)
plot_all_return_levels(return_periods, ncol_plot, nrow_plot, directory, 1:15, plot_lim = plot_lim)


#univariate simulated data
plot_lim <- c(100,300)

#1
directory <- "Simulated_Data_univariate/Station_1"
location_index <- 1
estimates <- find_univariate_data(directory,location_index)
names_param <- names(estimates$parameters)
names_param[2] <- "x"
names(estimates$parameters) <- names_param
estimates$parameters$location_indices <- rep(1, length(estimates$parameters$x))
plt_1 <- make_return_level_plot(return_periods, estimates, location_index, "1", directory, plot_lim = plot_lim)

#2
directory <- "Simulated_Data_univariate/Station_2"
location_index <- 1
estimates <- find_univariate_data(directory,location_index)
names(estimates$parameters) <- names_param
estimates$parameters$location_indices <- rep(1, length(estimates$parameters$x))
plt_2 <- make_return_level_plot(return_periods, estimates, location_index, "2", directory, plot_lim = plot_lim)

#3
directory <- "Simulated_Data_univariate/Station_3"
location_index <- 1
estimates <- find_univariate_data(directory,location_index)
names(estimates$parameters) <- names_param
estimates$parameters$location_indices <- rep(1, length(estimates$parameters$x))
plt_3 <- make_return_level_plot(return_periods, estimates, location_index, "3", directory, plot_lim = plot_lim)

#4
directory <- "Simulated_Data_univariate/Station_4"
location_index <- 1
estimates <- find_univariate_data(directory,location_index)
names(estimates$parameters) <- names_param
estimates$parameters$location_indices <- rep(1, length(estimates$parameters$x))
plt_4 <- make_return_level_plot(return_periods, estimates, location_index, "4", directory, plot_lim = plot_lim)

plt_1 + plt_2 + plt_3 + plt_4 +
  plot_layout(ncol=4, guides = "collect") &
  theme(legend.position = "bottom") &
  ggsave(paste0("/../Plots/return_level/univariate_simulated_data.pdf"),
         path=path, width = 10, height = 3+1/2, units = c("in"))



#Real Data
directory <- "Real_Data"
ncol_plot <- 4
nrow_plot <- 1
plot_lim <- c(50,350)
plot_all_return_levels(return_periods, ncol_plot, nrow_plot, directory, c("Helgeroa","Oscarsborg","Oslo","Viker"), plot_lim = plot_lim)

#Real Data common xi
directory <- "Real_Data_common_xi"
ncol_plot <- 4
nrow_plot <- 1
plot_lim <- c(50,350)
plot_all_return_levels(return_periods, ncol_plot, nrow_plot, directory, c("Helgeroa","Oscarsborg","Oslo","Viker"), plot_lim = plot_lim)


#univariate real data
# plot_lim <- c(0,350)

#Helgeroa
directory <- "Univariate_real_Data_Helgeroa"
location_index <- 1
estimates <- find_univariate_data(directory,location_index)
names_param <- names(estimates$parameters)
names_param[2] <- "x"
names(estimates$parameters) <- names_param
estimates$parameters$location_indices <- rep(1, length(estimates$parameters$x))
plt_Helgeroa <- make_return_level_plot(return_periods, estimates, location_index, "Helgeroa", directory, plot_lim = plot_lim)

#Oscarsborg
directory <- "Univariate_real_Data_Oscarsborg"
location_index <- 1
estimates <- find_univariate_data(directory,location_index)
names(estimates$parameters) <- names_param
estimates$parameters$location_indices <- rep(1, length(estimates$parameters$x))
plt_Oscarsborg <- make_return_level_plot(return_periods, estimates, location_index, "Oscarsborg", directory, plot_lim = plot_lim)

#Oslo
directory <- "Univariate_real_Data_Oslo"
location_index <- 1
estimates <- find_univariate_data(directory,location_index)
names(estimates$parameters) <- names_param
estimates$parameters$location_indices <- rep(1, length(estimates$parameters$x))
plt_Oslo <- make_return_level_plot(return_periods, estimates, location_index, "Oslo", directory, plot_lim = plot_lim)

#Viker
directory <- "Univariate_real_Data_Viker"
location_index <- 1
estimates <- find_univariate_data(directory,location_index)
names(estimates$parameters) <- names_param
estimates$parameters$location_indices <- rep(1, length(estimates$parameters$x))
plt_Viker <- make_return_level_plot(return_periods, estimates, location_index, "Viker", directory, plot_lim = plot_lim)

plt_Helgeroa + plt_Oscarsborg + plt_Oslo + plt_Viker +
  plot_layout(ncol=4, guides = "collect") &
  theme(legend.position = "bottom") &
  ggsave(paste0("/../Plots/return_level/univariate_real_data.pdf"),
         path=path, width = 10, height = 3+1/2, units = c("in"))
  