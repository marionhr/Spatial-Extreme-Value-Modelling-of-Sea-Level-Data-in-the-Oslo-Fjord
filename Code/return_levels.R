#!/usr/bin/env Rscript

set.seed(10)

#data functions:
apply_return_level_function_to_parameter_samples <- function(index_samples, m, location_index, estimates){
  ret <- return_levels(m, get(paste0("q.",location_index),estimates$data)[index_samples], get(paste0("lns.",location_index),estimates$data)[index_samples], estimates$data$xi[index_samples], estimates$parameters$alpha, estimates$parameters$beta)
  return(ret)
}

apply_return_periods <- function(m, estimates, location_index){
  ret <- lapply(1:nrow(estimates$data), apply_return_level_function_to_parameter_samples, m, location_index, estimates)
  ret <- as.data.frame(ret)
  names(ret) <- 1:nrow(estimates$data)
  ret <- t(ret)
  return(ret)
}

plot_return_levels <- function(return_levels, return_periods, location, directory, plot_lim, data, location_index, cv=0){
  plot_lim <- c(min(plot_lim[1],return_levels[1,3]),max(plot_lim[2],return_levels[3,]))
  if(cv==0 | location_index!=4){
    empirical_data <- data$parameters$x[data$parameters$location_indices==location_index]
  }else{
    empirical_data <- plot_lim[1]-10
  }
  
  plt <- ggplot()+
    geom_point(aes(x=length(empirical_data)/1:length(empirical_data), 
                   y=sort(empirical_data,decreasing = TRUE), 
                   colour="empirical \nreturn levels"))+
    geom_line(aes(x=return_periods, y=return_levels[2,], colour="median\n"))+
    geom_line(aes(x=return_periods, y=return_levels[1,], colour="95% credible \ninterval\n"),linetype = "dashed")+
    geom_line(aes(x=return_periods, y=return_levels[3,], colour="95% credible \ninterval\n"),linetype = "dashed")+
    xlab("Return period (years)")+
    ylab(paste0(location,": Return level (cm)"))+
    ylim(plot_lim)
    #guides(color=guide_legend(location))
    #ggsave(paste0("/../Plots/return_level/",directory,".",location,".pdf"),path=path, width = 5, height = 3, units = c("in"))
  return(plt)
}


make_return_level_plot <- function(m, estimates, location_index, locations, directory, plot_lim, cv=0){
  ret_lev <- mapply(apply_return_periods, m, MoreArgs=list("estimates"=estimates, "location_index"=location_index))
  
  median_conf_int <- mapply(median_and_confidence_interval, 1:length(m), MoreArgs = list("all_return_levels"=ret_lev))
  print(median_conf_int)
  plt <- plot_return_levels(median_conf_int, m, locations[location_index], directory, plot_lim, estimates, location_index, cv)
  return(plt)
}

plot_all_return_levels <- function(m, ncol_plot, nrow_plot, directory, locations, find_data_function = find_data, cv=0, plot_lim=c(0,1000)){
  estimates <- find_data_function(directory)
  nloc <- estimates$parameters$dim + cv
  plots <- map(1:nloc, ~make_return_level_plot(m, estimates, .x, locations=locations, directory=directory, plot_lim=plot_lim, cv=cv))
  ggarrange(plotlist=plots, ncol=ncol_plot, nrow = nrow_plot, common.legend = T, legend="bottom")+
    ggsave(paste0(path, "/../Plots/return_level/",directory,".pdf"),
           width = 10*ncol_plot/4, height = 12*nrow_plot/4+1/2, units = c("in"))
}

      
#return level plots
m <- c(1,2,5,10,20,30,40,50,100,150,200)

#Simulated Data 4 locations
directory <- "Simulated_Data_4_locations"
ncol_plot <- 4
nrow_plot <- 1
plot_lim <- c(100,300)
plot_all_return_levels(m, ncol_plot, nrow_plot, directory, 1:4, plot_lim = plot_lim)

#Simulated Data 4 locations common xi
directory <- "Simulated_Data_common_xi"
ncol_plot <- 4
nrow_plot <- 1
plot_lim <- c(100,300)
plot_all_return_levels(m, ncol_plot, nrow_plot, directory, 1:4, plot_lim = plot_lim)

#Simulated Data 4 locations different prior
directory <- "Simulated_Data_4_locations_different_prior"
ncol_plot <- 4
nrow_plot <- 1
plot_lim <- c(100,300)
plot_all_return_levels(m, ncol_plot, nrow_plot, directory, 1:15, plot_lim = plot_lim)

#Simulated Data 6 locations
directory <- "Simulated_Data_6_locations"
ncol_plot <- 4
nrow_plot <- 2
plot_lim <- c(100,300)
plot_all_return_levels(m, ncol_plot, nrow_plot, directory, 1:6, plot_lim = plot_lim)

#Simulated Data 8 locations
directory <- "Simulated_Data_8_locations"
ncol_plot <- 4
nrow_plot <- 2
plot_lim <- c(100,300)
plot_all_return_levels(m, ncol_plot, nrow_plot, directory, 1:8, plot_lim = plot_lim)

#Simulated Data 15 locations
directory <- "Simulated_Data_15_locations"
ncol_plot <- 4
nrow_plot <- 4
plot_lim <- c(100,300)
plot_all_return_levels(m, ncol_plot, nrow_plot, directory, 1:15, plot_lim = plot_lim)

#Simulated Data 15 locations different prior
directory <- "Simulated_Data_15_locations_different_prior"
ncol_plot <- 4
nrow_plot <- 4
plot_lim <- c(100,300)
plot_all_return_levels(m, ncol_plot, nrow_plot, directory, 1:15, plot_lim = plot_lim)

#Simulated Data 30 locations
directory <- "Simulated_Data_30_locations"
ncol_plot <- 5
nrow_plot <- 6
plot_lim <- c(100,300)
plot_all_return_levels(m, ncol_plot, nrow_plot, directory, 1:30, plot_lim = plot_lim)


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
plt_1 <- make_return_level_plot(m, estimates, location_index, "1", directory, plot_lim = plot_lim)

#2
directory <- "Simulated_Data_univariate/Station_2"
location_index <- 1
estimates <- find_univariate_data(directory,location_index)
names(estimates$parameters) <- names_param
estimates$parameters$location_indices <- rep(1, length(estimates$parameters$x))
plt_2 <- make_return_level_plot(m, estimates, location_index, "2", directory, plot_lim = plot_lim)

#3
directory <- "Simulated_Data_univariate/Station_3"
location_index <- 1
estimates <- find_univariate_data(directory,location_index)
names(estimates$parameters) <- names_param
estimates$parameters$location_indices <- rep(1, length(estimates$parameters$x))
plt_3 <- make_return_level_plot(m, estimates, location_index, "3", directory, plot_lim = plot_lim)

#4
directory <- "Simulated_Data_univariate/Station_4"
location_index <- 1
estimates <- find_univariate_data(directory,location_index)
names(estimates$parameters) <- names_param
estimates$parameters$location_indices <- rep(1, length(estimates$parameters$x))
plt_4 <- make_return_level_plot(m, estimates, location_index, "4", directory, plot_lim = plot_lim)

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
plot_all_return_levels(m, ncol_plot, nrow_plot, directory, c("Helgeroa","Oscarsborg","Oslo","Viker"), plot_lim = plot_lim)

#Real Data common xi
directory <- "Real_Data_common_xi"
ncol_plot <- 4
nrow_plot <- 1
plot_lim <- c(50,350)
plot_all_return_levels(m, ncol_plot, nrow_plot, directory, c("Helgeroa","Oscarsborg","Oslo","Viker"), plot_lim = plot_lim)

#CV real data
ncol_plot <- 4
nrow_plot <- 1
# plot_lim<- c(0,220)

directory <- "Real_Data_Helgeroa"
plot_all_return_levels(m, ncol_plot, nrow_plot, directory, c("Oscarsborg","Oslo","Viker","Helgeroa"), find_cv_data, 1, plot_lim = plot_lim)

directory <- "Real_Data_Oscarsborg"
plot_all_return_levels(m, ncol_plot, nrow_plot, directory, c("Helgeroa","Oslo","Viker","Oscarsborg"), find_cv_data, 1, plot_lim = plot_lim)

directory <- "Real_Data_Oslo"
plot_all_return_levels(m, ncol_plot, nrow_plot, directory, c("Helgeroa","Oscarsborg","Viker","Oslo"), find_cv_data, 1, plot_lim = plot_lim)

directory <- "Real_Data_Viker"
plot_all_return_levels(m, ncol_plot, nrow_plot, directory, c("Helgeroa","Oscarsborg","Oslo","Viker"), find_cv_data, 1, plot_lim = plot_lim)


#less data one location
# plot_lim <- c(0,400)

directory <- "Real_Data_Helgeroa_Less_Data"
plot_all_return_levels(m, ncol_plot, nrow_plot, directory, c("Helgeroa","Oscarsborg","Oslo","Viker"), plot_lim = plot_lim)

directory <- "Real_Data_Oscarsborg_Less_Data"
plot_all_return_levels(m, ncol_plot, nrow_plot, directory, c("Helgeroa","Oscarsborg","Oslo","Viker"), plot_lim = plot_lim)

directory <- "Real_Data_Oslo_Less_Data"
plot_all_return_levels(m, ncol_plot, nrow_plot, directory, c("Helgeroa","Oscarsborg","Oslo","Viker"), plot_lim = plot_lim)

directory <- "Real_Data_Viker_Less_Data"
plot_all_return_levels(m, ncol_plot, nrow_plot, directory, c("Helgeroa","Oscarsborg","Oslo","Viker"), plot_lim = plot_lim)


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
plt_Helgeroa <- make_return_level_plot(m, estimates, location_index, "Helgeroa", directory, plot_lim = plot_lim)

#Oscarsborg
directory <- "Univariate_real_Data_Oscarsborg"
location_index <- 1
estimates <- find_univariate_data(directory,location_index)
names(estimates$parameters) <- names_param
estimates$parameters$location_indices <- rep(1, length(estimates$parameters$x))
plt_Oscarsborg <- make_return_level_plot(m, estimates, location_index, "Oscarsborg", directory, plot_lim = plot_lim)

#Oslo
directory <- "Univariate_real_Data_Oslo"
location_index <- 1
estimates <- find_univariate_data(directory,location_index)
names(estimates$parameters) <- names_param
estimates$parameters$location_indices <- rep(1, length(estimates$parameters$x))
plt_Oslo <- make_return_level_plot(m, estimates, location_index, "Oslo", directory, plot_lim = plot_lim)

#Viker
directory <- "Univariate_real_Data_Viker"
location_index <- 1
estimates <- find_univariate_data(directory,location_index)
names(estimates$parameters) <- names_param
estimates$parameters$location_indices <- rep(1, length(estimates$parameters$x))
plt_Viker <- make_return_level_plot(m, estimates, location_index, "Viker", directory, plot_lim = plot_lim)

plt_Helgeroa + plt_Oscarsborg + plt_Oslo + plt_Viker +
  plot_layout(ncol=4, guides = "collect") &
  theme(legend.position = "bottom") &
  ggsave(paste0("/../Plots/return_level/univariate_real_data.pdf"),
         path=path, width = 10, height = 3+1/2, units = c("in"))
  