#!/usr/bin/env Rscript

#Maps with return level values for the Oslo fjord:


find_missing_location_parameters_map <- function(Covariate_matrix, data_map, distances, data, parameters, directory, dim_locations=3, dim=1){
  #function to find parameters q and lns for the locations where there is no data to estimate them
  #input:
  #Covariate_matrix: matrix of covariate values for all locations
  #data_map: data frame with latitude and longitude for all locations
  #distances: distances between locations in data_map
  #data: parameter estimates
  #parameters: list of parameter values
  #directory: where to save results
  #dim_locations: number of locations with data
  #dim: number of locations without data
  #output:
  #samples_q: samples for the parameter q
  #samples_lns: samples for the parameter lns
  ###########################################################################################
  
  
  conditional_normality <- function(Covariate_matrix, beta_vector, Correlation_matrix, dim_locations, dim_new_locations, parameter_estimate){
    #conditional normality of A given B, where A coincides with new locations and B coincides with the locations with data
    #input:
    #Covariate_matrix: matrix with covariate values for all locations
    #beta_vector: vector with parameter values for beta
    #Correlation_matrix: matrix with correlation values for all locations based on Matern correlation function
    #dim_locations: dimension of locations with data
    #dim_new_locations: dimension of locations without data
    #parameter_estimate: parameter estimates for q or lns based on data
    #output:
    #mean and variance for the new locations for either q or lns
    ###########################################################################################################
    
    #A|B
    Covariates_A <- Covariate_matrix[(dim_locations+1):(dim_locations+dim_new_locations),]
    Covariates_B <- Covariate_matrix[1:dim_locations,]
    
    Correlation_A <- Correlation_matrix[(dim_locations+1):(dim_locations+dim_new_locations),(dim_locations+1):(dim_locations+dim_new_locations)]
    Correlation_AB <- Correlation_matrix[1:dim_locations,(dim_locations+1):(dim_locations+dim_new_locations)]
    Correlation_B <- Correlation_matrix[1:dim_locations,1:dim_locations]
    
    mu_conditional <- Covariates_A %*% beta_vector +
      t(Correlation_AB) %*% solve(Correlation_B, parameter_estimate - Covariates_B %*% beta_vector)
    Sigma_conditional <- Correlation_A - 
      t(Correlation_AB) %*% solve(Correlation_B, Correlation_AB)
    
    return(list(mu=mu_conditional, Sigma=Sigma_conditional))
  }
  
  #find parameter estimates for q:
  
  apply_conditional_normality_to_datasamples_q <- function(index, data, parameters, Covariate_matrix, dim_locations, dim_new_locations){
    #input:
    #index: index for parameter samples
    #data: parameter estimates
    #parameters: list of parameter values
    #Covariate_matrix: matrix with covariate values for all locations
    #dim_locations: dimension of locations with data
    #dim_new_locations: dimension of locations without data
    #output:
    #samples for q
    ###################################################################
    Correlation_matrix_q <- data$sd_q[index]^2 * create_corr_mat(distances, data$range_q[index], nu=1)
    beta_vector_q <- data[index,(dim_locations*2+1+1):(dim_locations*2+1+parameters$ncovariates)] 
    beta_vector_q <- t(beta_vector_q)
    parameter_estimate_q <- data[index,1:dim_locations] 
    
    cond_norm <- conditional_normality(Covariate_matrix, 
                                       beta_vector_q, 
                                       Correlation_matrix_q, 
                                       dim_locations, 
                                       dim_new_locations, 
                                       parameter_estimate_q)
    std_norm <- rnorm(dim_new_locations)
    samples_q <- cond_norm$mu + t(chol(cond_norm$Sigma)) %*% std_norm 
    
    return(samples_q)
  }
  
  samples_q <- mapply(apply_conditional_normality_to_datasamples_q, 
                      1:nrow(data),
                      MoreArgs=list(data, 
                      parameters, 
                      Covariate_matrix, 
                      dim_locations, 
                      dim))
  samples_q <- as.data.frame(t(samples_q))

  #plot results:
  names(samples_q) <- paste0("q.",seq(1:dim))
  
  median_q_real_data <- colMedians(as.matrix(data[1:dim_locations]))
  median_q_new_locations <- colMedians(as.matrix(samples_q))
  
  ggplot() + 
    geom_polygon(data = fhidata::norway_map_counties, mapping = aes(x = long, y = lat, group = group, fill = hole), color = "black")+ 
    scale_fill_manual(values = c("white"))+
    theme_void()+ 
    coord_quickmap() + 
    coord_cartesian(xlim=c(9,12),ylim=c(58.8,60)) +
    geom_point(aes(x=data_map[1:dim_locations,]$longitude, y=data_map[1:dim_locations,]$latitude, colour=median_q_real_data, size = 0.8))+
    geom_point(aes(x=data_map[(dim_locations+1):(dim_locations+dim),]$longitude, y=data_map[(dim_locations+1):(dim_locations+dim),]$latitude, colour=median_q_new_locations, size = 0.7))+
    scale_size(range = c(2.5,3.5)) +
    scale_color_gradientn(colours=colorRamps::blue2green(10))+
    guides(fill=FALSE, size=FALSE)+
    labs(colour="q")+
    ggsave(paste0(path, "/../Plots/Map/",directory,"/map_all_q.pdf"),
           width = 5, height = 5, units = c("in"))
  
  #find parameter estimates for lns
  apply_conditional_normality_to_datasamples_lns <- function(index, data, parameters, Covariate_matrix, dim_locations, dim_new_locations){
    #input:
    #index: index for parameter samples
    #data: parameter estimates
    #parameters: list of parameter values
    #Covariate_matrix: matrix with covariate values for all locations
    #dim_locations: dimension of locations with data
    #dim_new_locations: dimension of locations without data
    #output:
    #samples for lns
    ###################################################################
    Correlation_matrix_lns <- data$sd_lns[index]^2 * create_corr_mat(distances, data$range_lns[index], nu=1)
    beta_vector_lns <- data[index,(dim_locations*2+3+1+parameters$ncovariates):(dim_locations*2+3+parameters$ncovariates*2)] 
    beta_vector_lns <- t(beta_vector_lns)
    parameter_estimate_lns <- data[index,(dim_locations+1):(dim_locations*2)] 
    
    cond_norm <- conditional_normality(Covariate_matrix, 
                                       beta_vector_lns, 
                                       Correlation_matrix_lns, 
                                       dim_locations, 
                                       dim_new_locations, 
                                       parameter_estimate_lns)
    std_norm <- rnorm(dim_new_locations)
    samples_lns <- cond_norm$mu + t(chol(cond_norm$Sigma)) %*% std_norm 
    
    return(samples_lns)
  }
  
  samples_lns <- mapply(apply_conditional_normality_to_datasamples_lns, 
                      1:nrow(data),
                      MoreArgs=list(data, 
                                    parameters, 
                                    Covariate_matrix, 
                                    dim_locations, 
                                    dim))
  samples_lns <- as.data.frame(t(samples_lns))
  
  #plot results
  names(samples_lns) <- paste0("lns.",seq(1:dim))
  
  median_lns_real_data <- colMedians(as.matrix(data[(dim_locations+1):(dim_locations*2)]))
  median_lns_new_locations <- colMedians(as.matrix(samples_lns))
  
  ggplot() + 
    geom_polygon(data = fhidata::norway_map_counties, mapping = aes(x = long, y = lat, group = group, fill = hole), color = "black")+ 
    scale_fill_manual(values = c("white"))+
    theme_void()+ 
    coord_quickmap() + 
    coord_cartesian(xlim=c(9,12),ylim=c(58.8,60)) +
    geom_point(aes(x=data_map[1:dim_locations,]$longitude, y=data_map[1:dim_locations,]$latitude, colour=median_lns_real_data, size = 0.8))+
    geom_point(aes(x=data_map[(dim_locations+1):(dim_locations+dim),]$longitude, y=data_map[(dim_locations+1):(dim_locations+dim),]$latitude, colour=median_lns_new_locations, size = 0.7))+
    scale_size(range = c(2.5,3.5)) +
    scale_color_gradientn(colours=colorRamps::blue2green(10))+
    guides(fill=FALSE, size=FALSE)+
    labs(colour="lns")+
    ggsave(paste0(path, "/../Plots/Map/",directory,"/map_all_lns.pdf"),
           width = 5, height = 5, units = c("in"))
  
  return(list("samples_q" = samples_q, "samples_lns" = samples_lns))
}



apply_return_levels_to_datasamples_map <- function(location_index, index_samples, samples_q, samples_lns, data, parameters, return_period, z_init=100){
  #input:
  #location_index: which location
  #index_samples: which sample
  #samples_q: samples for the parameter q
  #samples_lns: samples for the parameter lns
  #data: parameter estimates
  #parameters: list of parameter values
  #return_period: which return period
  #z_init: initial value for z in the nonlinear equation solver
  #output:
  #z: return level value for one sample and one location
  #############################################################
  z <- return_levels(return_period, samples_q[index_samples, location_index], samples_lns[index_samples, location_index], data$xi[index_samples], parameters$alpha, parameters$beta, z_init)
  return(z)
}

apply_return_levels_to_locations_map <- function(location_index, samples_q, samples_lns, data, parameters, return_period){
  #input:
  #location_index: which location
  #samples_q: samples for the parameter q
  #samples_lns: samples for the parameter lns
  #data: parameter estimates
  #parameters: list of parameter values
  #return_period: which return period
  #output:
  #return level values for all samples and one location
  #############################################################
  return(mapply(apply_return_levels_to_datasamples_map, index_samples=1:nrow(data), MoreArgs = list("location_index"=location_index,"samples_q"=samples_q, "samples_lns"=samples_lns, "data"=data, "parameters"=parameters, "return_period"=return_period)))
}


plot_map_results <- function(data_map, number_of_samples, median_and_confidence_interval_res, directory, dim_locations_estimated_param){
  #input:
  #data_map: data frame with latitude and longitude for all locations
  #number_of_samples: number of samples
  #median_and_confidence_interval_res: medians and confidence intervals for the additional locations
  #directory: where to save results
  #dim_locations_estimated_param: number of locations with data
  #output:
  #maps for the median and for confidence bounds
  #####################################################################################################
  
  ggplot() + 
    geom_polygon(data = fhidata::norway_map_counties, mapping = aes(x = long, y = lat, group = group, fill = hole), color = "black")+ 
    scale_fill_manual(values = c("white"))+
    theme_void()+ 
    coord_quickmap() + 
    coord_cartesian(xlim=c(9,12),ylim=c(58.8,60)) +
    geom_point(aes(x=data_map[(dim_locations_estimated_param+1):(dim_locations_estimated_param+number_of_samples),]$longitude, y=data_map[(dim_locations_estimated_param+1):(dim_locations_estimated_param+number_of_samples),]$latitude, colour=median_and_confidence_interval_res[1,], size = 0.7))+
    geom_point(aes(x=data_map[1:dim_locations_estimated_param,]$longitude, y=data_map[1:dim_locations_estimated_param,]$latitude, size = 0.8))+
    scale_size(range = c(2.5,3.5)) +
    scale_color_gradientn(colours=colorRamps::blue2green(10), limits=c(min(median_and_confidence_interval_res[1,]),max(median_and_confidence_interval_res[1,])))+
    guides(fill=FALSE, size=FALSE)+
    labs(colour="lower bound")+
    ggsave(paste0(path, "/../Plots/Map/",directory,"/map_return_level_lower_bound.pdf"),
           width = 4, height = 4, units = c("in"))
  
  
  ggplot() + 
    geom_polygon(data = fhidata::norway_map_counties, mapping = aes(x = long, y = lat, group = group, fill = hole), color = "black")+ 
    scale_fill_manual(values = c("white"))+
    theme_void()+ 
    coord_quickmap() + 
    coord_cartesian(xlim=c(9,12),ylim=c(58.8,60)) +
    geom_point(aes(x=data_map[(dim_locations_estimated_param+1):(dim_locations_estimated_param+number_of_samples),]$longitude, y=data_map[(dim_locations_estimated_param+1):(dim_locations_estimated_param+number_of_samples),]$latitude, colour=median_and_confidence_interval_res[2,], size = 0.7))+
    geom_point(aes(x=data_map[1:dim_locations_estimated_param,]$longitude, y=data_map[1:dim_locations_estimated_param,]$latitude, size = 0.8))+
    scale_size(range = c(2.5,3.5)) +
    scale_color_gradientn(colours=colorRamps::blue2green(10), limits=c(min(median_and_confidence_interval_res[2,]),max(median_and_confidence_interval_res[2,])))+
    guides(fill=FALSE, size=FALSE)+
    labs(colour="median")+
    ggsave(paste0(path, "/../Plots/Map/",directory,"/map_return_level_median.pdf"),
           width = 4, height = 4, units = c("in"))
  
  ggplot() + 
    geom_polygon(data = fhidata::norway_map_counties, mapping = aes(x = long, y = lat, group = group, fill = hole), color = "black")+ 
    scale_fill_manual(values = c("white"))+
    theme_void()+ 
    coord_quickmap() + 
    coord_cartesian(xlim=c(9,12),ylim=c(58.8,60)) +
    geom_point(aes(x=data_map[(dim_locations_estimated_param+1):(dim_locations_estimated_param+number_of_samples),]$longitude, y=data_map[(dim_locations_estimated_param+1):(dim_locations_estimated_param+number_of_samples),]$latitude, colour=median_and_confidence_interval_res[3,], size = 0.7))+
    geom_point(aes(x=data_map[1:dim_locations_estimated_param,]$longitude, y=data_map[1:dim_locations_estimated_param,]$latitude, size = 0.8))+
    scale_size(range = c(2.5,3.5)) +
    scale_color_gradientn(colours=colorRamps::blue2green(10), trans="log10", limits=c(min(median_and_confidence_interval_res[3,]),max(median_and_confidence_interval_res[3,])))+
    guides(fill=FALSE, size=FALSE)+
    labs(colour="upper bound")+
    ggsave(paste0(path, "/../Plots/Map/",directory,"/map_return_level_upper_bound.pdf"),
           width = 4, height = 4, units = c("in"))
  
}

find_locations_and_distances <- function(locations, nlocations){
  #input:
  #locations: data frame with longitude and latitude for locations with data
  #nlocations: number of additional locations in the maps
  #output:
  #list of data frame with all locations and distances between them
  ###########################################################################
  
  #read in additional locations. The locations are found so that they are as far as possible from existing locations
  norwegian_coast_data <- read.table(file=paste0(path,"/../Data/sorted_norwegian_south_coast_data.txt"))
  data_map <- norwegian_coast_data[1:nlocations,]
  data_map$location <- seq(1,nrow(data_map))
  data_map$lat <- as.numeric(as.character(data_map$lat))
  data_map$long <- as.numeric(as.character(data_map$long))
  names(data_map) <- c("latitude", "longitude","code","location")
  
  #remove duplicated locations
  data_map[which(is.element(data_map$latitude,locations$latitude) & is.element(data_map$longitude,locations$longitude)),] <- c(NA,NA,NA,NA)
  data_map <- data_map[complete.cases(data_map),]
  
  #covariate matrix
  first_cov_mat <- make_covariate_matrix(data_map)
  
  #remove locations that don't have covariate data
  data_map <- data_map[which(!is.na(first_cov_mat[,2])),]
  data_map <- data_map[which(!is.na(first_cov_mat[,3])),]
  data_map <- data_map[which(!is.na(first_cov_mat[,4])),]
  data_map <- data_map[which(!is.na(first_cov_mat[,5])),]
  
  #add in the locations with data
  data_map <- rbind(locations, data_map)
  data_map
  
  #find distance matrix
  distances <- find_distances(data_map)
  
  return(list("locations"=data_map, "distances"=distances))
}


make_map <- function(locations, nlocations, directory, dim_locations_estimated_param){
  #input:
  #locations: data frame with longitude and latitude
  #nlocations: number of additional locations in the maps
  #directory: where to save results
  #dim_locations_estimated_param: number of locations with estimated parameters based on data
  #output:
  #maps of return levels
  ##############################################################################################
  
  #locations and distances
  loc_and_dist <- find_locations_and_distances(locations,nlocations)
  data_map <- loc_and_dist$locations
  distances <- loc_and_dist$distances
  
  #find covariate data
  covariate_matrix <- make_covariate_matrix(data_map)
  
  #plot all the locations
  ggplot() + 
    geom_polygon(data = fhidata::norway_map_counties, mapping = aes(x = long, y = lat, group = group, fill = hole), color = "black")+ 
    scale_fill_manual(values = c("white"))+
    theme_void()+ 
    coord_quickmap() + 
    coord_cartesian(xlim=c(9,12),ylim=c(58.8,60)) +
    geom_point(aes(x=data_map$longitude, y=data_map$latitude), color="red")+
    guides(fill=FALSE, size=FALSE)+
    ggsave(paste0(path, "/../Plots/Map/",directory,"/data_map.pdf"), width = 10, height = 8, units = c("in"))
  
  #load parameter estimates based on data
  estimates <- find_data(directory)
  data <- estimates$data
  parameters <- estimates$parameters
  
  #find parameters q and lns for the locations that don't have data
  missing_locations_q_lns <- find_missing_location_parameters_map(covariate_matrix, data_map, distances, data, parameters, directory, dim_locations_estimated_param, nrow(data_map)-dim_locations_estimated_param)
  samples_q <- missing_locations_q_lns$samples_q
  samples_lns <- missing_locations_q_lns$samples_lns
  
  #return period
  return_period <- 200
  
  #find return levels
  return_levels_map <- mapply(apply_return_levels_to_locations_map, 1:ncol(samples_q), MoreArgs = list("samples_q"=samples_q, "samples_lns"=samples_lns, "data"=data, "parameters"=parameters, "return_period"=return_period))
  
  #find medians and confidence intervals
  median_and_confidence_interval_res <- mapply(median_and_confidence_interval, 1:ncol(samples_q), MoreArgs = list("all_return_levels"=return_levels_map))
  
  #make maps
  plot_map_results(data_map, ncol(samples_q), median_and_confidence_interval_res, directory, dim_locations_estimated_param)
}

############################################################################################################

#Real Data

make_map(locations, 1000, "Real_Data", 4)

#Simulated Data

dim <- 4
directory <- paste0("Simulated_Data_",dim,"_locations")
data_original_parameters <- readRDS(paste0(path, "/../Results/",directory,"/data_original_parameters.xml"))
names(data_original_parameters$locations) <- c("latitude","longitude","code","location")
make_map(data_original_parameters$locations, 1000, directory, nrow(data_original_parameters$locations))

dim <- 8
directory <- paste0("Simulated_Data_",dim,"_locations")
data_original_parameters <- readRDS(paste0(path, "/../Results/",directory,"/data_original_parameters.xml"))
names(data_original_parameters$locations) <- c("latitude","longitude","code","location")
make_map(data_original_parameters$locations, 1000, directory, nrow(data_original_parameters$locations))

dim <- 15
directory <- paste0("Simulated_Data_",dim,"_locations")
data_original_parameters <- readRDS(paste0(path, "/../Results/",directory,"/data_original_parameters.xml"))
names(data_original_parameters$locations) <- c("latitude","longitude","code","location")
make_map(data_original_parameters$locations, 1000, directory, nrow(data_original_parameters$locations))

# different prior
dim <- 15
directory <- paste0("Simulated_Data_",dim,"_locations_different_prior")
data_original_parameters <- readRDS(paste0(path, "/../Results/",directory,"/data_original_parameters.xml"))
names(data_original_parameters$locations) <- c("latitude","longitude","code","location")
make_map(data_original_parameters$locations, 1000, directory, nrow(data_original_parameters$locations))
