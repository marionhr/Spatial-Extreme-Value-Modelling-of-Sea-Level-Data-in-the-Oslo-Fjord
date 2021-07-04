#!/usr/bin/env Rscript

#Map:


find_missing_location_parameters_map <- function(Covariate_matrix, locations, distances, data, data_parameters, data_name, dim_locations=3, dim=1){
  
  conditional_normality <- function(Covariate_matrix, beta_vector, Correlation_matrix, dim_locations, dim_new_locations, parameter_estimate){
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
  
  #q
  apply_conditional_normality_to_datasamples_q <- function(index, data, data_parameters, Covariate_matrix, dim_locations, dim_new_locations){
    Correlation_matrix_q <- data$sd_q[index]^2 * create_corr_mat(distances, data$range_q[index], nu=1)
    beta_vector_q <- data[index,(dim_locations*2+1+1):(dim_locations*2+1+data_parameters$ncovariates)] 
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
                      data_parameters, 
                      Covariate_matrix, 
                      dim_locations, 
                      dim))
  samples_q <- as.data.frame(t(samples_q))

  #hist:
  names(samples_q) <- paste0("q.",seq(1:dim))
  # 
  # 
  # plots <- map(names(samples_q), ~hist_plot_no_real_value(samples_q, .x))
  # ggarrange(plotlist=plots)+
  #   ggsave(paste0(path, "/../Plots/Map/",data_name,"/hist_new_locations_q.pdf"),
  #          width = 10, height = 8, units = c("in"))
  
  median_q_real_data <- colMedians(as.matrix(data[1:dim_locations]))
  median_q_new_locations <- colMedians(as.matrix(samples_q))
  
  ggplot() + 
    geom_polygon(data = fhidata::norway_map_counties, mapping = aes(x = long, y = lat, group = group, fill = hole), color = "black")+ 
    scale_fill_manual(values = c("white"))+
    theme_void()+ 
    coord_quickmap() + 
    coord_cartesian(xlim=c(9,12),ylim=c(58.8,60)) +
    geom_point(aes(x=locations[1:dim_locations,]$longitude, y=locations[1:dim_locations,]$latitude, colour=median_q_real_data, size = 0.8))+
    geom_point(aes(x=locations[(dim_locations+1):(dim_locations+dim),]$longitude, y=locations[(dim_locations+1):(dim_locations+dim),]$latitude, colour=median_q_new_locations, size = 0.7))+
    scale_size(range = c(2.5,3.5)) +
    scale_color_gradientn(colours=colorRamps::blue2green(10))+
    guides(fill=FALSE, size=FALSE)+
    labs(colour="q")+
    ggsave(paste0(path, "/../Plots/Map/",data_name,"/map_all_q.pdf"),
           width = 5, height = 5, units = c("in"))
  
  #lns
  apply_conditional_normality_to_datasamples_lns <- function(index, data, data_parameters, Covariate_matrix, dim_locations, dim_new_locations){
    Correlation_matrix_lns <- data$sd_lns[index]^2 * create_corr_mat(distances, data$range_lns[index], nu=1)
    beta_vector_lns <- data[index,(dim_locations*2+3+1+data_parameters$ncovariates):(dim_locations*2+3+data_parameters$ncovariates*2)] 
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
                                    data_parameters, 
                                    Covariate_matrix, 
                                    dim_locations, 
                                    dim))
  samples_lns <- as.data.frame(t(samples_lns))
  
  #hist:
  names(samples_lns) <- paste0("lns.",seq(1:dim))
  # 
  # 
  # plots <- map(names(samples_lns), ~hist_plot_no_real_value(samples_lns, .x))
  # ggarrange(plotlist=plots)+
  #   ggsave(paste0(path, "/../Plots/Map/",data_name,"/hist_new_locations_lns.pdf"),
  #          width = 10, height = 8, units = c("in"))
  
  median_lns_real_data <- colMedians(as.matrix(data[(dim_locations+1):(dim_locations*2)]))
  median_lns_new_locations <- colMedians(as.matrix(samples_lns))
  
  ggplot() + 
    geom_polygon(data = fhidata::norway_map_counties, mapping = aes(x = long, y = lat, group = group, fill = hole), color = "black")+ 
    scale_fill_manual(values = c("white"))+
    theme_void()+ 
    coord_quickmap() + 
    coord_cartesian(xlim=c(9,12),ylim=c(58.8,60)) +
    geom_point(aes(x=locations[1:dim_locations,]$longitude, y=locations[1:dim_locations,]$latitude, colour=median_lns_real_data, size = 0.8))+
    geom_point(aes(x=locations[(dim_locations+1):(dim_locations+dim),]$longitude, y=locations[(dim_locations+1):(dim_locations+dim),]$latitude, colour=median_lns_new_locations, size = 0.7))+
    scale_size(range = c(2.5,3.5)) +
    scale_color_gradientn(colours=colorRamps::blue2green(10))+
    guides(fill=FALSE, size=FALSE)+
    labs(colour="lns")+
    ggsave(paste0(path, "/../Plots/Map/",data_name,"/map_all_lns.pdf"),
           width = 5, height = 5, units = c("in"))
  
  return(list("samples_q" = samples_q, "samples_lns" = samples_lns))
}



apply_return_levels_to_datasamples_map <- function(index_location, index_estimate, samples_q, samples_lns, data, data_parameters, m, z_init=100){
  z <- return_levels(m, samples_q[index_estimate, index_location], samples_lns[index_estimate, index_location], data$xi[index_estimate], data_parameters$alpha, data_parameters$beta, z_init)
  return(z)
}

apply_locations_map <- function(index_location, samples_q, samples_lns, data, data_parameters, m){
  return(mapply(apply_return_levels_to_datasamples_map, index_estimate=1:nrow(data), MoreArgs = list("index_location"=index_location,"samples_q"=samples_q, "samples_lns"=samples_lns, "data"=data, "data_parameters"=data_parameters, "m"=m)))
}


plot_map_results <- function(locations, data_map, samples_q, samples_lns, median_and_confidence_interval_res,data_name,dim_locations_estimated_param){
  
  ggplot() + 
    geom_polygon(data = fhidata::norway_map_counties, mapping = aes(x = long, y = lat, group = group, fill = hole), color = "black")+ 
    scale_fill_manual(values = c("white"))+
    theme_void()+ 
    coord_quickmap() + 
    coord_cartesian(xlim=c(9,12),ylim=c(58.8,60)) +
    geom_point(aes(x=data_map[(dim_locations_estimated_param+1):(dim_locations_estimated_param+ncol(samples_q)),]$longitude, y=data_map[(dim_locations_estimated_param+1):(dim_locations_estimated_param+ncol(samples_q)),]$latitude, colour=median_and_confidence_interval_res[1,], size = 0.7))+
    geom_point(aes(x=locations$longitude, y=locations$latitude, size = 0.8))+
    scale_size(range = c(2.5,3.5)) +
    scale_color_gradientn(colours=colorRamps::blue2green(10), limits=c(min(median_and_confidence_interval_res[1,]),max(median_and_confidence_interval_res[1,])))+
    guides(fill=FALSE, size=FALSE)+
    labs(colour="lower bound")+
    ggsave(paste0(path, "/../Plots/Map/",data_name,"/map_return_level_lower_bound.pdf"),
           width = 4, height = 4, units = c("in"))
  
  
  ggplot() + 
    geom_polygon(data = fhidata::norway_map_counties, mapping = aes(x = long, y = lat, group = group, fill = hole), color = "black")+ 
    scale_fill_manual(values = c("white"))+
    theme_void()+ 
    coord_quickmap() + 
    coord_cartesian(xlim=c(9,12),ylim=c(58.8,60)) +
    geom_point(aes(x=data_map[(dim_locations_estimated_param+1):(dim_locations_estimated_param+ncol(samples_q)),]$longitude, y=data_map[(dim_locations_estimated_param+1):(dim_locations_estimated_param+ncol(samples_q)),]$latitude, colour=median_and_confidence_interval_res[2,], size = 0.7))+
    geom_point(aes(x=locations$longitude, y=locations$latitude, size = 0.8))+
    scale_size(range = c(2.5,3.5)) +
    scale_color_gradientn(colours=colorRamps::blue2green(10), limits=c(min(median_and_confidence_interval_res[2,]),max(median_and_confidence_interval_res[2,])))+
    guides(fill=FALSE, size=FALSE)+
    labs(colour="median")+
    ggsave(paste0(path, "/../Plots/Map/",data_name,"/map_return_level_median.pdf"),
           width = 4, height = 4, units = c("in"))
  
  ggplot() + 
    geom_polygon(data = fhidata::norway_map_counties, mapping = aes(x = long, y = lat, group = group, fill = hole), color = "black")+ 
    scale_fill_manual(values = c("white"))+
    theme_void()+ 
    coord_quickmap() + 
    coord_cartesian(xlim=c(9,12),ylim=c(58.8,60)) +
    geom_point(aes(x=data_map[(dim_locations_estimated_param+1):(dim_locations_estimated_param+ncol(samples_q)),]$longitude, y=data_map[(dim_locations_estimated_param+1):(dim_locations_estimated_param+ncol(samples_q)),]$latitude, colour=median_and_confidence_interval_res[3,], size = 0.7))+
    geom_point(aes(x=locations$longitude, y=locations$latitude, size = 0.8))+
    scale_size(range = c(2.5,3.5)) +
    scale_color_gradientn(colours=colorRamps::blue2green(10), trans="log10", limits=c(min(median_and_confidence_interval_res[3,]),max(median_and_confidence_interval_res[3,])))+
    guides(fill=FALSE, size=FALSE)+
    labs(colour="upper bound")+
    ggsave(paste0(path, "/../Plots/Map/",data_name,"/map_return_level_upper_bound.pdf"),
           width = 4, height = 4, units = c("in"))
  
}

find_locations_and_distances <- function(locations, nlocations){
  norwegian_coast_data <- read.table(file=paste0(path,"/../Data/sorted_norwegian_south_coast_data.txt"))
  data_map <- norwegian_coast_data[1:nlocations,]
  data_map$location <- seq(1,nrow(data_map))
  data_map$lat <- as.numeric(data_map$lat)
  data_map$long <- as.numeric(data_map$long)
  names(data_map) <- c("latitude", "longitude","code","location")
  
  data_map[which(is.element(data_map$latitude,locations$latitude) & is.element(data_map$longitude,locations$longitude)),] <- c(NA,NA,NA,NA)
  data_map <- data_map[complete.cases(data_map),]
  
  first_cov_mat <- make_covariate_matrix(data_map)
  
  print(which(is.na(first_cov_mat[,2])))
  print(which(is.na(first_cov_mat[,3])))
  print(which(is.na(first_cov_mat[,4])))
  print(which(is.na(first_cov_mat[,5])))
  data_map <- data_map[which(!is.na(first_cov_mat[,2])),]
  data_map <- data_map[which(!is.na(first_cov_mat[,3])),]
  data_map <- data_map[which(!is.na(first_cov_mat[,4])),]
  data_map <- data_map[which(!is.na(first_cov_mat[,5])),]
  
  data_map <- rbind(locations, data_map)
  data_map
  
  plt <- ggplot() + 
    geom_polygon(data = fhidata::norway_map_counties, mapping = aes(x = long, y = lat, group = group, fill = hole), color = "black")+ 
    scale_fill_manual(values = c("white"))+
    theme_void()+ 
    coord_quickmap() + 
    coord_cartesian(xlim=c(9,12),ylim=c(58.8,60)) +
    geom_point(aes(x=data_map$longitude, y=data_map$latitude, colour="red"))
  print(plt)
  
  data_map.sf = st_as_sf(data_map,coords = c("longitude","latitude"),
                         crs="+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+towgs84=0,0,0")
  transformed.data_map.sf <- st_transform(data_map.sf, crs = CRS("+proj=utm +zone=32 +datum=WGS84"))
  distances = as.matrix(st_distance(transformed.data_map.sf))
  distances  # unit m
  distances <- distances/(10^3)  # unit km
  units(distances) = NULL
  
  return(list("locations"=data_map, "distances"=distances))
}


make_map <- function(locations, nlocations, data_name, dim_locations_estimated_param){
  #locations
  loc_and_dist <- find_locations_and_distances(locations,nlocations)
  data_map <- loc_and_dist$locations
  distances <- loc_and_dist$distances
  
  #extract the data:
  covariate_matrix <- make_covariate_matrix(data_map)
  
  ggplot() + 
    geom_polygon(data = fhidata::norway_map_counties, mapping = aes(x = long, y = lat, group = group, fill = hole), color = "black")+ 
    scale_fill_manual(values = c("white"))+
    theme_void()+ 
    coord_quickmap() + 
    coord_cartesian(xlim=c(9,12),ylim=c(58.8,60)) +
    geom_point(aes(x=data_map$longitude, y=data_map$latitude), color="red")+
    guides(fill=FALSE, size=FALSE)+
    ggsave(paste0(path, "/../Plots/Map/",data_name,"/data_map.pdf"), width = 10, height = 8, units = c("in"))
  
  #load results
  data_parameters <- find_data(data_name)
  data <- data_parameters$data
  
  missing_locations_q_lns <- find_missing_location_parameters_map(covariate_matrix, data_map, distances, data, data_parameters$parameters, data_name, dim_locations_estimated_param, nrow(data_map)-dim_locations_estimated_param)
  samples_q <- missing_locations_q_lns$samples_q
  samples_lns <- missing_locations_q_lns$samples_lns
  
  m <- 200
  return_levels_map <- mapply(apply_locations_map, 1:ncol(samples_q), MoreArgs = list("samples_q"=samples_q, "samples_lns"=samples_lns, "data"=data, "data_parameters"=data_parameters$parameters, "m"=m))
  median_and_confidence_interval_res <- mapply(median_and_confidence_interval, 1:ncol(samples_q), MoreArgs = list("all_return_levels"=return_levels_map))
  
  print(head(sort(median_and_confidence_interval_res[3,])))
  print(which(median_and_confidence_interval_res[3,]<200))
  
  plot_map_results(locations, data_map, samples_q, samples_lns, median_and_confidence_interval_res, data_name, dim_locations_estimated_param)
}

#Real Data

make_map(locations, 1000, "Real_Data", 4)
# 
# #Simulated Data
# 
# dim <- 4
# directory <- paste0("Simulated_Data_",dim,"_locations")
# data_original_parameters <- readRDS(paste0(path, "/../Results/",directory,"/data_original_parameters.xml"))
# names(data_original_parameters$locations) <- c("latitude","longitude","code","location")
# make_map(data_original_parameters$locations, 1000, directory, nrow(data_original_parameters$locations))

# different prior
# dim <- 4
# directory <- paste0("Simulated_Data_",dim,"_locations_different_prior")
# data_original_parameters <- readRDS(paste0(path, "/../Results/",directory,"/data_original_parameters.xml"))
# names(data_original_parameters$locations) <- c("latitude","longitude","code","location")
# make_map(data_original_parameters$locations, 1000, directory, nrow(data_original_parameters$locations))
#
# dim <- 6
# directory <- paste0("Simulated_Data_",dim,"_locations")
# data_original_parameters <- readRDS(paste0(path, "/../Results/",directory,"/data_original_parameters.xml"))
# names(data_original_parameters$locations) <- c("latitude","longitude","code","location")
# make_map(data_original_parameters$locations, 10, directory, nrow(data_original_parameters$locations))
# 
# dim <- 8
# directory <- paste0("Simulated_Data_",dim,"_locations")
# data_original_parameters <- readRDS(paste0(path, "/../Results/",directory,"/data_original_parameters.xml"))
# names(data_original_parameters$locations) <- c("latitude","longitude","code","location")
# make_map(data_original_parameters$locations, 1000, directory, nrow(data_original_parameters$locations))
# 
# dim <- 15
# directory <- paste0("Simulated_Data_",dim,"_locations")
# data_original_parameters <- readRDS(paste0(path, "/../Results/",directory,"/data_original_parameters.xml"))
# names(data_original_parameters$locations) <- c("latitude","longitude","code","location")
# make_map(data_original_parameters$locations, 1000, directory, nrow(data_original_parameters$locations))
# 
# different prior
# dim <- 15
# directory <- paste0("Simulated_Data_",dim,"_locations_different_prior")
# data_original_parameters <- readRDS(paste0(path, "/../Results/",directory,"/data_original_parameters.xml"))
# names(data_original_parameters$locations) <- c("latitude","longitude","code","location")
# make_map(data_original_parameters$locations, 1000, directory, nrow(data_original_parameters$locations))
# 
# dim <- 30
# directory <- paste0("Simulated_Data_",dim,"_locations")
# data_original_parameters <- readRDS(paste0(path, "/../Results/",directory,"/data_original_parameters.xml"))
# names(data_original_parameters$locations) <- c("latitude","longitude","code","location")
# make_map(data_original_parameters$locations, 1000, directory, nrow(data_original_parameters$locations))