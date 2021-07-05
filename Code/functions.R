#!/usr/bin/env Rscript

#####################################################################
#Reparametrization:
reparametrization <- function(q_alpha, lns_beta, xi, alpha, beta){
  #change the GEV parameters from quantile based parameters q and lns
  #to standard parameters mu and sigma
  #
  #input:
  #q_alpha: the alpha quantile of the GEV distribution
  #lns_beta: the difference between the 1-beta/2 quantile and the 
  #          beta/2 quantile of the GEV distribution
  #xi: shape parameter of the GEV distribution
  #
  #output:
  #mu: location parameter of the GEV distribution
  #sigma: scale parameter of the GEV distribution
  ###################################################################
  
  l_xi <- function(x, xi){
    return((-log(x))^(-xi))
  }
  
  l <- function(x){
    return(log(-log(x)))
  }
  
  if (xi!=0){
    mu <- q_alpha - exp(lns_beta)*(l_xi(alpha,xi)-1)/(l_xi(1-beta/2,xi)-l_xi(beta/2,xi))
    sigma <- exp(lns_beta)*xi/(l_xi(1-beta/2,xi)-l_xi(beta/2,xi))
  }
  else{
    mu <- q_alpha + exp(lns_beta)*(l(alpha)-1)/(l(beta/2)-l(1-beta/2))
    sigma <- exp(lns_beta)/(l(beta/2)-l(1-beta/2))
  }
  return (list("mu"=mu, "sigma"=sigma))
}


#####################################################################
#correlation matrix:
create_corr_mat <- function(distance_matrix, range, nu) {
  #Produces the Matern correlation matrix for given distance matrix
  #and parameters range and nu
  ###################################################################
  kappa <- (8*nu)^(1/2)/range
  N <- nrow(distance_matrix)
  corr <- matrix(NA, nrow=N, ncol=N)
  for (i in 1:N){
    corr[i, i] = 1 + 1e-10
  } 
  for (i in 2:N) {
    for (j in 1:(i - 1)) {
      corr[i,j] = besselK(distance_matrix[i, j] * kappa, nu) *
        (kappa * distance_matrix[i, j]) ^ nu * 2^(1 - nu) / gamma(nu)
      corr[j, i] = corr[i, j]
    }
  }
  return(corr)
}


#####################################################################
#plots of histograms:

#data: values to be used in the histograms
#value: real value to compare results to
#name: list of parameter names to make plots for
#prior_x + prior_y: the prior distributions to compare the histograms to

hist_plot <- function(data, value, name){
  plt <- ggplot()+
    geom_histogram(aes(x=get(name, data),y=..density..))+
    geom_vline(aes(xintercept = get(name, value)))+
    xlab(name)+
    xlim(c(max(min(get(name,data)),-10^4),min(max(get(name,data)),10^4)))
  return(plt)
}

hist_plot_no_real_value <- function(data, name){
  plt <- ggplot()+
    geom_histogram(aes(x=get(name, data),y=..density..))+
    xlab(name)+
    xlim(c(max(min(get(name,data)),-10^4),min(max(get(name,data)),10^4)))
  return(plt)
}

hist_plot_with_priors <- function(data, value, name, prior_x, prior_y){
  plt <- ggplot()+
    geom_histogram(aes(x=get(name, data),y=..density..))+
    geom_vline(aes(xintercept = get(name, value)))+
    geom_line(aes(x=get(name,prior_x), y=get(name, prior_y)))+
    xlab(name)
  return(plt)
}

hist_plot_no_real_value_with_priors <- function(data, name, prior_x, prior_y){
  plt <- ggplot()+
    geom_histogram(aes(x=get(name, data),y=..density..))+
    geom_line(aes(x=get(name,prior_x), y=get(name, prior_y)))+
    xlab(name)
  return(plt)
}

# find x values for the priors
find_x_values_prior <- function(data){
  return(list("xi"=seq(-10,10,0.1),
       "beta_q.1"=seq(min(data$beta_q.1),max(data$beta_q.1),0.1),
       "beta_q.2"=seq(min(data$beta_q.2),max(data$beta_q.2),0.1),
       "beta_q.3"=seq(min(data$beta_q.3),max(data$beta_q.3),0.1),
       "beta_q.4"=seq(min(data$beta_q.4),max(data$beta_q.4),0.1),
       "beta_q.5"=seq(min(data$beta_q.5),max(data$beta_q.5),0.1),
       "sd_q"=seq(min(data$sd_q),max(data$sd_q),0.1),
       "range_q"=seq(min(data$range_q),max(data$range_q),0.1),
       "beta_lns.1"=seq(min(data$beta_lns.1),max(data$beta_lns.1),0.1),
       "beta_lns.2"=seq(min(data$beta_lns.2),max(data$beta_lns.2),0.1),
       "beta_lns.3"=seq(min(data$beta_lns.3),max(data$beta_lns.3),0.1),
       "beta_lns.4"=seq(min(data$beta_lns.4),max(data$beta_lns.4),0.1),
       "beta_lns.5"=seq(min(data$beta_lns.5),max(data$beta_lns.5),0.1),
       "sd_lns"=seq(min(data$sd_lns),max(data$sd_lns),0.1),
       "range_lns"=seq(min(data$range_lns),max(data$range_q),0.1)))
}

#find y values of the priors based on parameter names
prior_plot_func <- function(name, x_values_prior){
  if (grepl("xi", name, fixed = TRUE)){
    return(dnorm(get(name,x_values_prior),0,10))
  }
  if (grepl("beta", name, fixed = TRUE)){
    return(dnorm(get(name,x_values_prior),0,100))
  }
  if (grepl("sd", name, fixed = TRUE)){
    return(dgamma(get(name,x_values_prior),5/4,1/4))
  }
  if (grepl("range", name, fixed = TRUE)){
    return(dgamma(get(name,x_values_prior),5,1/12))
  }
}


#####################################################################
#normalize data
normalize_data <- function(data){
  return((data-min(data, na.rm = TRUE))/(max(data, na.rm = TRUE)-min(data, na.rm = TRUE)))
}


#####################################################################
#return level functions:

return_level_function <- function(z, m, mu, sigma, xi){
  #equation for the return levels and return periods for the GEV 
  #distribution with standard parameters
  ###################################################################
  ret <- pgev(z,mu,sigma,xi)-(1-1/m)
  return(ret)
}

return_levels <- function(m, q, lns, xi, alpha, beta, z_init=100){
  #find the return levels z of the reparametrized GEV distribution 
  #given return period m
  ###################################################################
  reparameterized <- reparametrization(q, lns, xi, alpha, beta)
  z <- nleqslv(x=z_init, fn=return_level_function, jac=NULL, m, reparameterized$mu, reparameterized$sigma, xi)
  return(z$x)
}

#median and confidence bounds of return levels
median_and_confidence_interval <- function(index, all_return_levels){
  #find median and confidence bounds given estimates of the return 
  #levels and an index indicating which column the return level 
  #estimates are located in
  ###################################################################
  return_levels <- all_return_levels[,index]
  sort_return_levels <- sort(return_levels)
  return(c(sort_return_levels[0.025*length(sort_return_levels)], 
           sort_return_levels[0.5*length(sort_return_levels)], 
           sort_return_levels[0.975*length(sort_return_levels)]))
}


#####################################################################
#find data:
find_data <- function(directory){
  #extract the parameter estimates from files in the given directory
  ###################################################################
  data_q_alpha <- read.table(file=paste0(path, "/../Results/",directory,"/data_q_alpha.xml"))
  data_lns_beta <- read.table(file=paste0(path, "/../Results/",directory,"/data_lns_beta.xml"))
  data_xi <- read.table(file=paste0(path, "/../Results/",directory,"/data_xi.xml"))
  data_beta_q <- read.table(file=paste0(path, "/../Results/",directory,"/data_beta_q.xml"))
  data_sd_q <- read.table(file=paste0(path, "/../Results/",directory,"/data_sd_q.xml"))
  data_range_q <- read.table(file=paste0(path, "/../Results/",directory,"/data_range_q.xml"))
  data_beta_lns <- read.table(file=paste0(path, "/../Results/",directory,"/data_beta_lns.xml"))
  data_sd_lns <- read.table(file=paste0(path, "/../Results/",directory,"/data_sd_lns.xml"))
  data_range_lns <- read.table(file=paste0(path, "/../Results/",directory,"/data_range_lns.xml"))
  data_parameters <- readRDS(paste0(path, "/../Results/",directory,"/data_parameters.xml"))
  
  data <- c(data_q_alpha, data_lns_beta, data_xi,data_beta_q,data_sd_q,data_range_q,data_beta_lns,data_sd_lns,data_range_lns)
  data <- as.data.frame(data)
  return(list("data"=data, "parameters"=data_parameters))
}

find_univariate_data <- function(directory, location_index){
  #extract the parameter estimates of the univariate model
  ###################################################################
  data_q_alpha <- read.table(file=paste0(path, "/../Results/",directory,"/data_q_alpha.xml"))
  data_lns_beta <- read.table(file=paste0(path, "/../Results/",directory,"/data_lns_beta.xml"))
  data_xi <- read.table(file=paste0(path, "/../Results/",directory,"/data_xi.xml"))
  data_parameters <- readRDS(paste0(path, "/../Results/",directory,"/data_parameters.xml"))
  data <- c(data_q_alpha, data_lns_beta, data_xi)
  data <- as.data.frame(data)
  names(data) <- c(paste0("q.",location_index),paste0("lns.",location_index),"xi")
  
  return(list("data"=data, "parameters"=data_parameters))
}


#####################################################################
#covariate matrix
make_covariate_matrix <- function(locations){
  #extract covariate data for locations of interest
  ###################################################################
  u10_raster <- raster(paste0(path,"/../Data/u10.nc"))
  v10_raster <- raster(paste0(path,"/../Data/v10.nc"))
  msl_raster <- raster(paste0(path,"/../Data/msl.nc"))
  tide_raster <- raster(paste0(path,"/../Data/tide.nc"))
  
  location_data_u10 <- raster::extract(u10_raster, cbind(locations$longitude,
                                                         locations$latitude))
  location_data_v10 <- raster::extract(v10_raster, cbind(locations$longitude,
                                                         locations$latitude))
  location_data_msl <- raster::extract(msl_raster, cbind(locations$longitude,
                                                         locations$latitude))
  location_data_tide <- raster::extract(tide_raster, cbind(locations$longitude,
                                                           locations$latitude))
  
  #covariate matrix
  dim <- nrow(locations)
  ncovariates <- 5
  covariate_matrix <- matrix(NA, ncol=ncovariates, nrow=dim)
  covariate_matrix[,1] <- rep(1,dim)
  covariate_matrix[,2] <- t(location_data_u10)
  covariate_matrix[,3] <- t(location_data_v10)
  covariate_matrix[,4] <- t(location_data_msl)
  covariate_matrix[,5] <- t(location_data_tide)
  covariate_matrix
  return(covariate_matrix)
}

#####################################################################

find_distances <- function(locations){
  #find distances between locations
  #input: data frame with longitude and latitude as columns
  #output: distance matrix
  ###################################################################
  names(locations)
  names(simulated_locations)
  latitude_name <- {if(is.element("lat",names(locations))) "lat" else "latitude"}
  longitude_name <- {if(is.element("long",names(locations))) "long" else "longitude"}
  
  locations.sf = st_as_sf(locations,coords = c(longitude_name,latitude_name),
                          crs="+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+towgs84=0,0,0")
  transformed.locations.sf <- st_transform(locations.sf, crs = CRS("+proj=utm +zone=32 +datum=WGS84"))
  distances <- as.matrix(st_distance(transformed.locations.sf))
  distances #units m
  distances <- distances/(10^3)  #unit km
  units(distances) <- NULL
  return(distances)
}

#########################################################################

date_to_year <- function(data){
  #input: data frame with columns for: year, month, day and hour
  #output: list of years with decimals representing the month, day and hour
  ###########################################################################
  return(data$year+1/(as.numeric(strftime(paste0(data$year,"-",12,"-",31),format="%j")))*(as.numeric(strftime(paste0(data$year, "-",data$month,"-",data$day),format="%j"))+data$hour*(1/24)))
}


