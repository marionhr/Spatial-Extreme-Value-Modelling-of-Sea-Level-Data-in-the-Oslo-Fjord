#!/usr/bin/env Rscript

set.seed(10)

#Generate random data:
dim <- 15  #4, 8, 15,   if choose 8 remember to change s as well
nu <- 1

#locations:
sorted_norwegian_coast_data <- read.table(file=paste0(path,"/../Data/sorted_norwegian_south_coast_data.txt"))
if (dim>4){
  simulated_locations <- as.data.frame(rbind(as.matrix(locations[,c(3,4,2,1)]),as.matrix(sorted_norwegian_coast_data[1:(dim-4),])))
} else {
  simulated_locations <- locations[,c(3,4,2,1)]
}

names(simulated_locations) <- names(sorted_norwegian_coast_data)
simulated_locations$location <- seq(1,dim)
simulated_locations$lat <- as.numeric(as.character(simulated_locations$lat))
simulated_locations$long <- as.numeric(as.character(simulated_locations$long))

simulated_distances <- find_distances(simulated_locations)

#parameter values:
xi <- -0.05

#covariate matrix
u10_raster <- raster(paste0(path,"/../Data/u10.nc"))
v10_raster <- raster(paste0(path,"/../Data/v10.nc"))
msl_raster <- raster(paste0(path,"/../Data/msl.nc"))
tide_raster <- raster(paste0(path,"/../Data/tide.nc"))

location_data_u10 <- raster::extract(u10_raster, cbind(simulated_locations$long,
                                                       simulated_locations$lat))
location_data_v10 <- raster::extract(v10_raster, cbind(simulated_locations$long,
                                                       simulated_locations$lat))
location_data_msl <- raster::extract(msl_raster, cbind(simulated_locations$long,
                                                       simulated_locations$lat))
location_data_tide <- raster::extract(tide_raster, cbind(simulated_locations$long,
                                                         simulated_locations$lat))
ncovariates <- 5
covariate_matrix <- matrix(NA, ncol=ncovariates, nrow=dim)
covariate_matrix[,1] <- rep(1,dim)
covariate_matrix[,2] <- t(location_data_u10)
covariate_matrix[,3] <- t(location_data_v10)
covariate_matrix[,4] <- t(location_data_msl)
covariate_matrix[,5] <- t(location_data_tide)
covariate_matrix

s <- round(runif(dim, 50, 100))
# s[5:dim] <- round(runif(dim, 5, 15))
# s[5:6] <- c(1,1)
s

alpha <- 0.5
beta <- 0.05

beta_q <- c(90,9,12,-18,15)
covariate_matrix %*% beta_q

sd_q <- 0.3
range_q <- 100

beta_lns <- c(4.3,0.2,0.1,-0.6,0.2)
covariate_matrix %*% beta_lns

sd_lns <- 0.2
range_lns <- 80

max(covariate_matrix %*% beta_q)
min(covariate_matrix %*% beta_q)

max(covariate_matrix %*% beta_lns)
min(covariate_matrix %*% beta_lns)

#Multinormal distribution for q and lns:

mean_q <- covariate_matrix %*% beta_q
Sigma_q <- sd_q^2*create_corr_mat(simulated_distances, range_q, nu)

mean_lns <- covariate_matrix %*% beta_lns
Sigma_lns <- sd_lns^2*create_corr_mat(simulated_distances, range_lns, nu)

q_alpha <- mvrnorm(1, mean_q, Sigma_q)
lns_beta <- mvrnorm(1, mean_lns, Sigma_lns) 


#generate data from the gev distribution

location_indices <- rep(c(1:dim), s)
location_indices

reparametrized <- reparametrization(q_alpha, lns_beta, xi, alpha, beta)
X_vec <- rgev(length(location_indices), reparametrized$mu[location_indices], reparametrized$sigma[location_indices], xi)
print(X_vec)

min(X_vec)
max(X_vec)

which(X_vec>134.68154)

#Plot locations
ggplot() + 
  geom_polygon(data = fhidata::norway_map_counties, mapping = aes(x = long, y = lat, group = group, fill = hole), color = "black")+ 
  scale_fill_manual(values = c("white"))+
  theme_void()+ 
  coord_quickmap() + 
  coord_cartesian(xlim=c(9,12),ylim=c(58.8,60)) +
  geom_point(aes(x=simulated_locations$long, y=simulated_locations$lat), color="red")+
  #geom_text(aes(x=simulated_locations$long, y=simulated_locations$lat, label=simulated_locations$location),position=position_jitter(width=0.01,height=0.01), color="red") + 
  geom_label_repel(aes(x=simulated_locations$long, y=simulated_locations$lat, label=simulated_locations$location), color="red", max.overlaps = 15) + 
  guides(fill=FALSE, size=FALSE)+
  ggsave(paste0(path, "/../Plots/Extra_Plots/",dim,"_Simulated_Locations.pdf"), width = 5, height = 4, units = c("in"))

