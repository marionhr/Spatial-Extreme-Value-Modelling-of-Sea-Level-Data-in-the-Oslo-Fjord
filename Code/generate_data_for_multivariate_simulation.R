#!/usr/bin/env Rscript

set.seed(10)

#Generate random data used in the simulation study

#dimentions
dim <- 4  #4, 8, 15


#locations:
sorted_norwegian_coast_data <- read.table(file=paste0(path,"/../Data/sorted_norwegian_south_coast_data.txt"))

if (dim>4){ #the first four locations are the same as for the real data
  simulated_locations <- as.data.frame(rbind(as.matrix(locations[,c(3,4,2,1)]),as.matrix(sorted_norwegian_coast_data[1:(dim-4),])))
} else {
  simulated_locations <- locations[,c(3,4,2,1)]
}

names(simulated_locations) <- names(sorted_norwegian_coast_data)
simulated_locations$location <- seq(1,dim)
simulated_locations$lat <- as.numeric(as.character(simulated_locations$lat))
simulated_locations$long <- as.numeric(as.character(simulated_locations$long))

#distances:
simulated_distances <- find_distances(simulated_locations)

#covariate matrix
make_covariate_matrix(simulated_locations)

#number of data per location
s <- round(runif(dim, 50, 100))
if (dim == 8){
  s[5:dim] <- round(runif(dim, 5, 15))
}

#parameters:
xi <- -0.05
alpha <- 0.5
beta <- 0.05
nu <- 1

beta_q <- c(90,9,12,-18,15)
sd_q <- 0.3
range_q <- 100

beta_lns <- c(4.3,0.2,0.1,-0.6,0.2)
sd_lns <- 0.2
range_lns <- 80

#generate data from the multinormal distribution for q and lns:

mean_q <- covariate_matrix %*% beta_q
Sigma_q <- sd_q^2*create_corr_mat(simulated_distances, range_q, nu)

mean_lns <- covariate_matrix %*% beta_lns
Sigma_lns <- sd_lns^2*create_corr_mat(simulated_distances, range_lns, nu)

q_alpha <- mvrnorm(1, mean_q, Sigma_q)
lns_beta <- mvrnorm(1, mean_lns, Sigma_lns) 


#generate data from the gev distribution

location_indices <- rep(c(1:dim), s)

reparametrized <- reparametrization(q_alpha, lns_beta, xi, alpha, beta)
X_vec <- rgev(length(location_indices), reparametrized$mu[location_indices], reparametrized$sigma[location_indices], xi)


#Plot locations
ggplot() + 
  geom_polygon(data = fhidata::norway_map_counties, mapping = aes(x = long, y = lat, group = group, fill = hole), color = "black")+ 
  scale_fill_manual(values = c("white"))+
  theme_void()+ 
  coord_quickmap() + 
  coord_cartesian(xlim=c(9,12),ylim=c(58.8,60)) +
  geom_point(aes(x=simulated_locations$long, y=simulated_locations$lat), color="red")+
  geom_label_repel(aes(x=simulated_locations$long, y=simulated_locations$lat, label=simulated_locations$location), color="red", max.overlaps = 15) + 
  guides(fill=FALSE, size=FALSE)+
  ggsave(paste0(path, "/../Plots/Extra_Plots/",dim,"_Simulated_Locations.pdf"), width = 5, height = 4, units = c("in"))

