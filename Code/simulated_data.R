#!/usr/bin/env Rscript

set.seed(10)


#Find parameter estimates for simulated data using the multivariate GEV model and STAN

#if choose different prior:
different_prior <- FALSE   #remember to also change the prior in the STAN file: reparametrized_GEV.stan

save_in_repository <- paste0("Simulated_Data_",dim,"_locations")
if (different_prior==TRUE){
  save_in_repository <- paste0("Simulated_Data_",dim,"_locations_different_prior")
}


#prepare for STAN
GEV_dat <- list(dim = dim,
                N = sum(s),
                x = X_vec,
                location_indices = location_indices,
                covariates = covariate_matrix,
                distance_matrix = simulated_distances,
                alpha = alpha,
                beta = beta,
                nu = nu,
                ncovariates = ncovariates,
                xi_limits = c(-0.5,0.5))


#fit STAN
fit <- stan(file = paste0(path, "/reparametrized_GEV.stan"), 
            data = GEV_dat,
            iter = 2*10^3,
            control = list(adapt_delta = 0.99,
                           max_treedepth = 20))


#save results
data_q_alpha <- rstan::extract(fit, pars = c("q"))
data_lns_beta <- rstan::extract(fit, pars = c("lns"))
data_xi <- rstan::extract(fit, pars = c("xi"))
data_beta_q <- rstan::extract(fit, pars = c("beta_q"))
data_sd_q <- rstan::extract(fit, pars = c("sd_q"))
data_range_q <- rstan::extract(fit, pars = c("range_q"))
data_beta_lns <- rstan::extract(fit, pars = c("beta_lns"))
data_sd_lns <- rstan::extract(fit, pars = c("sd_lns"))
data_range_lns <- rstan::extract(fit, pars = c("range_lns"))
data_parameters <- GEV_dat
data_original_parameters <- list(q_alpha = q_alpha,
                                 lns_beta = lns_beta,
                                 xi = xi,
                                 beta_q = beta_q,
                                 sd_q = sd_q,
                                 range_q = range_q,
                                 beta_lns = beta_lns,
                                 sd_lns = sd_lns,
                                 range_lns = range_lns,
                                 locations = simulated_locations)

write.table(data_q_alpha, file=paste0(path, "/../Results/",save_in_repository,"/data_q_alpha.xml"))
write.table(data_lns_beta, file=paste0(path, "/../Results/",save_in_repository,"/data_lns_beta.xml"))
write.table(data_xi, file=paste0(path, "/../Results/",save_in_repository,"/data_xi.xml"))
write.table(data_beta_q, file=paste0(path, "/../Results/",save_in_repository,"/data_beta_q.xml"))
write.table(data_sd_q, file=paste0(path, "/../Results/",save_in_repository,"/data_sd_q.xml"))
write.table(data_range_q, file=paste0(path, "/../Results/",save_in_repository,"/data_range_q.xml"))
write.table(data_beta_lns, file=paste0(path, "/../Results/",save_in_repository,"/data_beta_lns.xml"))
write.table(data_sd_lns, file=paste0(path, "/../Results/",save_in_repository,"/data_sd_lns.xml"))
write.table(data_range_lns, file=paste0(path, "/../Results/",save_in_repository,"/data_range_lns.xml"))
saveRDS(data_parameters, file=paste0(path, "/../Results/",save_in_repository,"/data_parameters.xml"))
saveRDS(data_original_parameters, file=paste0(path, "/../Results/",save_in_repository,"/data_original_parameters.xml"))


#plots
stan_trace(fit, inc_warmup = TRUE, nrow=3, 
           pars = c("xi",
                    "beta_q",
                    "sd_q",
                    "range_q",
                    "beta_lns",
                    "sd_lns",
                    "range_lns"))+
  ggsave(paste0(path, "/../Plots/",save_in_repository,"/trace.pdf"),
         width = 10, height = 10, units = c("in"))

stan_hist(fit, nrow=3, 
          pars = c("xi",
                   "beta_q",
                   "sd_q",
                   "range_q",
                   "beta_lns",
                   "sd_lns",
                   "range_lns"))+
  ggsave(paste0(path, "/../Plots/",save_in_repository,"/hist.pdf"),
         width = 10, height = 10, units = c("in"))


pdf(paste0(path, "/../Plots/",save_in_repository,"/pairs.pdf"))
pairs(fit, pars = c("xi",
                    "beta_q",
                    "sd_q",
                    "range_q",
                    "beta_lns",
                    "sd_lns",
                    "range_lns",
                    "lp__"))
dev.off()

################################################################################
################################################################################


#load results
results <- find_data(save_in_repository)
data <- results$data
data_parameters <- results$parameters
data_original_parameters <- readRDS(paste0(path, "/../Results/",save_in_repository,"/data_original_parameters.xml"))


#prepare data and parameters:
parameter_names <- names(data[(data_parameters$dim*2+1):ncol(data)])
parameter_set <-  setNames(c(data_original_parameters$xi, 
                             data_original_parameters$beta_q, 
                             data_original_parameters$sd_q, 
                             data_original_parameters$range_q, 
                             data_original_parameters$beta_lns, 
                             data_original_parameters$sd_lns,
                             data_original_parameters$range_lns),
                           parameter_names)
data_set <- data[(data_parameters$dim*2+1):ncol(data)]

#priors:
x_values_prior <- find_x_values_prior(data)
y_values_prior <- mapply(prior_plot_func, parameter_names, MoreArgs = list(x_values_prior))

if (different_prior==TRUE){
  y_values_prior$range_q <- dgamma(get("range_q",x_values_prior),1.5,1/40)
  y_values_prior$range_lns <- dgamma(get("range_lns",x_values_prior),1.5,1/40)
}

x_values_prior$xi <- data_parameters$xi_limits[1]+(data_parameters$xi_limits[2]-data_parameters$xi_limits[1])*exp(x_values_prior$xi)/(1+exp(x_values_prior$xi))

#plot histograms of the parameters
plots <- map(parameter_names, ~hist_plot_with_priors(data_set, as.data.frame(t(parameter_set)), .x, x_values_prior, y_values_prior))
ggarrange(plotlist=plots[1], ncol=1, nrow = 1)+
  ggsave(paste0(path, "/../Plots/",save_in_repository,"/hist_with_real_xi.pdf"),
         width = 10/4, height = 8*1/4, units = c("in"))
ggarrange(plotlist=plots[2:8], ncol=4, nrow = 2)+
  ggsave(paste0(path, "/../Plots/",save_in_repository,"/hist_with_real_q.pdf"),
         width = 10, height = 8*2/4, units = c("in"))
ggarrange(plotlist=plots[9:15], ncol=4, nrow = 2)+
  ggsave(paste0(path, "/../Plots/",save_in_repository,"/hist_with_real_lns.pdf"),
         width = 10, height = 8*2/4, units = c("in"))

#plot histograms of q
parameter_names <- names(data[1:(data_parameters$dim)])
data_set <- data[1:(data_parameters$dim)]
parameter_set <- as.data.frame(data_original_parameters$q_alpha)
names(parameter_set) <- parameter_names

plots <- map(parameter_names, ~hist_plot(data_set, parameter_set, .x))
ggarrange(plotlist=plots, ncol=4, nrow = ceiling(length(plots)/4))+
  ggsave(paste0(path, "/../Plots/",save_in_repository,"/hist_parameters_q.pdf"),
         width = 10, height = 8*ceiling(length(plots)/4)/4, units = c("in"))

#plot histograms of lns
parameter_names <- names(data[(data_parameters$dim+1):(data_parameters$dim*2)])
data_set <- data[(data_parameters$dim+1):(data_parameters$dim*2)]
parameter_set <- as.data.frame(data_original_parameters$lns_beta)
names(parameter_set) <- parameter_names

plots <- map(parameter_names, ~hist_plot(data_set, parameter_set, .x))
ggarrange(plotlist=plots, ncol=4, nrow = ceiling(length(plots)/4))+
  ggsave(paste0(path, "/../Plots/",save_in_repository,"/hist_parameters_lns.pdf"),
         width = 10, height = 8*ceiling(length(plots)/4)/4, units = c("in"))
