#!/usr/bin/env Rscript

set.seed(10)

#Find parameter estimates for real data using the multivariate GEV model and STAN


#parameters:
dim <- ncol(distances)

alpha <- 0.5
beta <- 0.05
nu <- 1

N_data <- c(nrow(Helgeroa_yearly_max_data),nrow(Oscarsborg_yearly_max_data),nrow(Oslo_yearly_max_data),nrow(Viker_yearly_max_data))
X_vec <- c(Helgeroa_yearly_max_data$yearly_max,Oscarsborg_yearly_max_data$yearly_max,Oslo_yearly_max_data$yearly_max,Viker_yearly_max_data$yearly_max)
location_indices <- rep(seq(1,4),N_data)

covariate_matrix <- read.table(paste0(path,"/../Data/covariate_matrix.txt"))
ncovariates <- ncol(covariate_matrix)

#prepare for STAN:
GEV_dat <- list(dim = dim,
                N = sum(N_data),
                x = X_vec,
                location_indices = location_indices,
                covariates = covariate_matrix,
                distance_matrix = distances,
                alpha = alpha,
                beta = beta,
                nu = nu,
                xi_limits = c(-0.5,0.5),
                ncovariates = ncovariates)

#fit STAN
fit <- stan(file = paste0(path, "/reparametrized_GEV.stan"), 
            data = GEV_dat,
            iter = 2*10^3,
            control = list(adapt_delta = 0.995,
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

write.table(data_q_alpha, file=paste0(path, "/../Results/Real_Data/data_q_alpha.xml"))
write.table(data_lns_beta, file=paste0(path, "/../Results/Real_Data/data_lns_beta.xml"))
write.table(data_xi, file=paste0(path, "/../Results/Real_Data/data_xi.xml"))
write.table(data_beta_q, file=paste0(path, "/../Results/Real_Data/data_beta_q.xml"))
write.table(data_sd_q, file=paste0(path, "/../Results/Real_Data/data_sd_q.xml"))
write.table(data_range_q, file=paste0(path, "/../Results/Real_Data/data_range_q.xml"))
write.table(data_beta_lns, file=paste0(path, "/../Results/Real_Data/data_beta_lns.xml"))
write.table(data_sd_lns, file=paste0(path, "/../Results/Real_Data/data_sd_lns.xml"))
write.table(data_range_lns, file=paste0(path, "/../Results/Real_Data/data_range_lns.xml"))
saveRDS(data_parameters, file=paste0(path, "/../Results/Real_Data/data_parameters.xml"))


#plot results
stan_trace(fit, inc_warmup = TRUE, nrow=3, 
           pars = c("xi",
                    "beta_q",
                    "sd_q",
                    "range_q",
                    "beta_lns",
                    "sd_lns",
                    "range_lns"))+
  ggsave(paste0(path, "/../Plots/Real_Data/trace.pdf"),
         width = 10, height = 10, units = c("in"))

stan_hist(fit, nrow=3, 
          pars = c("xi",
                   "beta_q",
                   "sd_q",
                   "range_q",
                   "beta_lns",
                   "sd_lns",
                   "range_lns"))+
  ggsave(paste0(path, "/../Plots/Real_Data/hist.pdf"),
         width = 10, height = 10, units = c("in"))

pdf(paste0(path, "/../Plots/Real_Data/pairs.pdf"))
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
results <- find_data("Real_Data")
data <- results$data
data_parameters <- results$parameters


#prepare data and parameters:
data_set <- data[c((data_parameters$dim*2+1):(data_parameters$dim*2+5+2*data_parameters$ncovariates))]
data_set <- as.data.frame(data_set)
parameter_names <- names(data_set)

#priors
x_values_prior <- find_x_values_prior(data)
y_values_prior <- mapply(prior_plot_func, parameter_names, MoreArgs = list(x_values_prior))
x_values_prior$xi <- data_parameters$xi_limits[1]+(data_parameters$xi_limits[2]-data_parameters$xi_limits[1])*exp(x_values_prior$xi)/(1+exp(x_values_prior$xi))

#plot histograms of the parameters
plots <- map(parameter_names, ~histogram_of_results(data_set, .x, prior_x=x_values_prior, prior_y=y_values_prior))
ggarrange(plotlist=plots[1], ncol=1, nrow = 1)+
  ggsave(paste0(path, "/../Plots/Real_Data/hist_with_real_xi.pdf"),
         width = 10*1/4, height = 8*1/4, units = c("in"))
ggarrange(plotlist=plots[2:8], ncol=4, nrow = 2)+
  ggsave(paste0(path, "/../Plots/Real_Data/hist_with_real_q.pdf"),
         width = 10, height = 8*2/4, units = c("in"))
ggarrange(plotlist=plots[9:15], ncol=4, nrow = 2)+
  ggsave(paste0(path, "/../Plots/Real_Data/hist_with_real_lns.pdf"),
         width = 10, height = 8*2/4, units = c("in"))

#plot histograms of q
data_set <- data[c(1:data_parameters$dim)]
data_set <- as.data.frame(data_set)
names(data_set) <- paste0("q.",locations$location)

plots <- map(names(data_set), ~histogram_of_results(data_set, .x))
ggarrange(plotlist=plots, ncol=4, nrow = ceiling(length(plots)/4))+
  ggsave(paste0(path, "/../Plots/Real_Data/hist_q.pdf"),
         width = 10, height = 8*ceiling(length(plots)/4)/4, units = c("in"))

#plot histograms of lns
data_set <- data[c((1+data_parameters$dim):(data_parameters$dim*2))]
data_set <- as.data.frame(data_set)
names(data_set) <- paste0("lns.",locations$location)

plots <- map(names(data_set), ~histogram_of_results(data_set, .x))
ggarrange(plotlist=plots, ncol=4, nrow = ceiling(length(plots)/4))+
  ggsave(paste0(path, "/../Plots/Real_Data/hist_lns.pdf"),
         width = 10, height = 8*ceiling(length(plots)/4)/4, units = c("in"))
