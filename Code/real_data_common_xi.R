#!/usr/bin/env Rscript

set.seed(10)


#parameters:
dim <- ncol(distances)
alpha <- 0.5
beta <- 0.05

N_data <- c(nrow(Helgeroa_yearly_max_data),nrow(Oscarsborg_yearly_max_data),nrow(Oslo_yearly_max_data),nrow(Viker_yearly_max_data))
X_vec <- c(Helgeroa_yearly_max_data$yearly_max,Oscarsborg_yearly_max_data$yearly_max,Oslo_yearly_max_data$yearly_max,Viker_yearly_max_data$yearly_max)
location_indices <- rep(seq(1,4),N_data)


#STAN:
GEV_dat <- list(dim = dim,
                N = sum(N_data),
                x = X_vec,
                location_indices = location_indices,
                alpha = alpha,
                beta = beta,
                xi_limits = c(-0.5,0.5))


fit <- stan(file = paste0(path, "/reparametrized_GEV_common_xi.stan"), 
            data = GEV_dat,
            iter = 2*10^3,
            control = list(adapt_delta = 0.9,#0.995,
                           max_treedepth = 10))#20))

data_q_alpha <- rstan::extract(fit, pars = c("q"))
data_lns_beta <- rstan::extract(fit, pars = c("lns"))
data_xi <- rstan::extract(fit, pars = c("xi"))
data_parameters <- GEV_dat

write.table(data_q_alpha, file=paste0(path, "/../Results/Real_Data_common_xi/data_q_alpha.xml"))
write.table(data_lns_beta, file=paste0(path, "/../Results/Real_Data_common_xi/data_lns_beta.xml"))
write.table(data_xi, file=paste0(path, "/../Results/Real_Data_common_xi/data_xi.xml"))
saveRDS(data_parameters, file=paste0(path, "/../Results/Real_Data_common_xi/data_parameters.xml"))

print(fit, probs=c(0.025, 0.5, 0.975))

stan_trace(fit, inc_warmup = TRUE, nrow=3, 
           pars = c("xi",
                    "q",
                    "lns"))+
  ggsave(paste0(path, "/../Plots/Real_Data_common_xi/trace.pdf"),
         width = 10, height = 10, units = c("in"))

stan_hist(fit, nrow=3, 
          pars = c("xi",
                   "q",
                   "lns"))+
  ggsave(paste0(path, "/../Plots/Real_Data_common_xi/hist.pdf"),
         width = 10, height = 10, units = c("in"))

pdf(paste0(path, "/../Plots/Real_Data_common_xi/pairs.pdf"))
pairs(fit, pars = c("xi",
                    "q",
                    "lns",
                    "lp__"))
dev.off()



#load results
data_q_alpha <- read.table(file=paste0(path, "/../Results/Real_Data_common_xi/data_q_alpha.xml"))
data_lns_beta <- read.table(file=paste0(path, "/../Results/Real_Data_common_xi/data_lns_beta.xml"))
data_xi <- read.table(file=paste0(path, "/../Results/Real_Data_common_xi/data_xi.xml"))
data_parameters <- readRDS(paste0(path, "/../Results/Real_Data_common_xi/data_parameters.xml"))

data <- c(data_q_alpha, data_lns_beta, data_xi)
data <- as.data.frame(data)

print("results")
for (name in names(data)){
  cat(name, mean(get(name, data)),"\n")
}

#Plot

#parameters:
data_set <- data
data_set <- as.data.frame(data_set)

prior_x <- list("q.1"=seq(min(data$q.1),max(data$q.1),0.1),
                "q.2"=seq(min(data$q.2),max(data$q.2),0.1),
                "q.3"=seq(min(data$q.3),max(data$q.3),0.1),
                "q.4"=seq(min(data$q.4),max(data$q.4),0.1),
                "lns.1"=seq(min(data$lns.1),max(data$lns.1),0.01),
                "lns.2"=seq(min(data$lns.2),max(data$lns.2),0.01),
                "lns.3"=seq(min(data$lns.3),max(data$lns.3),0.01),
                "lns.4"=seq(min(data$lns.4),max(data$lns.4),0.01),
                "xi"=seq(-10,10,0.01))

prior_y <- list("q.1"=dnorm(prior_x$q.1, 100, 50),
                "q.2"=dnorm(prior_x$q.2, 100, 50),
                "q.3"=dnorm(prior_x$q.3, 100, 50),
                "q.4"=dnorm(prior_x$q.4, 100, 50),
                "lns.1"=dnorm(prior_x$lns.1, 0, 10),
                "lns.2"=dnorm(prior_x$lns.2, 0, 10),
                "lns.3"=dnorm(prior_x$lns.3, 0, 10),
                "lns.4"=dnorm(prior_x$lns.4, 0, 10),
                "xi"=dnorm(prior_x$xi, 0, 10))

prior_x$xi <- data_parameters$xi_limits[1]+(data_parameters$xi_limits[2]-data_parameters$xi_limits[1])*exp(prior_x$xi)/(1+exp(prior_x$xi))
#y_values_prior <- as.data.frame(y_values_prior)

parameter_names <- c(paste0("q.",c("Helgeroa","Oscarsborg","Oslo","Viker")),paste0("lns.",c("Helgeroa","Oscarsborg","Oslo","Viker")),"xi")
names(data_set) <- parameter_names
names(prior_x) <- parameter_names
names(prior_y) <- parameter_names

plots <- map(parameter_names, ~hist_plot_no_real_value_with_priors(data_set, .x, prior_x, prior_y))
ggarrange(plotlist=plots[9], ncol=1, nrow = 1)+
  ggsave(paste0(path, "/../Plots/Real_Data_common_xi/hist_with_real_xi.pdf"),
         width = 10*1/4, height = 8*1/4, units = c("in"))
ggarrange(plotlist=plots[1:4], ncol=4, nrow = 1)+
  ggsave(paste0(path, "/../Plots/Real_Data_common_xi/hist_with_real_q.pdf"),
         width = 10, height = 8*1/4, units = c("in"))
ggarrange(plotlist=plots[5:8], ncol=4, nrow = 1)+
  ggsave(paste0(path, "/../Plots/Real_Data_common_xi/hist_with_real_lns.pdf"),
         width = 10, height = 8*1/4, units = c("in"))

