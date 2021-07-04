#!/usr/bin/env Rscript

set.seed(10)


univariate_real_data <- function(y, location_name){
  
  N <- length(y)
  
  alpha <- 0.5
  beta <- 0.05
  
  #STAN:
  GEV_dat <- list(N = N,
                  y = y,
                  alpha = alpha,
                  beta = beta,
                  xi_limits = c(-0.5,0.5))
  if (location_name == "Viker"){
    GEV_dat$xi_limits <- c(-0.95,0.95)
  }
  
  fit <- stan(file = paste0(path, "/reparametrized_GEV_univariate.stan"), 
              data = GEV_dat,
              iter = 2*10^3,
              control = list(adapt_delta = 0.999999,
                             max_treedepth = 15))
  
  
  #save results
  data_q_alpha <- rstan::extract(fit, pars = c("q"))
  data_lns_beta <- rstan::extract(fit, pars = c("lns"))
  data_xi <- rstan::extract(fit, pars = c("xi"))
  data_parameters <- GEV_dat
  
  write.table(data_q_alpha, file=paste0(path, "/../Results/Univariate_real_Data_",location_name,"/data_q_alpha.xml"))
  write.table(data_lns_beta, file=paste0(path, "/../Results/Univariate_real_Data_",location_name,"/data_lns_beta.xml"))
  write.table(data_xi, file=paste0(path, "/../Results/Univariate_real_Data_",location_name,"/data_xi.xml"))
  saveRDS(data_parameters, file=paste0(path, "/../Results/Univariate_real_Data_",location_name,"/data_parameters.xml"))
  
  print(fit, probs=c(0.025, 0.5, 0.975))
  
  #Plots:
  stan_trace(fit, inc_warmup = TRUE, nrow=3, 
             pars = c("xi",
                      "q",
                      "lns"))+
    ggsave(paste0(path, "/../Plots/Univariate_real_Data_",location_name,"/trace.pdf"),
           width = 5, height = 5, units = c("in"))
  
  
  stan_hist(fit, nrow=3, 
            pars = c("xi",
                     "q",
                     "lns"))+
    ggsave(paste0(path, "/../Plots/Univariate_real_Data_",location_name,"/hist.pdf"),
           width = 5, height = 5, units = c("in"))
  
  
  pdf(paste0(path, "/../Plots/Univariate_real_Data_",location_name,"/pairs.pdf"))
  pairs(fit, pars = c("xi",
                      "q",
                      "lns",
                      "lp__"))
  dev.off()
}
  
load_real_results_function <- function(location_name){
  #load results
  data_q_alpha <- read.table(file=paste0(path, "/../Results/Univariate_real_Data_",location_name,"/data_q_alpha.xml"))
  data_lns_beta <- read.table(file=paste0(path, "/../Results/Univariate_real_Data_",location_name,"/data_lns_beta.xml"))
  data_xi <- read.table(file=paste0(path, "/../Results/Univariate_real_Data_",location_name,"/data_xi.xml"))
  data_parameters <- readRDS(paste0(path, "/../Results/Univariate_real_Data_",location_name,"/data_parameters.xml"))
  
  data <- c(data_q_alpha, data_lns_beta, data_xi)
  data <- as.data.frame(data)
  
  print("results")
  for (name in names(data)){
    cat(name, mean(get(name, data)),"\n")
  }
  
  
  #parameters:
  data_set <- data
  data_set <- as.data.frame(data_set)
  parameter_names <- paste0(location_name,".",names(data_set))
  
  prior_x <- list("q"=seq(min(data$q),max(data$q),0.1),
                  "lns"=seq(min(data$lns),max(data$lns),0.01),
                  "xi"=seq(-10,10,0.001))
  
  prior_y <- list("q"=dnorm(prior_x$q, 100, 50),
                  "lns"=dnorm(prior_x$lns, 0, 10),
                  "xi"=dnorm(prior_x$xi, 0, 10))
  
  prior_x$xi <- data_parameters$xi_limits[1]+(data_parameters$xi_limits[2]-data_parameters$xi_limits[1])*exp(prior_x$xi)/(1+exp(prior_x$xi))
  names(data_set) <- parameter_names
  names(prior_x) <- parameter_names
  names(prior_y) <- parameter_names
  
  plots <- map(parameter_names, ~hist_plot_no_real_value_with_priors(data_set, .x, prior_x, prior_y))
  ggarrange(plotlist=plots)+
    ggsave(paste0(path, "/../Plots/Univariate_real_Data_",location_name,"/hist_with_real.pdf"),
           width = 5, height = 4, units = c("in"))
  return(plots)
}

univariate_real_data(Helgeroa_yearly_max_data$yearly_max, "Helgeroa")
univariate_real_data(Oslo_yearly_max_data$yearly_max, "Oslo")
univariate_real_data(Oscarsborg_yearly_max_data$yearly_max, "Oscarsborg")
univariate_real_data(Viker_yearly_max_data$yearly_max, "Viker")

plots_Helgeroa <- load_real_results_function("Helgeroa")
plots_Oscarsborg <- load_real_results_function("Oscarsborg")
plots_Oslo <- load_real_results_function("Oslo")
plots_Viker <- load_real_results_function("Viker")

plots <- c(plots_Helgeroa, plots_Oscarsborg, plots_Oslo, plots_Viker)
ggarrange(plotlist=plots[c(3,6,9,12)],ncol=4,nrow=1)+
  ggsave(paste0(path, "/../Plots/Univariate_real_Data/hist_with_real_xi.pdf"),
         width = 10, height = 8*1/4, units = c("in"))
ggarrange(plotlist=plots[c(1,4,7,10)],ncol=4,nrow=1)+
  ggsave(paste0(path, "/../Plots/Univariate_real_Data/hist_with_real_q.pdf"),
         width = 10, height = 8*1/4, units = c("in"))
ggarrange(plotlist=plots[c(2,5,8,11)],ncol=4,nrow=1)+
  ggsave(paste0(path, "/../Plots/Univariate_real_Data/hist_with_real_lns.pdf"),
         width = 10, height = 8*1/4, units = c("in"))
