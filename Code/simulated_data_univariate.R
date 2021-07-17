#!/usr/bin/env Rscript

set.seed(10)

#Find parameter estimates for the simulated data using the univariate GEV model and STAN

#Data

univariate_simulated_data <- function(station_number){
  #Parameter estimation for the univariate model for the simulated data
  #input: 
  #station_number: the number of the location
  #################################################################
  y <- X_vec[location_indices==station_number]
  N <- length(y)
  
  
  #prepare for STAN:
  GEV_dat <- list(N = N,
                  y = y,
                  alpha = alpha,
                  beta = beta,
                  xi_limits = c(-0.5,0.6))
  
  #fit STAN
  fit <- stan(file = paste0(path, "/reparametrized_GEV_univariate.stan"), 
              data = GEV_dat,
              iter = 2*10^4,
              control = list(adapt_delta = 0.99,
                             max_treedepth = 15))
  
  
  #save results
  data_q_alpha <- rstan::extract(fit, pars = c("q"))
  data_lns_beta <- rstan::extract(fit, pars = c("lns"))
  data_xi <- rstan::extract(fit, pars = c("xi"))
  data_parameters <- GEV_dat
  data_original_parameters <- list(xi = xi,
                                   q = q_alpha[station_number],
                                   lns = lns_beta[station_number])
  
  write.table(data_q_alpha, file=paste0(path, "/../Results/Simulated_Data_univariate/Station_",station_number,"/data_q_alpha.xml"))
  write.table(data_lns_beta, file=paste0(path, "/../Results/Simulated_Data_univariate/Station_",station_number,"/data_lns_beta.xml"))
  write.table(data_xi, file=paste0(path, "/../Results/Simulated_Data_univariate/Station_",station_number,"/data_xi.xml"))
  saveRDS(data_parameters, file=paste0(path, "/../Results/Simulated_Data_univariate/Station_",station_number,"/data_parameters.xml"))
  saveRDS(data_original_parameters, file=paste0(path, "/../Results/Simulated_Data_univariate/Station_",station_number,"/data_original_parameters.xml"))

  
  #Plots:
  stan_trace(fit, inc_warmup = TRUE, nrow=3, 
             pars = c("xi",
                      "q",
                      "lns"))+
    ggsave(paste0(path, "/../Plots/Simulated_Data_univariate/Station_",station_number,"/trace.pdf"),
           width = 5, height = 5, units = c("in"))
  
  
  stan_hist(fit, nrow=3, 
            pars = c("xi",
                     "q",
                     "lns"))+
    ggsave(paste0(path, "/../Plots/Simulated_Data_univariate/Station_",station_number,"/hist.pdf"),
           width = 5, height = 5, units = c("in"))
  
  
  pdf(paste0(path, "/../Plots/Simulated_Data_univariate/Station_",station_number,"/pairs.pdf"))
  pairs(fit, pars = c("xi",
                      "q",
                      "lns",
                      "lp__"))
  dev.off()
  
}



load_simulated_results_function <- function(station_number){
  #load the results and make histogram plots for the parameter estimates
  #station_number: 1,2,3 or 4
  ########################################################################
  
  #load results
  data_q_alpha <- read.table(file=paste0(path, "/../Results/Simulated_Data_univariate/Station_",station_number,"/data_q_alpha.xml"))
  data_lns_beta <- read.table(file=paste0(path, "/../Results/Simulated_Data_univariate/Station_",station_number,"/data_lns_beta.xml"))
  data_xi <- read.table(file=paste0(path, "/../Results/Simulated_Data_univariate/Station_",station_number,"/data_xi.xml"))
  data_parameters <- readRDS(paste0(path, "/../Results/Simulated_Data_univariate/Station_",station_number,"/data_parameters.xml"))
  data_original_parameters <- readRDS(paste0(path, "/../Results/Simulated_Data_univariate/Station_",station_number,"/data_original_parameters.xml"))
  
  data <- c(data_xi, data_q_alpha, data_lns_beta)
  data <- as.data.frame(data)


  #priors
  prior_x <- list("xi"=seq(min(data$xi),max(data$xi),0.01),
                  "q"=seq(min(data$q),max(data$q),0.1),
                  "lns"=seq(min(data$lns),max(data$lns),0.01))
  
  prior_y <- list("xi"=dnorm(prior_x$xi, 0, 10),
                  "q"=dnorm(prior_x$q, 100, 50),
                  "lns"=dnorm(prior_x$lns, 0, 10))
  
  #parameter names for the plots
  parameter_names <- paste0(names(data),".",station_number)  
  
  #rename data frames and lists
  names(data) <- parameter_names
  names(data_original_parameters) <- parameter_names
  names(prior_x) <- parameter_names
  names(prior_y) <- parameter_names
  
  #names for x label in the plots
  plot_names <- list(bquote(xi ~"."~ .(station_number)),paste0("q.",station_number),paste0("lns.",station_number))
  names(plot_names) <- parameter_names
  
  
  #plot histograms
  plots <- map(parameter_names, ~hist_plot_with_priors(data, data_original_parameters, .x, plot_names, prior_x, prior_y))
  return(plots)
}

#find parameter estimates
univariate_simulated_data(1)
univariate_simulated_data(2)
univariate_simulated_data(3)
univariate_simulated_data(4)

#make plots
plots1 <- load_simulated_results_function(1)
plots2 <- load_simulated_results_function(2)
plots3 <- load_simulated_results_function(3)
plots4 <- load_simulated_results_function(4)

plots <- c(plots1,plots2,plots3,plots4)
ggarrange(plotlist=plots[c(1,4,7,10)],ncol=4,nrow=1)+
  ggsave(paste0(path, "/../Plots/Simulated_Data_univariate/hist_with_real_xi.pdf"),
         width = 10, height = 8*1/4, units = c("in"))
ggarrange(plotlist=plots[c(2,5,8,11)],ncol=4,nrow=1)+
  ggsave(paste0(path, "/../Plots/Simulated_Data_univariate/hist_with_real_q.pdf"),
         width = 10, height = 8*1/4, units = c("in"))
ggarrange(plotlist=plots[c(3,6,9,12)],ncol=4,nrow=1)+
  ggsave(paste0(path, "/../Plots/Simulated_Data_univariate/hist_with_real_lns.pdf"),
         width = 10, height = 8*1/4, units = c("in"))

