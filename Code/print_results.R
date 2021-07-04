Results1 <- find_data("Simulated_Data_15_locations")
Results2 <- find_data("Simulated_Data_15_locations_different_prior")
original_parameters <- readRDS(paste0(path, "/../Results/Simulated_Data_15_locations/data_original_parameters.xml"))
original_parameters$beta_q <- t(original_parameters$beta_q)
original_parameters$beta_lns <- t(original_parameters$beta_lns)
original_parameters <- as.data.frame(original_parameters[1:9])

param <- names(Results1$data)
names(original_parameters) <- param

print_confidence_interval_and_real_values <- function(name, Results1, Results2, original_parameters){
  vec1 <- get(name, Results1$data)
  vec2 <- get(name, Results2$data)
  orig <- get(name, original_parameters)
  return(c(orig, sort(vec1)[0.025*length(vec1)], sort(vec1)[0.975*length(vec1)], sort(vec2)[0.025*length(vec2)], sort(vec2)[0.975*length(vec2)]))
}

original_parameters


res <- mapply(print_confidence_interval_and_real_values, param, MoreArgs= list(Results1, Results2, original_parameters))
t(res)

library(xtable)

tab<-xtable(t(res), caption= "summary statistics of air pollution data", 
            align=c("|c","|c","|c","|c","|c", "|c|"))       
tab


diff1 <- res[2,]-res[4,]
diff2 <- res[3,]-res[5,]

max(abs(diff1[-c(38,45)]))
max(abs(diff2[-c(38,45)]))

max(abs(diff1[1:15]))
max(abs(diff2[1:15]))

max(abs(diff1[16:30]))
max(abs(diff2[16:30]))

max(abs(diff1[c(31:37,39,44)]))
max(abs(diff2[c(31:37,39,44)]))

max(abs(diff1[c(38,45)]))
max(abs(diff2[c(38,45)]))
