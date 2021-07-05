#!/usr/bin/env Rscript

#the path:
path <- "C:/Users/mario/Git/Spatial-Extreme-Value-Modelling-of-Sea-Level-Data-in-the-Oslo-Fjord/Code"
path <- "/home/shomeb/m/marionhr/Spatial-Extreme-Value-Modelling-of-Sea-Level-Data-in-the-Oslo-Fjord/Code"

#libraries:
library(tidyverse)
library(ggplot2)
library(forecast)
library(ggpubr)
library(mev)
library(ismev)
library(extRemes)
library(lubridate)
library(nleqslv)
library(latex2exp)
library(scales)
library(matlib)
library(evd)
library(hexbin)
library(spatstat)
library(patchwork)
library(mvtnorm)
library(copula)
library(evd)
library(dplyr)
library(geosphere)
library(geoR)
library(sf)
library(sp)
library(scoringRules)
library(fhidata)
library(robustbase)
library(raster)
library(rgdal)
library(purrr)
library(gridExtra)
library(patchwork)
library(geoR)
library(Rfast)
library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
library(rgdal) # package for geospatial analysis
library(ggrepel)
library(fields)
library(colorRamps)

#Stan libraries
library("rstan") # observe startup messages
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


#Code:

#Functions that are reused in multiple files
source(paste0(path,"/functions.R"))

#Prepare the sea level data
source(paste0(path,"/data.R"))
source(paste0(path,"/remove_trend.R"))
source(paste0(path,"/missing_data.R"))
source(paste0(path,"/yearly_max_data.R"))

#Parameter estimation for the Real Data
source(paste0(path,"/univariate_real_data.R"))  #Univariate
source(paste0(path,"/real_data_common_xi.R"))   #Common xi
source(paste0(path,"/real_data.R"))             #Multivariate

#Generate Simulated data
source(paste0(path,"/generate_data_for_multivariate_simulation.R"))

#Parameter estimation for the Simulated Data
source(paste0(path,"/simulated_data_univariate.R"))   #Univariate
source(paste0(path,"/simulated_data_common_xi.R"))    #Common xi
source(paste0(path,"/simulated_data.R"))              #Multivariate

#Return level plots
source(paste0(path,"/return_levels.R"))

#Spatial maps of return levels
source(paste0(path,"/map.R"))


#Stan files:
#reparametrized_GEV.stan                    #Parameter estimation for the multivariate/spatial model
#reparametrized_GEV_common_xi.stan          #Parameter estimation for the model where the xi is the only GEV parameter with a spatial representation
#reparametrized_GEV_univariate.stan         #Parameter estimation for the univariate model


#Other files:
#source(paste0(path,"/covariate_data.R"))   #Preparation of the covariate data for analysis
#source(paste0(path,"/norwegian_coast.R"))  #Preparation of GPS coordinates for the Oslo Fjord
#source(paste0(path,"/plots.R"))            #Plots
#source(paste0(path,"/print_results.R"))    #Prints out some useful information used in the thesis



#Install Packages:
#install.packages(c("tidyverse","ggplot2","forecast","ggpubr","mev","ismev","extRemes","lubridate","nleqslv","latex2exp","scales","matlib","evd","hexbin","spatstat","patchwork","mvtnorm","copula","evd","dplyr","geosphere","geoR","sf","sp","scoringRules","fhidata","robustbase","raster","rgdal","purrr","gridExtra","patchwork","geoR","Rfast","ncdf4","raster","rgdal ","ggrepel","fields"))

#Stan:
#Sys.setenv(DOWNLOAD_STATIC_LIBV8 = 1) # only necessary for Linux without the nodejs library / headers
#install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)

#fhidata
#install.packages("devtools")
#require(devtools)
#install_version("fhidata", version = "2019.8.27", repos = "http://cran.us.r-project.org")
