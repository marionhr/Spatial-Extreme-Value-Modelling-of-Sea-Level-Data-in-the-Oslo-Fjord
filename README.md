# Spatial Extreme Value Modelling of Sea Level Data in the Oslo Fjord

Code for master thesis by Marion Helen Røed

To run the code open the main.R file.\
The main file contains all the necessary code files.\
The path needs to be changed.\
The files functions.R and data.R should always be run first.\
The files labeled #Prepare the sea level data or #Generate Simulated data need to be run before the parameter estimation for real data and simulated data respectivelly.\
Parameter estimation and spatial map construction can be quite time consuming.

To run simulation studies the parameter dim in generate_data_for_multivariate_simulation.R needs to be set to 4, 8 or 15.\
In the study using a different prior for 15 locations, the parameter different_prior in simulated_data.R needs to be set to TRUE and the priors of the range needs to be changed in reparametrized_GEV.stan on lines 175/176 and 179/180.

