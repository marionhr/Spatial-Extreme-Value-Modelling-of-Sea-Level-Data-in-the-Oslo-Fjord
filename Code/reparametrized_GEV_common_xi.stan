functions {
  // Define gev density function

  real gev_lpdf(real y, real mu, real lnsigma, real xi) {
    //input:
    //y: data
    //mu, lnsigma, xi: parameters of GEV
    //output:
    //log of the GEV probability distribution function
    
    if (1+xi*(y-mu)/exp(lnsigma) <= 0)
      reject("the gev function does not exist for the parameters");
    if (xi!=0)
      return log(1/exp(lnsigma)*(1+xi*(y-mu)/exp(lnsigma))^(-1/xi-1)*exp(-(1+xi*(y-mu)/exp(lnsigma))^(-1/xi)));
    else
      return log(1/exp(lnsigma)*exp(-(y-mu)/exp(lnsigma))^(xi+1)*exp(-(exp(-(y-mu)/exp(lnsigma)))));
  }
  
  real gev_vec_lpdf(vector y, vector mu, vector lnsigma, real xi, int N){
    //input:
    //y: data
    //mu, lnsigma, xi: parameters of GEV
    //N: length of data vector
    //output:
    //sum over the data of the log of the GEV probability distribution function
    //equals: log of the product over data of the GEV probability distribution function
    
    vector[N] ret;
    for (i in 1:N){
      ret[i]=gev_lpdf(y[i] | mu[i], lnsigma[i], xi);
    }
    return sum(ret);
  }
  
  real l_xi(real x, real xi){
    //used in reparametrization
    return (-log(x))^(-xi);
  }
  
  real l(real x){
    //used in reparametrization
    return log(-log(x));
  }
  
  matrix reparametrization(vector q_alpha, vector lns_beta, real xi, real alpha, real beta, int dim){
    //input:
    //q_alpha, lns_beta: quantile based parameters of the GEV
    //xi: parameter of the GEV
    //alpha: quantile related to q
    //beta: quantile related to lns
    //dim: length of parameter vectors
    //output:
    //GEV parameters mu and sigma
    
    vector[dim] mu;
    vector[dim] lnsigma;
    matrix[dim, 2] ret;
    if (xi!=0){
      mu = q_alpha - exp(lns_beta)*(l_xi(alpha,xi)-1)/(l_xi(1-beta/2,xi)-l_xi(beta/2,xi));
      lnsigma = log(exp(lns_beta)*xi/(l_xi(1-beta/2,xi)-l_xi(beta/2,xi)));
    }
    else{
      mu = q_alpha + exp(lns_beta)*(l(alpha)-1)/(l(beta/2)-l(1-beta/2));
      lnsigma = log(exp(lns_beta)/(l(beta/2)-l(1-beta/2)));
    }
    ret[,1]=mu;
    ret[,2]=lnsigma;
    return ret;
  }
}



data {
  int<lower=0> dim;                    //number of locations
  int<lower=0> N;                      //number of observations
  vector[N] x;                         //observations
  int location_indices[N];             //which observations belong to which locations
  
  // Quantiles for location and spread in GEV-reparametrisation
  real<lower=0,upper=1> alpha;
  real<lower=0,upper=0.5> beta;
  
  vector[2] xi_limits;                 // limits for the parameter xi
}


parameters {
  //hyperparameters
  real xi_param;
  // parameters for GEV
  vector[dim] q;
  vector[dim] lns;
}


transformed parameters {

  // Parameter for GEV
  real xi;
  
  // specify parameters for GEV
  xi = xi_limits[1] + (xi_limits[2]-xi_limits[1])*exp(xi_param)/(1+exp(xi_param));
}

model {
  // matrix for mu and sigma for the GEV distribution
  matrix[dim,2] parameters_standard_gev;

  // Priors for the hyperparameters
  q ~ normal(100,50);
  lns ~ normal(0,10);
  xi_param ~ normal(0,10);
  
  // Likelihood
  parameters_standard_gev = reparametrization(q, lns, xi, alpha, beta, dim);
  x ~ gev_vec(parameters_standard_gev[location_indices,1], parameters_standard_gev[location_indices,2], xi, N);
}





