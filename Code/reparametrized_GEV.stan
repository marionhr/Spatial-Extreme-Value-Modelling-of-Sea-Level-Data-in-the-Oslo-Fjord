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
  
  matrix create_corr_mat(matrix dist_mat, real range, int nu) {
    //Matern correlation function
    //input:
    //dist_mat: distance matrix
    //range, nu: parameters for Matern
    
    real kappa = (8*nu)^(1.0/2.0)/range;
    int N = rows(dist_mat);
    matrix[N, N] corr;
    for (i in 1:N) corr[i, i] = 1 + 1e-10;
    for (i in 2:N) {
      for (j in 1:(i - 1)) {
        corr[i,j] = modified_bessel_second_kind(nu, dist_mat[i, j] * kappa) *
          (kappa * dist_mat[i, j]) ^ nu * 2^(1 - nu) / tgamma(nu); 
        corr[j, i] = corr[i, j];
      }
    }
    return corr;
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
  
  int nu;                              //parameter for Matern
  int ncovariates;                     //number of covariates
  
  // Quantiles for location and spread in GEV-reparametrisation
  real<lower=0,upper=1> alpha;
  real<lower=0,upper=0.5> beta;
  
  
  matrix[dim,ncovariates] covariates;  // the covariate matrix
  matrix[dim,dim] distance_matrix;     // the distance matrix
  
  vector[2] xi_limits;                 // limits for the parameter xi
}


parameters {
  //hyperparameters
  real xi_param;
    
  vector[ncovariates] beta_q;
  vector[ncovariates] beta_lns;
  
  real<lower=0> sd_q;
  real<lower=0> range_q;
  real<lower=0> sd_lns;
  real<lower=0> range_lns;

  // standard normal distributions used generate multivariate normal estimates
  vector[dim] std_norm_q;
  vector[dim] std_norm_lns;
}


transformed parameters {
  // parameters for the Gaussian field
  vector[dim] mean_q;
  vector[dim] mean_lns;
  matrix[dim,dim] Sigma_q;
  matrix[dim,dim] Sigma_lns;
  
  // cholesky decomposed correlation matrix used generate multivariate normal estimates
  matrix[dim, dim] L_q;
  matrix[dim, dim] L_lns;
  
  // parameters for GEV
  vector[dim] q;
  vector[dim] lns;
  real xi;
  
  // specify parameters for the Gaussian field
  mean_q = covariates*beta_q;
  mean_lns = covariates*beta_lns;
  Sigma_q = create_corr_mat(distance_matrix, (range_q), nu);
  Sigma_lns = create_corr_mat(distance_matrix, (range_lns), nu);
  
  // cholesky decomposition
  L_q = cholesky_decompose(Sigma_q);
  L_lns = cholesky_decompose(Sigma_lns);
  
  // specify parameters for GEV
  q = mean_q + sd_q * (L_q * std_norm_q);
  lns = mean_lns + sd_lns* (L_lns * std_norm_lns);
  xi = xi_limits[1] + (xi_limits[2]-xi_limits[1])*exp(xi_param)/(1+exp(xi_param));
}

model {
  // matrix for mu and sigma for the GEV distribution
  matrix[dim,2] parameters_standard_gev;

  // Priors for the hyperparameters
  xi_param ~ normal(0,10);
  
  beta_q ~ normal(0, 100);
  sd_q ~ gamma(5.0/4.0, 1.0/4.0);
  range_q ~ gamma(5, 1.0/12.0);
  //range_q ~ gamma(1.5, 1.0/40.0); //different prior
  
  beta_lns ~ normal(0,100);
  sd_lns ~ gamma(5.0/4.0, 1.0/4.0);
  range_lns ~ gamma(5, 1.0/12.0);
  //range_lns ~ gamma(1.5, 1.0/40.0); //different prior
  
  //gaussian field
  std_norm_q ~ std_normal();
  std_norm_lns ~ std_normal();
  
  // Likelihood
  parameters_standard_gev = reparametrization(q, lns, xi, alpha, beta, dim);
  x ~ gev_vec(parameters_standard_gev[location_indices,1], parameters_standard_gev[location_indices,2], xi, N);
}



