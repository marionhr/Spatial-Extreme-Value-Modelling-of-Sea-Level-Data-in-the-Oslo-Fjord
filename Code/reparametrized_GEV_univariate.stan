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
  
  real l_xi(real x, real xi){
    //used in reparametrization
    return (-log(x))^(-xi);
  }
  
  real l(real x){
    //used in reparametrization
    return log(-log(x));
  }
  
  vector reparametrization(real q_alpha, real lns_beta, real xi, real alpha, real beta){
    //input:
    //q_alpha, lns_beta: quantile based parameters of the GEV
    //xi: parameter of the GEV
    //alpha: quantile related to q
    //beta: quantile related to lns
    //output:
    //GEV parameters mu and sigma
    real mu;
    real lnsigma;
    vector[2] ret;
    if (xi!=0){
      mu = q_alpha - exp(lns_beta)*(l_xi(alpha,xi)-1)/(l_xi(1-beta/2,xi)-l_xi(beta/2,xi));
      lnsigma = log(exp(lns_beta)*xi/(l_xi(1-beta/2,xi)-l_xi(beta/2,xi)));
    }
    else{
      mu = q_alpha + exp(lns_beta)*(l(alpha)-1)/(l(beta/2)-l(1-beta/2));
      lnsigma = log(exp(lns_beta)/(l(beta/2)-l(1-beta/2)));
    }
    ret[1]=mu;
    ret[2]=lnsigma;
    return ret;
  }
}


data {
  int<lower=0> N; //Number of data points
  vector[N] y;    //The data
  
  // Quantiles for location and spread in GEV-reparametrisation
  real<lower=0,upper=1> alpha;
  real<lower=0,upper=0.5> beta;
  
  vector[2] xi_limits;
}

parameters {
  real q;
  real lns;
  real xi_param;
}

transformed parameters {
  real xi;
  
  xi = xi_limits[1] + (xi_limits[2]-xi_limits[1])*exp(xi_param)/(1+exp(xi_param));
}

model {
  vector[2] ret;
  
  // Priors
  q ~ normal(100,50);
  lns ~ normal(0,10);
  xi_param ~ normal(0,10);
  
  // Likelihood
  ret = reparametrization(q, lns, xi, alpha, beta);  
  
  for(n in 1:N){
    target +=  gev_lpdf(y[n] | ret[1], ret[2], xi);
  }
}
