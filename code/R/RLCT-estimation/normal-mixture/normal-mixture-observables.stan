functions {
  real dirichletd1_lpdf(real rho, real alpha){
    return (rho*(1-rho))^(alpha-1);
  }
  
  real mixture_lpdf(real x, real rho, real b) {
    return rho*normal_lpdf(x| 0, 1)+(1-rho)*normal_lpdf(x| b, 1);
  }
}
data {
  int<lower=0> n;        // number of samples 
  vector[n] x;           // observations
  real<lower=0> beta;    // inverse temperature
  real<lower=0> sigma;
  real<lower=0> alpha;
}
parameters {
  real<lower=0, upper=1> rho;
  real b;
}
model {
  // priors
  rho ~ dirichletd1(alpha);
  b ~ normal(0, sigma);
  
  // likelihood
  for (i in 1:n) {
    target += beta*mixture_lpdf(x[i]| rho, b);
  }
}
generated quantities {
  real log_prob_data = 0;    //compute the log probability of the data using posterior samples
  for (i in 1:n) {
    log_prob_data += mixture_lpdf(x[i] | rho, b);
  }
}
