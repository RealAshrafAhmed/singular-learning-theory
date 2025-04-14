functions {
  real mixture_lpdf(real x, real rho, vector mu) {
    return log_sum_exp(log(rho)+normal_lpdf(x | mu[1], 1),
                       log(1-rho)+normal_lpdf(x | mu[2], 1));
  }
}
data {
  int<lower=0> n;      // number of samples 
  vector[n] x;         // observations
  real<lower=0> beta;  // inverse temperature
  real<lower=0> alpha; // generalized parameter for the prior distribution of the mixing parameter
}
parameters {
  real rho;
  vector[2] mu;
}
model {
  // priors
  rho ~ uniform(0,1);
  mu ~ multi_normal([0, 0], [[1, 0], [0, 1]]);
  
  // likelihood
  for (i in 1:n) {
    target += beta*mixture_lpdf(x[i]| rho, mu);
  }
}
generated quantities {
  real log_prob_data = 0; //compute the log probability of the data using posterior samples
  for (i in 1:n) {
    log_prob_data += mixture_lpdf(x[i] | rho, mu);
  }
}
