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
}
parameters {
  real rho_raw;
  vector[2] mu;
}
transformed parameters {
 real rho = inv_logit(rho_raw);
}
model {
  // priors
  rho_raw ~ normal(0, 2);
  mu ~ normal(0, 2);
  
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
