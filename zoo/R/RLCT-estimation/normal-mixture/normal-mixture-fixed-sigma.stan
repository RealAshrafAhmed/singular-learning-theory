functions {
  real mixture_lpdf(real x, vector rho, vector mu, real sigma) {
    int k = num_elements(rho);
    vector[k] log_weights;
    for(j in 1:k) {
      log_weights[j] = log(rho[j]) + normal_lpdf(x | mu[j], sigma);
    }
    return log_sum_exp(log_weights);
  }
}
data {
  int<lower=0> n;     // number of samples 
  vector[n] x;        // observations
  real<lower=0> beta; // inverse temperature
  int<lower=1> k;     // number of components to mix
  real<lower=0> sigma;
}
parameters {
  simplex[k] rho;
  vector[k] mu;
}
model {
  // priors
  rho ~ dirichlet(rep_vector(1, k));
  for(j in 1:k) {
    mu[j] ~ normal(0, 2);
  }
  
  // likelihood
  for (i in 1:n) {
    target += beta*mixture_lpdf(x[i]| rho, mu, sigma);
  }
}
generated quantities {
  real log_prob_data = 0; //compute the log probability of the data using posterior samples
  for (i in 1:n) {
    log_prob_data += mixture_lpdf(x[i] | rho, mu, sigma);
  }
}
