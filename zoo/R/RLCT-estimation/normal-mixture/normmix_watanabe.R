rm(list=ls())
library(rstan)
library(data.table)
library(ggplot2)

options(mc.cores=parallel::detectCores()) #init number of cores

model_code_str = "
functions {
  real mixture_lpdf(vector x, real rho, vector mu, real sigma1, real sigma2) {
    return log_sum_exp(log(rho)+normal_lpdf(x | mu[1], sigma1),
                       log(1-rho)+normal_lpdf(x | mu[2], sigma2));
  }
}
data {
  int<lower=0> n;               // number of samples 
  real<lower=0> beta;  //inverse temperature
  vector[n] x;                  // observations
}
parameters {
  real<lower=0, upper=1> rho;
  vector[2] mu;
}
model {
  target += beta*mixture_lpdf(x | rho, mu, 1, 1);
  
  // priors
  rho ~ uniform(0,1);
  mu ~ multi_normal([0,0], [[2, 0], [0, 2]]);
}
generated quantities {
  real log_prob_data; //compute the log probability of the data using posterior samples
  log_prob_data = mixture_lpdf(x | rho, mu, 1, 1);
}
"

fitstan = function(data, beta=1, size, model, chains=1, warmup=1000) {
  # data: iid observations from some unknown distribution
  # beta: inverse temperature
  # size: number of samples generates from mc
  # warmup: number of samples used to estimate the covariance used by HMC
  # chains: number of chains generated
  stan_data = list(n=length(data), beta=beta, x=data)
  fit = sampling(model, 
                 data = stan_data,
                 chains = chains, # instead of repeating the fit, just let stan generate all the chains we need for that n_s and beta
                 iter = size+warmup, # always number of samples needed + warmup
                 warmup = warmup,
                 verbose = FALSE,
                 control = list(adapt_delta = 0.99, max_treedepth=12),
                 refresh=1000)
  return(fit)
}

n=50
c_factors = c(1/100, 1/0, 1, 5, 10, 20, 50, 100)*1/log(n)
d=10
Delta_beta = d/log(n)

total_sims = 10
RLCT_estimates = matrix(data=NA,nrow=total_sims, ncol=length(c_factors))
model = stan_model(model_code=model_code_str)

pb = txtProgressBar(min = 0, max = total_sims, initial = 0) 
data = rnorm(n, 0, 1)
for(k in 1:length(c_factors)) {
  c = c_factors[k]
  beta=c/log(n)
  for(i in 1:total_sims) {
    fit1= fitstan(data=data,
                  beta=beta,
                  model=model,
                  size=2000,
                  warmup=2000,
                  chains=4)
    
    log_prob_data1 = extract(fit1, par=c("log_prob_data"))$log_prob_data
    
    # beta2 = beta1+Delta_beta
    # fit2= fitstan(data=data,
    #              beta=beta2,
    #              model=model,
    #              size=5000,
    #              warmup=500,
    #              chains=4)
    # 
    # log_prob_data2 = extract(fit2, par=c("log_prob_data"))$log_prob_data
    #watanabe estimator
    # RLCT_estimates[i,1] = -(mean(log_prob_data2)-mean(log_prob_data1))/(1/beta2-1/beta1)
    
    #Toru estimator
    RLCT_estimates[i,k] = beta1^2*(mean(log_prob_data1^2)-mean(log_prob_data1)^2)
    setTxtProgressBar(pb,i)
  }
}

RLCT_estdf = data.frame(RLCT=RLCT_estimates)
# ggplot(RLCT_estdf, aes(x=RLCT.1)) + geom_histogram() + 
#   geom_vline(xintercept = 3/4, color="red") +
#   geom_vline(xintercept = mean(RLCT_estimates), color="lightblue") +
#   ggtitle(sprintf("n=%s, c=%s, d=%s, simulations=%s", n, c, d, total_sims))

ggplot(RLCT_estdf, aes(x=RLCT.2)) + geom_histogram() + 
  geom_vline(xintercept = 3/4, color="red") +
  geom_vline(xintercept = mean(RLCT_estimates[,2]), color="lightblue") +
  ggtitle(sprintf("n=%s, c=%s, d=%s, simulations=%s", n, c, d, total_sims))


pairs(fit1, pars = c("mu[1]", "mu[2]", "rho"), las = 1)
traceplot(fit1, pars = c("mu[1]", "mu[2]", "rho"), inc_warmup = TRUE, nrow = 3)
traceplot(fit2, pars = c("mu[1]", "mu[2]", "rho"), inc_warmup = TRUE, nrow = 3)

stan_par(fit1, par=c("mu[1]"), chain = 1)

print(fit1)
stan_diag(fit1, "divergence")
