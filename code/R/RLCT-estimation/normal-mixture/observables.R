rm(list=ls())
library(rstan)
library(data.table)
# library(bayesplot)
# library(threejs)
# source("lib/stan_utility.R")

set.seed(333) # ensure reproducible results.

# set the parallelism based on # of cores of the machine
options(mc.cores=parallel::detectCores())

fitstan = function(data, testdata, alpha, beta=1, sigma=1, size, model, chains=1, warmup=1000) {
  # function to fit a finite mixture normal model.
  #   the algorithm used is the default NUTS (an adaptive HMC)
  #
  # vector<real> data : iid observations from some unknown distribution
  # real beta         : inverse temperature
  # int k             : number of components to mix 
  # int size          : number of posterior samples post warmup
  # stan_model model  : stan model binary to use
  # int chains        : number of chains to generate in parallel
  # int warmup        : number of samples used to during warmup phase of the chain
  # return            : returns a fitted stan model (i.e. posterior samples and results of the generated quantities block)
  
  # input data to stan, if you are going to change the stan file data block
  # make sure you change this one too.
  stan_data = list(n=length(data), beta=beta, x=data, alpha=alpha, sigma=sigma)
  
  # posterior sampling time
  fit = sampling(model, 
                 data = stan_data,
                 chains = chains, # instead of repeating the fit, just let stan generate all the chains we need for that n_s and beta
                 iter = size+warmup, # always number of samples needed + warmup
                 warmup = warmup,
                 verbose = FALSE,
                 control = list(adapt_delta = .95), # dont use anything lower than .9, might need tunning if the observations size is large
                 refresh=0) # how often to print, this can be too noise
  return(fit)
}

# Compile the stan model once for all experiments.
# since the model isn't changing at all during the experiments, there is no need to recompile it
# we could turn off stan recompilation but this is more explicit
model = stan_model(file="normal-mixture-observables.stan")

# model specs, change this only if you want to change the model (2 vs 3 normal mixture)
# and the number of chains used. we kept the number of chains at this level since it is more 
# complicated to and time consuming to run different chain conditions in a single run
k=2 # number of components
m=2 # number of parallel chains to compute the estimator \hat{\lambda}^m

# *****************************************
# Running a single chain to spot check the model and run some diagnostics
# *****************************************
# data needed, this doesn't impact the experiments themselves, this is only for spot checking
n=50
data = rnorm(n, 0, 1)
testdata = rnorm(100, 0, 1) # test data to evaluate posterior expectations

fit = fitstan(data=data,
              testdata,
              beta=1,
              alpha=1,
              sigma=10,
              model=model,
              size=4000, # this is set based on chain sim study
              warmup=2000, # be careful of lower values, they impact chain congerence a lot
              chains=m)

extract(fit, )
rstan::expose_stan_functions(fit)

w = extract(fit, par=c("rho", "b"), permuted=TRUE)
predictions = matrix(data=NA, nrow=length(testdata), ncol=4)
predictions[,1] = testdata
for(i in 1:length(testdata)) {
  x_tilde = predictions[i,1]
  predictions[i,2] = dnorm(x_tilde, 0, 1, log=TRUE)
  predictive_result = 0
  for(j in 1:length(w$rho)) {
    predictive_result = predictive_result + exp(mixture_lpdf(x_tilde, w$rho[j], w$b[j]))
  }
  predictions[i,3] = predictive_result
  predictions[i,4] = exp(predictions[i,2])*(predictions[i,2]-predictions[i,3])
}

# Let visualize the chain and run some diagnostics
# dynamically building the params based on number of components
mu_par = c("mu")
rho_par = c("rho")
if(k > 1) {
  for(i in 1:k) {
    mu_par[i] = sprintf("mu[%s]", i)
  }
  
  for(i in 1:(k-1)) {
    rho_par[i] = sprintf("rho[%s]", i)
  }
}

model_par = c(mu_par, rho_par)
fit_np <- nuts_params(fit)

# just to make the graph tidy, if we have a lot of params pair, break them on different figures
# if(length(model_par > 5)) {
#   mcmc_pairs(
#     as.array(fit),
#     np = fit_np,
#     pars = c(mu_par, rho_par),
#     off_diag_args = list(size = 0.75)
#   )
# } else {
#   mcmc_pairs(
#     as.array(fit),
#     np = fit_np,
#     pars = mu_par,
#     off_diag_args = list(size = 0.75)
#   )
#   
#   mcmc_pairs(
#     as.array(fit),
#     np = fit_np,
#     pars = rho_par,
#     off_diag_args = list(size = 0.75)
#   )
# }


# fit_draws <- extract(fit, permuted=FALSE)

# static scatter 3d plot of the chain
# s3d <- scatterplot3d(x = draws.posterior$mu[,1],
#                      y = draws.posterior$mu[,2],
#                      z = draws.posterior$rho,
#                      xlab = "mu[1]",
#                      ylab = "mu[2]",
#                      zlab = "rho",
#                      main = "Markov Chain Draws from Posterior (3D)",
#                      pch = 16,    # Solid circles
#                      color = "blue",
#                      grid = TRUE,
#                      box = TRUE)

# scatterplot3js(x = fit_draws[,1,3],
#           y = fit_draws[,1,4],
#           z = fit_draws[,1,1],
#           color = "blue",
#           size = 0.1,
#           xlab = "mu[1]",
#           ylab = "mu[2]",
#           zlab = "rho[1]",
#           main = sprintf("MCMC draws from the posterior of a %s normal mixture, chain index %s", k, 1))

# general stan diagnostics
# check_all_diagnostics(fit)

# diagnosing chain mixing
stan_ac(fit, par=c("mu[1]","mu[2]", "rho[1]"), lags=200, nrow=3)
traceplot(fit,
          pars = c("mu[1]", "mu[2]", "rho[1]"),
          inc_warmup = TRUE,
          nrow = 3,
          window=c(200, 4000)) # let's only look at the post warmup

# checking for divergences
# fit_np <- nuts_params(fit)
# div_style <- parcoord_style_np(div_color = "green", div_size = 0.05, div_alpha = 0.4)
# 
# mcmc_parcoord(
#   as.array(fit),
#   np = fit_np,
#   pars = c(mu_par, rho_par),
#   size = 0.25,
#   alpha = 0.1,
#   np_style = div_style
# )


# mcmc_trace(as.array(fit), pars = "mu[1]", np = fit_np) + xlab("Post-warmup iteration")
# mcmc_trace(as.array(fit), pars = "mu[2]", np = fit_np) + xlab("Post-warmup iteration")
# mcmc_trace(as.array(fit), pars = "rho", np = fit_np) + xlab("Post-warmup iteration")

# checking for divergence with equal units
# mcmc_parcoord(
#   draws,
#   transform = function(x) {(x - mean(x)) / sd(x)},
#   size = 0.25,
#   alpha = 0.1,
#   np = np,
#   np_style = div_style
# )

# tunning diagnostics
# stan_par(fit, par=c("mu[1]"), chain = 1)

# ******************************************
# Time for generating some data for analysis
# ******************************************

# observations sample size
n_factors = c(10, 50, 100, 250, 500, 1000, 2000)
# n_factors = c(50)

# inverse temperature factor, 1 is optimal as per the paper
# c_factors = c(1/10, 1, 1.5, 2, 5, 10)
c_factors = c(1, 1.5, 2, 5, 10) 

# different markov chain size
chain_sizes = seq(2, 10, 2)*1000
total_sims = 5000

# We are going to keep appending data to an output file so we don't have
# to repeat everything from scratch every the process fails.
RLCT_estimates <- data.table(
  trial=integer(),
  RLCT=numeric(),
  beta=numeric(),
  c=numeric(),
  n=integer(),
  chain_size=integer()
)

datafile <- sprintf("k%s-normal-mixture-fixed-sigma-m%s.csv", k, m)
if (file.exists(datafile)) {
  # instead of replacing the file every time, let's resume where we were
  print(sprintf("File %s already exists, let's try to resume from where we ended", datafile))
  # file.remove(RLCT_estimates_file)
  data_sofar <- as.matrix(read.table(datafile, sep= ",",header=TRUE))
  last_n = max(data_sofar[,5])
  last_chain_size = max(data_sofar[,6])
  last_c = max(data_sofar[,4])
  last_trial = max(data_sofar[,1])
  print(sprintf("resume from n=%s, c=%s, chain_size=%s, trial=%s", 
                last_n, 
                last_c, 
                last_chain_size, 
                last_trial))
} else {
  fwrite(RLCT_estimates, 
         file = datafile, 
         col.names=TRUE, row.names = FALSE)
}

print(sprintf("Creatted file %s to store RLCT estimates using a Toru estimator.", datafile))

pb = txtProgressBar(min = 0, max = total_sims*length(c_factors), initial = 0) 

for(n in n_factors) {
  data = rnorm(n, 0, 1) # a new sample per factor and keeping it fixed for the other conditions
  for(j in 1:length(c_factors)) { # iterator per temperature factor
    c = c_factors[j]
    beta=c/log(n)
    for(chain_size in chain_sizes) { # iteration per chain size
      for(i in 1:total_sims) { # repeat for total_sims conditions to approx estimator variance
        if(n <= last_trial & c <= last_c & chain_size <= last_chain_size & i <= last_trial) {
          next; # skip to the latest simulation
        }
        fit = fitstan(data=data,
                      beta=beta,
                      k=k,
                      model=model,
                      size=chain_size,
                      warmup=2000, # be careful when changing this, consult chain size analysis
                      chains=m)    # number of chains to compute \hat{\lambda}^m
        
        draws = extract(fit, 
                        par=c("log_prob_data"), 
                        permuted=FALSE) # must be FALSE so we can extract the individual chains
        
        # if m > 1, we need to compute the estimator per chain and average them later
        # if m == 1, we are just averaging over a single value (i.e. it is just the single value)
        per_chain_estimate = matrix(data=NA, nrow=m, ncol=1)
        for(s in 1:m) { 
          per_chain_estimate[s] = beta^2*(mean(draws[,s,]^2)-mean(draws[,s,])^2)
        }
        
        more_estimates = data.table(trial = i, 
                                    RLCT = mean(per_chain_estimate), 
                                    beta=c/log(n),
                                    c=c,
                                    n=n,
                                    chain_size=chain_size)
        
        fwrite(more_estimates, 
               file = datafile, 
               append = TRUE, col.names = FALSE
        )
        
        #print progress so far
        print(
          sprintf("Updated file with simulation %s/%s for beta=%.3f, chains=%s, c=%s, chain size=%s, and n=%s",
                  i, total_sims, beta, m, c, chain_size, n)
        )
      }
    }
  }
}
print(sprintf("All done with %s", datafile))

# to analyze, open analysis.R

