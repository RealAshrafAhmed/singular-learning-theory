rm(list=ls())
library(rstan)
library(data.table)
# library(bayesplot)
source("./two-components-known-sigma/globals.R")
source(paste0(basedir, "/fit.R"))

# *****************************************
# Running a single chain to spot check the model and run some diagnostics
# *****************************************
# data needed, this doesn't impact the experiments themselves, this is only for spot checking
n=25
data = rnorm(n, 0, 1)

fit = fitstan(data=data,
              beta=1,
              model=model,
              size=4000, # this is set based on chain sim study
              warmup=4000, # be careful of lower values, they impact chain congerence a lot
              chains=2)

# Let visualize the chain and run some diagnostics
# general stan dia  gnostics
check_all_diagnostics(fit$fit)

# diagnosing chain mixing
stan_ac(fit$fit, par=c("mu[1]","mu[2]", "rho"), lags=50, nrow=3)
traceplot(fit$fit,
          pars = c("mu[1]", "mu[2]", "rho"),
          inc_warmup = TRUE,
          nrow = 3,
          window=c(4000, 8000)) # let's only look at the post warmup

model_par = c("mu[1]", "mu[2]", "rho")
fit_np <- nuts_params(fit$fit)

mcmc_pairs(
  as.array(fit),
  np = fit_np,
  pars = model_par,
  off_diag_args = list(size = 0.75)
)

draws.posterior <- extract(fit$fit, permuted=TRUE)

# checking for divergences
div_style <- parcoord_style_np(div_size = 0.05, div_alpha = 0.4)
mcmc_parcoord(
  as.array(fit),
  np = fit_np,
  pars = c("mu[1]", "mu[2]", "rho"),
  size = 0.25,
  alpha = 0.1,
  np_style = div_style
)


# mcmc_trace(as.array(fit), pars = c("mu[1]", "mu[2]", "rho"), np = fit_np, nrow=3)

# checking for divergence with equal units
mcmc_parcoord(
  as.array(fit),
  transform = function(x) {(x - mean(x)) / sd(x)},
  pars = c("mu[1]", "mu[2]", "rho"),
  size = 0.25,
  alpha = 0.1,
  np = fit_np,
  np_style = div_style
)

# tunning diagnostics
stan_par(fit, par=c("mu[1]"), chain = 1)

# static scatter 3d plot of the chain
library(scatterplot3d)
s3d <- scatterplot3d(x = draws.posterior$mu[,1],
                     y = draws.posterior$mu[,2],
                     z = draws.posterior$rho,
                     xlab = "mu[1]",
                     ylab = "mu[2]",
                     zlab = "rho",
                     main = "Posterior MCMC samples",
                     pch = 16,    # Solid circles
                     color = "blue",
                     # angle=65,
                     cex.symbols=.2,
                     grid = TRUE,
                     box = TRUE)