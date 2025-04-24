rm(list=ls())
library(rstan)
library(data.table)
source("./two-components-known-sigma/globals.R")
source(paste0(basedir, "/fit.R"))

# ******************************************
# Time for generating some data for analysis
# ******************************************

# observations sample size
n_factors = c(10, 100, 250, 1000)
total_sims = 1000
n_test = 1000

# generate the samples and save them. that way we can rerun the simulations without worry about 
# data variability
for(n in n_factors) {
  datafile <- paste0(basedir, "/data/observations-n", n, "-times", total_sims,".csv")
  if (file.exists(datafile)) { # only create data if they don't already exist
    print(sprintf("Skipping file %s, it already exists!", datafile))
  } else {
    data = matrix(data=rnorm(n*total_sims, mean=0, sd=1), nrow=n, ncol=total_sims)
    fwrite(data.frame(data),
           file = datafile, 
           col.names=TRUE, row.names = FALSE)
  }
}

# generate samples to estimate the generalization error expectation
testfile <- paste0(basedir, "/data/test-n", n_test, ".csv")
if (file.exists(testfile)) { # only create data if they don't already exist
  print(sprintf("Skipping file %s, it already exists!", testfile))
} else {
  fwrite(data.frame(x=rnorm(n=n_test, mean=0, sd=1)),
         file = testfile, 
         col.names=TRUE, row.names = FALSE)
}

predictive_lpdf = function(w_samples) {
  # evaluate the likelihood of x under the parameters
  
  w_fun = function(x) {
    return(function(row) {
      exp(mixture_lpdf(x, rho=row["rho"], mu=c(row["mu[1]"], row["mu[2]"])))
    })
  }
  
  ell = function(x) { # compute the predictive
    return(log(mean(apply(w_samples, 1, FUN=w_fun(x)))))
  }
  
  return(ell)
}

# We are going to keep appending data to an output file so we don't have
# to repeat everything from scratch every the process fails.
estimates <- data.table(
  ctrial=integer(),
  ge=numeric(),
  cn=integer(),
  cchain_size=integer(),
  cname=character()
)

datafile <- paste0(basedir, "/data/ge-estimates.csv")
if (file.exists(datafile)) {
  # since we already have a a file, let's load the data it has and use it to check
  # later so we can skip them
  print(sprintf("File %s already exists, let's try to resume from where we ended", datafile))
  data_sofar <- as.data.table(read.table(datafile, sep= ",",header=TRUE))
} else {
  fwrite(estimates, 
         file = datafile, 
         col.names=TRUE, row.names = FALSE)
  data_sofar <- as.data.table(read.table(datafile, sep= ",",header=TRUE))
}

print(sprintf("Created file %s to store free energy estimates", datafile))


chain_size = 6000
estimators = c("Gn", "Tn", "Cn")

pb = txtProgressBar(min = 0, max = total_sims, initial = 0)

for(i in 1:total_sims) { # repeat for total_sims conditions to approx estimator variance
  for(n in n_factors) {
    for(estimator in estimators) {
      observationsfile <- paste0(basedir, "/data/observations-n", n, "-times", total_sims, ".csv")
      data = as.matrix(read.table(observationsfile, sep= ",",header=TRUE))
      testfile <- paste0(basedir, "/data/test-n", n_test, ".csv")
      test_data = as.matrix(read.table(testfile, sep= ",",header=TRUE))[,1]

      # check if the data has already been generated and saved in the file
      match = data_sofar[ctrial==i & cn==n & cchain_size==chain_size & cname==estimator]
      # print(match)
      if(dim(match)[1]==1){
        print(sprintf("skipping estimator=%s, n=%s, chain_size=%s, trial=%s", estimator, n, chain_size, i))
        next;
      } else if(dim(match)[1] > 1) {
        print(sprintf("%s duplicate entries for n=%s, chain_size=%s, trial=%s", dim(match)[1], n, chain_size, i))
        print("potential data corruption, aborting")
        stop()
      }
      
      more_estimates = NA
      
      if(estimator == "Gn") {
        model_fit = fitstan(data=data[,i],
                            beta=1,
                            model=model,
                            size=chain_size,
                            warmup=4000, # be careful when changing this, consult chain size analysis
                            chains=2)    # number of chains to compute \hat{\lambda}^m
        
        draws = as.data.table(extract(model_fit$fit,
                                      par=c("rho", "mu[1]", "mu[2]"),
                                      permuted=TRUE))
        
        pred_lpdf = predictive_lpdf(draws)
        test_data_lpdf = sapply(test_data, FUN=pred_lpdf)
        true_data_lpdf = dnorm(test_data, mean=0, sd=1, log=TRUE)
        ge = mean(true_data_lpdf)-mean(test_data_lpdf)
        
        # compute the generalization error estimate
        more_estimates = data.table(ctrial = i,
                                    ge = ge,
                                    cn=n,
                                    cchain_size=chain_size,
                                    cname="Gn")
      } else if(estimator == "Tn") {
        model_fit = fitstan(data=data[,i],
                            beta=1,
                            model=model,
                            size=chain_size,
                            warmup=4000, # be careful when changing this, consult chain size analysis
                            chains=2)    # number of chains to compute \hat{\lambda}^m
        
        draws = as.data.table(extract(model_fit$fit,
                                      par=c("rho", "mu[1]", "mu[2]"),
                                      permuted=TRUE))
        pred_lpdf = predictive_lpdf(draws)
        test_data_lpdf = sapply(data[,i], FUN=pred_lpdf)
        
        # compute the generalization error estimate
        more_estimates = data.table(ctrial = i,
                                    ge = mean(-1*test_data_lpdf),
                                    cn=n,
                                    cchain_size=chain_size,
                                    cname="Tn")
      } else if(estimator == "Cn") {
        # leave one out estimate
        # we need to drop a value and refit and take the average
        estimate = rep(0, n)
        for(j in 1:n) {
          model_fit = fitstan(data=data[,i][-j],
                              beta=1,
                              model=model,
                              size=chain_size,
                              warmup=4000, # be careful when changing this, consult chain size analysis
                              chains=1)    # number of chains to compute \hat{\lambda}^m
          
          draws = as.data.table(extract(model_fit$fit,
                                        par=c("rho", "mu[1]", "mu[2]"),
                                        permuted=TRUE))
          pred_lpdf = predictive_lpdf(draws)
          test_data_lpdf = sapply(data[,i], FUN=pred_lpdf)
          estimate[j] = mean(-1*test_data_lpdf)
          print(paste0("estimating cross validation ", n, "\\", j ," dropout out"))
        }
        # compute the generalization error estimate
        more_estimates = data.table(ctrial = i,
                                    ge = mean(estimate),
                                    cn=n,
                                    cchain_size=chain_size,
                                    cname="Cn")
      }
  
      fwrite(more_estimates, 
             file = datafile, 
             append = TRUE, col.names = FALSE
      )
      
      #print progress so far
      print(
        sprintf("Updated file with simulation %s/%s for %s chain size=%s, and n=%s",
                i, total_sims, estimator, chain_size, n)
      )
    }
  }
}
print(sprintf("All done with %s", datafile))

# to analyze, open analysis.R

