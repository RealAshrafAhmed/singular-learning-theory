rm(list=ls())
library(rstan)
library(data.table)
source("./two-components-known-sigma/globals.R")
source(paste0(basedir, "/fit.R"))

# ******************************************
# Time for generating some data for analysis
# ******************************************

# observations sample size
n_factors = c(10, 100, 1000)

# generate the samples and save them. that way we can rerun the simulations without worry about 
# data variability
for(n in n_factors) {
  datafile <- paste0(basedir, "/data/observations-n", n, ".csv")
  if (file.exists(datafile)) { # only create data if they don't already exist
    print(sprintf("Skipping file %s, it already exists!", datafile))
  } else {
    fwrite(data.frame(x=rnorm(n=n, mean=0, sd=1)),
           file = datafile, 
           col.names=TRUE, row.names = FALSE)
  }
}

# We are going to keep appending data to an output file so we don't have
# to repeat everything from scratch every the process fails.
estimates <- data.table(
  ctrial=integer(),
  cFE=numeric(),
  cbeta=numeric(),
  cm=numeric(),
  cn=integer(),
  cchain_size=integer(),
  cname=character()
)

datafile <- paste0(basedir, "/data/fe-estimates.csv")
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


m_factors = c(-1, 0, 1, 10) # -1 for AFE, 0 for WBIC, >1 for WsBIC
c=1
chain_sizes = c(6000, 10000)
total_sims = 100

pb = txtProgressBar(min = 0, max = total_sims, initial = 0)

for(n in n_factors) {
  observationsfile <- paste0(basedir, "/data/observations-n", n, ".csv")
  data = as.matrix(read.table(observationsfile, sep= ",",header=TRUE))[,1]
  beta=c/log(n)
  for(chain_size in chain_sizes) {
    for(i in 1:total_sims) { # repeat for total_sims conditions to approx estimator variance
      for(m in m_factors) {
        more_estimates = NA
        # check if the data has already been generated and saved in the file
        match = data_sofar[ctrial==i & cm==m & cn==n & cchain_size==chain_size]
        # print(match)
        if(dim(match)[1]==1){
          print(sprintf("skipping m=%s, n=%s, c=%s, chain_size=%s, trial=%s", m, n, c, chain_size, i))
          next;
        } else if(dim(match)[1] > 1) {
          print(sprintf("%s duplicate entries for n=%s, c=%s, chain_size=%s, trial=%s", dim(match)[1], n, c, chain_size, i))
          print("potential data corruption, aborting")
          stop()
        }
        
        chains=1
        if(m > 1) { # need m chains for WsBIC with Lambda^m
          chains=m
        }
        
        model_fit = fitstan(data=data,
                      beta=beta,
                      model=model,
                      size=chain_size,
                      warmup=2000, # be careful when changing this, consult chain size analysis
                      chains=chains)    # number of chains to compute \hat{\lambda}^m

        if(m == -1) { # AFE
          # we don't need draws just the MAP estimator, ideally we compute both AFE(MAP) and AFE(MLE)
          # rstan optimizing api doesn't offer the MLE
          modeL_optim = optimizing(model, model_fit$data)
          AFE = n*empirical_nll(data=data, par=modeL_optim$par) + 3/4*log(n)
          
          more_estimates = data.table(ctrial = i, 
                                      cFE = AFE,
                                      cbeta=beta,
                                      cm=m,
                                      cn=n,
                                      cchain_size=chain_size,
                                      cname="AFE")
        } else if(m == 0) { # WBIC
          draws = extract(model_fit$fit, 
                          par=c("log_prob_data"), 
                          permuted=FALSE) # must be FALSE so we can extract the individual chains
          
          more_estimates = rbind(data.table(ctrial = i,
                                            cFE = -1*mean(draws[,1,]),
                                            cbeta=beta,
                                            cm=m,
                                            cn=n,
                                            cchain_size=chain_size,
                                            cname="WBIC"))
        } else { # m > 0, Naiive BIC
          # compute RLCT estimate
          draws = extract(model_fit$fit, 
                          par=c("log_prob_data"), 
                          permuted=FALSE) # must be FALSE so we can extract the individual chains
          
          per_chain_estimate = matrix(data=NA, nrow=1, ncol=1)
          for(s in 1:m) { 
            per_chain_estimate[s] = beta^2*(mean(draws[,s,]^2)-mean(draws[,s,])^2)
          }
          
          RLCT_hat = mean(per_chain_estimate)
          
          #find MAP estimate and compute the WsBIC
          modeL_optim = optimizing(model, model_fit$data)
          naiveBIC = n*empirical_nll(data=data, par=modeL_optim$par) + RLCT_hat*log(n)
          more_estimates = data.table(ctrial = i, 
                                      cFE = naiveBIC,
                                      cbeta=c/log(n),
                                      cm=m,
                                      cn=n,
                                      cchain_size=chain_size,
                                      cname=sprintf("naiveBIC(%s)", m))
        }

        fwrite(more_estimates, 
               file = datafile, 
               append = TRUE, col.names = FALSE
        )
        
        #print progress so far
        print(
          sprintf("Updated file with simulation %s/%s for beta=%.3f, c=%s, chain size=%s, m=%s, and n=%s",
                  i, total_sims, beta, c, chain_size, m, n)
        )
      }
    }
  }
}
print(sprintf("All done with %s", datafile))

# to analyze, open analysis.R

