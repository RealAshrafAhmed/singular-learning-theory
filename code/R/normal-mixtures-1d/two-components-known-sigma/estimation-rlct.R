rm(list=ls())
library(parallel)
library(rstan)
library(data.table)
source("./two-components-known-sigma/globals.R")
source(paste0(basedir, "/fit.R"))
library(tictoc) # to time things


# ******************************************
# Time for generating some data for analysis
# ******************************************

# observations sample size
# n_factors = c(10, 50, 100, 250, 500, 1000, 2000)
n_factors = c(10, 50, 100, 150, 200, 250, 300)#, 400, 500, 600, 700, 800, 900, 1000)
# n_factors = c(400, 500, 600, 700, 800, 900, 1000)#, 400, 500, 600, 700, 800, 900, 1000)
total_sims = 100
total_chains_per_sim = 50

# generate the samples and save them. that way we can rerun the simulations without worry about 
# data variability
for(n in n_factors) {
  datafile <- paste0(basedir, "/data/observations/n", n, "-times", total_sims,".csv")
  if (file.exists(datafile)) { # only create data if they don't already exist
    print(sprintf("Skipping file %s, it already exists!", datafile))
  } else {
    data = matrix(data=rnorm(n*total_sims, mean=0, sd=1), nrow=n, ncol=total_sims)
    fwrite(data,
           file = datafile, 
           col.names=TRUE, row.names = FALSE)
  }
}


# inverse temperature factor, 1 is optimal as per the paper
# c_factors = c(1/10, 1, 1.5, 2, 5, 10)
c_factors <- c(.5, .75, 1, 1.25, 1.5, 1.75, 2, 2.5)#, 5, 10)

# different markov chain size
chain_sizes = c(8000)
# m=1

# We are going to keep appending data to an output file so we don't have
# to repeat everything from scratch every the process fails.
RLCT_estimates <- data.table(
  trial=integer(),
  RLCT=numeric(),
  beta=numeric(),
  c=numeric(),
  n=integer(),
  chain_index=integer(),
  chain_size=integer()
)

datafile <- paste0(basedir, "/data/estimates/rlct-m", total_chains_per_sim ,".csv")
if (file.exists(datafile)) {
  # since we already have a a file, let's load the data it has and use it to check
  # later so we can skip them
  print(sprintf("File %s already exists, let's try to resume from where we ended", datafile))
  data_sofar <- as.matrix(read.table(datafile, sep= ",",header=TRUE))
} else {
  fwrite(RLCT_estimates, 
         file = datafile, 
         col.names=TRUE, row.names = FALSE)
  data_sofar <- as.matrix(read.table(datafile, sep= ",",header=TRUE))
}

print(sprintf("Created file %s to store RLCT estimates using a Toru estimator.", datafile))

pb = txtProgressBar(min = 0, max = total_sims*length(c_factors), initial = 0) 

for(i in 1:total_sims) { # repeat for total_sims conditions to approx estimator variance
  tic(paste0("trial=", i))
  for(j in 1:length(c_factors)) { # iterator per temperature factor
    tic(paste0("c", c_factors[j]))
    # mclapply(n_factors, function(n) {
      # print(paste0("handling n", n))
    for(n in n_factors) {
      observationsfile <- paste0(basedir, "/data/observations/n", n, "-times", total_sims,".csv")
      data = as.matrix(read.table(observationsfile, sep= ",",header=TRUE))[,i]
      c = c_factors[j]
      beta=c/log(n)
      for(chain_size in chain_sizes) { # iteration per chain size
        # check if the data has already been generated and saved in the file
        # print(paste0("checking if we've already collected data before trial=", i, " n=", n))
        match = data_sofar[data_sofar[,1]==i
                           & data_sofar[,4]==c
                           & data_sofar[,5]==n
                           & data_sofar[,6]==total_chains_per_sim
                           & data_sofar[,7]==chain_size]
        if(length(match)==7){
          print(sprintf("skipping n=%s, c=%s, chain_size=%s, trial=%s", n, c, chain_size, i))
          next;
        } else if(length(match) > 7) {
          print(sprintf("%s duplicate entries for n=%s, c=%s, chain_size=%s, trial=%s", length(match), n, c, chain_size, i))
          print("potential data corruption, aborting")
          stop()
        }
        # print("No match, let's do it!.")
        tic(paste0("n", n))
        
        model_fit = fitstan(data=data,
                      beta=beta,
                      model=model,
                      size=chain_size,
                      warmup=4000, # be careful when changing this, consult chain size analysis
                      chains=total_chains_per_sim)    # number of chains to compute \hat{\lambda}^m
        
        draws = extract(model_fit$fit, 
                            par=c("log_prob_data"), 
                            permuted=FALSE) # must be FALSE so we can extract the individual chains
        
        # if m > 1, we need to compute the estimator per chain and average them later
        # if m == 1, we are just averaging over a single value (i.e. it is just the single value)
        per_chain_estimate = matrix(data=NA, nrow=1, ncol=1)
        for(s in 1:total_chains_per_sim) { 
          per_chain_estimate = beta^2*(mean(draws[,s,]^2)-mean(draws[,s,])^2)
          more_estimates = data.table(trial = i, 
                                      RLCT = per_chain_estimate,
                                      beta=c/log(n),
                                      c=c,
                                      n=n,
                                      chain_index=s,
                                      chain_size=chain_size)
          
          fwrite(more_estimates, 
                 file = datafile, 
                 append = TRUE, col.names = FALSE
          )
        }
        
        #print progress so far
        print(
          sprintf("Updated file with simulation %s/%s for beta=%.5f, c=%s, total_chains_per_sim=%s, chain size=%s, and n=%s",
                  i, total_sims, beta, c, total_chains_per_sim, chain_size, n)
        )
        toc()
      }
    }
      # , mc.cores = detectCores()-2)
    toc()
  }
  toc()
}
print(sprintf("All done with %s", datafile))

# to analyze, open analysis.R

