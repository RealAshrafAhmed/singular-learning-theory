# Compile the stan model once for all experiments.
# since the model isn't changing at all during the experiments, there is no need to recompile it
# we could turn off stan recompilation but this is more explicit
model = stan_model(file=paste0(basedir, "/model.stan"))  

fitstan = function(data, beta=1, sigma=1, size, model, chains=1, warmup=1000) {
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
  stan_data = list(n=length(data), x=data, beta=beta)
  
  expose_stan_functions(model)
  
  # posterior sampling time
  fit = sampling(model, 
                 data = stan_data,
                 chains = chains, # instead of repeating the fit, just let stan generate all the chains we need for that n_s and beta
                 iter = size+warmup, # always number of samples needed + warmup
                 warmup = warmup,
                 verbose = FALSE,
                 control = list(adapt_delta = .95), # dont use anything lower than .9, might need tunning if the observations size is large
                 refresh=0) # how often to print, this can be too noise
  return(list(fit=fit, data=stan_data))
}

empirical_nll = function(data, par) {
  n = length(data)
  result = rep(0, n)
  for(i in 1:n) {
    result[i] = mixture_lpdf(x=data[i], rho=par["rho"], mu=c(par["mu[1]"], par["mu[2]"]))
  }
  
  return(-1*mean(result))
}