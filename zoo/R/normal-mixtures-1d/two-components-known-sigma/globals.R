basedir <- './two-components-known-sigma'
source("lib/stan_utilities.R")

# set.seed(333) # for reproducability between machines

# set the parallelism based on # of cores of the machine
options(mc.cores=parallel::detectCores())
