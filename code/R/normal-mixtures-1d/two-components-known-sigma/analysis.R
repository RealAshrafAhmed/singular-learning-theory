rm(list=ls())
library(ggplot2)
library(data.table)
source("./two-components-known-sigma/globals.R")

datafile <- paste0(basedir, "/data/k2-normal-mixture-fixed-sigma-m1.csv")
RLCT_data <- as.matrix(read.table(datafile, sep= ",",header=TRUE))

ggplot(RLCT_data, aes(x=factor(c), y=RLCT, color=factor(chain_size)))+
  geom_boxplot()+
  geom_hline(yintercept = 3/4, color="red")+
  facet_wrap(~factor(n))

datafile <- paste0(basedir, "/data/gfe-estimates.csv")
fe_data <- as.data.table(read.table(datafile, 
                                sep= ",",
                                header=TRUE,
                                # colClasses = c(integer(), numeric(), numeric(), integer(), integer(), integer(), character()
                                ))
n_factors = unique(fe_data$n)

# compute the true leading orders of FE for redline plotting
fe_true = rep(0, length(n_factors))
for(i in 1:length(n_factors)) {
  n = n_factors[i]
  print(n)
  observationsfile <- paste0(basedir, "/data/observations-n", n, ".csv")
  data = as.matrix(read.table(observationsfile, sep= ",",header=TRUE))[,1]
  Sn = -1*mean(dnorm(data, mean=0, sd=1, log=TRUE))
  fe_true[i] = n*Sn + 3/4*log(n)
}


ggplot(fe_data, aes(x=factor(chain_size), y=FE, color=name))+geom_boxplot()+
  # geom_hline(yintercept = 3/4, color="red")+
  geom_hline(data = data.frame(yintercept=fe_true, n=n_factors), aes(yintercept = yintercept), color = "red", linetype = "dashed") +
  facet_wrap(~factor(n))
