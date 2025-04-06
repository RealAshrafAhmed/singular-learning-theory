rm(list=ls())

RLCT_estimates_file <- "k2-normal-mixture-fixed-sigma-m2-d2025-04-06.csv"
RLCT_estimates <- as.matrix(read.table(RLCT_estimates_file, sep= ",",header=TRUE))

ggplot(RLCT_estimates, aes(x=factor(chain_size), y=RLCT, color=factor(chain_size)))+
  geom_boxplot()+
  geom_hline(yintercept = 3/4, color="red")+
  facet_wrap(~n)