rm(list=ls())

datafile <- "k2-normal-mixture-fixed-sigma-m2-d2025-04-06.csv"
RLCT_data <- as.matrix(read.table(datafile, sep= ",",header=TRUE))

ggplot(RLCT_data, aes(x=factor(chain_size), y=RLCT, color=factor(chain_size)))+
  geom_boxplot()+
  geom_hline(yintercept = 3/4, color="red")+
  facet_wrap(~n)