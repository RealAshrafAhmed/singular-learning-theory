rm(list=ls())
library(ggplot2)
library(data.table)
source("./two-components-known-sigma/globals.R")

datafile <- paste0(basedir, "/k2-normal-mixture-fixed-sigma-m1.csv")
RLCT_data <- as.matrix(read.table(datafile, sep= ",",header=TRUE))

RLCT_data = RLCT_data[RLCT_data[,5]<500,] #filter out unfinished samples
ggplot(RLCT_data, aes(x=factor(c), y=RLCT, color=factor(chain_size)))+
  geom_boxplot()+
  geom_hline(yintercept = 3/4, color="red")+
  facet_wrap(~factor(n))

RLCT_data_100 = RLCT_data[RLCT_data[,5]==100,] #filter out unfinished samples
ggplot(RLCT_data_100, aes(x=factor(), y=RLCT, color=factor(chain_size)))+
  geom_boxplot()+
  geom_hline(yintercept = 3/4, color="red")+
  labs(title = "Toru RLCT estimates of a two normal mixture (n=100)")

