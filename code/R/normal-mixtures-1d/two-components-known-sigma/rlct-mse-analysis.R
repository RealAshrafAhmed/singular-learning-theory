rm(list=ls())
library(ggplot2)
# library(ggh4x)
library(data.table)
library(ggthemes)
source("./two-components-known-sigma/globals.R")

datafile <- paste0(basedir, "/data/estimates/rlct-m1.csv")
rlctdf <- read.table(datafile, sep= ",",header=TRUE)

cfactors = unique(rlctdf$c)
# cfactors = c(1., 2.)
nfactors = unique(rlctdf$n)
max_chainsize = max(rlctdf$chain_size)

mse_data <- data.table(
  c=numeric(),
  n=integer(),
  bias=numeric(),
  variance=numeric(),
  mse=numeric()
)

for(c in cfactors) {
  for(n in nfactors) {
    match = rlctdf[rlctdf$n==n & rlctdf$c == c & rlctdf$chain_size==max_chainsize,]
    estimates = match$RLCT
    bias = mean(estimates)-0.75
    variance = mean(estimates^2) - (mean(estimates))^2
    mse = bias + variance
    mse_data=rbind(mse_data, data.table(
      c=c,
      n=n,
      bias=bias,
      variance=variance,
      mse=mse
    ))
  }
}


ggplot(mse_data, aes(x = factor(n), y = factor(c), fill = bias)) +
  geom_tile() + 
  theme_wsj() + scale_fill_gradient2(low = "blue", mid = "white", high = "red")


ggplot(mse_data, aes(x=factor(n), y=mse, group=factor(c), color=factor(c)))+
  geom_line() +scale_y_log10() + theme_hc()

ggplot(mse_data, aes(x=factor(n), y=variance, group=factor(c), color=factor(c)))+
  geom_line() + theme_hc()

ggplot(mse_data, aes(x=factor(n), y=bias, group=factor(c), color=factor(c)))+
  geom_line() + theme_hc()
