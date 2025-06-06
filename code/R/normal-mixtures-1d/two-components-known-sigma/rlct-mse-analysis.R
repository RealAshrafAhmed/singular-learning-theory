rm(list=ls())
library(ggplot2)
library(ggpubr)
# library(ggh4x)
library(data.table)
library(ggthemes)
library(dplyr)
source("./two-components-known-sigma/globals.R")

datafile <- paste0(basedir, "/data/estimates/rlct-m50.csv")
rlctdt <- read.table(datafile, sep= ",",header=TRUE)

rlctdt <- rlctdt[rlctdt$c <= 2.5 & rlctdt$n < 600, ]

cfactors = unique(rlctdt$c)
# cfactors = c(1., 2.)
nfactors = unique(rlctdt$n)
max_chainsize = max(rlctdt$chain_size)


mse_data <- data.table(
  c=numeric(),
  n=integer(),
  bias=numeric(),
  bias_sd=numeric(),
  variance=numeric(),
  # variance_sd=numeric(),
  mse=numeric()
  # mse_sd=numeric()
)

for(c in cfactors) {
  for(n in nfactors) {
    match = rlctdt[rlctdt$n==n & rlctdt$c == c & rlctdt$chain_size==max_chainsize,]
    estimates = match$RLCT
    # abort()
    # bias = mean(estimates-0.75)
    bias = estimates-0.75
    # bias_sd = sd(bias)
    variance = mean(estimates^2) - (mean(estimates))^2
    mse = mean(bias)^2 + variance
    mse_data=rbind(mse_data, data.table(
      c=c,
      n=n,
      bias=mean(bias),
      bias_sd=sd(bias),
      variance=variance,
      mse=mse
    ))
  }
}


ggplot(mse_data, aes(x = factor(n), y = factor(c), fill = bias)) +
  geom_tile() + 
  theme_wsj() + scale_fill_gradient2(low = "blue", mid = "white", high = "red")


# mse <- ggplot(mse_data, aes(x=factor(n), y=mse, group=factor(c), color=factor(c)))+
#   geom_line() +scale_y_log10() + theme_hc()

vplot <- ggplot(mse_data, aes(x=factor(n), y=variance, group=factor(c), color=factor(c)))+
  geom_line() + theme_hc() +
  ggtitle("variance")

bplot <- ggplot(mse_data, aes(x=factor(n), y=bias, group=factor(c), color=factor(c)))+
  geom_line() + theme_hc() +
  ggtitle("bias") + ylim(-0.08, 0.02)

mplot <- ggplot(mse_data, aes(x=factor(n), y=mse, group=factor(c), color=factor(c)))+
  geom_line() + theme_hc() +
  ggtitle("mse")

splot <- ggplot(rlctdt[rlctdt$n==100,], aes(x=factor(c), y=RLCT, color=factor(c), kind=factor(trial)))+
  geom_boxplot() + theme(legend.position="none")

together <- ggarrange(
  bplot, vplot, mplot,
  # labels=c("bias", "variance", "mse"),
  common.legend = TRUE, legend = "right", ncol=3
)

ggarrange(
  splot, 
  together,
  nrow=2
)


# plot<- ggarrange(ba,mi,fa, ncol=3, nrow=1, common.legend = TRUE,legend="bottom")

# +
#   ggtitle(
#     "Estimated RLCT MSE for different sample size and inverse temperature factors",
#     subtitle="true=standard normal, model=two component normal mixture with known var, samples=20, number of chains per trial=50.")
# annotate_figure(together, top = text_grob("Dive depths (m)", 
#                                       color = "red", face = "bold", size = 14))

ggplot(rlctdt, aes(x=factor(c), y=RLCT, color=factor(c), kind=factor(trial)))+
  geom_boxplot() + theme(legend.position="none") +
  facet_wrap(~factor(n))



rlctdf = setDT(rlctdt)
grouped_df_dplyr <- rlctdf[, 
                           .(count_other = .N), 
                           by = .(n, c), 
                           .SDcols = c("trial", "n", "c")
]