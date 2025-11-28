rm(list=ls())
library(data.table)
library(ggplot2)

# datafile <- "k2-normal-mixture-fixed-sigma-m2-d2025-04-06.csv"
# RLCT_data <- as.matrix(read.table(datafile, sep= ",",header=TRUE))
# 
# ggplot(RLCT_data, aes(x=factor(chain_size), y=RLCT, color=factor(chain_size)))+
#   geom_boxplot()+
#   geom_hline(yintercept = 3/4, color="red")+
#   facet_wrap(~n)


alphas = seq(0, 5, 0.5)

new_drichilet = function(alpha) {
  fun = function(x) {
    return((x*(1-x))^(1-alpha))
  }
  
  return(fun)
}

raw_data = list()
data_table = data.table(a=numeric(),
                        alpha=numeric(),
                        p_a=numeric())

xticks = seq(0, 1, 0.05)
for(i in 1:length(alphas)) {
  alpha = alphas[i]
  ddrichilet = new_drichilet(alpha)
  for(j in 1:length(xticks)) {
    a = xticks[j]
    data_table = rbind(data_table, data.table(a = a,
                                              p_a = ddrichilet(a),
                                              alpha = alpha))
  }
}

ggplot(data_table, aes(x=factor(a), y=p_a, group=factor(alpha), color=factor(alpha)))+geom_line()
