rm(list=ls())

# datafile <- "k2-normal-mixture-fixed-sigma-m2-d2025-04-06.csv"
# RLCT_data <- as.matrix(read.table(datafile, sep= ",",header=TRUE))
# 
# ggplot(RLCT_data, aes(x=factor(chain_size), y=RLCT, color=factor(chain_size)))+
#   geom_boxplot()+
#   geom_hline(yintercept = 3/4, color="red")+
#   facet_wrap(~n)


alpha = seq(0, 2, 0.5)
new_drichilet = function(alpha) {
  fun = function(x) {
    return((x*(1-x))^(1-alpha))
  }
}

xticks = seq(0, 1, 0.1)
data = matrix(data=NA, ncol=length(alpha)+1, nrow=length(xticks))
data[,1] = xticks
for(i in 1:length(alpha)) {
  data[,i+1] = new_drichilet(alpha[i])(xticks) 
}


