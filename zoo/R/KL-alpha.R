
library(data.table)
t = seq(0, 5, 0.01)
# plot(t, -log(t)+t-1)

s_alpha = function(alpha) {
  fun = function(t) {
    return (1-t^alpha)/alpha
  }
  
  return(fun)
}

s_data = data.table(
  alpha=numeric(),
  t=numeric(),
  s=numeric()
)
alpha = seq(0, 1, 0.1) #difference values should give us different divergences
for(a in alpha) {
  if(a == 0) {
    s = function(t) {
      return(-log(t)+t-1)
    }
  } else {
    s = s_alpha(a)
  }
  for(i in t) {
    more_data = data.table(alpha=a,
                           t=i,
                           s=s(i))
    s_data = rbind(s_data, more_data)
  }
}

library(ggplot2)

ggplot(s_data, aes(x=t, y=s, color=factor(alpha)))+
  geom_line()