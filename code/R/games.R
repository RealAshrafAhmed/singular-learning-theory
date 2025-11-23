rm(list=ls())

t = seq(-10,10,0.05)
bump_function = function(t) {
  if(t <= 0) {
    return(0)
  }
  
  return(exp(-1/t))
}

plot(t, 1/t)