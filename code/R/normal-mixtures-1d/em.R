rm(list=ls())

rho = sample(c(0,1), 10, replace=TRUE, prob=c(0.8, 0.2))
prop1 = rnorm(10, 0, 1)
prop2 = rnorm(10, 2, 1)
data = rep(0, 10)
for(i in 1:10) {
  if(rho[i] == 0) {
    data[i] = prop1[i]
  } else {
    data[i] = prop2[i]
  }
}

rho_hat = c(0.5)
for(i in 1:30){
  pi_hat_next = pi_hat[i]*dnorm
}