library(MASS) 
data(galaxies)

X = galaxies/1000

library(mclust, quietly=TRUE)
fit = Mclust(X, G=4, Model="V")

plot(fit, what="density", main="", xlab="Velocity (Mm/s)") 
rug(X)

library(sBIC)
gMix = GaussianMixtures(maxNumComponents=10, phi=1, restarts=100)

set.seed(1234)
m = sBIC(X, gMix) 
print(m)

matplot( cbind(m$BIC-m$BIC[1],m$sBIC-m$sBIC[1]), 
         pch=c(1,3), 
         col="black", 
         xlab="Numberofcomponents", ylab=expression(BIC-BIC(M[1])), las=1,xaxt="n" ) axis(1,at=1:10) legend("topleft", c(expression(BIC), expression(bar(sBIC)[1])), pch=c(1,3), y.intersp=1.2)

