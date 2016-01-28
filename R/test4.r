# this is a test program 

if( !('isr3.check' %in% ls()) ) { 
  print("loading isr3.so")
  dyn.load('~/src/isr3/src/covarTree.so')
  isr3.check <- T
}

source('~/src/isr3/R/mnp.R')

set.seed(100)
mu <- rnorm(100)
sigma <- rWishart(1,104,diag(100))[,,1]

set.seed(100)
runTime <- proc.time()
for(i in 1:100) r1 <- RMVN(mu,sigma)
print(proc.time() - runTime)

set.seed(100)
runTime <- proc.time()
for(i in 1:100) r2 <- RMVN2(mu,sigma)
print(proc.time() - runTime)

library(mvtnorm)
set.seed(100)
runTime <- proc.time()
for(i in 1:100) r2 <- mvtnorm::rmvnorm(1,mu,sigma,method="chol")
print(proc.time() - runTime)

print(max(abs(r1-r2)))


