# this is a test program 
library(ISR3)
library(MNP)
library(mvtnorm)

set.seed(100)
mu <- rnorm(100)
sigma <- rWishart(1,104,diag(100))[,,1]
iter <- 1

set.seed(100)
runTime <- proc.time()
for(i in 1:iter) r1 <- RMVN(mu,sigma)
print(proc.time() - runTime)

set.seed(100)
runTime <- proc.time()
for(i in 1:iter) r2 <- RMVN2(mu,sigma)
print(proc.time() - runTime)

library(mvtnorm)
set.seed(100)
runTime <- proc.time()
for(i in 1:iter) r3 <- mvtnorm::rmvnorm(1,mu,sigma,method="chol")
print(proc.time() - runTime)

print(max(abs(r1-r2)))


