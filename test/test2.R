library(ISR3)

######################################################################
# Test Two : rebuild speed test 
######################################################################






# simulation parameters
set.seed(100)
n <- 80 
p <- 5 
missing <- 20

# these are dependent variables that we will record model parameters for
varList <-  sprintf("Var%d",1:p)

# variable relationsips

S <- rWishart(1,p+ceiling(sqrt(p)),diag(p))[,,1]

# cholesky
U <- chol(S)

X <- matrix(rnorm(n*p), nrow=n) %*% U
colnames(X) <- varList

###############################################
################ start here ###################
###############################################

fitMatrix <- matrix(1,p,p)
fitMatrix[upper.tri(fitMatrix,T)] <- 0
fitMatrix <- fitMatrix==1
colnames(fitMatrix) <- varList 
rownames(fitMatrix) <- varList 

X[sample(1:(n*p),size=missing)] <- NA


isr.out <- isr(X, fitMatrix, mcmcIter=1,sampleRate=1,intercept=T)
















