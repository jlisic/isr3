source('sweepTree.R')
source('rebuildCovar.R')
source('getCovs.R')

######################################################################
# Test Two : rebuild speed test 
######################################################################






# simulation parameters
set.seed(100)
n <- 100000 
p <- 200

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

XX <- t(X) %*% X

fitMatrix <- matrix(1,p,p)
fitMatrix[upper.tri(fitMatrix,T)] <- 0
fitMatrix <- fitMatrix==1
colnames(fitMatrix) <- varList 
rownames(fitMatrix) <- varList 



# create tree
sweep.time <- proc.time()
E<-  sweepTree(XX,fitMatrix)
print(proc.time() - sweep.time)



###############################################
################ start here ###################
###############################################
## prep
X.df <- as.data.frame(X)


E2 <- E

fit.time <- proc.time()
for( i in 1:nrow(fitMatrix) ) {
  
  a.string <- paste( colnames(fitMatrix)[fitMatrix[i,]], collapse=" + ") 

  if( a.string != "" ) {
    a <- sprintf("%s ~ -1  + %s",
    rownames(fitMatrix)[i],
    paste( colnames(fitMatrix)[fitMatrix[i,]], collapse=" + ") 
    )
  } else {
    a <- sprintf("%s ~ -1 ", rownames(fitMatrix)[i])
  }

  fit <- lm( as.formula(a) , data=X.df)

  E2[ rownames(fitMatrix)[i], names(fit$coefficients) ] <- fit$coefficients
  E2[ rownames(fitMatrix)[i], ncol(E) ] <- sum(fit$residuals^2) 

}
print(proc.time() - fit.time)


print(max(abs(E - E2)))


###############################################
################  ###################
###############################################

betas <- E[,rownames(E)]
sigmas <- E[,ncol(E)]

cat("*********************REBUILD COVAR********************\n")
fit.time <- proc.time()
CISR3 <- rebuildCovar( E )
print(proc.time() - fit.time)

fit.time <- proc.time()
CISR <- getCovs(betas,sigmas) 
print(proc.time() - fit.time)



print(max(abs(CISR3 - CISR)))





