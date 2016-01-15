source('sweepTree.R')

######################################################################
# The basic idea here is to take a set of variables and covariates
# Example One
######################################################################

# this is a list of all dependent and independent variables
covarList <- c("CoVar1", "CoVar2", "Var1", "Var2", "Var3")

# these are dependent variables that we will record model parameters for
varList <-  c("Var1", "Var2", "Var3")


# simulation parameters
set.seed(100)
n <- 100000 


# variable relationsips
covarMatrix <- c( 
#          CoVar1  CoVar2  Var1  Var2  Var3
# 1. CoVar1
                2.00,      
# 1. CoVar2
                0.33,      2.00,
# 1. Var1
                0.50,      0.00,    2.00,
# 2. Var2
                0.25,      0.50,    0.75,    2.00,
# 3. Var3
                0.60,      0.10,    0.10,    0.20,    2.00
)

# simulation of variables under the variable relationships

#upper
p <- floor(sqrt(2*length(covarMatrix))) 
S <- matrix(0,p,p)
S[upper.tri(S,T)] <- covarMatrix
S <- t(S)
S[upper.tri(S,T)] <- covarMatrix

# cholesky
#U <- chol(S)

#X <- matrix(rnorm(n*p), nrow=n) %*% U

#XX <- t(X) %*% X
XX <- S

colnames(XX) <- covarList 
rownames(XX) <- covarList 


fitMatrix <- matrix( c( 
#          CoVar1  CoVar2  Var1  Var2  Var3
# 1. Var1
                T,      T,    F,    F,    F,
# 2. Var2
                T,      F,    T,    F,    F,
# 3. Var3
                T,      T,    T,    T,    F 
),nrow=3,byrow=T)

colnames(fitMatrix) <- covarList 
rownames(fitMatrix) <- varList 



# create tree
#myTree <- 
E<-  sweepTree(XX,fitMatrix)



