

######################################################################
# The basic idea here is to take a set of variables and covariates
# Example One
######################################################################

covarList <- c("CoVar1", "CoVar2", "Var1", "Var2", "Var3")
varList <-  c("Var1", "Var2", "Var3")

set.seed(100)

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

#upper
n <- floor(sqrt(2*length(covarMatrix))) 
S <- matrix(0,n,n)
S[upper.tri(S,T)] <- covarMatrix
S <- t(S)
S[upper.tri(S,T)] <- covarMatrix

# cholesky
U <- chol(S)


covarMatrix <- matrix( c( 
#          CoVar1  CoVar2  Var1  Var2  Var3
# 1. Var1
                1,      1,    0,    0,    0,
# 2. Var2
                1,      0,    1,    0,    0,
# 3. Var3
                1,      1,    1,    1,    0
),nrow=3,byrow=T)

colnames(covarMatrix) <- covarList 
rownames(covarMatrix) <- varList 

# create tree
#myTree <- covarTree(covarList, varList)



