# this is a test program 
library(isrnass)
#source('~/src/isr3/R/mnp.R')

if( !('isr3.check' %in% ls()) ) { 
  print("loading isr3.so")
  dyn.load('~/src/isr3/src/covarTree.so')
  isr3.check <- T
}



#rmvnorm <- function(count, mu, sigma) { 
#  return(RMVN2(mu,sigma))
#}




# input:
# X, matrix of missing with NA's
# M, matrix of correlation relationships
#    If a variable has a missing variable, an entry must exist in M.

isr <- function(X, M, mcmcIter=10, Xobserved, intercept=T) {
  
  # get missing values in X observed
  if( missing(Xobserved) ) Xobserved <- !is.na(X)  

  if(intercept) {
    X.names <- colnames(X)
    X <- cbind(1,X)
    colnames(X) <- c("B_0",X.names)
    
    M.names <- colnames(M)
    M <- cbind(1,M)
    colnames(M) <- c("B_0",M.names)

    # add col to Xobserved
    Xobserved <- cbind(1,Xobserved)
  }

  # get dimensions
  n <- nrow(X)
  p <- ncol(X)
  b <- nrow(M)

  Sigma <- matrix(0,p,p) 
  est <- matrix(0,b,p+1) 

  if( p != ncol(M) ) stop(sprintf('incompatable dimensions between X and M'))
  
  # get rownames from M 
  regVars <- rownames(M)
  if(is.null(regVars)) stop(sprintf('Row names required for M'))

  # get positions in M of rownames
  regIndex <- as.integer(sapply(regVars, FUN=match, colnames(M) ))
  
  # create a reverse mapping
  mapIndex <- rep(0,p)
  mapIndex[regIndex] <- 1:length(regIndex) 
  
  # Useful example 
  #
  # M:
  #         B_0 CoVar1  CoVar2  Var1  Var2  Var3
  # 0. Var1 T,  T,      T,      F,    F,    F
  # 1. Var2 T,  T,      F,      T,    F,    F
  # 2. Var3 T,  T,      T,      T,    T,    F 
  #
  # regIndex <-  3 4 5 
  # mapIndex <- -1, -1, -1, 0, 1 2
  
  r.result <- .C("Risr",
    as.double(c(X)),
    as.integer(c(Xobserved)),
    as.integer(c(t(M))),
    as.integer(regIndex -1),
    as.integer(mapIndex -1),
    as.integer(n),
    as.integer(p),
    as.integer(b),
    as.double(c(Sigma)),
    as.double(c(est)),
    as.integer(mcmcIter)
  )

  E <- matrix( r.result[[1]], nrow=n)
  colnames(E) <- colnames(M) 

  return(E)

}






set.seed(100)
mcmciter <- 1000 


# this is a list of all dependent and independent variables
covarList <- c("CoVar1", "CoVar2", "Var1", "Var2", "Var3")

# these are dependent variables that we will record model parameters for
varList <-  c("Var1", "Var2", "Var3")


# simulation parameters
set.seed(100)
n <- 50000
p <- 5

# generate some data N_5(0,v)
covarMatrix <- rWishart(1,p+ceiling(sqrt(p)),diag(p))[,,1]
x <- matrix(rnorm(p*n),ncol=p)


# simulation of variables under the variable relationships
U <- chol(covarMatrix)

X <- matrix(rnorm(n*p), nrow=n) %*% U
colnames(X) <- covarList

###############################################
################ ISR 3      ###################
###############################################

XX <- t(X) %*% X

fitMatrix <- matrix( c( 
#          CoVar1  CoVar2  Var1  Var2  Var3
# 1. Var1
                T,      T,    F,    F,    F,
# 2. Var2
                T,      T,    T,    F,    F,
# 3. Var3
                T,      T,    T,    T,    F 
),nrow=3,byrow=T)

colnames(fitMatrix) <- covarList 
rownames(fitMatrix) <- varList 

dat.tmp <- X
obs.tmp <- (X < Inf) * 1
obs.tmp[7:10,varList] <- 0 


# perform ISR 
set.seed(100)
sweep.time <- proc.time()
E<-  isr(dat.tmp,fitMatrix,mcmcIter=mcmciter,obs.tmp)
print(proc.time() - sweep.time)




###############################################
################ ISR 2      ###################
###############################################

# perform pstep from old isr

#    column  missingness  variable  tranformation  category  Group
#    number  indicator    type      type                     Indicator
index.tmp <- matrix(
  c( 1,       0,           3,        1,             0,        0,
     2,       0,           3,        1,             0,        0,
     3,       2,           3,        1,            -1,        0, 
     4,       2,           3,        1,            -1,        0, 
     5,       2,           3,        1,            -1,        0  ), nrow=5, byrow=TRUE)
 
rownames(index.tmp) <- c('CoVar1','CoVar2','Var1','Var2','Var3')
colnames(index.tmp) <- c( 'Column Number', 'Missingness Indicator', 'Variable Type', 'Transformation Type', 'Category', 'Group Indicator')


source('~/Documents/ISR/isrnass/R/PStep.R')
source('~/Documents/ISR/isrnass/R/fixVariableNames.R')
source('~/Documents/ISR/isrnass/R/myswp.R')
source('~/Documents/ISR/isrnass/R/myrswp.R')
source('~/Documents/ISR/isrnass/R/inOrder.R')
library(mvtnorm)
ord <-  c(3,4,5)
badcovar <- c()
intersects <- c()
grps <- c()


set.seed(100)
sweep.time <- proc.time()
initParams=PStep2(
  dat=dat.tmp,
  index=index.tmp,
  ord=ord,
  badcovar=badcovar,
  intersects=intersects,
  grps=grps,
  method="chol"
)

### ISR.ARMS ###
dat.new=ISR.ARMS(
  data=dat.tmp,
  obs=obs.tmp,
  initParams=initParams,
  index=index.tmp,
  mis=ord,               
  badcovar=badcovar,
  intersects=intersects,
  grps=grps,
  mcmciter=mcmciter,
  imp=1,
  where=c(1)
)
print(proc.time() - sweep.time)


print("max diff")
print( max(abs(cbind(1,dat.new) - E)))


