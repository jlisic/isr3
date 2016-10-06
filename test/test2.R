library(ISR3)
#library(isrnass)
######################################################################
# Test Two : rebuild speed test 
######################################################################



# simulation parameters
set.seed(100)
n <- 20 
p <- 5 
missing <- 20 
iter <- 5 

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

######## BUILD INPUTS FOR ISR3 ################
fitMatrix <- matrix(1,p,p)
fitMatrix[upper.tri(fitMatrix,T)] <- 0
fitMatrix <- fitMatrix==1
colnames(fitMatrix) <- varList 
rownames(fitMatrix) <- varList 

X[sample(1:(n*p),size=missing)] <- NA


######### BUILD INPUTS FOR NASSISR ############
#obs <- !is.na(X)
#dat.tmp <- X
#Xbar <- colMeans(X,na.rm=T) 
#dat.tmp[!obs] <- matrix( rep(Xbar,nrow(X)), byrow=T,ncol=length(Xbar))[!obs]  
#
## get the ordering etc...
#ord <- 1:5 
#badcovar <- c()
#intersects <- c()
#grps <- c()
#
## create an index
#index <- cbind(
#               1:5,          # orig. col. index
#               rep(2,5),    # missingness
#               rep(1,5),     # variable type
#               rep(1,5),          # transformation type
#               rep(-1,5),    # category
#               rep(0,5)      # group
#              )
#
#rownames( index ) <- varList 
#colnames(index) <- c(
#  'Column Number',
#  'Missingness Indicator',
#  'Variable Type',
#  'Transformation Type',
#  'Category',
#  'Group Indicator'
#)



#### ISR ####
set.seed(100)
start.time <- proc.time()
isr.out <- isr(X, fitMatrix,mi=iter, burnIn=10,thinning=10,intercept=T)
print( proc.time() - start.time )


#### PStep2 ###
#set.seed(100)
#start.time <- proc.time()
#for( i in 1:iter) {
#params=PStep2(
#  dat=dat.tmp,
#  index=index,
#  ord=ord,
#  badcovar=badcovar,
#  intersects=intersects,
#  grps=grps,
#  method="chol"
#)
#
#### IStep ###
#imputes=IStep(
#  dat.tmp,
#  obs,
#  params,
#  n=nrow(dat.tmp),
#  p1=ncol(dat.tmp),
#  mis=ord,
#  p=length(ord)
#  )
#
#  dat.tmp <- imputes[[1]]
#}
#print( proc.time() - start.time )





#set.seed(100)
#mice.out <- mice(X,m=5,maxit=50, method='pmm') 
#print(isr.out$imputed[,,iter] - imputes[[1]])


