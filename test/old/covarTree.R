source('~/src/rSWP/swp.R')


treeNode <- list(covarName=c(),yesSub=c(),noSub=c(),varList=c(),split=F)
#treeNode <- data.frame(covarName=list(),yesSub=list(),noSub=list(),varList=list())

covarTree <- function( covarList, varList ) {

  if( length(covarList) > 0 ) {
    newNode <- treeNode # create node
    newNode$covarName <- covarList[1] # label node

    node <- newNode # assign node
             
    # find variables in varList that include covarList[1] 
    varListYes <- intersect( 
            varList ,
            rownames( covarMatrix[ covarMatrix[,covarList[1]]==1,,drop=F] )
          )
    # find the ones that do not
    varListNo <- setdiff( varList, varListYes) 
  
    if(length(varListYes)>0) {
      if(length(covarList) == 1) {
        newNode$varList <- varListYes 
      } else {
        newNode$yesSub <- covarTree(covarList[-1], varListYes)
      }
    }
    if(length(varListNo)>0) {
      if(length(covarList) == 1) {
        newNode$varList <- varListNo 
      } else {
        newNode$noSub <- covarTree(covarList[-1], varListNo) 
      }
    }

    if( !is.null(newNode$yesSub) & !is.null(newNode$noSub) ) newNode$split = T

    return( newNode )
  }

  # return items in list
  return( varList )
}


printTree <- function( x, sweepVar=c() ) {

  if( !is.null(x) ) {
    if( is.null(x$varList) ) {
#      if( x$split ) {
#        matrixNameCache[length(matrixCache) +1] <- list(paste(sweepVar,collapse=' '))
#        matrixNameCache[length(matrixCache) +1] <- list(currentMatrix)
#        assign('matrixCache', matrixCache, envir=.GlobalEnv)
#        assign('matrixNameCache', matrixNameCache, envir=.GlobalEnv)
#      }
      printTree(x$yesSub, c(sweepVar, x$covarName) )
      printTree(x$noSub, sweepVar)
    } else {
      # write out error var and beta
      for(y in x$varList){
      c(sweepVar,y, )
      }

    }
  }
}



matrixCache <- list()
betaMatrix  <- covarMatrix
sigmaMatrix <- covarMatrix[,varList]
printTree(myTree)

######################################################################
# The basic idea here is to take a set of variables and covariates
# Example One
######################################################################

covarList <- c("CoVar1", "CoVar2", "Var1", "Var2", "Var3")
varList <-  c("Var1", "Var2", "Var3")

Sigma <- rWishart(1, df=10, diag(5) )[,,1]
print(eigen(Sigma)$values)
rownames(Sigma) <- covarList
colnames(Sigma) <- covarList


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
myTree <- covarTree(covarList, varList)

print("My Tree")
printTree(myTree)




######################################################################
# The basic idea here is to take a set of variables and covariates
# Example Two 
######################################################################


if( F ) { 

# number of covar
p <- 100 
m <- 200


covarMatrix <- matrix( runif(m*(m+p)) ,nrow=m,,byrow=T)
covarMatrix <- (covarMatrix > .5) * 1

varList <-  sprintf("Var%d",1:m)
covarList <- c(sprintf("Covar%d",1:p),varList)


colnames(covarMatrix) <- covarList 
rownames(covarMatrix) <- varList 

# create tree
myTree <- covarTree(covarList, varList)


print("My Tree")
runTime <- proc.time() 
printTree(myTree)
print( proc.time() - runTime)

}






