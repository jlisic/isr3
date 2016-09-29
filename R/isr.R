#' Iterative Sequential Regression 
#' 
#' \code{isr} performs imputation of missing values based on an optionally
#' specified model.  Missingness is assumed to be missing at random (MAR). 
#'
#' @param X A matrix of points to be imputed or used for covariates for the
#'   by isr.  \code{NA} values are considered missing.  Distinct column names 
#'   are required for each variable. 
#' @param M A matrix specifying the relationships between each of the variables.
#'   Each column name of M must match to the column names of X.  Each row of 
#'   M identifies a variable that will be imputed for, all other variables are
#'   treated as covariates for all variables.  Only the lower diagonal of this
#'   matrix is used to identify relationships, relationships between the 
#'   variable and itself are also ignored.  If missing, dependence is assumed
#'   between all variables.
#' @param mcmcIter A scalar indicating the number of iterations to run.
#' @param sampleRate  A scalar that provides thinning for the MCMC routine.
#' @param intercept A logical value identifying if the imputation model should 
#'   have an intercept.
#'
#' @examples 
#' set.seed(100)
#' mcmciter <- 100 
#' 
#' # these are dependent variables that we will record model parameters for
#' varList <-  c("Var1", "Var2", "Var3")
#' 
#' # simulation parameters
#' set.seed(100)
#' n <- 30
#' p <- 5
#' 
#' # generate a covar matrix
#' covarMatrix <- rWishart(1,p+1,diag(p))[,,1]
#' 
#' # simulation of variables under the variable relationships
#' U <- chol(covarMatrix)
#' 
#' X <- matrix(rnorm(n*p), nrow=n) %*% U
#' colnames(X) <- covarList
#' 
#' # specify relationships
#' fitMatrix <- matrix( c( 
#' #             Var1  Var2  Var3
#' # 1. Var1
#'                 F,    F,    F,
#' # 2. Var2
#'                 T,    F,    F,
#' # 3. Var3
#'                 T,    T,    F 
#' ),nrow=3,byrow=T)
#' 
#' colnames(fitMatrix) <- covarList 
#' rownames(fitMatrix) <- varList 
#' 
#' E<-  isr(X,fitMatrix,mcmcIter=mcmciter)
#'
#' @useDynLib ISR3
#' @export
#'
#' @references
#' Robbins, M. W., & White, T. K. (2011). Farm commodity payments and imputation in the Agricultural Resource Management Survey. American journal of agricultural economics, DOI: 10.1093/ajae/aaq166.
#'

isr <- function(X, M, mcmcIter=100, sampleRate=20, intercept=T) {
  
  # get missing values in X observed
  Xobserved <- !is.na(X)

  # initial impute via column means
  Xbar <- colMeans(X,na.rm=T) 
  X[!Xobserved] <- matrix( rep(Xbar,nrow(X)), byrow=T,ncol=length(Xbar))[!Xobserved]  

  # check if X is a matrix
  if( !is.matrix(X) ) stop("X is not a matrix.")
  
  # create colnames if they don't exist
  if( is.null(colnames(X)) )  colnames(X) <- sprintf("X%d",1:ncol(X))

  # handle missing M
  if(missing(M)) {
    M <- matrix(T,nrow=ncol(X),ncol=ncol(X)) 
    M[upper.tri(M,T)] <- 0
    colnames(M) <- colnames(X) 
    rownames(M) <- colnames(X) 
  } 
  
  # check if M is a matrix
  if( !is.matrix(M) ) stop("M is not a matrix.")

  # check if rownames for M are in X
  if( sum(rownames(M) %in% colnames(X)) != nrow(M) ) stop("Some row names in M are not column names in X.")

  # check if colnames in X and M match 
  if( sum(colnames(X) %in% colnames(M)) != ncol(X) ) stop(sprintf("Column names between X and M do not match, %d != %d.", sum(colnames(X) %in% colnames(M), ncol(X) ))) 

  # ensure the diagonal is false
  diag(M) <- F


  # handle intercept
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

  # saved results
  est <- matrix(0,b,p+1) 
  S <- rep( 0, n*b* ceiling(mcmcIter/sampleRate) )  

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
    as.double( c(X) ),        #  1
    as.integer(c(Xobserved)), #  2
    as.integer(c(t(M))),      #  3
    as.integer(regIndex -1),  #  4
    as.integer(mapIndex -1),  #  5
    as.integer(n),            #  6
    as.integer(p),            #  7
    as.integer(b),            #  8
    as.double(c(S)),          #  9
    as.double(c(est)),        # 10
    as.integer(mcmcIter),     # 11
    as.integer(sampleRate)    # 12
  )

  E <- matrix( r.result[[1]], nrow=n)
  colnames(E) <- colnames(M) 
  S <- matrix( r.result[[9]], nrow=n*ceiling(mcmcIter/sampleRate) )  

  return( list( X=E, imputed=S) )
}






