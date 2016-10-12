#' Iterative Sequential Regression 
#' 
#' \code{isr} performs imputation of missing values based on an optionally
#' specified model.  Missingness is assumed to be missing at random (MAR). 
#'
#' @param X A matrix of points to be imputed or used for covariates by isr.  
#'   \code{NA} values are considered missing.  If column names are used, 
#'   duplicate column names are not allowed. 
#' @param M A boolean valued optional matrix specifying the factorized pdf of the joint multivariate normal distribution of the variables requiring imputation.
#'   A description of the factorized pdf is provided in the details.
#'   The column names of \code{M} must match the column names of \code{X}, and the rows names of \code{M} must be a subset of the column names in \code{X}, in the same order as in \code{X}. 
#'   Variables requiring imputation are each associated with a row in \code{M}; the conditional relationship to variables in \code{X} is indicated by the boolean valued elements of each row vector.
#'   A value of \code{TRUE} indicates conditional dependence, likewise a value of \code{FALSE} indicates conditional independence.  
#'   Because this is a factorized pdf, the variable in the first row of \code{M} cannot specify a conditional dependence with a variable in a later row of \code{M}. 
#'   If \code{M} is missing, dependence is assumed between all variables being imputed.  No missing values are allowed.
#' @param Xinit An optional matrix with the same dimensions of \code{X}, with no missing values.
#'   All values of \code{Xinit} should match those of \code{X}, with the exception of missing
#'   values.  Values of \code{Xinit} that share an index with a missing value in \code{X} are
#'   treated as initial imputations.  If Xinit is not specified, variable means are used as initial imputations.
#' @param mi A scalar indicating the number of imputations to return 
#' @param burnIn A scalar indicating the number of iterations to burn in before
#'   returning imputations.  Note, that burnIn is the total number of iterations, no thinning is performed until multiple imputation generation starts. 
#' @param thinning  A scalar that represents the amount of thinning for the MCMC routine.  A value of one implies no thinning.
#' @param intercept A logical value identifying if the imputation model should 
#'   have an intercept.
#' @return This function returns a list with two elements: \code{param} a three dimensional array  
#'   of parameter estimates of the factored pdf.  The last dimension is an index for the multiple imputations.
#'   \code{imputed} a three dimensional array of \code{X} with imputed values, the last dimension is an 
#'   index for the multiple imputations. 
#' @details  The ISR algorithm performs Bayesian multivariate normal imputation.  This imputation follows two steps, an imputation step and a prediction step.
#'   In the imputation step, the missing values are imputed from a Normal-Inverse-Wishart model with non-informative priors.
#'   In the prediction step, the parameters are estimated using both the observed and imputed values.
#'    
#'   Imputation of parameters are done through the conditional factoring of the joint pdf.
#'   A conditional factoring is an expansion of the joint pdf of all
#'   the dependent variables in \code{X}.  e.g. Pr(X|Z) = Pr(X1,X2,X3|Z) = Pr(X1,Z) Pr(X2|X1,Z) Pr(X3|X1,X2,Z),
#'   where the right hand side is the fully conditional specification for the dependent variables X1-X3 and independent variable Z.
#'
#' @examples 

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
#' 
#' # specify relationships
#' fitMatrix <- matrix( c( 
#'   #  Covar2    CoVar1   Var1     Var2     Var3
#'      # 1. Var1
#'        TRUE,    TRUE,   FALSE,   FALSE,   FALSE,
#'      # 2. Var2
#'        TRUE,    TRUE,   FALSE,    FALSE,   FALSE,
#'      # 3. Var3
#'        TRUE,    TRUE,   TRUE,    TRUE,    FALSE 
#'  ),nrow=3,byrow=TRUE)
#' 
#' covarList <- c('Covar2', 'CoVar1', 'Var1', 'Var2','Var3')
#' 
#' # setup names
#' colnames(fitMatrix) <- covarList 
#' rownames(fitMatrix) <- covarList[-1:-2] 
#' colnames(X) <- covarList
#' 
#' XImputed <- isr(X,fitMatrix)
#'
#' @useDynLib ISR3
#' @export
#'
#' @references
#' Robbins, M. W., & White, T. K. (2011). Farm commodity payments and imputation in the Agricultural Resource Management Survey. American journal of agricultural economics, DOI: 10.1093/ajae/aaq166.
#'

isr <- function(X, M, Xinit, mi=1, burnIn=100, thinning=20, intercept=T) {
  
  # get missing values in X observed
  Xobserved <- !is.na(X)

  # ensure mi is a positive integer
  mi <- as.integer(mi)
  if( mi < 1 ) stop("mi must be larger than 0.")
  
  # ensure thinning is a positive integer
  thinning <- as.integer(thinning)
  if( thinning < 1 ) stop("thinning must be larger than 0.")
  
  # ensure burnIn is an integer
  burnIn <- as.integer(burnIn)
  if( burnIn < 0 ) stop("burnIn must be larger than or equal to 0.")


  # initial impute via column means
  if( missing(Xinit) ) {
    Xbar <- colMeans(X,na.rm=T) 
    X[!Xobserved] <- matrix( rep(Xbar,nrow(X)), byrow=T,ncol=length(Xbar))[!Xobserved]  
  } else {
    X <- Xinit
  }

  # check if X is a matrix
  if( !is.matrix(X) ) stop("X is not a matrix.")
  
  # create colnames if they don't exist
  if( is.null(colnames(X)) )  colnames(X) <- sprintf("X%d",1:ncol(X))

  # handle missing M
  if(missing(M)) {
    M <- matrix(TRUE,nrow=ncol(X),ncol=ncol(X)) 
    M[upper.tri(M,TRUE)] <- 0
    colnames(M) <- colnames(X) 
    rownames(M) <- colnames(X) 
  } 
  
  # check if M is a matrix
  if( !is.matrix(M) ) stop("M is not a matrix.")

  # check if rownames for M are in X
  if( sum(rownames(M) %in% colnames(X)) != nrow(M) ) stop("Some row names in M are not column names in X.")

  # check if colnames in X and M match 
  if( sum(colnames(X) == colnames(M)) != ncol(X) ) stop(sprintf("Column names between X and M do not match, %d != %d.", sum(colnames(X) %in% colnames(M), ncol(X) ))) 

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
  est <- rep(0,b*(p+1)*mi)
  S <- rep( 0, n*b*mi ) 

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
    as.double(S),          #  9
    as.double(est),        # 10
    as.integer(burnIn),       # 11
    as.integer(thinning),     # 12
    as.integer(mi)            # 13
  )


  E <- r.result[[10]]
  
  S <- matrix( r.result[[9]], nrow=n*mi )  

  E <- array(E,dim=c(ncol(M)+1,nrow(M),mi),
            list(c(colnames(M),'S2'),rownames(M),c()) ) 

  # need to do something to fix this for with and without intercept
  if( intercept ) {
    print(c(n,NCOL(S),mi))
    print(colnames(M)[-1])
    S <- array( S,dim=c(n,NCOL(S),mi), list(c(), rownames(M),c() ) ) 
  } else {
    S <- array( S,dim=c(n,NCOL(S),mi), list(c(), rownames(M),c() ) ) 
  }

  return( list( param=E, imputed=S) )
}






