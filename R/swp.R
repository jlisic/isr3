#' Sweep Function
#' 
#' \code{SWP} performs the sweep operator.
#'
#' @param V A symmetric matrix to be swept; this matrix cannot contain missing
#'   data or infinite values.
#' @param b An array of integers or column names to sweep.
#' @return The swept matrix \code{V}.  Sweeping will not occur if the column
#'   being swept has a zero-valued diagonal element.
#' @details This program applies the sweep operator as defined in (Dempster 1969).  
#' @examples 
#' set.seed(100)
#' # generate a symmetric positive definite matrix
#' Sigma <- rWishart(1,4,diag(3))[,,1]
#' # sweep all the columns to produce the inverse
#' Sigma.inv <- SWP(Sigma,1:3)
#' @useDynLib ISR3
#' @export
#' @references
#' Dempster, A.P. (1969). \emph{Elements of continuous multivariate analysis}. Reading, MA: Addison-Wesley.
SWP <- function( V, b ) {

  if( is.character(b)[1] ) b <- sapply(b, match,colnames(V) )
  if(length(b) < 1) stop(sprintf('no column selected')) 

  p <- dim(V)
  if( is.null(p[1]) )    stop(sprintf('not a matrix'))
  if( p[1] != p[2] )     stop(sprintf('not a square matrix'))
  if( sum(b > p[1]) > 0) stop(sprintf('number of columns is less than the col/row selected'))

  r.result <- .C("RVSWP",
    as.double(V[lower.tri(V,T)]), 
    as.integer(b-1), 
    as.integer(p),
    as.integer(length(b))
  )
      
  V[lower.tri(V,T)] <- r.result[[1]] 
  V <- t(V)
  V[lower.tri(V,T)] <- r.result[[1]] 

  return(V) 
}


#' Reverse Sweep Function
#' 
#' \code{RSWP} performs the reverse sweep operator.
#'
#' @param V A symmetric matrix to be reverse swept; this matrix cannot contain missing 
#'   data or infinite values.
#' @param b An array of integers or column names to reverse sweep.
#' @return The reverse swept matrix \code{V}.  Reverse sweeping will not occur if the 
#'   column being swept has a zero-valued diagonal element.
#' @details This program applies the reverse sweep operator as defined in 
#'   (Dempster 1969).  
#' @examples 
#' set.seed(100)
#' # generate symmetric positive definite matrix
#' Sigma <- rWishart(1,4,diag(3))[,,1]
#' # sweep all the columns to produce the inverse
#' # and then reverse sweep them all back to Sigma
#' Sigma2 <- RSWP(SWP(Sigma,1:3),1:3)
#' @useDynLib ISR3
#' @export
#' @references
#' Dempster, A.P. (1969). \emph{Elements of continuous multivariate analysis}. Reading, MA: Addison-Wesley.
RSWP <- function( V, b ) {

  if( is.character(b)[1] ) b <- sapply(b, match,colnames(V) )
  if(length(b) < 1) stop(sprintf('no column selected')) 

  p <- dim(V)
  if( is.null(p[1]) )    stop(sprintf('not a matrix'))
  if( p[1] != p[2] )     stop(sprintf('not a square matrix'))
  if( sum(b > p[1]) > 0) stop(sprintf('number of columns is less than the col/row selected'))

  r.result <- .C("RVRevSWP",
    as.double(V[lower.tri(V,T)]), 
    as.integer(b -1), 
    as.integer(p),
    as.integer(length(b))
  )
           
  V[lower.tri(V,T)] <- r.result[[1]] 
  V <- t(V)
  V[lower.tri(V,T)] <- r.result[[1]] 

  return(V) 
}






