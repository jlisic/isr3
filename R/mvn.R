
RMVN2 <- function( mu, sigma) { 

  n <- NROW(sigma)
  x <- rep(0,n) 
  if( length(mu) != n) stop("mu does not have correct dimensions")
  
  result <- .C("RMVN2",  
               as.double(x),
               as.double(mu), 
               as.double(sigma[lower.tri(sigma,T)]), 
               as.integer(n)
               )
  return(result[[1]])
}




#RWISH <- function( S, v) { 
#
#  n <- NROW(S)
#  x <- rep(0,n^2) 
#  if( v <= n) stop("invalid value of df")
#  
#  result <- .C("RWISH",  
#               as.double(x),
#               as.double(S), 
#               as.integer(v),
#               as.integer(n)
#               )
#
#  return(matrix(result[[1]],nrow=n))
#}
#

