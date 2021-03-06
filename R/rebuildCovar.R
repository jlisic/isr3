# test function


rebuildCovar <- function( E ) {


  # check M properties
  if( !is.matrix(E) )    stop(sprintf('E is not a matrix'))
  E <- cbind(E[,rownames(E)],E[,ncol(E)])
  n <- nrow(E) 


  r.result <- .C("RRebuildCovar",
    as.double(rep(0,n*(n+1)/2)), 
    as.double(c(t(E))),
    as.integer(n)
  )

  E <- E[,-1]
  E[lower.tri(E,T)] <- r.result[[1]]  
  E <- t(E)
  E[lower.tri(E,T)] <- r.result[[1]]  

  return(E) 
}







