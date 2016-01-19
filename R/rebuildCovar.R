# this is a test program 

if( !('isr3.check' %in% ls()) ) { 
  print("loading isr3.so")
  dyn.load('~/src/isr3/src/covarTree.so')
  isr3.check <- T
}

#SWP
rebuildCovar <- function( E ) {


  # check M properties
  if( !is.matrix(E) )    stop(sprintf('E is not a matrix'))
  m <- nrow(E)
  n <- ncol(E) - 1

  r.result <- .C("RRebuildCovar",
    as.double(rep(0,m*(m+1)/2)), 
    as.double(c(t(E))),
    as.integer(0:(m-1)),
    as.integer(n),
    as.integer(m)
  )

  E <- r.result[[1]]  
  #E <- matrix( r.result[[1]] ,ncol=p[1]+1,byrow=T) 
  #E <- E[regIndex,]
  #E <- E * cbind(M,1)

  return(E) 
}







