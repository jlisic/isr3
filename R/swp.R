if( !('isr3.check' %in% ls()) ) { 
  print("loading isr3.so")
  dyn.load('~/src/isr3/src/covarTree.so')
  isr3.check <- T
}


#SWP
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
    as.integer(p)
  )
      
  V[lower.tri(V,T)] <- r.result[[1]] 
  V <- t(V)
  V[lower.tri(V,T)] <- r.result[[1]] 

  return(V) 
}


#RSWP
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
    as.integer(p)
  )
           
  V[lower.tri(V,T)] <- r.result[[1]] 
  V <- t(V)
  V[lower.tri(V,T)] <- r.result[[1]] 

  return(V) 
}






