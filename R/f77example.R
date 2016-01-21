# this is a test program 


if( !('isr3.check' %in% ls()) ) { 
  print("loading isr3.so")
  dyn.load('~/src/isr3/src/covarTree.so')
  isr3.check <- T
}



# fortran through C example
RcholInv <- function( A, B ) {

  if( !is.matrix(A) )    stop(sprintf('A is not a matrix'))
  if( !is.matrix(B) )    stop(sprintf('B is not a matrix'))

  if( nrow(A) != ncol(A) ) stop( sprintf('X is not square'))

  r.result <- .C("RcholInv",
    as.double(A),
    as.double(B),
    as.integer(nrow(B)),
    as.integer(ncol(B))
  )

  X <- matrix(r.result[[2]],ncol=ncol(A))

  return(X) 
}







set.seed(100)

# generate some data N_5(0,v)
v <- rWishart(1,10,diag(5))[,,1]
x <- matrix(rnorm(50),ncol=5)

# y:= x * v **T
v.L <- chol2inv(chol(v)) * lower.tri(v.L,T)
y <- x %*% v.L 


out <- RcholInv( t(chol(v)), x ) 

#print(out)

#print( ex.inv - chol2inv( chol( ex) ) )







