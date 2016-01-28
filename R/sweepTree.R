# this is a test program 

if( !('isr3.check' %in% ls()) ) { 
  print("loading isr3.so")
  dyn.load('~/src/isr3/src/covarTree.so')
  isr3.check <- T
}

#SWP
print("This test needs to be rewritten to use X, instead of V")
sweepTree <- function( V, M,n=0 ) {

  # check V properties
  p <- dim(V)
  if( is.null(p[1]) )    stop(sprintf('V is not a matrix'))
  if( p[1] != p[2] )     stop(sprintf('V is not a square matrix'))

  # check that M is logical
  if(!is.logical(M)) stop(sprintf('M is not logical (boolean)'))

  # check M properties
  m <- dim(M)
  if( !is.matrix(M) )    stop(sprintf('M is not a matrix'))

  # get rownames from M 
  regVars <- rownames(M)
  if(is.null(regVars)) stop(sprintf('Row names required for V and M'))
  #cat(sprintf("Performing Regression on:\n%s\n", paste(regVars,collapse=',')))

  # get positions in M of rownames
  regIndex <- as.integer(sapply(regVars, FUN=match, rownames(V) ))
  if(is.na(sum(regIndex))) stop(sprintf('Row names differ between V and M'))
  

  # create a reverse mapping
  mapIndex <- rep(0,p[1])
  mapIndex[regIndex] <- 1:length(regIndex) 
  
  
  print("regIndex") 
  print(regIndex-1) 
  
  print("mapIndex") 
  print(mapIndex-1) 

  r.result <- .C("RSweepTree",
    as.double(V[lower.tri(V,T)]), 
    as.integer(c(t(M))),
    as.integer(regIndex -1),
    as.double(rep(0,m[1]*(p[1]+1))),
    as.integer(mapIndex-1), 
    as.integer(p[1]),
    as.integer(m[1]),
    as.integer(n)
  )

  E <- matrix( r.result[[4]] ,ncol=p[1]+1,byrow=T) 
  #E <- E[regIndex,]
  E <- E * cbind(M,1)

  return(E) 
}







