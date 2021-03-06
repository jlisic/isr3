

sweepTree <- function( V, M,df=NCOL(V) ) {

  # check V properties
  p <- dim(V)
  if( is.null(p[1]) )    stop(sprintf('V is not a matrix'))
  if( p[1] != p[2] )     stop(sprintf('V is not a square matrix'))

  # if missing create B
  if(missing(M)) {M <- as.logical(matrix(1,p,p)); print('missing M') }
  
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
  
  
#  print("regIndex") 
#  print(regIndex-1) 
#  
#  print("mapIndex") 
#  print(mapIndex-1) 
#
#  print(M)
  
  r.result <- .C("RSweepTree",
    as.double(V[lower.tri(V,TRUE)]), 
    as.integer(c(t(M))),
    as.integer(regIndex -1),
    as.double(rep(0,m[1]*(p[1]+1))), # estimate
    as.integer(mapIndex-1), 
    as.integer(p[1]),
    as.integer(m[1]),
    as.integer(df)
  )

  E <- matrix( r.result[[4]] ,ncol=p[1]+1,byrow=TRUE) 
  E <- E * cbind(M,1)

  return(E) 
}







