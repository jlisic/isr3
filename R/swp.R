# this is a test program 

if( !('sweep.check' %in% ls()) ) { 
  print("loading swp.so")
  dyn.load('~/src/rSWP/swp.so')
  sweep.check <- T
}

#SWP
SWP <- function( V, b ) {

  if( is.character(b)[1] ) b <- sapply(b, match,colnames(V) )
  if(length(b) < 1) stop(sprintf('no column selected')) 

  p <- dim(V)
  if( is.null(p[1]) )    stop(sprintf('not a matrix'))
  if( p[1] != p[2] )     stop(sprintf('not a square matrix'))
  if( sum(b > p[1]) > 0) stop(sprintf('number of columns is less than the col/row selected'))

  r.result <- .C("RSWP",
    as.double(c(V)), 
    as.integer(b-1), 
    as.integer(length(b)), 
    as.integer(p)
  )
         
  r.result <- matrix(r.result[[1]],nrow=p,byrow=T) 
  colnames(r.result) <- colnames(V)
  rownames(r.result) <- rownames(V)

  return(r.result) 
}



x <- matrix( c(
10.31288,  2.448485, -2.412443, -1.393328, -6.486046,
2.448485,  7.749054, -2.442433, 0.5318202,  1.395485,
-2.412443, -2.442433,  11.15893,  1.965747,  6.863963,
-1.393328, 0.5318202,  1.965747,  10.67714,  1.181446,
-6.486046,  1.395485,  6.863963,  1.181446,  12.50144
), nrow=5)






