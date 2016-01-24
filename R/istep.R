# this is a test program 
source('~/Documents/ISR/isrnass/R/myrswp.R')

if( !('isr3.check' %in% ls()) ) { 
  print("loading isr3.so")
  dyn.load('~/src/isr3/src/covarTree.so')
  isr3.check <- T
}



#  void RIStep(
#     double * X,   // observations matrix of dim n by p
#     double * S,   // variance of dim p by p
#     double * Beta, //Beta matrix
#     int * MIndex,  //matrix of missing value indexes m by p
#     int * nPtr,
#     int * pPtr
#     ) {    

RIStep <- function(X, S, B, M) {

  n <- nrow(X)
  p <- nrow(S)
  m <- nrow(M)
 
  # handle vector valued B 
  if( is.null(nrow(B)) ) {
    b = 1 
  } else {
    b = ncol(B)
  }

  if( p != ncol(X) ) stop(sprintf('incompatable dimensions between X and S'))
  if( p != ncol(X) ) stop(sprintf('incompatable dimensions between M and X'))
  if( n != nrow(X) ) stop(sprintf('incompatable dimensions between M and X'))
  
  r.result <- .C("RIStep",
    as.double(c(X)),
    as.double(c(S)),
    as.double(c(B)),
    as.integer(c(M)),
    as.integer(n),
    as.integer(p),
    as.integer(b)
  )

  return(r.result)

}







set.seed(100)

# generate some data N_5(0,v)
v <- rWishart(1,10,diag(5))[,,1]
x <- matrix(rnorm(50),ncol=5)

# y:= x * v **T
v.L <- chol2inv(chol(v)) * lower.tri(v,T)
y <- x %*% v.L 


out <- RcholInv( t(chol(v)), x ) 

#print(out)

#print( ex.inv - chol2inv( chol( ex) ) )



IStep.out <- RIStep( x, v, v , x < Inf) 



y <- (x - x %*% v)
S.inv <- chol2inv(chol(v))

z <- x 

i <- 3
  RS.inv <- myrswp(S.inv,i,method='chol')
  a <- y[,-i] %*% diag(RS.inv[-i,i])


for( i in 1:ncol(v) ) {
  RS.inv <- myrswp(S.inv,i,method='chol')
  z[,i] <- y[,-i] %*% RS.inv[-i,i,drop=F]
}

#v <- v[1:2,1:2]
#
#
#
#i <- 2 
#vii <- v[i,i] - v[i,-i] %*% a[-i,-i] %*% v[i,-i]
#
#a.swp <- myrswp(-a,i,method="chol")
#
#vii2 <- a.swp[i,i]
#
#v[2,2] - v[1,2]^2/v[1,1] 
#
#
#
#print(vii)
#print(vii2)
#


