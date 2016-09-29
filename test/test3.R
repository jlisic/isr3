# test for swp and rswp
source('swp.R')
#source('myswp.R')
#source('myrswp.R')


maxERR <- 2e-8

set.seed(100)

S <- rWishart(1,12,diag(3))[,,1]

print(solve(S))


#o1 <- myswp(S,1,method='chol')
#o2 <- myswp(o1,2,method='chol')
#o3 <- myswp(o2,3,method='chol')
#r2 <- myrswp(o3,3,method='chol')
#r1 <- myrswp(o2,2,method='chol')
#r0 <- myrswp(o1,1,method='chol')


a1 <- SWP(S,1)
a2 <- SWP(a1,2)
a3 <- SWP(a2,3)
b2 <- RSWP(o3,3)
b1 <- RSWP(o2,2)
b0 <- RSWP(o1,1)


max.diff <- max( c(
#max(abs(a1 - o1)),
#max(abs(a2 - o2)),
#max(abs(a3 - o3)),
#
#max(abs(r2 - b2)),
#max(abs(r1 - b1)),
#max(abs(r0 - b0)),
#
#max(abs(r0 - S)),

                   max(abs(a3 + chol2inv(chol(S)))),
                   max(abs(b0 - S))
)
)

if( max.diff > maxERR) stop(sprintf("Sweep has unsuitable error levels %f",max.diff))


