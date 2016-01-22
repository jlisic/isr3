# test for swp and rswp
source('swp.R')
source('myswp.R')
source('myrswp.R')

set.seed(100)

S <- rWishart(1,12,diag(3))[,,1]

print(solve(S))

a1 <- SWP(S,1)
a2 <- SWP(a1,2)
a3 <- SWP(a2,3)

o1 <- myswp(S,1,method='chol')
o2 <- myswp(o1,2,method='chol')
o3 <- myswp(o2,3,method='chol')

print(a1 - o1)
print(a2 - o2)
print(a3 - o3)

#r2 <- RSWP(-a3,3)
#r1 <- RSWP(r2,2)
#r0 <- RSWP(r1,1)

