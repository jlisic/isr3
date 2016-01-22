myswp <-
function (V, b, method="") {
    p <- ncol(V)
 
    #a is the set of variables not in b
    u <- is.na(match(1:p, b))
    a <- (1:p)[u]

    # create output data set
    out <- 0 * V
    dimnames(out) <- dimnames(V)

    # if all vars are in p then just solve
    if (length(a) == 0) { 
      return(-solve(V))
    # if b is of length(0) then return the input matrix
    } else if (length(a) == p) { 
      return(V)
    # finally if we have work to do, we start sweeping!
    } else {

#  .   a b 
#  a Saa Sab
#  b     Sbb 
        Saa <- V[a, a, drop = FALSE]
        Sab <- V[a, b, drop = FALSE]
        Sbb <- V[b, b, drop = FALSE]
#  .   a b 
#  a Saa Sab
#  b     Sbb^-1 
       if( method == "chol" ) { 
         # this is to allow the cholesky decomp to run if Sbb is Neg. Def. 
         if( Sbb[1] < 0) {
           Sbb=chol2inv(chol(-1 * Sbb))
           Sbb= -1 * Sbb
         } else {
           Sbb=chol2inv(chol(Sbb))
         } 
       } else {
         Sbb <- solve(Sbb)
       }

        
#  .   a b 
#  a Saa Sab
#  b     Sbb^-1 
        B <- Sab %*% Sbb
        out[a, a] <- Saa - B %*% t(Sab)
        out[a, b] <- B
        out[b, a] <- t(B)
        out[b, b] <- -1 * Sbb

        return(out)
    }
}
