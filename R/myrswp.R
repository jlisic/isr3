myrswp <-
function (V, b, method="") {

    p <- ncol(V)
    u <- is.na(match(1:p, b))
    a <- (1:p)[u]
    out <- 0 * V

    dimnames(out) <- dimnames(V)

    if (length(a) == 0) 
        return(-solve(V))
    else if (length(a) == p) 
        return(V)
    else {
#  .   a b 
#  a Saa Sab
#  b     Sbb 
        Saa <- V[a, a, drop = FALSE]
        Sab <- V[a, b, drop = FALSE]
        Sbb <- V[b, b, drop = FALSE]
    
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

        
        B <- Sab %*% (Sbb)
        out[a, a] <- Saa - B %*% t(Sab)
        out[a, b] <- -B
        out[b, a] <- -t(B)
        out[b, b] <- -Sbb
        return(out)
    }
}
