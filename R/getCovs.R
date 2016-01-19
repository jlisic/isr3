getCovs <-
function (betas,sigmas,p=length(sigmas)) {
  ### betas is a p x p matrix with 0's in the diagonal and the upper triangle.
  #   beta example
  #   0   0   0 
  #   b21 0   0 
  #   b31 b32 0
 
   # covars example 
   #  0 0 0
   #  0 0 0
   #  0 0 0
  covars=matrix(0,p,p)

  for(j in 1:p) {
    if(j==1){
      #  sigma[1]  0 0
      #  0         0 0
      #  0         0 0
      covars[1,1]=sigmas[1]

    # we handle the case of j = 2 separately to avoid extra matrix casts
    } else if(j == 2) {
      covars[2,1] = betas[2,1] * covars[1,1] 
      covars[1,2] = covars[2,1]
      covars[2,2] = sigmas[2] + betas[2,1] *  sigmas[1]  * betas[2,1] 
      
      print(covars)
    
    # This is the general case for j > 2 
    } else if(j > 2) {

    # create range up to j - 1
      j.range <- 1:(j-1)

    # example:
    #   coefs = [ b31 b32 ] 
      coefs=matrix( betas[j,j.range] ,nrow=1,ncol=j-1)

      print(coefs)
      print(covars[j.range,j.range])
      print(coefs %*% covars[j.range,j.range])
       
    # example:
    #   covars[j,j.range] = [ b31 b32 ] [ sigma1 b21                          ]
    #                                   [ b21    sigma2 +  b21 * sigma1 * b21 ]
      covars[j      ,j.range] = coefs %*% covars[j.range,j.range]
      covars[j.range,j      ] = covars[j, j.range]
    
    # example:
    #   covars[j,j] = sigma[j] + [ b31 b32 ] [ sigma1 b21                          ]  [b31]
    #                                        [ b21    sigma2 +  b21 * sigma1 * b21 ]  [b32]
      covars[j      ,j      ] = sigmas[j] + coefs %*% covars[j.range,j.range] %*% t(coefs)
    }
  }

  return(covars)
}
