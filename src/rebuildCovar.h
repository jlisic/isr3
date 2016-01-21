/* header file for intQueue */

#ifndef REBUILDCOVAR_HEADER

#define REBUILDCOVAR_HEADER

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "sweepTree.h"


/* 
 * A function to rebuild the covariance matrix from
 * a linear model specification.
 *
 * Input:
 *
 * x: is an m by (p+1) (double) matrix of covariates and variances (col major)
 * y: is an m*(m+1)/2 double array
 *
 * Output:
 *
 * y: is an m by m covariance matrix in triangular + diag form  
 */
void rebuildCovar( 
  double * x, 
  double * y, 
  int p
);


/* R interface */
void RRebuildCovar( 
  double * x,          // upper (lower in R) triangular matrix including diag
  double * est,        // p by p matrix of parameter estimates
  int  *   nPtr        // number of rows/cols in x
); 

#endif




