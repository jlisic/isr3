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
void rebuildTree( 
  double * x, 
  double * y, 
  int * index, 
  int p, 
  int m
);


/* R interface */
void RRebuildTree( 
  double * x,          // upper (lower in R) triangular matrix including diag
  int *   M,          // m by p matrix of model parameter inclusions 
  int  * regIndex,   // variables (row indexes) that will be regressed
  double * est,        // p by p matrix of parameter estimates
  int  *   index,      // identify index with row number e.g. (-1,-1,-1,0,1,2)
  int  *   pPtr,       // number of rows/cols in x
  int  *   mPtr        // number of rows in M 
); 

#endif




