/* header file for intQueue */

#ifndef SWEEPTREE_HEADER

#define SWEEPTREE_HEADER

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <R_ext/Utils.h>
#include <R.h>
#include <Rmath.h>
#include "covarTree.h"

/* sweep function */
void VSWP(
  double * v,
  int i,
  int n 
);


/* sweep test function */
void VSWP2(
  double * x,
  int i,
  int n 
);

/* reverse sweep function */
void VRevSWP(
  double * v,
  int i,
  int n 
); 

/* R sweep function interface */
void RVSWP(
  double * v,
  int * i,
  int * n 
);

/* R interface for RMVN2 */
void RMVN2(
    double * sample,
    double * mean,
    double * var,
    int * sizePtr
    ); 

/* R reverse sweep function interface*/
void RVRevSWP(
  double * v,
  int * i,
  int * n 
); 

void printFullMatrix ( double * x, int n, int m );

void printCovarMatrix ( double * x, int k );

void copyCovarMatrix ( double * x, double * y, int k ); 

/* function to copy the contents from a triangular array X to a  matrix Y */
void copyMatrixFromLowerTriangularArray(double * X, double * Y, int n); 

/* function to copy the contents from a matrix X to a triangular array Y */
void copyMatrixToLowerTriangularArray(double * X, double * Y, int n); 

/* sweepTree */
void sweepTree( 
    covarTreePtr x, 
    double * V, 
    int k, 
    double ** matrixCache, 
    int  *   index,      // identify index with row number e.g. (-1,-1,-1,0,1,2)
    double * estimates ,
    bool * M,
    int n
    ); 

/* save parameters */
void saveParameterEstimates( double * V, int k, int i, int * index, double * estimates, bool * M, int df ); 

/* R interface */
void RSweepTree( 
  double * x,          // upper (lower in R) triangular matrix including diag
  int *   M,          // m by p matrix of model parameter inclusions 
  int  * regIndex,   // variables (row indexes) that will be regressed
  double * est,        // p by p matrix of parameter estimates
  int  *   index,      // identify index with row number e.g. (-1,-1,-1,0,1,2)
  int  *   pPtr,       // number of rows/cols in x
  int  *   mPtr,       // number of rows in M 
  int *n                // number of observations used to generate x
); 

/* Gibbs sampling approach to generate a deviate from a MVN dist with mean=mean and var=var */
void ArMVN(                      
	  double *sample,         // sample array of length size
	  double *mean,           // mean array of length size 
	  double *var,            // the lower triangular array of the variance 
	  int size                // length of mean and dim of var 
 ); 

#endif




