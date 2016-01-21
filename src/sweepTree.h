/* header file for intQueue */

#ifndef SWEEPTREE_HEADER

#define SWEEPTREE_HEADER

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "covarTree.h"


/* sweep function */
void VSWP(
  double * v,
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

/* R reverse sweep function interface*/
void RVRevSWP(
  double * v,
  int * i,
  int * n 
); 

void printFullMatrix ( double * x, int n, int m );

void printCovarMatrix ( double * x, int k );

void copyCovarMatrix ( double * x, double * y, int k ); 

/* sweepTree */
void sweepTree( 
    covarTreePtr x, 
    double * V, 
    int k, 
    double ** matrixCache, 
    int  *   index,      // identify index with row number e.g. (-1,-1,-1,0,1,2)
    double * estimates 
    ); 

/* save parameters */
void saveParameterEstimates( double * V, int k, int i, int * index, double * estimates ); 

/* R interface */
void RSweepTree( 
  double * x,          // upper (lower in R) triangular matrix including diag
  int *   M,          // m by p matrix of model parameter inclusions 
  int  * regIndex,   // variables (row indexes) that will be regressed
  double * est,        // p by p matrix of parameter estimates
  int  *   index,      // identify index with row number e.g. (-1,-1,-1,0,1,2)
  int  *   pPtr,       // number of rows/cols in x
  int  *   mPtr        // number of rows in M 
); 

#endif




