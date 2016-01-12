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

void printFullMatrix ( double * x, int n, int m );

void printCovarMatrix ( double * x, int k );

void copyCovarMatrix ( double * x, double * y, int k ); 

/* sweepTree */
void sweepTree( covarTreePtr x, double * V, int k, double ** matrixCache, double * estimates ); 

/* save parameters */
void saveParameterEstimates( double * V, int k, int i, double * estimates ); 




#endif




