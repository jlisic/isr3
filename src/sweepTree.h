/* header file for intQueue */

#ifndef SWEEPTREE_HEADER

#define SWEEPTREE_HEADER

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "covarTree.h"


/* sweep function */
void SWP(
  double **X,             /* The Matrix to work on */
  int k,                  /* The row to sweep */
  int size);              /* The dim. of X */

/* sweep function wrapper for array based matricies */
/* R interface */
void RSWP(
  double * x,
  int * kPtr,
  int * kSizePtr,
  int * sizePtr
); 

/* R interface */
double * ASWP(
  double * x,
  int * kPtr,
  int  kSizePtr,
  int  sizePtr
);

void printCovarMatrix ( double * x, int k );
void copyCovarMatrix ( double * x, double * y, int k ); 

/* sweepTree */
void sweepTree( covarTreePtr x, double * V, int k, double ** matrixCache ); 


#endif




