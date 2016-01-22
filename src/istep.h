/* header file for intQueue */

#ifndef ISTEP_HEADER

#define ISTEP_HEADER

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <Rmath.h>
#include <R.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>


// print out matrix in column major continuous array form 
void RprintMatrixDouble( double *x , int n, int m ); 


void RIStep(
  double * X,   // observations matrix of dim n by p
  double * S,   // variance of dim p by p 
  double * Beta, //Beta matrix
  int * MIndex, // matrix of missing value indexes m by p 
  int * nPtr,
  int * pPtr,
  int * bPtr
  ); 


void RcholInv(
  double * v,
  double * x,
  int * m,
  int * n 
  ); 

#endif



