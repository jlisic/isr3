#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Rmath.h>
#include <R.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
  

void RprintMatrix( double *x , int n, int m ) {

  int i,j ;
  
  Rprintf("\t");
  for( i = 0; i < n; i++) Rprintf("%d\t",i); 
  Rprintf("\n");

  for( i = 0; i < n; i++) {
    Rprintf("%d\t",i); 
    for( j = 0; j < n; j++) Rprintf("%5.3f, ", x[j*n + i]); 
    Rprintf("\n");
    }


}



/*
 * 0. Prep work
 *
 * Create Sweep Tree (Function Done)
 * Subset X by missing items 
 *
 * *** P Step ***
 *
 * 1. Run Sweep Tree (Function Done)
 *
 * 2. Generate deviates per Sweep Tree 
 *
 * *** I Step ***
 *
 * 3. Rebuild Covar matrix (Function Done)
 *
 * 4. Calculate XB of missing
 *
 * do in fortran
 *
 * 5. Calculate cholesky decomp and inverse
 *
 * do in fortran
 *
 * 6. Reverse sweep to get a part of the conditonal mean, and the cond var
 *
 * 7. Do some matrix multiplication to get conditional means
 *
 * 8. Generate deviates 
 *
 * Go back to 1 until done
 *
 * 0 Clean up
 *
*/

void RIStep(
  double * X,  // observations matrix of dim n by p
  double * S,  // variance of dim p by p 
  int * MIndex   // matrix of missing value indexes m by p 
  int * n,
  int * p,
  int * m 
  ) {

  int errorCode;
  int i,j;
  int n = *nPtr;
  int p = *pPtr;
  int m = *mPtr;
  int rowIndexI, rowIndexJ;

  int * M;
  double * XB;

  //double alpha = 1.0;

  Rprintf("Input:\n");
  Rprintf("X:\n");
  RprintMatrix(X, n, p);
  Rprintf("S:\n");
  RprintMatrix(S, p, p);

  // copy over missing
  M = calloc( sizeof(double), m * p);
  for( i = 0; i < m; i++)
  {
    rowIndexI = MIndex[i];
    for( j = 0; j < p; j++) M[rowIndex*p + j];
  } 

  // missing XBeta values 
  XB = calloc( sizeof(double), m * p);
  for( i = 0; i < m; i++) {
    rowIndexI = i*p;
    for( j = 0; j < p; j++) { 
      rowIndexJ = j*p;
      for( k = 0; k < p; k++) XB[rowIndexI+k] += M[rowIndexI + k] * Beta[rowIndexJ + k];
    }
  }
  
  Rprintf("XB:\n");
  RprintMatrix(XB, m, p);

  free(M);
  free(XB);

/*
  // calculate inverse
  F77_CALL(dpotri)(
      "L",        // UPLO
      n,         // N 
      v,          // A 
      n,         // LDA
      &errorCode  // error code
      );

  if(errorCode != 0) Rprintf("LAPACK dpotrs failed with error code = %d\n", errorCode);

  Rprintf("CholInv output:\n");
  RprintMatrix(v, *n, *n);
   
  // multiply cholesky inverse by matrix 
  // v := (alpha) v * x
  F77_CALL(dtrmm)(
      "R",        // SIDE 
      "L",        // UPLO
      "N",        // TRANSA
      "N",        // DIAG (Unit Traingular, i.e. diag = 1) 
       m ,        // M rows of B
       n ,        // N cols fo B
       &alpha ,   // ALPHA
       v,         // A
       n ,        // LDA rows of A (in this context)
       x ,        // B
       m          // LDB 
      );
*/
  return;
}


void RcholInv(
  double * v,
  double * x,
  int * m,
  int * n 
  ) {

  int errorCode;
  int i,j;
  double alpha = 1.0;

  Rprintf("Input:\n");
  RprintMatrix(v, *n, *n);



  // calculate inverse
  F77_CALL(dpotri)(
      "L",        // UPLO
      n,         // N 
      v,          // A 
      n,         // LDA
      &errorCode  // error code
      );

  if(errorCode != 0) Rprintf("LAPACK dpotrs failed with error code = %d\n", errorCode);

  Rprintf("CholInv output:\n");
  RprintMatrix(v, *n, *n);
   
  // multiply cholesky inverse by matrix 
  // v := (alpha) v * x
  F77_CALL(dtrmm)(
      "R",        // SIDE 
      "L",        // UPLO
      "N",        // TRANSA
      "N",        // DIAG (Unit Traingular, i.e. diag = 1) 
       m ,        // M rows of B
       n ,        // N cols fo B
       &alpha ,   // ALPHA
       v,         // A
       n ,        // LDA rows of A (in this context)
       x ,        // B
       m          // LDB 
      );
  return;
}

