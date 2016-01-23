#include "istep.h"



void RprintMatrixDouble( double *x , int n, int m ) {

  int i,j ;
  
  Rprintf("\t");
  for( i = 0; i < m; i++) Rprintf("%d\t",i); 
  Rprintf("\n");

  for( i = 0; i < n; i++) {
    Rprintf("%d\t",i); 
    for( j = 0; j < m; j++) Rprintf("%5.2f, ", x[j*n + i]); 
    Rprintf("\n");
    }

  return;
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
 * 6. Get conditional means and variances 
 *
 * 8. Generate deviates 
 *
 * Go back to 1 until done
 *
 * 0 Clean up
 *
*/

void RIStep(
  double * X,   // observations matrix of dim n by p
  double * S,   // variance of dim p by p 
  double * B, //Beta matrix
  int * MIndicator, // matrix of missing value indicators n by p 
  int * nPtr,
  int * pPtr,
  int * bPtr
  ) {

  int errorCode;
  int i,j,k;
  int n = *nPtr;
  int p = *pPtr;
  int b = *bPtr;

  int rowIndexI, rowIndexJ;
  bool * observed;  // vector of length n, T if observed F if not observed 

  double alpha = 1.0;
  double beta  = 0.0;

  double * XB;
  double * SA; // S in tri array mode
  double * Y;  // final storage

  int c; // place holder for index in multiplication of an imputation step

  Rprintf("Input:\n");
  Rprintf("X:\n");
  RprintMatrixDouble(X, n, p);
  Rprintf("S:\n");
  RprintMatrixDouble(S, p, p);


  /****************** Step 0 *******************/


  // identify completely observered rows in X (startup cost )
  observed = calloc( sizeof(bool), n );
  
  // allocate our final data 
  Y = calloc( sizeof(double), n*p );
  
  // pre allocate SA
  SA = calloc( sizeof(double), (p*(p+1))/2);

  // find the number of missing values and their index (startup cost)
  // We set the value to -1 if Mindex[i] 
  for( i = 0; i < n; i++)
  {
    for( j = 0; j < p; j++) 
      if( MIndicator[i*p + j] == 1) break;

    if( j < p ) observed[i] = true;
  } 

  for( i = 0; i < n; i++) Rprintf("(%d) %s\n", i, observed[i] ? "True" : "False");
 

  /****************** Step 4 *******************/
  /* Calculate XB = X * Beta                   */
  /****************** Step 4 *******************/

  // missing XBeta values 
  XB = calloc( sizeof(double), n * b);
  
  
  // multiply matrix by matrix 
  // C := (alpha) A * B + (beta) * C 
  F77_CALL(dgemm)(
      "N",        // TRANSA 
      "N",        // TRANSB 
       &n,        // M (rows of A) 
       &b,        // N (cols of B) 
       &p,        // K (cols of A)
       &alpha,    // ALPHA
       X,         // A
       &n,        // LDA rows of A (in this context)
       B,         // B
       &b,        // LDB 
       &beta,     // BETA 
       XB,        // C 
       &n         // LDC 
      );
  
  Rprintf("XB:\n");
  RprintMatrixDouble(XB, n, p);


  /****************** Step 5 *******************/
  /* Calculate Inverse of Matrix               */
  /* via Cholesky Factorization                */
  /****************** Step 5 *******************/
  
  // calculate cholesky decomposition 
  F77_CALL(dpotrf)(
      "L",        // UPLO
      &p,         // N (dim of A)
      S,          // A 
      &p,         // LDA
      &errorCode  // error code
      );

  if(errorCode != 0) Rprintf("LAPACK dpotrf failed with error code = %d\n", errorCode);
  
  Rprintf("S.chol:\n");
  RprintMatrixDouble(S, p, p);

  // calculate inverse
  F77_CALL(dpotri)(
      "L",        // UPLO
      &p,         // N (dim of A)
      S,          // A 
      &p,         // LDA
      &errorCode  // error code
      );

  if(errorCode != 0) Rprintf("LAPACK dpotrs failed with error code = %d\n", errorCode);

  Rprintf("S.inv:\n");
  RprintMatrixDouble(S, p, p);
 

  /****************** Step 6 *******************/
  /* Get conditional means and variances       */ 
  /****************** Step 6 *******************/

  // (x - XB)(Sigma.inv_{ii} - Sigma.inv_{-1,-1} 

  for( i = 0; i < n*p; i++) Y[i] = 0;

  for( i = 0; i < p; i++) {
    copyMatrixToLowerTriangularArray(S, SA, p); 
  //Rprintf("SA:\n");
  //printCovarMatrix(SA,p);  


    VRevSWP( SA, i, p); 
  
    /*********************************************
     * 0 
     * 1  5 
     * 2  6  9
     * 3  7 10 12
     * 4  8 11 13 14
     ********************************************/
    Rprintf("X: %d\n", i);
    RprintMatrixDouble(X, n, p);
    Rprintf("XB: %d \n", i);
    RprintMatrixDouble(XB, n, b);
    Rprintf("SA %d:\n", i);
    printCovarMatrix(SA,p);  
  
    c = i;
    for(j=0; j <i; j++, c += p - j) {
      for(k=0; k <n; k++) {
        printf("j = %d, i = %d, Y[%d] = X[%d] - XB[%d]) * SA[%d] = %f\n",  
            j, i, i*n+k, j*n+k, j*n+k, c, (X[j*n +k] - XB[j*n +k]) * SA[c]);  

        Y[i*n+k] += (X[j*n +k] - XB[j*n +k]) * SA[c];  
      }
    }
    printf("c = %d\n", c);
    for(j=i+1; j <p; j++) { 
      c++; 
      for(k=0; k <n; k++) {
        printf("j = %d, i = %d, Y[%d] = X[%d] - XB[%d]) * SA[%d] = %f\n",  
            j, i, i*n+k, j*n+k, j*n+k, c, (X[j*n +k] - XB[j*n +k]) * SA[c]);  
        Y[i*n+k] += (X[j*n +k] - XB[j*n +k]) * SA[c];  
      }
    }
  }
  
  RprintMatrixDouble(Y, n, p);

  //for(i=0;i<p;i++) {
  //VRevSWP( SA, i, p); 
  //
  //
  // //store new missing in X
  //for(j=0; j <n; j++) { 
  //
  //}
  //}


  free(SA);
  free(observed);
  free(XB);


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
  RprintMatrixDouble(v, *n, *n);



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
  RprintMatrixDouble(v, *n, *n);
   
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

