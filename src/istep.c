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


void RprintMatrixInt( int *x , int n, int m ) {

  int i,j ;
  
  Rprintf("\t");
  for( i = 0; i < m; i++) Rprintf("%d\t",i); 
  Rprintf("\n");

  for( i = 0; i < n; i++) {
    Rprintf("%d\t",i); 
    for( j = 0; j < m; j++) Rprintf("%d, ", x[j*n + i]); 
    Rprintf("\n");
    }

  return;
}


void RprintMatrixBool( bool *x , int n, int m ) {

  int i,j ;
  
  Rprintf("\t");
  for( i = 0; i < m; i++) Rprintf("%d\t",i); 
  Rprintf("\n");

  for( i = 0; i < n; i++) {
    Rprintf("%d\t",i); 
    for( j = 0; j < m; j++) Rprintf("%s, ", x[j*n + i] ? "T" : "F" ); 
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
 * *** Assumptions ***
 * X is fully observed for all variables before 
 * column maxObsIndex.
 *
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

void Risr(
  double * X,       // observations matrix of dim n by p
  int * MIndicator, // matrix of missing value indicators n by p 
  int * M,          // size of est (minus one column) but binary (except it's in integer due to R)
  int * regIndex,
  int * index,
  int * nPtr,
  int * pPtr,
  int * bPtr,       // the length of beta
  double * S,       // variance of dim p by p [ only returned ] 
  double * est      // estimate               [ only returned ]
  ) {

  int errorCode;  // lapack error code
  int i,j,k;
  int n = *nPtr;
  int p = *pPtr;
  int b = *bPtr;
  int maxObsIndex = p - b; // maximum observed index

  bool * MBool;      // converts M into boolean for covartree
  double ** cache;   // an array of covariance matricies used in the sweepTree.

  double * estRebuild; // a reduced est for the rebuild covar function

  int rowIndexI, rowIndexJ;
  bool * observed;  // vector of length n, T if observed F if not observed 

  double alpha = 1.0; // used for fortran BLAS calls
  double beta  = 0.0; // used for fortran BLAS calls

  double * Beta; // Beta coefficients 
  double * XB;   // X times Beta
  double * SA;   // S in tri array mode
  double * Y;    // final storage
  double * XX;   // X**T X

  int c; // place holder for index in multiplication of an imputation step
  
  int cacheSize = 0; // used to determine cache size for sweep tree

  /*
  Rprintf("Input:\n");
  Rprintf("X:\n");
  RprintMatrixDouble(X, n, p);
  Rprintf("MIndicator:\n");
  RprintMatrixInt(MIndicator, n, p);
  Rprintf("M:\n");
  RprintMatrixInt(M, p, b);
  */

  /****************** Step 0 *******************/

  GetRNGstate();

  /************* P Step Setup *********/
  
  // sweep tree for P step 
  covarTreePtr myTree = NULL;

  // create MBool for the conversion
  MBool = calloc( sizeof( bool *), p*b);
  for(j=0;j<p*b;j++) MBool[j] = (bool) M[j];

  estRebuild = calloc( sizeof(double), (b+1)*b);

  // create tree
  //myTree = createCovarTree( NULL, MBool, p, regIndex[0], 0, &cacheSize); 
  for(i=0;i<b;i++) {
    myTree = createCovarTree(myTree, &(MBool[i*p]), p, regIndex[i], 0, &cacheSize); 
  }

  Rprintf("M:\n");
  RprintMatrixInt(M, p, b);
  Rprintf("MBool:\n");
  RprintMatrixBool(MBool, p, b);
  

  // allocate space for the cache
  cache = calloc( sizeof(double *) , cacheSize+1 );

  /************* I Step Setup *********/

  // identify completely observered rows in X (startup cost )
  observed = calloc( sizeof(bool), n );
  
  // allocate our final data 
  Y = calloc( sizeof(double), n*p );
  
  // pre allocate SA
  SA = calloc( sizeof(double), (p*(p+1))/2);
  
  // pre allocate Beta 
  Beta = calloc( sizeof(double), p*b );

  // pre allocate XB
  XB = calloc( sizeof(double), n * b);

  // pre allocate X**T X
  XX = calloc( sizeof(double), p*p);  

  // find the number of missing values and their index (startup cost)
  // We set the value to -1 if Mindex[i] 
  for( i = 0; i < n; i++)
  {
    for( j = 0; j < p; j++) 
      if( MIndicator[i*p + j] == 1) break;

    if( j < p ) observed[i] = true;
  } 

  //for( i = 0; i < n; i++) Rprintf("(%d) %s\n", i, observed[i] ? "True" : "False");
 
  /****************** Step 1 *******************/
  /* Calculate Beta and sigma of the           */
  /* conditional distributions.                */
  /****************** Step 1 *******************/
               /* and */ 
  /****************** Step 2 *******************/
  /* Generate deviates per the sweep tree      */
  /****************** Step 2 *******************/
  
  // multiply matrix by matrix 
  // C := (alpha) op(A) * op(B) + (beta) * C 
  F77_CALL(dgemm)(
      "T",        // TRANSA 
      "N",        // TRANSB 
       &p,        // M (rows of op(A) 
       &p,        // N (cols of B) 
       &n,        // K (cols of A)
       &alpha,    // ALPHA
       X,         // A
       &n,        // LDA rows of A (in this context)
       X,         // B
       &n,        // LDB 
       &beta,     // Beta (scalar for BLAS) 
       XX,        // C 
       &p         // LDC 
      );
  
//Rprintf("XX:\n");
//RprintMatrixDouble(XX, p, p);
  copyMatrixToLowerTriangularArray(XX, SA, p); 
  
  sweepTree(myTree, SA, p, cache, index, est, MBool,n);
  
  Rprintf("Est:\n");
  RprintMatrixDouble(est, p+1, b);
  
 


  /****************** Step 3 *******************/
  /* Rebuild covariance                        */
  /****************** Step 3 *******************/

  // subset est by the last b variables
  /*  est example:
   *  b = 3
   *  p = 6
   *
   *       B_1 CoVar1 CoVar2 Var1 Var2 Var3 Sigma
   * Var1  0   1      2      3    4    5    6
   * Var2  7   8      9     10   11   12   13
   * Var3 14  15     16     17   18   19   20
  */ 

  for(i=0; i < b; i++) 
    for(j=b; j < p +1; j++) estRebuild[i*(b+1)+j-b ] = est[i*(p+1) + j]; 
  
  Rprintf("Est Rebuilt:\n");
  RprintMatrixDouble(estRebuild, b+1, b);

  rebuildCovar( estRebuild, SA, b);
  
  Rprintf("Sigma Rebuilt:\n");
  printCovarMatrix(SA,b);  


  /****************** Step 4 *******************/
  /* Calculate XB = X * Beta                   */
  /****************** Step 4 *******************/
  
  for(i=0; i < b; i++) 
    for(j=0; j < p ; j++) Beta[i*p+j] = est[i*(p+1) + j]; 
  
  Rprintf("X:\n");
  RprintMatrixDouble(X, n, p);
  
  Rprintf("B:\n");
  RprintMatrixDouble(Beta, p, b);

  // missing XBeta values 
  
  
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
       Beta,         // B
       &p,        // LDB 
       &beta,     // BETA 
       XB,        // C 
       &n         // LDC 
      );
  
  Rprintf("XB:\n");
  RprintMatrixDouble(XB, n, b);


  /****************** Step 5 *******************/
  /* Calculate Inverse of Matrix               */
  /* via Cholesky Factorization                */
  /****************** Step 5 *******************/
  
//  // calculate cholesky decomposition 
//  F77_CALL(dpotrf)(
//      "L",        // UPLO
//      &p,         // N (dim of A)
//      S,          // A 
//      &p,         // LDA
//      &errorCode  // error code
//      );
//
//  if(errorCode != 0) Rprintf("LAPACK dpotrf failed with error code = %d\n", errorCode);
//  
//  Rprintf("S.chol:\n");
//  RprintMatrixDouble(S, p, p);
//
//  // calculate inverse
//  f77_CALL(dpotri)(
//      "L",        // UPLO
//      &p,         // N (dim of A)
//      S,          // A 
//      &p,         // LDA
//      &errorCode  // error code
//      );
//
//  if(errorCode != 0) Rprintf("LAPACK dpotrs failed with error code = %d\n", errorCode);
//
//  Rprintf("S.inv:\n");
//  RprintMatrixDouble(S, p, p);
// 
//
//  /****************** Step 6 *******************/
//  /* Get conditional means and variances       */ 
//  /****************** Step 6 *******************/
//
//  // (x - XB)(Sigma.inv_{ii} - Sigma.inv_{-1,-1} 
//
//  for( i = 0; i < n*p; i++) Y[i] = 0;
//
//  for( i = 0; i < p; i++) {
//    copyMatrixToLowerTriangularArray(S, SA, p); 
//  //Rprintf("SA:\n");
//  //printCovarMatrix(SA,p);  
//
//
//    VRevSWP( SA, i, p); 
//  
//    /*********************************************
//     * 0 
//     * 1  5 
//     * 2  6  9
//     * 3  7 10 12
//     * 4  8 11 13 14
//     ********************************************/
//    Rprintf("X: %d\n", i);
//    RprintMatrixDouble(X, n, p);
//    Rprintf("XB: %d \n", i);
//    RprintMatrixDouble(XB, n, b);
//    Rprintf("SA %d:\n", i);
//    printCovarMatrix(SA,p);  
//  
//    c = i;
//    for(j=0; j <i; j++, c += p - j) {
//      for(k=0; k <n; k++) {
//        printf("j = %d, i = %d, Y[%d] = X[%d] - XB[%d]) * SA[%d] = %f\n",  
//            j, i, i*n+k, j*n+k, j*n+k, c, (X[j*n +k] - XB[j*n +k]) * SA[c]);  
//
//        Y[i*n+k] += (X[j*n +k] - XB[j*n +k]) * SA[c];  
//      }
//    }
//    printf("c = %d\n", c);
//    for(j=i+1; j <p; j++) { 
//      c++; 
//      for(k=0; k <n; k++) {
//        printf("j = %d, i = %d, Y[%d] = X[%d] - XB[%d]) * SA[%d] = %f\n",  
//            j, i, i*n+k, j*n+k, j*n+k, c, (X[j*n +k] - XB[j*n +k]) * SA[c]);  
//        Y[i*n+k] += (X[j*n +k] - XB[j*n +k]) * SA[c];  
//      }
//    }
//  }
//  
//  RprintMatrixDouble(Y, n, p);

  PutRNGstate();

  /* P Step */
  free(MBool);
  free(cache);
  deleteCovarTree(myTree);
 
  /* I Step */ 
  free(Y);
  free(SA);
  free(XB);
  free(Beta);
  free(observed);

  /* General */
  free(XX);

  return;
}

