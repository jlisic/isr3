#include "istep.h"


void RprintMatrixDouble( double *x , int n, int m ) {

  int i,j ;
  
  Rprintf("\t");
  for( i = 0; i < m; i++) Rprintf("%d\t",i); 
  Rprintf("\n");

  for( i = 0; i < n; i++) {
    Rprintf("%d\t",i); 
    for( j = 0; j < m; j++) Rprintf("%10.7f, ", x[j*n + i]); 
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


void RprintMatrixBool( int *x , int n, int m ) {

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
  double * est,     // estimate               [ only returned ]
  int * maxIter,    // number of iters
  int * sampleRate, // how often to sample
  int * miPtr       // number of multiple imputations
  ) {

  int errorCode;  // lapack error code
  int i,j,k,m;
  int n = *nPtr;  // number of records (observations, observed or not)
  int p = *pPtr;  // number of variables (including intercept)
  int b = *bPtr;  // number of variables being imputed for
  int mi = *miPtr;

  double ** cache;   // an array of covariance matricies used in the sweepTree.

  double * estRebuild; // a reduced est for the rebuild covar function

  double alpha = 1.0; // used for fortran BLAS calls
  double beta  = 0.0; // used for fortran BLAS calls

  int * MIndicator_tmp;   // pointer offset for MIndicator 
  
  double * Beta; // Beta coefficients 
  double * X_tmp;   // X times Beta
  double * XB;   // X times Beta
  double * XB_tmp;   // X times Beta
  double * SA;   // S in tri array mode
  double * SA_tmp;   // S in tri array mode
  double * XX;   // X**T X
  double * Y;

  double s2; //used to store the variance
  int c; // place holder for index in multiplication of an imputation step
  int cacheSize = 0; // used to determine cache size for sweep tree
  int one = 1; // used for dgemv call
  int iter; // MCMC iter


  /****************** Step 0 *******************/

  GetRNGstate();

  /************* P Step Setup *********/
  
  // sweep tree for P step 
  covarTreePtr myTree = NULL;

  // create MBool for the conversion
  estRebuild = calloc( sizeof(double), (b+1)*b);

  // create tree
  //myTree = createCovarTree( NULL, MBool, p, regIndex[0], 0, &cacheSize); 
  for(i=0;i<b;i++) {
    myTree = createCovarTree(myTree, &(M[i*p]), p, regIndex[i], 0, &cacheSize); 
  }

  // allocate space for the cache
  cache = calloc( sizeof(double *) , cacheSize+1 );

  /************* I Step Setup *********/
  // identify completely observered rows in X (startup cost )
  
  // allocate our final data 
  Y = calloc( sizeof(double), n );
  
  // pre allocate SA
  SA = calloc( sizeof(double), (p*(p+1))/2);
  SA_tmp = calloc( sizeof(double), (p*(p+1))/2);
  
  // pre allocate Beta 
  Beta = calloc( sizeof(double), p*b );

  // pre allocate XB and copy over X
  XB = calloc( sizeof(double), n * p);
  for(i=0; i<n*p; i++) XB[i] = X[i];

  // pre allocate X**T X
  XX = calloc( sizeof(double), p*p);  


  for( iter=0; iter < ((*maxIter) + mi * (*sampleRate) ); iter++) {
 
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
    
    copyMatrixToLowerTriangularArray(XX, SA, p); 
    
    sweepTree(myTree, SA, p, cache, index, est, M,n);

  
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
      for(j=0; j <= b; j++) estRebuild[i*(b+1)+j ] = est[i*(p+1)+p-b+j]; 
  
    rebuildCovar( estRebuild, SA, b);
    
    /****************** Step 4 *******************/
    /* Calculate XB = X * Beta                   */
    /****************** Step 4 *******************/
    
    for(i=0; i < b; i++) 
      for(j=0; j < p ; j++) Beta[i*p+j] = est[i*(p+1) + j]; 
  
    // missing XBeta values 
    
    // multiply matrix by vector, updating the vector each step 
    // y := (alpha) A * x + (beta) * y 
    for(i=0,m=p-b; i<b; i++,m++) {
    F77_CALL(dgemv)(
        "N",        // TRANS 
         &n,        // M (rows of A) 
         &p,        // N (cols of A) 
         &alpha,    // ALPHA
         XB,         // A
         &n,        // LDA rows of A (in this context)
         &(Beta[p*i]),         // B
         &one,        // INCX 
         &beta,     // BETA 
         Y,     // Y 
         &one         // INCY 
        );
      for(j=0; j < n; j++) XB[m*n +j] = Y[j];
    } 
  
  
    /****************** Step 5 *******************/
    /* Calculate Inverse of Matrix               */
    /* via Cholesky Factorization                */
    /****************** Step 5 *******************/
    
    // calculate cholesky decomposition in packed form
    F77_CALL(dpptrf)(
        "L",        // UPLO
        &b,         // N (dim of A)
        SA,          // A 
        &errorCode  // error code
        );
  
    if(errorCode != 0) Rprintf("LAPACK dpotrf failed with error code = %d\n", errorCode);
  
    // calculate inverse
    F77_CALL(dpptri)(
        "L",        // UPLO
        &b,         // N (dim of A)
        SA,          // A 
        &errorCode  // error code
        );
  
    if(errorCode != 0) Rprintf("LAPACK dpotrs failed with error code = %d\n", errorCode);
 

    /****************** Step 6 *******************/
    /* Get conditional means and variances       */ 
    /****************** Step 6 *******************/
    // (x - XB)(Sigma.inv_{ii} - Sigma.inv_{-1,-1} 

  
    // create offsets 
    X_tmp = &(X[n*(p-b)]);
    XB_tmp = &(XB[n*(p-b)]);
    MIndicator_tmp = &(MIndicator[n*(p-b)]);
    
   
    for( i = 0; i < b; i++) {
      // make a copy of SA
      for( j = 0; j < (b*(b+1))/2; j++) SA_tmp[j] = SA[j];
      
      // calculate variance for conditional distribution 
      VRevSWP( SA_tmp, i, b); 
  
      /*********************************************
       * example indexing for SA_tmp
       * 0 
       * 1  5 
       * 2  6  9
       * 3  7 10 12
       * 4  8 11 13 14
       ********************************************/
  
      /* there is probably a better way to rewrite this */ 
      s2 = sqrt( -1 * SA_tmp[i*b - (i*(i-1))/2]);
 
      // write XB + e to column i of X_tmp for missing values. 
      for(k=0; k <n; k++) 
        if(MIndicator_tmp[i*n+k] == 0) 
          X_tmp[i*n+k] = XB_tmp[i*n+k] + norm_rand() * s2;
        
     
      // add (X - XB) * invSigma[notX,x] to X_tmp 
      c = i;
      for(j=0; j <i; j++, c += b - j) {
        s2 = SA_tmp[c];
        for(k=0; k <n; k++) 
          if(MIndicator_tmp[i*n+k] == 0) 
            X_tmp[i*n+k] += (X_tmp[j*n +k] - XB_tmp[j*n +k]) * s2; 
          
      }
      for(j=i+1; j <b; j++) { 
        c++; 
        s2 = SA_tmp[c];
        for(k=0; k <n; k++) 
          if(MIndicator_tmp[i*n+k] == 0) 
            X_tmp[i*n+k] += (X_tmp[j*n +k] - XB_tmp[j*n +k]) * s2;  
          
      }
    }

//    Rprintf("X_tmp %d\n",iter);
//    RprintMatrixDouble( X_tmp , n, b ); 

    // save result
    if( iter >= *maxIter) { 
      if( (iter - *maxIter) % (*sampleRate) == 0) {
        printf("iter = %d\n", iter);
        for(i=0;i<n*b;i++) S[i] = X_tmp[i];
        S = &(S[n*b]);
      }
    }
  }
    
  X_tmp = NULL;
  XB_tmp = NULL;
  MIndicator_tmp = NULL;

  PutRNGstate();

  /* P Step */
  free(cache);
  deleteCovarTree(myTree);
 
  /* I Step */ 
  free(SA);
  free(SA_tmp);
  free(XB);
  free(Beta);
  free(Y);

  /* General */
  free(XX);

  return;
}

