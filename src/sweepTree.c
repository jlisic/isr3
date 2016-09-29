#include "sweepTree.h"




/* copy covar matrix function */
void copyCovarMatrix ( double * x, double * y, int k ) {
  int i;
  for(i=0; i < (k*(k+1))/2; i++) x[i] = y[i];
  return;
}




/* copy matrix function */
void printFullMatrix ( double * x, int n, int m ) {
  int i, j;

  for( i = 0; i < m; i++) printf("%d\t", i); 
  printf("\n");
  
  for( i = 0; i < n; i++) {
    for( j = 0; j < m; j++) {
      printf("%5.3f, ", x[m*i + j]);
    }
    printf("\n");
  }
  return;
}




/* print covar matrix function */
void printCovarMatrix ( double * x,  int n ) {
  int i, j, m;
  m = 0;
  for( i = 0; i < n; i++) {
    for( j = 0; j < n; j++) {
      if( j < i) {
        //printf("            ");
        printf("\t ");
      } else {
        printf("%5.4f,", x[m]);
        m++;
      }
    }
    printf("\n");
  }
  return;
}




/* perform sweep */
void VSWP(
  double * v,
  int i,
  int n 
) {
  
  int j,k,m,N;
  double * x = NULL;
  double * y = NULL;
  double * z = NULL;
  double xii;
  int * row;
 
  N = (n*(n+1))/2;

  if( i == 0 ) { 
    xii = v[0];
  } else {
    xii = v[i*n - (i*(i-1))/2];
  }
  
  if( xii < 10e-20 ) {
    printf("SWP: singular matrix.\n");
    return;
  }

  // allocate if needed
  if( x == NULL) {
    x = calloc( sizeof(double *), N );
    y = calloc( sizeof(double *), N );
    z = calloc( sizeof(double *), n);
    row = calloc(sizeof( int * ), n);

    // reduce row calculations
    row[0] = 0;
    for( j = 1; j < n; j++) row[j] += row[j-1] + n - j + 1;
  }

  xii = 1/xii;
 
  /* create easy to access index */ 
  for(j=0; j<i; j++) {
    z[j] = v[row[j] - j +i]; 
  }
  for(   ; j<n; j++) {
    z[j] = v[row[i] - i + j]; 
  }

  /* fill the rest of the buffer */ 
  m = 0;
  for( j = 0; j < n; j++ ) {
    for( k = j; k < n; k++ ) {
    x[m] = z[k] * xii; 
    m++; 
    }
  }
  
  /* create a second buffer and multiply by xii */
  m = 0;
  for( j = 0; j < n; j++ ) {
    for( k = j; k < n; k++ ) {
    y[m] = z[j]; 
    m++; 
    }
  }

  /* perform my math */
  for(j=0; j<N; j++)  {
    v[j] = v[j] - y[j] * x[j];
  };

  /* for the case of i == j replace x[j]*xii by xii */
  for(j=0; j<i; j++) {
    v[row[j] - j + i] = x[j];
  }
  for(   ; j<n; j++) {
    v[row[i] -i + j] = x[j];
  }
  v[ row[i] ] = -1 * xii;
    
  // free needed 
  free(z); 
  free(x);
  free(y);
  free(row);
  
  return;
}




/* perform reverse sweep */
void VRevSWP(
  double * v,
  int i,
  int n 
) {
  
  int j,k,m,N;
  double * x = NULL;
  double * y = NULL;
  double * z = NULL;
  double xii;
  int * row;
 
  N = (n*(n+1))/2;

  if( i == 0 ) { 
    xii = v[0];
  } else {
    xii = v[i*n - (i*(i-1))/2];
  }
  
  if( fabs(xii) < 10e-20 ) {
    printf("SWP: singular matrix.\n");
    return;
  }

  // allocate if needed
  if( x == NULL) {
    x = calloc( sizeof(double *), N );
    y = calloc( sizeof(double *), N );
    z = calloc( sizeof(double *), n);
    row = calloc(sizeof( int * ), n);

    // reduce row calculations
    row[0] = 0;
    for( j = 1; j < n; j++) row[j] += row[j-1] + n - j + 1;
  }

  xii = 1/xii;
 
  /* create easy to access index */ 
  for(j=0; j<i; j++) {
    z[j] = v[row[j] - j +i]; 
  }
  for(   ; j<n; j++) {
    z[j] = v[row[i] - i + j]; 
  }

  /* fill the rest of the buffer */ 
  m = 0;
  for( j = 0; j < n; j++ ) {
    for( k = j; k < n; k++ ) {
    x[m] = z[k] * xii; 
    m++; 
    }
  }
  
  /* create a second buffer and multiply by xii */
  m = 0;
  for( j = 0; j < n; j++ ) {
    for( k = j; k < n; k++ ) {
    y[m] = z[j]; 
    m++; 
    }
  }


  /* perform my math */
  for(j=0; j<N; j++)  v[j] = v[j] - y[j] * x[j];

  /* for the case of i == j replace x[j]*xii by xii */
  for(j=0; j<i; j++) {
    v[row[j] - j + i] = -1 * x[j];
  }
  for(   ; j<n; j++) {
    v[row[i] -i + j] = - 1 * x[j];
  }
  v[ row[i] ] = -1 * xii;
    
  // free needed 
  free(z); 
  free(x);
  free(y);
  free(row);
  
  return;
}




/* R interface to perform reverse sweep */
void RVSWP(
  double * v,
  int * i,
  int * n 
) {
  VSWP(v, *i, *n); 
  return;
}




/* R interface to perform reverse sweep */
void RVRevSWP(
  double * v,
  int * i,
  int * n 
) {
  VRevSWP(v, *i, *n); 
  return;
}




/* copy matrix to lower triangular array */
/* X is the matrix and Y is the lower triangular array */
void copyMatrixToLowerTriangularArray(double * X, double * Y, int n) {

  int i,j,k,indexI;

  k=0;

  for(i=0;i<n;i++) {
    indexI=i*n;
    for(j=i;j<n;j++) Y[k++] = X[indexI +j];
  }

  return;
}




/* copy matrix from lower triangular array */
/* X is the triangular array and Y is the matrix */
void copyMatrixFromLowerTriangularArray(double * X, double * Y, int n) {

  int i,j,k,indexI;
  k = 0;
  for(i=0;i<n;i++) {
    indexI=i*n;
    for(j=i;j<n;j++) Y[indexI+j] = X[k++];
  }

  return;
}


/* get the index from an upper triangular square matrix form a row and columns 
 * r - row (starting at 0)
 * c - column (starting at 0)
 * n - dim of square matrix
 * */
static inline int rc2ut( int row, int col, int n) {

  int tmp;

  if( col < row ){
    tmp =  col;
    col = row;
    row = tmp;
  } 

  return( n * row - row*(row - 1)/ 2 + col - row );   
}



/* save Parameters */
void saveParameterEstimates( 
    double * V, 
    int k,             // variabls + parameters
    int i,             // index from index to save estimate 
    int * index, 
    double * estimates,
    int * M,
    int df 
    ) {

  int j,o,p,u;
  int v = 0;
  int m = 0;
  int l = 0;
  
  //int n = (k*(k+1))/2; // dim for lower diagonal including diagonal
  double * sample; 
  double * mean;
  double * var;
  double conditionalVar;
  int * usedColumn; 

  if( M != NULL ) { 

    usedColumn = calloc( i+1, sizeof( int ) );
    
    // figure out how many covariates there are 
    for(j=0; j < i; j++) {
      if( M[k*index[i]+j] ) {
        usedColumn[m] = j;  
        m++;
      }
    }

    var = calloc( sizeof( double ), (m*(m+1))/2 );
    mean  = calloc( sizeof( double ), m);
    sample = calloc( sizeof( double ), m);

    //conditionalVar = -1.0 * V[i*k - (i)*(i-1)/2]/rchisq(df - m); 
    conditionalVar = -1.0 * V[i*k - (i)*(i-1)/2]/(double)(df - m); 

    for(j=0; j < m; j++){
      o = usedColumn[j];
      mean[l++] = V[ rc2ut(o,i,k) ];
      for(p=j; p < m; p++) {
        u = usedColumn[p];
        var[v++] = conditionalVar * V[ rc2ut(o,u,k) ]; 
      }
    }

    // generate a deviate 
    //ArMVN( sample, mean, var, j);
    for(j=0;j<m;j++) sample[j] = mean[j];
 
    // write results to estimates 
    for(j=0,m=0; j<i; j++) if( M[k*index[i] +j] ) estimates[(k+1)*index[i] + j] = sample[m++]; 
    estimates[(k+1)*index[i] +k] = -1.0 * conditionalVar; 

    free(usedColumn);
    free(sample);
    free(mean);
    free(var);
  }

  return;
}



/* r interface for ArMVN */
void RMVN2(
    double * sample,
    double * mean,
    double * var,
    int * sizePtr
    ) {

  GetRNGstate();
  ArMVN(sample, mean, var , *sizePtr); 
  PutRNGstate();

  return;
}



/* Gibbs sampling approach to generate a deviate from a MVN dist with mean=mean and var=var */
void ArMVN(                      
	  double *sample,         // sample array of length size
	  double *mean,           // mean array of length size 
	  double *var,            // the lower triangular array of the variance 
	  int size                // length of mean and dim of var 
 ) {
  int i,j,k,l;
  int n = ((size+1)*(size+2))/2;
  double * S = calloc( sizeof(double), n);
  double conditionalMean;
    
  /* draw from mult. normal using SWP */
  /* S:
   *   0   1   2   3
   * 0 -1    
   * 1 M0  V00 
   * 2 M1  V10 V11 
   * 3 M2  V20 V21 V22 
   *
   *    0
   *    1  6  
   *    2  7 11 
   *    3  8 12 15
   *    4  9 13 16 18
   *    5 10 14 17 19 20
   *
   *    move along row i:
   *    row i, col 0:  c =  i 
   *    row i, col 1:  c += size   
   *    row i, col 2:  c += size - 1
   *    row i, col 3:  c += size - 2 
   *    row i, col 4:  c += size - 3 
   *    row i, col 5:  c += size - 4 
   *    ...
   *    row i, col j:  c += size - j +1  
   *
   *    move along diag
   *    diag 0: 0
   *    diag 1: c+= size + 1
   *    diag 2: c+= size
   *    diag 3: c+= size -1 
   *    diag 4: c+= size -2
   *    diag 5: c+= size -3
   *
   *    diag j = c+= size -j +2
   */
  S[0] = -1;
  for(i=1;i<=size;i++) S[i]=mean[i-1];
  for(i=size+1,j=0;i<n;i++,j++) S[i]=var[j];

  // sample 0 = Z * sd00 + M0
  sample[0]= norm_rand() * sqrt(S[size+1]) + S[1];

  // init index for diagonal
  l = size+1;
  for(i=2;i<=size;i++){

    l += size - i + 2; 

    // perform sweep
    VSWP(S,i-1,size+1);
   
    // get the conditional mean 
    conditionalMean=S[i];
  
    // multiply the conditional mean times beta 
    k=i;
    for(j=1;j<i;j++){
      k += size - j + 1; 
      conditionalMean+=sample[j-1]*S[k];
    }
  
    sample[i-1]= norm_rand() * sqrt(S[l]) + conditionalMean;
  }
  
  free(S);
  return;
}




/* sweepTree */
void sweepTree( 
    covarTreePtr x, 
    double * V, 
    int k, 
    double ** matrixCache, 
    int * index,
    double * estimates,
    int * M,
    int n 
  ) {

  int i;

  // check for null trees
  if( x == NULL ) return; 
 
  // this is where we write it out 
  if( x->varList != NULL) {
    for( i = 0; i < x->varListLength; i++) {
      //printf("Variable = %d\n", x->varList[i] );
      saveParameterEstimates(V, k, x->varList[i], index, estimates,M,n);
    }
    return;
  }

  // to add, if both are not null, then cache
  /*    
          |
        1. cache, 3. use cache
       /     \
  2. cache   4. no cache  
     /  \      \
  */
  if( (x->yes != NULL) & (x->no != NULL) ) {
    matrixCache[x->cacheIndex] = calloc( sizeof(double), (k*(k+1))/2 );
    copyCovarMatrix(matrixCache[x->cacheIndex],V,k);
    //printCovarMatrix(matrixCache[x->cacheIndex],k);
  }

  // sweep to the left (yes)
  if( x->yes != NULL) {
    VSWP(V,x->index,k);
    //printCovarMatrix(V,k);
    // move to the yes
    sweepTree(x->yes,V,k,matrixCache,index,estimates,M,n);
  }

  // don't sweep to the right (no)
  if( x->no != NULL) {
    // get matrixCache if x->yes is not null
    if( x->yes != NULL) copyCovarMatrix(V,matrixCache[x->cacheIndex],k);

    // move to the no 
    sweepTree(x->no,V,k,matrixCache,index,estimates,M,n);

    // free the matrixCache 
    if( matrixCache[x->cacheIndex] != NULL ) {
      free( matrixCache[x->cacheIndex] );
      matrixCache[x->cacheIndex] = NULL;
    }
  }

  return;
}



/* R interface */
void RSweepTree( 
  double * x,          // upper (lower in R) triangular matrix including diag
  int *   M,           // m by p matrix of model parameter inclusions 
  int  * regIndex,     // variables (row indexes) that will be regressed
  double * est,        // m by p +2 matrix of parameter estimates and degrees of freedom
                       // index p is for var, index p+1 is for df
  int  *   index,      // identify index with row number e.g. (-1,-1,-1,0,1,2)
  int  *   pPtr,       // number of rows/cols in x
  int  *   mPtr,       // number of rows in M 
  int * n              // number of obs used to calc x
) {


  covarTreePtr myTree = NULL;
  int i,j;
  int cacheSize;
  double ** cache;

  int p = * pPtr;
  int m = * mPtr;
  cacheSize = 0;


  // create tree
  myTree = createCovarTree( NULL, M, p, regIndex[0], 0, &cacheSize); 
  for(i=1;i<m;i++) {
    myTree = createCovarTree(myTree, &(M[i*p]), p, regIndex[i], 0, &cacheSize); 
  }

  // allocate space for the cache
  cache = calloc( sizeof(double *) , cacheSize+1 );

  // estimate parameters for model through tree
  sweepTree(myTree, x, p, cache, index, est, M,*n);

  free(cache);
  deleteCovarTree(myTree);

  return;
}




/* test function */
#ifdef TEST_SWEEPTREE 
int main( void ) {

  covarTreePtr myTree = NULL;
  int * covarList = calloc(5, sizeof(int));
  int i;
  int cacheSize;
  double ** cache;
  double * est;
  int n = 5;

  // symmetric matrix
  double x[] = {  
    10.31288,  2.448485, -2.412443, -1.393328, -6.486046, 
               7.749054, -2.442433, 0.5318202,  1.395485,
                          11.15893,  1.965747,  6.863963,
                                     10.67714,  1.181446,
                                                12.50144 
  };

  //             0  1  2 [3  4]
  int index[] = {0, 0, 0, 0, 1};  
  int m = 2;
  //printf("\n Matrix = %p \n", (void *) x);
  //printCovarMatrix(x,n);
  cacheSize = 0;


  for( i=0; i < 3; i++) covarList[i] = i % 2 == 0 ? 1 : 0; 
  covarList[3] = 0;
  covarList[4] = 0;
  //for( i=0; i < 5; i++) printf("covarList[%d] = %d\n", i , (int) covarList[i]);

  myTree = createCovarTree( NULL, covarList, n, 3,  0, &cacheSize); 


  for( i=0; i < 4; i++) covarList[i] = 1; 
  covarList[4] = 0;

  //for( i=0; i < 5; i++) printf("covarList[%d] = %d\n", i , (int) covarList[i]);

  myTree = createCovarTree( myTree, covarList, n, 4 , 0,  &cacheSize); 


  //printf("Printing Tree\n");
  //printCovarTree(myTree);

  est = calloc( sizeof(double), m * (n+1) );

  //printf("Cache Size = %d\n",cacheSize);
  cache = calloc( sizeof(double *) , cacheSize+1 );
  sweepTree(myTree, x, 5, cache, index, est,NULL,0);


  //printFullMatrix( est, m, n+1);
  free(cache);
  free(est);

  //printf("Deleteing Tree\n");
  deleteCovarTree(myTree);

  free( covarList );

  return(0);
}
#endif



#ifdef TEST_SWEEPTREE2 
int main() {

  int n = 5;
  int i;

  double * v = calloc(sizeof(double), (n*(n+1))/2);
  double * x = calloc(sizeof(double), n*n);


  for( i=0; i < n*n; i++) x[i] = i+1; 


  copyMatrixToLowerTriangularArray(x,v,n);
  printCovarMatrix(v,n);  
  
  for( i=0; i < n*n; i++) x[i] = 0; 
  
  copyMatrixFromLowerTriangularArray(v,x,n);
 
  for( i=0; i < (n*(n+1))/2; i++) v[i] = 0; 
  
  copyMatrixToLowerTriangularArray(x,v,n);
  printCovarMatrix(v,n);  


  free(x);
  free(v);

  return 0;
}
#endif 
