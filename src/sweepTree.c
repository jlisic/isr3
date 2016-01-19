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




/* copy covar matrix function */
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
 
  N = n*(n+1)/2;

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




/* save Parameters */
void saveParameterEstimates( 
    double * V, 
    int k, 
    int i, 
    int * index, 
    double * estimates 
    ) {
  int j;
    
  int * row = calloc(sizeof( int * ), k);

  // reduce row calculations
  row[0] = 0;
  for( j = 1; j < k; j++) row[j] += row[j-1] + k - j + 1;
  
  // horizontal across
  for(j=0; j<i; j++) estimates[(k+1)*index[i] + j] = V[row[j] - j +i]; 
  
  // save sigma
  estimates[(k+1)*index[i] +k] = V[row[i]]; 

  free(row);

  return;
}




/* function to print a covarTree */
void sweepTree( 
    covarTreePtr x, 
    double * V, 
    int k, 
    double ** matrixCache, 
    int * index,
    double * estimates 
  ) {

  int i;

  if( x == NULL ) return; 
 
  // this is where we write it out 
  if( x->varList != NULL) {
    for( i = 0; i < x->varListLength; i++) {
      //printf("Variable = %d\n", x->varList[i] );
      saveParameterEstimates(V, k, x->varList[i], index, estimates);
    }
    return;
  }

  // to add, if both are not null, then cache
  if( (x->yes != NULL) & (x->no != NULL) ) {
    matrixCache[x->cacheIndex] = calloc( sizeof(double), (k*(k+1))/2 );
    copyCovarMatrix(matrixCache[x->cacheIndex],V,k);
    //printCovarMatrix(matrixCache[x->cacheIndex],k);
  }

  if( x->yes != NULL) {
    VSWP(V,x->index,k);
    //printCovarMatrix(V,k);
    sweepTree(x->yes,V,k,matrixCache,index,estimates);
  }
  if( x->no != NULL) {
    // get matrixCache if x->yes is not null
    if( x->yes != NULL) copyCovarMatrix(V,matrixCache[x->cacheIndex],k);

    // sweep in matrixCache 
    sweepTree(x->no,V,k,matrixCache,index,estimates);

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
  int *   M,          // m by p matrix of model parameter inclusions 
  int  * regIndex,   // variables (row indexes) that will be regressed
  double * est,        // p by p matrix of parameter estimates
  int  *   index,      // identify index with row number e.g. (-1,-1,-1,0,1,2)
  int  *   pPtr,       // number of rows/cols in x
  int  *   mPtr        // number of rows in M 
) {


  covarTreePtr myTree = NULL;
  int i,j;
  int cacheSize;
  double ** cache;
  bool * MBool;        //array to convert int from R to boolean

  int p = * pPtr;
  int m = * mPtr;
  cacheSize = 0;

//debug
//printCovarMatrix(x,p);

  // fix an annoying R interface issue, e.g. .C does not support logical to boolean
  MBool = calloc( sizeof( bool *), p);
  for(j=0;j<p;j++) MBool[j] = (bool) M[j];

  // create tree
  myTree = createCovarTree( NULL, MBool, p, regIndex[0], 0, &cacheSize); 
  for(i=1;i<m;i++) {
    for(j=0;j<p;j++) MBool[j] = (bool) M[i*p + j];
    
    myTree = createCovarTree(myTree, MBool, p, regIndex[i], 0, &cacheSize); 
  }

  // allocate space for the cache
  cache = calloc( sizeof(double *) , cacheSize+1 );

  // estimate parameters for model through tree
  sweepTree(myTree, x, p, cache, index, est);

  free(MBool);
  free(cache);
  deleteCovarTree(myTree);

  return;
}




/* test function */
#ifdef TEST_SWEEPTREE 
int main( void ) {

  covarTreePtr myTree = NULL;
  bool * covarList = calloc(sizeof(bool), 5 );
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


  for( i=0; i < 3; i++) covarList[i] = i % 2 == 0 ? true : false; 
  covarList[3] = false;
  covarList[4] = false;
  //for( i=0; i < 5; i++) printf("covarList[%d] = %d\n", i , (int) covarList[i]);

  myTree = createCovarTree( NULL, covarList, n, 3,  0, &cacheSize); 


  for( i=0; i < 4; i++) covarList[i] = true; 
  covarList[4] = false;

  //for( i=0; i < 5; i++) printf("covarList[%d] = %d\n", i , (int) covarList[i]);

  myTree = createCovarTree( myTree, covarList, n, 4 , 0,  &cacheSize); 


  //printf("Printing Tree\n");
  //printCovarTree(myTree);

  est = calloc( sizeof(double), m * (n+1) );

  //printf("Cache Size = %d\n",cacheSize);
  cache = calloc( sizeof(double *) , cacheSize+1 );
  sweepTree(myTree, x, 5, cache, index, est);


  //printFullMatrix( est, m, n+1);
  free(cache);
  free(est);

  //printf("Deleteing Tree\n");
  deleteCovarTree(myTree);

  free( covarList );

  return(0);
}
#endif


