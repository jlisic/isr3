#include "sweepTree.h"



/******************************************************************
 *   This function is a part of MNP: R Package for Estimating the
 *     Multinomial Probit Models by Kosuke Imai and David A. van Dyk.
 *       Copyright: GPL version 2 or later.
 ********************************************************************/
void SWP(
  double **X,             /* The Matrix to work on */
  int k,                  /* The row to sweep */
  int size)               /* The dim. of X */
{
  int i,j;
  
  if(X[k][k] < 10e-20) {
    printf("SWP: singular matrix.\n");
    return;
  }
  else X[k][k]=-1/X[k][k];

  for(i=0;i<size;i++)
    if(i!=k){
      X[i][k]=-X[i][k]*X[k][k];
      X[k][i]=X[i][k];
    }

  for(i=0;i<size;i++) 
  {
    for(j=0;j<size;j++)
    { 
      if(i!=k && j!=k) X[i][j]=X[i][j]+X[i][k]*X[k][j]/X[k][k];
    }
  }

  return;  
}




/* R interface */
void RSWP(
  double * x,
  int * kPtr,
  int * kSizePtr,
  int * sizePtr
) {
  
  size_t i;
  double ** X;
  
  X = calloc( *sizePtr, sizeof(double *) );
  for(i = 0; i < *sizePtr; i++) X[i] = &x[ *sizePtr * i ];
 
  for(i=0; i < *kSizePtr; i++) SWP(X, kPtr[i], *sizePtr);
  
  free(X);
  return;
}



/* Array interface */
double * ASWP(
  double * x,
  int * kPtr,
  int kSizePtr,
  int sizePtr
) {
  
  size_t i;
  double ** X;
  
  X = calloc( sizePtr, sizeof(double *) );
  for(i = 0; i < sizePtr; i++) X[i] = &x[ sizePtr * i ];
 
  for(i=0; i < kSizePtr; i++) SWP(X, kPtr[i], sizePtr);
  
  free(X);
  return(x);
}




/* print covar matrix function */
void printCovarMatrix ( double * x, int k ) {

  int i,j;

  printf("\n\t");
  for(i=0; i < k; i++) printf("%d\t", i);
  printf("\n");

  for(i=0; i < k; i++){
    printf("%d\t", i);
    for(j=0; j < k; j++) printf("%5.3f\t", x[i*k + j]);
    printf("\n");
  }

  return;
}




/* copy covar matrix function */
void copyCovarMatrix ( double * x, double * y, int k ) {
  int i;
  for(i=0; i < k*k; i++) x[i] = y[i];
  return;
}




/* save Parameters */
void saveParameterEstimates( V, k, i, estimates ) {
  for(int j = 0; j <= i; j++) estimates[k*i + j] = V[k*i + j]; // copy the row
  return;
}




//// need to take care of cache now
// need to create a function to get an index from cache index
// possibly create a function that creates the base tree, and create
// a separate function to add to the tree

/* function to print a covarTree */
void sweepTree( 
    covarTreePtr x, 
    double * V, 
    int k, 
    double ** matrixCache, 
    double * sigma, 
    double * beta
  ) {

  if( x == NULL ) return; 
 
  // this is where we write it out 
  if( x->varList != NULL) {
    for( i = 1; i < x->varListLength; i++) saveParameterEstimates(V, k, x->varList[i], estimates);
    return;
  }

  // to add, if both are not null, then cache
  if( (x->yes != NULL) & (x->no != NULL) ) {
    matrixCache[x->cacheIndex] = calloc( sizeof(double), k*k);
    copyCovarMatrix(matrixCache[x->cacheIndex],V,k);
  }

  if( x->yes != NULL) {
    V=ASWP(V,&(x->index),1,k);
    //printCovarMatrix(V,k);
    sweepTree(x->yes,V,k, matrixCache);
  }
  if( x->no != NULL) {
    // get matrixCache if x->yes is not null
    if( x->yes != NULL) copyCovarMatrix(V,matrixCache[x->cacheIndex],k);

    // sweep in matrixCache 
    sweepTree(x->no,V,k,matrixCache);

    // free the matrixCache 
    if( matrixCache[x->cacheIndex] != NULL ) {
      free( matrixCache[x->cacheIndex] );
      matrixCache[x->cacheIndex] == NULL;
    }
  }

  return;
}









/* test function */
int main( void ) {

  covarTreePtr myTree = NULL;
  bool * covarList = calloc(sizeof(bool), 5 );
  int i;
  int cacheSize;
  double ** cache;

  // symmetric matrix
  double x[] = {  
    10.31288,  2.448485, -2.412443, -1.393328, -6.486046, 
    2.448485,  7.749054, -2.442433, 0.5318202,  1.395485,
   -2.412443, -2.442433,  11.15893,  1.965747,  6.863963,
   -1.393328, 0.5318202,  1.965747,  10.67714,  1.181446,
   -6.486046,  1.395485,  6.863963,  1.181446,  12.50144 };

  //printCovarMatrix(x,5);
  cacheSize = 0;


  for( i=0; i < 5; i++) covarList[i] = i % 2 == 0 ? true : false; 
  for( i=0; i < 5; i++) printf("covarList[%d] = %d\n", i , (int) covarList[i]);

  myTree = createCovarTree( NULL, covarList, 5, 2, 0, &cacheSize); 


  for( i=0; i < 5; i++) covarList[i] = true; 
  for( i=0; i < 5; i++) printf("covarList[%d] = %d\n", i , (int) covarList[i]);

  myTree = createCovarTree( myTree, covarList, 5, 3, 0,  &cacheSize); 


  printf("Printing Tree\n");
  //printCovarTree(myTree);
  
  printf("Cache Size = %d\n",cacheSize);
  cache = calloc( sizeof(double *) , cacheSize+1 );
  sweepTree(myTree, x, 5, cache);
  for( i=0; i <= cacheSize; i++) {
    printf("%d: %p\n", i, (void *) cache[i] );
    if( cache[i] != NULL ) free( cache[i] );
  }
  free(cache);


  printf("Deleteing Tree\n");
  deleteCovarTree(myTree);

  free( covarList );

  return(0);
}



