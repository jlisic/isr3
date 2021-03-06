#include "rebuildCovar.h"






/* 
 * A function to rebuild the covariance matrix from
 * a linear model specification.
 *
 * Input:
 *
 * x: is an m by (p+1) (double) matrix of covariates and variances (col major)
 * y: is an m*(m+1)/2 double array
 *
 * Output:
 *
 * y: is an m by m covariance matrix in triangular + diag form  
 */
void rebuildCovar( double * x, double * y, int p) {

  int i,j,k;
  double tmp;
  double * beta = NULL;
  int * row = calloc(p, sizeof(int));

  /*
   * row: is an index array for assisting in indexing
   *
   *  0
   *  1  6
   *  2  7 11
   *  3  8 12 15
   *  4  9 13 16 18
   *  5 10 14 17 19 20
  */
    
  row[0] = 0;
  for( i = 1; i < p; i++) row[i] += row[i-1] + p - i + 1;

  // m < p
  for(i = 0; i < p; i++) {
    //printf("y:i= %d\n",i );  
    //printCovarMatrix(y,p);

    // beta vector
    beta = &(x[i*(p+1)]);
   
    // update variance   
    y[ row[i] ] = beta[p];


    /* example p = 3
     *
     * [0 .]
     * [1 3] 
     *  2 4 5
     *
     * [0 .] beta[1]
     * [1 3] beta[2] 
     *
     *
     */
    
    /* example p = 4 
     *
     *  0 . . .
     *  1 4 . . 
     *  2 5 7 .
     *  3 6 8 9
     *
     * [0 1 2] beta[1]
     * [. 4 5] beta[2] 
     * [. . 7] beta[3]
     *
     * for i = 3, j = 0 , iterate by k for each element
     * [0 1 2] beta[1]
     *         beta[2] 
     *         beta[3]
     *
     * for i = 3, j = 1 , iterate by k for each element
     *         beta[1]
     * [1 4 5] beta[2] 
     *         beta[3]
     *
     * for i = 3, j = 2 , iterate by k for each element
     *         beta[1]
     *         beta[2] 
     * [2 5 7] beta[3]
     */

    // iterate over row of y
    for(j=0;j<i;j++) {
      tmp=0;  // tmp variable

      
      // multiply beta times the covar matrix to update covar matrix for row i
      // unalligned memory multiply 
      //printf("*** i=%d j= %d***\n",i,j);
      for(k=0; k<j; k++) {
        //printf("%d: (%d) %5.4f * %5.4f (unalligned)\n",k, row[k]-k+j, y[row[k]-k+j], beta[k]);
        tmp += y[row[k]-k+j] * beta[k];
      }
      // alligned memory multiply 
      for(   ; k<i; k++) {
        //printf("%d: (%d) %5.4f * %5.4f (alligned)\n",k, row[j]-j+k, y[row[j]-j+k], beta[k]);
        tmp += y[row[j]-j+k] * beta[k]; 
      }
      
      //printf("%d %d tmp = %f\n",i, row[j] +i, tmp);
      // update covariance 
      y[ row[j] -j +i ] = tmp;
      
      // update marginal variance 
      //printf("%d %d tmp = %f + %f * %f \n",i, row[i] +i, y[row[i]], tmp, beta[j] );
      y[ row[i] ] += tmp * beta[j];
    }

  }

  beta = NULL;
  free(row);
  return;
}




/* R interface */
void RRebuildCovar( 
  double * x,          // upper (lower in R) triangular matrix including diag
  double * est,        // p by p matrix of parameter estimates
  int  *   nPtr        // number of rows/cols in x
) {

  int n = * nPtr;

//debug
//printFullMatrix( est, n, n+1);

  rebuildCovar( est, x, n);

//debug 
//printf("Rebuilt\n");
//printCovarMatrix(x,n);

  return;
}




/* test function */
#ifdef TEST_REBUILDCOVAR


int main( void ) {


  printf("\n\nRebuild Tree\n\n");

  covarTreePtr myTree = NULL;
  bool * covarList = calloc(5, sizeof(bool));
  int i;
  int cacheSize;
  double ** cache;
  double * est;
  double * y;
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
  int vars[] = {3,4};  
  int m = 2;
  
  
  //printf("\n Matrix = %p \n", (void *) x);
  //printCovarMatrix(x,n);
  cacheSize = 0;


  for( i=0; i < 3; i++) covarList[i] = i % 2 == 0 ? true : false; 
  covarList[3] = false;
  covarList[4] = false;
  //for( i=0; i < 5; i++) printf("covarList[%d] = %d\n", i , (int) covarList[i]);

  myTree = createCovarTree( NULL, covarList, n, 3, 0, &cacheSize); 


  for( i=0; i < 4; i++) covarList[i] = true; 
  covarList[4] = false;

  //for( i=0; i < 5; i++) printf("covarList[%d] = %d\n", i , (int) covarList[i]);

  myTree = createCovarTree( myTree, covarList, n, 4, 0,  &cacheSize); 


  //printf("Printing Tree\n");
  //printCovarTree(myTree);

  est = calloc( n * (n+1), sizeof(double) );

  //printf("Cache Size = %d\n",cacheSize);
  cache = calloc( cacheSize +1, sizeof(double *) );
  sweepTree(myTree, x, 5, cache, index, est);


  //printFullMatrix( est, m, n+1);

  y = calloc( m*(m+1)/2, sizeof(double) );
  rebuildCovar( est, y, vars, n, m);

  //printf("Rebuilt\n");
  //printCovarMatrix(y,m);


  free(cache);
  free(est);
  free(y);

  //printf("Deleteing Tree\n");
  deleteCovarTree(myTree);

  free( covarList );

  return(0);
}
#endif


