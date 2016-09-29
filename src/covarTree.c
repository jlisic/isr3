#include "covarTree.h"



/* function to create a new covarTree */
/*
 * The covar tree goes to the left when a covar is included, and to the
 * right when it is not.  Consider the case of regressing X3, X4, and X5 
 * on X = {X1|...|X5}, where X3 and X5 are near-singular.
 *
 *    X1 X2 X3 X4 X5
 * X3  1  1  0  0 0
 * X4  1  1  1  0 0
 * X5  1  1  0  1 0
 *
 * Tree Structure:
 *                          X1=1
 *            X2=1
 *      X3=1         X3=0                  
 *      X4=0      X4=1  X4=0  
 *      X5=0      X5=0  X5=0   
 *
 *   Cache points are created when both children are populated.
 *   In this example there is a cache point at X2=1 and X3=0|X2=1.
 *   When we create the covar tree, we don't set these points yet,
 *   instead we create a mapping to an array of matrix pointers
 *   that will be populated at a later time.
 */
covarTreePtr createCovarTree( 
    covarTreePtr x,      // root node of tree
    int * covarList,    // array of inclusion and excusion of covariates
    int covarListLength, // number of covariates 
    int varIndex,        // varIndex to store
    int covarIndex,      // covar starting point usually 0
    int * cacheIndex     // index of cache
    ) {
  
  if( covarListLength > 0 ) {  // not the last covar

    // allocate and initialize 
    if( x == NULL ) {
      //printf("Allocating %d\n", covarListLength);
      x = malloc( sizeof( covarTree ) );
      x->index = covarIndex;        // covariate index 
      x->cacheIndex = ++(*cacheIndex);              // location in cache array 
      x->varList = NULL;              // list of variables  
      x->varListLength = 0;           // length of the list of variables
      x->yes = NULL;                  // yes, the covariate index is swept
      x->no = NULL;                   // no, the covariate index is not swept
    }

    // if the covarlist is true, set x->yes
    if( covarList[0] ) {
      //if( x->no != NULL) x->cacheIndex = ++(*cacheIndex); 
      x->yes = createCovarTree( x->yes, covarList+1, covarListLength-1, varIndex, covarIndex+1,cacheIndex );
    // otherwise set x->no
    } else {
      //if( x->yes != NULL) x->cacheIndex = ++(*cacheIndex); 
      x->no = createCovarTree( x->no, covarList+1, covarListLength-1, varIndex, covarIndex+1, cacheIndex );
    }

  // if covar list length = 0
  } else { 
  // add on variables 
    
    // allocate and initialize 
    if( x == NULL ) {
      //printf("Allocating %d\n", covarListLength);
      x = malloc( sizeof( covarTree ) );
      x->index = -1;                  // covariate index 
      x->cacheIndex = 0;              // location in cache array 
      x->varList = NULL;              // list of variables  
      x->varListLength = 0;           // length of the list of variables
      x->yes = NULL;                  // yes, the covariate index is swept
      x->no = NULL;                   // no, the covariate index is not swept
    }
    //printf("At the End %d\n", varIndex);

    // this is a bit messy, but works for low duplicate sets of covariates
    int * varListNew = malloc( sizeof(int) * (x->varListLength + 1) );      // create new array of size + 1
    int i;
    for( i=0; i < x->varListLength; i++) varListNew[i] = x->varList[i]; // copy the old data into the new array
    free(x->varList);                                                       // delete old array
    x->varList = varListNew;                                             // assign new array to node
    x->varList[x->varListLength] = varIndex;                                   // add the new varIndex
    x->varListLength++;
  }

  return(x);
} 




/* function to delete a covarTree */
void deleteCovarTree( covarTreePtr x ) {
  
  if( x == NULL ) return; 

  deleteCovarTree(x->yes);
  deleteCovarTree(x->no);

  //printf("Deleting %d\n", x->index);

  if( x->varList != NULL) free(x->varList);
  free(x);
  x = NULL;

  return;
}




/* function to print a covarTree */
void printCovarTree( covarTreePtr x ) {

  if( x == NULL ) return; 
  
  if( x->varList != NULL) {
    int i;
    printf("(%d", x->varList[0]);
    for( i = 1; i < x->varListLength; i++) printf(", %d", x->varList[i] ); 
    printf(")\n");
    return;
  }

  if( x->yes != NULL) {
    printf("node=%d,%d Yes: ", (int) x->index, (int) x->cacheIndex );
    printCovarTree(x->yes);
  }
  if( x->no != NULL) {
    printf("node=%d,%d No : ", (int) x->index, (int) x->cacheIndex );
    printCovarTree(x->no);
  }

  return;
}




/* need to update */
/* test function */
/*
int main( void ) {

  covarTreePtr myTree;

  bool * covarList = calloc(sizeof(bool), 5 );
  int i;

  for( i=0; i < 5; i++) covarList[i] = i % 2 == 0 ? true : false; 
  for( i=0; i < 5; i++) printf("covarList[%d] = %d\n", i , (int) covarList[i]);

  myTree = createCovarTree( NULL, covarList, 5, 2); 
  
  for( i=0; i < 5; i++) covarList[i] = true; 
  for( i=0; i < 5; i++) printf("covarList[%d] = %d\n", i , (int) covarList[i]);

  myTree = createCovarTree( myTree, covarList, 5, 3); 
  myTree = createCovarTree( myTree, covarList, 5, 4); 


  printf("Printing Tree\n");
  printCovarTree(myTree);

  printf("Deleteing Tree\n");
  deleteCovarTree(myTree);

  free( covarList );

  return(0);
}
*/
