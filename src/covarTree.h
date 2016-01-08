/* header file for intQueue */

#ifndef COVARTREE_HEADER

#define COVARTREE_HEADER

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

/* create a queue ADT */
struct covarTree { 
  int index;
  int * varList;
  int varListLength;
  int cacheIndex;
  struct covarTree * yes;
  struct covarTree * no;
};

typedef struct covarTree covarTree;
typedef struct covarTree * covarTreePtr;

/* function to create a new covarTree */
covarTreePtr createCovarTree( covarTreePtr, bool * , int, int, int, int * );

/* function to delete a covarTree */
void deleteCovarTree( covarTreePtr x ); 

/* function to delete a covarTree */
void printCovarTree( covarTreePtr x ); 

#endif




