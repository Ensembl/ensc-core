#ifndef __VECTOR_H__
#define __VECTOR_H__

#include "EnsC.h"
#include <stdio.h>
#include <stdlib.h>

typedef struct VectorStruct {
  void **elements;
  int nElement;
  int isSpecial;
} Vector;

#ifdef __MAIN_C__
  Vector emptyVectorData = {NULL,0,1};
  Vector *emptyVector = &emptyVectorData;
#else
  extern Vector *emptyVector;
#endif

Vector *Vector_new();
void *Vector_addElement(Vector *vector, void *elem);
#define Vector_getNumElement(v) (v)->nElement
void *Vector_getElementAt(Vector *v, int ind);
void Vector_free(Vector *vector, int freeFunc());

void Vector_append(Vector *dest, Vector *src);
void *Vector_setElementAt(Vector *v, int ind, void *elem);
void Vector_sort(Vector *v, SortCompFunc sortFunc);




#endif
