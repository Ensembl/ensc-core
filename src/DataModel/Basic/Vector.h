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

#ifdef __VECTOR_MAIN__
  Vector  emptyVectorData = {NULL,0,1};
  Vector *emptyVector = &emptyVectorData;
  void   *singleEntryVectorArray[1] = { NULL };
  Vector  singleEntryVectorData = {singleEntryVectorArray,1,1};
  Vector *singleEntryVector = &singleEntryVectorData;
#else
  extern Vector *emptyVector;
  extern Vector *singleEntryVector;
#endif

Vector *Vector_new();
void *Vector_addElement(Vector *vector, void *elem);
#define Vector_getNumElement(v) (v)->nElement
void *Vector_getElementAt(Vector *v, int ind);
void Vector_free(Vector *vector, int freeFunc());

void Vector_append(Vector *dest, Vector *src);
void *Vector_setElementAt(Vector *v, int ind, void *elem);
void Vector_sort(Vector *v, SortCompFunc sortFunc);
void *Vector_getLastElement(Vector *v);

#endif
