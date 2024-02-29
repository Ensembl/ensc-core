/*
 * Copyright [1999-2024] EMBL-European Bioinformatics Institute
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *      http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef __VECTOR_H__
#define __VECTOR_H__

#include "EnsC.h"
#include "Class.h"
#include "Object.h"
#include <stdio.h>
#include <stdlib.h>

typedef void   (*Vector_ElementFreeFunc)();

typedef struct VectorStruct Vector;

OBJECTFUNC_TYPES(Vector)

typedef struct VectorFuncsStruct {
  OBJECTFUNCS_DATA(Vector)
} VectorFuncs;

#define FUNCSTRUCTTYPE VectorFuncs
struct VectorStruct {
  OBJECT_DATA
  void **elements;
  int nElement;
  int nAlloced;
  int batchSize;
  int isSpecial;
  Vector_ElementFreeFunc freeElement;  
};
#undef FUNCSTRUCTTYPE

#define VECTOR_DEFAULTBATCHSIZE 10

Vector *Vector_new();
void *Vector_addElement(Vector *vector, void *elem);
#define Vector_getNumElement(v) (v)->nElement
//void *Vector_getElementAt(Vector *v, int ind);
#define Vector_getElementAt(v, ind) (v)->elements[(ind)]

void Vector_free(Vector *vector);
void Vector_setFreeFunc(Vector *vector, void freeFunc());
void Vector_setBatchSize(Vector *vector, int batchSize);
Vector *Vector_newFromArray(void **array, int nInArray);

void Vector_append(Vector *dest, Vector *src);
void Vector_reverse(Vector *v);
void *Vector_setElementAt(Vector *v, int ind, void *elem);
void Vector_sort(Vector *v, SortCompFunc sortFunc);
void *Vector_getLastElement(Vector *v);
void Vector_setNumElement(Vector *v, int nElem);
void *Vector_removeElementAt(Vector *v, int ind);
void Vector_removeAll(Vector *v);
Vector *Vector_copy(Vector *v);
void **Vector_toArray(Vector *v);
void *Vector_insertElementAt(Vector *v, int ind, void *elem);


#ifdef __VECTOR_MAIN__
  VectorFuncs
    vectorFuncs = {
                    Vector_free,
                    NULL, // shallowCopy
                    NULL  // deepCopy
                   };
#else
  extern VectorFuncs vectorFuncs;
#endif

#ifdef __VECTOR_MAIN__
  Vector  emptyVectorData = {CLASS_VECTOR,-1, &vectorFuncs, NULL,0,0,0,1,NULL};
  Vector *emptyVector = &emptyVectorData;
  void   *singleEntryVectorArray[1] = { NULL };
  Vector  singleEntryVectorData = {CLASS_VECTOR,-1, &vectorFuncs, singleEntryVectorArray,1,1,0,1,NULL};
  Vector *singleEntryVector = &singleEntryVectorData;
#else
  extern Vector *emptyVector;
  extern Vector *singleEntryVector;
#endif




#endif
