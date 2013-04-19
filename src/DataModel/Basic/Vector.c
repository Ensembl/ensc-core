#define __VECTOR_MAIN__
#include "Vector.h"
#undef __VECTOR_MAIN__
#include "EnsC.h"

Vector *Vector_new() {
  Vector *vector;

  if ((vector = (Vector *)calloc(1,sizeof(Vector))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for Vector\n");
    return NULL;
  }

  vector->objectType = CLASS_VECTOR;

  Object_incRefCount(vector);

  vector->funcs = &vectorFuncs;

  return vector;
}

Vector *Vector_newFromArray(void **array, int nInArray) {
  Vector *vector = Vector_new();

  int i;
  for (i=0;i<nInArray;i++) {
    Vector_addElement(vector, array[i]);
  }

  return vector;
}

void Vector_setFreeFunc(Vector *v, void freeElement()) {
  v->freeElement = freeElement;
}

void *Vector_getElementAt(Vector *v, int ind) {
  if (ind < 0 || ind >= v->nElement) {
    fprintf(stderr,"ERROR: Invalid element index %d in Vector_getElementAt\n",ind);
    exit(1);
  }
  return v->elements[ind];
}

void *Vector_getLastElement(Vector *v) {
  if (!v->nElement) {
    fprintf(stderr,"ERROR: No elements in vector in Vector_getLastElement\n");
    exit(1);
  }
  return v->elements[v->nElement-1];
}

void *Vector_setElementAt(Vector *v, int ind, void *elem) {
  if (ind < 0) {
    fprintf(stderr,"ERROR: Invalid element index %d in Vector_setElementAt\n",ind);
    exit(1);
  } else if (ind >= v->nElement) {
    Vector_setNumElement(v, ind+1);
  }
/* NIY free old one
*/
  
  v->elements[ind] = elem;

  return v->elements[ind];
}

void *Vector_removeElementAt(Vector *v, int ind) {
  void *removed;
  int i;

  if (ind < 0) {
    fprintf(stderr,"ERROR: Invalid element index %d in Vector_removeElementAt\n",ind);
    exit(1);
  } else if (ind >= v->nElement) {
    fprintf(stderr,"ERROR: Invalid element index %d in Vector_removeElementAt\n",ind);
    exit(1);
  }
/* NIY free old one
*/
  
  removed = v->elements[ind];
  
  for (i=ind+1; i<v->nElement; i++) {
    v->elements[i-1] = v->elements[i];
  }

  //Vector_setNumElement(v, v->nElement-1);
  v->nElement--;

  return removed;
}

void *Vector_insertElementAt(Vector *v, int ind, void *elem) {
  void *removed;
  int i;

  if (ind < 0) {
    fprintf(stderr,"ERROR: Invalid element index %d in Vector_insertElementAt\n",ind);
    exit(1);
  } else if (ind > v->nElement) {
    fprintf(stderr,"ERROR: Invalid element index %d Vector_insertElementAt\n",ind);
    exit(1);
  }
  
  Vector_setNumElement(v, v->nElement+1);

  for (i=v->nElement-2; i>=ind; i--) {
    v->elements[i+1] = v->elements[i];
  }

  return Vector_setElementAt(v, ind, elem);
}

void Vector_reverse(Vector *v) {
  int up =0;
  int down = v->nElement-1;
  void *tmp;

  for (; up<down; up++, down--) {
    tmp = v->elements[up];
    v->elements[up] = v->elements[down];
    v->elements[down] = tmp;
  }
}

void Vector_append(Vector *dest, Vector *src) {
  int i;

  for (i=0; i<Vector_getNumElement(src); i++) {
    Vector_addElement(dest, Vector_getElementAt(src,i));
  }

  return;
}

void Vector_sort(Vector *v, SortCompFunc sortFunc) {
  qsort(v->elements,v->nElement,sizeof(void *),sortFunc);
}

void *Vector_addElement(Vector *v, void *elem) {
  if (elem == NULL) {
    fprintf(stderr, "WARNING: Element null in Vector_addElement call\n");
  }

  if (!v->nElement) v->elements = NULL;

  Vector_setNumElement(v,v->nElement+1);

  v->elements[v->nElement-1] = elem;

  return elem;
}

void Vector_setNumElement(Vector *v, int nElem) {
  int i;

  if ((v->elements = (void **)realloc(v->elements,nElem*sizeof(void *))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for elem array\n");
    return;
  }

  for (i=v->nElement; i<nElem; i++) {
    v->elements[i] = NULL;
  }

  v->nElement = nElem;

  return;
}

Vector *Vector_copy(Vector *v) {
  int i;
  Vector *newV = Vector_new();


  Vector_setNumElement(newV, v->nElement);

  memcpy(newV->elements,v->elements,v->nElement*sizeof(void *));

  return newV;
}

void Vector_free(Vector *v) {
  int i;

  //printf("Vector free called\n");
  if (v->isSpecial) {
    printf(" - special vector so returning\n");
    return;
  }

  Object_decRefCount(v);

  if (Object_getRefCount(v) > 0) {
    return;
  } else if (Object_getRefCount(v) < 0) {
    fprintf(stderr,"Error: Negative reference count for Vector\n"
                   "       Freeing it anyway\n");
  }


  if (v->freeElement) {
    for (i=0;i<v->nElement;i++) {
      if (v->elements[i]) {
        v->freeElement(v->elements[i]);
      }
    }
  }
  free(v->elements);
  free(v);
  //printf(" - done freeing vector\n");
}
