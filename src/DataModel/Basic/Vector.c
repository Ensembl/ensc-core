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

  return vector;
}

void *Vector_getElementAt(Vector *v, int ind) {
  if (ind < 0 || ind >= v->nElement) {
    fprintf(stderr,"ERROR: Invalid element index %d\n",ind);
    exit(1);
  }
  return v->elements[ind];
}

void *Vector_getLastElement(Vector *v) {
  if (!v->nElement) {
    fprintf(stderr,"ERROR: No elements in vector\n");
    exit(1);
  }
  return v->elements[v->nElement-1];
}

void *Vector_setElementAt(Vector *v, int ind, void *elem) {
  if (ind < 0 || ind >= v->nElement) {
    fprintf(stderr,"ERROR: Invalid element index %d\n",ind);
    exit(1);
  }
/* NIY free old one
*/
  
  v->elements[ind] = elem;

  return v->elements[ind];
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

  v->nElement++;
  if ((v->elements = (void **)realloc(v->elements,
               v->nElement*sizeof(void *))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for elem array\n");
    return NULL;
  }

  v->elements[v->nElement-1] = elem;

  return elem;
}

void Vector_free(Vector *v, int freeFunc()) {
  int i;

  if (v->isSpecial) {
    return;
  }

  if (freeFunc) {
    for (i=0;i<v->nElement;i++) {
      if (v->elements[i]) {
        freeFunc(v->elements[i]);
      }
    }
  }
  free(v->elements);
  free(v);
}
