#define __MAIN_C__
#include "Vector.h"
#undef __MAIN_C__

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

void *Vector_addElement(Vector *v, void *elem) {
  if (elem == NULL) {
    fprintf(stderr, "ERROR: Element null in Vector_addElement call\n");
    return NULL;
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
      freeFunc(v->elements[i]);
    }
  }
  free(v->elements);
  free(v);
}
