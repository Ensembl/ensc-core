#define __MAIN_C__
#include "Set.h"
#undef __MAIN_C__

Set *Set_new() {
  Set *set;
  if ((set = (Set *)calloc(1,sizeof(Set))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for Set\n");
    return NULL;
  }

  return set;
}

void *Set_getElementAt(Set *s, int ind) {
  if (ind < 0 || ind >= s->nElement) {
    fprintf(stderr,"ERROR: Invalid element index %d\n",ind);
    exit(1);
  }
  return s->elements[ind];
}

void *Set_addElement(Set *s, void *elem) {
  if (elem == NULL) {
    fprintf(stderr, "ERROR: Element null in Set_addElement call\n");
    return NULL;
  }

  if (!s->nElement) s->elements = NULL;

  s->nElement++;
  if ((s->elements = (void **)realloc(s->elements,
               s->nElement*sizeof(void *))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for elem array\n");
    return NULL;
  }

  s->elements[s->nElement-1] = elem;

  return elem;
}

void Set_free(Set *s, int freeFunc()) {
  int i;

  if (s->isSpecial) {
    return;
  }

  if (freeFunc) {
    for (i=0;i<s->nElement;i++) {
      freeFunc(s->elements[i]);
    }
  }
  free(s->elements);
  free(s);
}
