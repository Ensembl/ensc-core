#ifndef __SET_H__
#define __SET_H__

#include <stdio.h>
#include <stdlib.h>

typedef struct SetStruct {
  void **elements;
  int nElement;
  int isSpecial;
} Set;

#ifdef __MAIN_C__
Set emptySetData = {NULL,0,1};
Set *emptySet = &emptySetData;

#else

extern Set *emptySet;

#endif

Set *Set_new();
void *Set_addElement(Set *set, void *elem);
#define Set_getNumElement(s) (s)->nElement
void *Set_getElementAt(Set *s, int ind);
void Set_free(Set *set, int freeFunc());

#endif
