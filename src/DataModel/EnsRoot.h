#ifndef __ENSROOT_H__
#define __ENSROOT_H__

#include "Object.h"

/* FUNCSTRUCTTYPE must be #defined where ENSROOT_DATA is used */
typedef struct NoTypeFuncsStruct {
  int i; // Get compiler complaints if its empty
} NoTypeFuncs;

#define ENSROOT_DATA \
  OBJECT_DATA \
  int referenceCount; \
  FUNCSTRUCTTYPE *funcs;


#endif
