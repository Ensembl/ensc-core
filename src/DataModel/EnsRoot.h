#ifndef __ENSROOT_H__
#define __ENSROOT_H__

#include "Object.h"

/* FUNCSTRUCTTYPE must be #defined where ENSROOT_DATA is used */

#define ENSROOT_DATA \
  OBJECT_DATA \
  int referenceCount; \
  FUNCSTRUCTTYPE *funcs;

/*
struct EnsRootStruct {
  ENSROOT_DATA
};
*/

#endif
