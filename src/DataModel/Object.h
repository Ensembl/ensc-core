#ifndef __OBJECT_H__
#define __OBJECT_H__

#include "Class.h"

typedef struct ObjectStruct Object;

#define OBJECT_DATA \
  ClassType objectType; \
  int referenceCount;

struct ObjectStruct {
  OBJECT_DATA
};

#define Object_incRefCount(obj) (obj)->referenceCount++
#define Object_decRefCount(obj) (obj)->referenceCount--

#define Object_getRefCount(obj) (obj)->referenceCount

#endif
