#ifndef __OBJECT_H__
#define __OBJECT_H__

#include "Class.h"

typedef struct ObjectStruct Object;
#define OBJECTFUNC_TYPES(CLASSTYPE) \
typedef void (*CLASSTYPE ## _FreeFunc)(CLASSTYPE *);

#define OBJECTFUNCS_DATA(CLASSTYPE) \
  CLASSTYPE ## _FreeFunc free;

OBJECTFUNC_TYPES(Object)

typedef struct ObjectFuncsStruct {
  OBJECTFUNCS_DATA(Object)
} ObjectFuncs;




/* FUNCSTRUCTTYPE must be #defined where OBJECT_DATA is used */
#define OBJECT_DATA \
  ClassType objectType; \
  int referenceCount; \
  FUNCSTRUCTTYPE *funcs;


#define FUNCSTRUCTTYPE Object
struct ObjectStruct {
  OBJECT_DATA
};
#undef FUNCSTRUCTTYPE

#define Object_incRefCount(obj) (obj)->referenceCount++
#define Object_decRefCount(obj) (obj)->referenceCount--

#define Object_getRefCount(obj) (obj)->referenceCount

#endif
