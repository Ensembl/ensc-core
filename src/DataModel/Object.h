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


#define FUNCSTRUCTTYPE ObjectFuncs
struct ObjectStruct {
  OBJECT_DATA
};
#undef FUNCSTRUCTTYPE

void Object_freeImpl(Object *obj);

#define Object_incRefCount(obj) (obj)->referenceCount++
#define Object_decRefCount(obj) (obj)->referenceCount--

#define Object_getRefCount(obj) (obj)->referenceCount


#define Object_free(obj) \
      ((obj)->funcs->free == NULL ? \
         (fprintf(stderr,"Error: Null pointer for free func - bye\n"),  exit(1), (void *)NULL) : \
         ((obj)->funcs->free((obj))))

#endif
