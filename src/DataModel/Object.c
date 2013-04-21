#define __OBJECT_MAIN__
#include <stdio.h>
#include <stdlib.h>
#include "Object.h"
#include "Class.h"
#undef __OBJECT_MAIN__

void Object_freeImpl(Object *obj) {
   (obj)->funcs->free == NULL ? \
        (fprintf(stderr,"Error: Null pointer for free func - bye\n"),  exit(1), (void *)NULL) : \
        ((obj)->funcs->free((obj)), (void *)NULL);
}

void Object_errorUnimplementedMethod(Object *obj, char *methodName) {
  fprintf(stderr, "ERROR: Unimplemented method %s in object of type %s\n", methodName, Class_findByType(obj->objectType)->name);
  exit(1);
}
