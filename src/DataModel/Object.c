#define __OBJECT_MAIN__
#include <stdio.h>
#include "Object.h"
#undef __OBJECT_MAIN__

void Object_freeImpl(Object *obj) {
   (obj)->funcs->free == NULL ? \
        (fprintf(stderr,"Error: Null pointer for free func - bye\n"),  exit(1), (void *)NULL) : \
        ((obj)->funcs->free((obj)));
}
