/*
 * Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
 * Copyright [2016-2017] EMBL-European Bioinformatics Institute
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *      http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef __OBJECT_H__
#define __OBJECT_H__

#include "Class.h"

typedef struct ObjectStruct Object;
#define OBJECTFUNC_TYPES(CLASSTYPE) \
typedef void (*CLASSTYPE ## _FreeFunc)(CLASSTYPE *); \
typedef CLASSTYPE *(*CLASSTYPE ## _ShallowCopyFunc)(CLASSTYPE *); \
typedef CLASSTYPE *(*CLASSTYPE ## _DeepCopyFunc)(CLASSTYPE *);

#define OBJECTFUNCS_DATA(CLASSTYPE) \
  CLASSTYPE ## _FreeFunc free; \
  CLASSTYPE ## _ShallowCopyFunc shallowCopy; \
  CLASSTYPE ## _DeepCopyFunc deepCopy;

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
void Object_errorUnimplementedMethod(Object *obj, char *methodName);

#define Object_incRefCount(obj) (obj)->referenceCount++
//void Object_incRefCount(Object *obj);
#define Object_decRefCount(obj) (obj)->referenceCount--
//void Object_decRefCount(Object *obj);

#define Object_getRefCount(obj) (obj)->referenceCount

// Comment out to reduce warnings void Object_errorUnimplementedMethod(Object *obj, char *methodName);

#define Object_free(obj) \
      ((obj)->funcs->free == NULL ? \
         (fprintf(stderr,"Error: Null pointer for free func - bye\n"),  exit(1), (void *)NULL) : \
         ((obj)->funcs->free((obj)), (void *)NULL))

#define Object_shallowCopy(obj) \
      ((obj)->funcs->shallowCopy == NULL ? \
         (fprintf(stderr,"Error: Null pointer for shallowCopy func - bye\n"),  exit(1), (void *)NULL) : \
         ((obj)->funcs->shallowCopy((obj))))

#endif
