/*
 * Copyright [1999-2017] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

#define __CLONE_MAIN__
#include "Clone.h"
#undef __CLONE_MAIN__
#include "CloneAdaptor.h"

Clone *Clone_new() {
  Clone *cl;

  if ((cl = (Clone *)calloc(1,sizeof(Clone))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for cl\n");
    return NULL;
  }

  cl->objectType = CLASS_CLONE;

  cl->funcs = &cloneFuncs;

  Object_incRefCount(cl);

  return cl;
}

char *Clone_setName(Clone *cl, char *name) {
  if ((cl->name = (char *)malloc(strlen(name)+1)) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for clone name\n");
    return NULL;
  }

  strcpy(cl->name,name);

  return cl->name;
}

char *Clone_setEmblAcc(Clone *cl, char *acc) {
  if ((cl->emblAcc = (char *)malloc(strlen(acc)+1)) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for clone acc\n");
    return NULL;
  }

  strcpy(cl->emblAcc,acc);

  return cl->emblAcc;
}

void Clone_free(Clone *clone) {
  Object_decRefCount(clone);

  if (Object_getRefCount(clone) > 0) {
    return;
  } else if (Object_getRefCount(clone) < 0) {
    fprintf(stderr,"Error: Negative reference count for Clone\n"
                   "       Freeing it anyway\n");
  }

  if (clone->name) EcoString_freeStr(ecoSTable, clone->name);
  if (clone->emblAcc) EcoString_freeStr(ecoSTable, clone->emblAcc);

  free(clone);
}


