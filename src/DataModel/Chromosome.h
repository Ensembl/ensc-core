/*
 * Copyright [1999-2024] EMBL-European Bioinformatics Institute
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

#ifndef __CHROMOSOME_H__
#define __CHROMOSOME_H__

#include "DataModelTypes.h"
#include "Storable.h"
#include "Object.h"

OBJECTFUNC_TYPES(Chromosome)

typedef struct ChromosomeFuncsStruct {
  OBJECTFUNCS_DATA(Chromosome)
} ChromosomeFuncs;

#define FUNCSTRUCTTYPE ChromosomeFuncs
struct ChromosomeStruct {
  OBJECT_DATA
  Storable st;
  int length;
  char *name;
};
#undef FUNCSTRUCTTYPE


Chromosome *Chromosome_new(void);

#define Chromosome_setDbID(c,dbID) Storable_setDbID(&((c)->st),(dbID))
#define Chromosome_getDbID(c) Storable_getDbID(&((c)->st))

#define Chromosome_setAdaptor(c,ad) Storable_setAdaptor(&((c)->st),(ad))
#define Chromosome_getAdaptor(c) Storable_getAdaptor(&((c)->st))

#define Chromosome_setLength(c,len) (c)->length = (len)
#define Chromosome_getLength(c) (c)->length

ECOSTRING Chromosome_setName(Chromosome *c,char *name);
#define Chromosome_getName(c) (c)->name

void Chromosome_free(Chromosome *chr);

#ifdef __CHROMOSOME_MAIN__
  ChromosomeFuncs
    chromosomeFuncs = {
                    Chromosome_free,
                    NULL, // shallowCopy
                    NULL  // deepCopy
                   };
#else
  extern ChromosomeFuncs chromosomeFuncs;
#endif


#endif
