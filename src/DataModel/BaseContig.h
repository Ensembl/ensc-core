/*
 * Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

#ifndef __BASECONTIG_H__
#define __BASECONTIG_H__

#include "DataModelTypes.h"
#include "Storable.h"
#include "EnsRoot.h"
#include "Sequence.h"

#define BASECONTIGFUNC_TYPES(CLASSTYPE) \
  SEQUENCEFUNC_TYPES(CLASSTYPE)

#define BASECONTIGFUNCS_DATA(CLASSTYPE) \
  SEQUENCEFUNCS_DATA(CLASSTYPE)
  

  
BASECONTIGFUNC_TYPES(BaseContig)

typedef struct BaseContigFuncsStruct {
  BASECONTIGFUNCS_DATA(BaseContig)
} BaseContigFuncs;

#define BASECONTIG_DATA \
  SEQUENCE_DATA \
  Storable st; \
  long start; \
  long end;

#undef FUNCSTRUCTTYPE
#define FUNCSTRUCTTYPE BaseContigFuncs
struct BaseContigStruct {
  BASECONTIG_DATA
};
#undef FUNCSTRUCTTYPE

#define BaseContig_setDbID(bc,dbID) Storable_setDbID(&((bc)->st),dbID)
#define BaseContig_getDbID(bc) Storable_getDbID(&((bc)->st))

#define BaseContig_setAdaptor(bc,ad) Storable_setAdaptor(&((bc)->st),ad)
#define BaseContig_getAdaptor(bc) Storable_getAdaptor(&((bc)->st))

#define BaseContig_getObjectType(bc) (bc)->objectType

#define BaseContig_getSubSeq(seq,start,end,strand) \
      ((seq)->funcs->getSubSeq == NULL ? \
         (fprintf(stderr,"Error: Null pointer for getSubSeq - bye\n"),  exit(1), (char *)NULL) : \
         ((seq)->funcs->getSubSeq((seq),(start),(end),(strand))))

#define BaseContig_getSeq(seq) ((seq)->funcs->getSeq == NULL ? ((seq)->seq) : \
                                                               ((seq)->funcs->getSeq((seq))))

#define BaseContig_getName(seq) ((seq)->funcs->getName == NULL ? ((seq)->name) : \
                                                                 ((seq)->funcs->getName((seq))))


#endif
