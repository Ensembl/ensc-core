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

#ifndef __RAWCONTIG_H__
#define __RAWCONTIG_H__

#include "DataModelTypes.h"
#include "Storable.h"
#include "BaseContig.h"
#include "Vector.h"
#include "SequenceAdaptor.h"

BASECONTIGFUNC_TYPES(RawContig)

typedef struct RawContigFuncsStruct {
  BASECONTIGFUNCS_DATA(RawContig)
} RawContigFuncs;

#define FUNCSTRUCTTYPE RawContigFuncs
struct RawContigStruct {
  BASECONTIG_DATA
  int emblOffset;
  IDType cloneId;
};
#undef FUNCSTRUCTTYPE

RawContig *RawContig_new(void);

#define RawContig_setDbID(rc,id) BaseContig_setDbID((rc),(id))
#define RawContig_getDbID(rc) BaseContig_getDbID((rc))

#define RawContig_setAdaptor(rc,ad) BaseContig_setAdaptor((rc),(ad))
#define RawContig_getAdaptor(rc) BaseContig_getAdaptor((rc))

ECOSTRING RawContig_setName(RawContig *rc, char *name);
ECOSTRING RawContig_getName(RawContig *rc);

#define RawContig_setCloneID(rc,cid) (rc)->cloneId = (cid)
long RawContig_getCloneID(RawContig *rc);

#define RawContig_setEMBLOffset(rc,eo) (rc)->emblOffset = (eo)
int RawContig_getEMBLOffset(RawContig *rc);

#define RawContig_setLength(rc,l) (rc)->length = (l)
int RawContig_getLength(RawContig *rc);

Vector *RawContig_getAllSimpleFeatures(RawContig *rc, char *logicName, double *scoreP);
Vector *RawContig_getAllPredictionTranscripts(RawContig *rc, char *logicName);
Vector *RawContig_getAllRepeatFeatures(RawContig *rc, char *logicName);
Vector *RawContig_getAllDNAAlignFeatures(RawContig *rc, char *logicName, double *scoreP);
Vector *RawContig_getAllProteinAlignFeatures(RawContig *rc, char *logicName, double *scoreP);
char *RawContig_getSubSeq(RawContig *contig, int start, int end, int strand);
char *RawContig_getSeq(RawContig *contig);



#ifdef __RAWCONTIG_MAIN__
  RawContigFuncs rawContigFuncs = {
                           NULL, // free
                           NULL, // shallowCopy
                           NULL, // deepCopy
                           RawContig_getName,
                           RawContig_getSeq,
                           RawContig_getSubSeq
                          };
#else
  extern RawContigFuncs rawContigFuncs;
#endif


#endif
