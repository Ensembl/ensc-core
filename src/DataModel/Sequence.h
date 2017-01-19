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

#ifndef __SEQUENCE_H__
#define __SEQUENCE_H__

#include "DataModelTypes.h"
#include "EnsRoot.h"
#include "EnsC.h"

#define SEQUENCEFUNC_TYPES(CLASSTYPE) \
  OBJECTFUNC_TYPES(CLASSTYPE) \
  typedef ECOSTRING (*CLASSTYPE ## _GetNameFunc)(CLASSTYPE *seq); \
  typedef char * (*CLASSTYPE ## _GetSeqFunc)(CLASSTYPE *seq); \
  typedef char * (*CLASSTYPE ## _GetSubSeqFunc)(CLASSTYPE *seq, int start, int end, int strand);


#define SEQUENCEFUNCS_DATA(CLASSTYPE) \
  OBJECTFUNCS_DATA(CLASSTYPE) \
  CLASSTYPE ## _GetNameFunc getName; \
  CLASSTYPE ## _GetSeqFunc getSeq; \
  CLASSTYPE ## _GetSubSeqFunc getSubSeq;

SEQUENCEFUNC_TYPES(Sequence)

typedef struct SequenceFuncsStruct {
  SEQUENCEFUNCS_DATA(Sequence)
} SequenceFuncs;

#define SEQUENCE_DATA \
  ENSROOT_DATA \
  char *seq; \
  ECOSTRING name; \
  int length;

#define FUNCSTRUCTTYPE SequenceFuncs
struct SequenceStruct {
  SEQUENCE_DATA
};
#undef FUNCSTRUCTTYPE

void Sequence_freePtrs(Sequence *seq);
Sequence *Sequence_new();

ECOSTRING Sequence_setName(Sequence *seq, char *name);
char *Sequence_setSeq(Sequence *seq, char *seqstr);

#ifdef __SEQUENCE_MAIN__
 SequenceFuncs sequenceFuncs;
#else
 extern SequenceFuncs sequenceFuncs;
#endif


#endif
