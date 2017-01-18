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

#ifndef __PREDICTIONEXON_H__
#define __PREDICTIONEXON_H__

#include "EnsC.h"
#include "DataModelTypes.h"
#include "SeqFeature.h"
#include "Slice.h"

SEQFEATUREFUNC_TYPES(PredictionExon)

typedef struct PredictionExonFuncsStruct {
  SEQFEATUREFUNCS_DATA(PredictionExon)
} PredictionExonFuncs;



#define PREDICTIONEXON_DATA \
  SEQFEATURE_DATA \
  char *displayLabel; \
  char *seqCacheString;

#define FUNCSTRUCTTYPE PredictionExonFuncs
struct PredictionExonStruct {
  PREDICTIONEXON_DATA
};
#undef FUNCSTRUCTTYPE

#define PredictionExon_setSeqCacheString(exon,seq) (exon)->seqCacheString = (seq)
#define PredictionExon_getSeqCacheString(exon) (exon)->seqCacheString

#define PredictionExon_setStart(predictionExon,start) SeqFeature_setStart((predictionExon),start)
#define PredictionExon_getStart(predictionExon) SeqFeature_getStart((predictionExon))

#define PredictionExon_setEnd(predictionExon,end) SeqFeature_setEnd((predictionExon),end)
#define PredictionExon_getEnd(predictionExon) SeqFeature_getEnd((predictionExon))

#define PredictionExon_setScore(predictionExon,score) SeqFeature_setScore((predictionExon),(score))
#define PredictionExon_getScore(predictionExon) SeqFeature_getScore((predictionExon))

#define PredictionExon_setpValue(predictionExon,pv) SeqFeature_setpValue((predictionExon),(pv))
#define PredictionExon_getpValue(predictionExon) SeqFeature_getpValue((predictionExon))

#define PredictionExon_setPhase(predictionExon,p) SeqFeature_setPhase((predictionExon),(p))
#define PredictionExon_getPhase(predictionExon) SeqFeature_getPhase((predictionExon))

#define PredictionExon_setStrand(predictionExon,strand) SeqFeature_setStrand((predictionExon),(strand))
#define PredictionExon_getStrand(predictionExon) SeqFeature_getStrand((predictionExon))

#define PredictionExon_setDbID(predictionExon,dbID) SeqFeature_setDbID((predictionExon),(dbID))
#define PredictionExon_getDbID(predictionExon) SeqFeature_getDbID((predictionExon))

#define PredictionExon_setAdaptor(predictionExon,ad) SeqFeature_setAdaptor((predictionExon),(ad))
#define PredictionExon_getAdaptor(predictionExon) SeqFeature_getAdaptor((predictionExon))

#define PredictionExon_setAnalysis(predictionExon,ana) SeqFeature_setAnalysis((predictionExon),(ana))
#define PredictionExon_getAnalysis(predictionExon) SeqFeature_getAnalysis((predictionExon))

#define PredictionExon_getLength(predictionExon) SeqFeature_getLength((predictionExon))

#define PredictionExon_setSlice(predictionExon,slice) SeqFeature_setSlice((predictionExon),(slice))
#define PredictionExon_getSlice(predictionExon) SeqFeature_getSlice((predictionExon))

#define PredictionExon_free(predictionExon) SeqFeature_free((predictionExon))

PredictionExon *PredictionExon_new();

ECOSTRING PredictionExon_setDisplayLabel(PredictionExon *sf, char *label);
#define PredictionExon_getDisplayLabel(predictionExon) (predictionExon)->displayLabel

#define PredictionExon_shallowCopy(pe) SeqFeature_shallowCopy((pe))

#define PredictionExon_transfer(pe, slice) SeqFeature_transfer((pe), (slice))

char  *PredictionExon_getSeqString(PredictionExon *exon);

void PredictionExon_freeImpl(PredictionExon *sf);
PredictionExon *PredictionExon_shallowCopyImpl(PredictionExon *exon);

void PredictionExon_loadGenomicMapper(Exon *exon, Mapper *mapper, IDType id, int start);

#ifdef __PREDICTIONEXON_MAIN__
 PredictionExonFuncs
   predictionExonFuncs = {
                      PredictionExon_freeImpl,
                      PredictionExon_shallowCopyImpl, // shallowCopy
                      NULL, // deepCopy
                      NULL, // getStart
                      NULL, // setStart
                      NULL, // getEnd
                      NULL, // setEnd
                      NULL, // getStrand
                      NULL, // setStrand
                      NULL, // getSeq
                      NULL, // setSeq
                      NULL, // getLength
                      NULL, // reverseComplement
                     };
#else
 extern PredictionExonFuncs predictionExonFuncs;
#endif


#endif
