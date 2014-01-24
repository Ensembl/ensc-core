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

#ifndef __GENOMICALIGN_H__
#define __GENOMICALIGN_H__

#include "ComparaDataModelTypes.h"
#include "EnsRoot.h"
#include "Storable.h"
#include "DNAFrag.h"

OBJECTFUNC_TYPES(GenomicAlign)

typedef struct GenomicAlignFuncsStruct {
  OBJECTFUNCS_DATA(GenomicAlign)
} GenomicAlignFuncs;

#define FUNCSTRUCTTYPE GenomicAlignFuncs
struct GenomicAlignStruct {
  ENSROOT_DATA
  Storable st;
  char *cigarString;
  DNAFrag *consensusDNAFrag;
  int consensusStart;
  int consensusEnd;
  DNAFrag *queryDNAFrag;
  int queryStart;
  int queryEnd;
  signed char queryStrand;
  ECOSTRING alignmentType;
  double percentId;
  double score;
};
#undef FUNCSTRUCTTYPE


#define GenomicAlign_setDbID(ga, id) Storable_setDbID(&((ga)->st), (id))
#define GenomicAlign_getDbID(ga) Storable_getDbID(&((ga)->st)) 

#define GenomicAlign_setAdaptor(ga, a) Storable_setAdaptor(&((ga)->st), (a))
#define GenomicAlign_getAdaptor(ga) Storable_getAdaptor(&((ga)->st)) 

char *GenomicAlign_setCigarString(GenomicAlign *ga, char *cs);
#define GenomicAlign_getCigarString(ga) (ga)->cigarString

char *GenomicAlign_setAlignmentType(GenomicAlign *ga, char *at);
#define GenomicAlign_getAlignmentType(ga) (ga)->alignmentType

#define GenomicAlign_setConsensusStart(ga, cs) (ga)->consensusStart = (cs)
#define GenomicAlign_getConsensusStart(ga) (ga)->consensusStart

#define GenomicAlign_setConsensusEnd(ga, ce) (ga)->consensusEnd = (ce)
#define GenomicAlign_getConsensusEnd(ga) (ga)->consensusEnd

#define GenomicAlign_setQueryStart(ga, qs) (ga)->queryStart = (qs)
#define GenomicAlign_getQueryStart(ga) (ga)->queryStart

#define GenomicAlign_setQueryEnd(ga, qe) (ga)->queryEnd = (qe)
#define GenomicAlign_getQueryEnd(ga) (ga)->queryEnd

#define GenomicAlign_setQueryStrand(ga, qs) (ga)->queryStrand = (qs)
#define GenomicAlign_getQueryStrand(ga) (ga)->queryStrand

#define GenomicAlign_setPercentId(ga, pid) (ga)->percentId = (pid)
#define GenomicAlign_getPercentId(ga) (ga)->percentId

#define GenomicAlign_setScore(ga, s) (ga)->score = (s)
#define GenomicAlign_getScore(ga) (ga)->score

#define GenomicAlign_setConsensusDNAFrag(ga, df) (ga)->consensusDNAFrag = (df)
#define GenomicAlign_getConsensusDNAFrag(ga) (ga)->consensusDNAFrag

#define GenomicAlign_setQueryDNAFrag(ga, df) (ga)->queryDNAFrag = (df)
#define GenomicAlign_getQueryDNAFrag(ga) (ga)->queryDNAFrag

GenomicAlign *GenomicAlign_new();

void GenomicAlign_free(GenomicAlign *ga);

#define GENOMICALIGN_CONSENSUS    1<<1
#define GENOMICALIGN_ALIGNSLICES  1<<2
#define GENOMICALIGN_FIXCONSENSUS 1<<3
#define GENOMICALIGN_FIXQUERY     1<<4

#ifdef __GENOMICALIGN_MAIN__
  GenomicAlignFuncs
    genomicAlignFuncs = {
                    GenomicAlign_free
                   };
#else
  extern GenomicAlignFuncs genomicAlignFuncs;
#endif

#endif
