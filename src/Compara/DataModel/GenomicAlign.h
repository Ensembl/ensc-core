#ifndef __GENOMICALIGN_H__
#define __GENOMICALIGN_H__

#include "ComparaDataModelTypes.h"
#include "EnsRoot.h"
#include "Storable.h"
#include "DNAFrag.h"

#define FUNCSTRUCTTYPE NoTypeFuncs 
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
  char *alignmentType;
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

#define GENOMICALIGN_CONSENSUS    1<<1
#define GENOMICALIGN_ALIGNSLICES  1<<2
#define GENOMICALIGN_FIXCONSENSUS 1<<3
#define GENOMICALIGN_FIXQUERY     1<<4
#endif
