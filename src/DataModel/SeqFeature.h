#ifndef __SEQFEATURE_H__
#define __SEQFEATURE_H__

#include "EnsC.h"
#include "Storable.h"
#include "Sequence.h"
#include "Analysis.h"
#include "BaseContig.h"

typedef struct SeqFeatureStruct {
  int         start;
  int         end;
  signed char phase;
  signed char endPhase;
  signed char frame;
  signed char strand;
  char *      seqName;
  Storable    st;
  Analysis *  analysis;
  double      score;
  double      pValue;
  double      percentId;
  BaseContig *contig;
} SeqFeature; 

#define SeqFeature_setStart(sf,s) (sf)->start = (s)
#define SeqFeature_getStart(sf) (sf)->start

#define SeqFeature_setEnd(sf,e) (sf)->end = (e)
#define SeqFeature_getEnd(sf) (sf)->end

#define SeqFeature_setPhase(sf,p) (sf)->phase = (p)
#define SeqFeature_getPhase(sf) (sf)->phase

#define SeqFeature_setEndPhase(sf,ep) (sf)->endPhase = (ep)
#define SeqFeature_getEndPhase(sf) (sf)->endPhase

#define SeqFeature_setStrand(sf,strnd) (sf)->strand = strnd
#define SeqFeature_getStrand(sf) (sf)->strand

#define SeqFeature_setAnalysis(sf,ana) (sf)->analysis = ana
#define SeqFeature_getAnalysis(sf) (sf)->analysis

char *SeqFeature_setStableId(SeqFeature *sf, char *stableId);
#define SeqFeature_getStableId(sf) (sf)->stableId

#define SeqFeature_setDbID(sf,dbID) Storable_setDbID(&((sf)->st),dbID)
#define SeqFeature_getDbID(sf) Storable_getDbID(&((sf)->st))

#define SeqFeature_setAdaptor(sf,ad) Storable_setAdaptor(&((sf)->st),ad)
#define SeqFeature_getAdaptor(sf) Storable_getAdaptor(&((sf)->st))

#define SeqFeature_setContig(sf,c) (sf)->contig = (BaseContig *)(c)
#define SeqFeature_getContig(sf) (sf)->contig

#define SeqFeature_getLength(sf) ((sf)->end - (sf)->start + 1)

#endif
