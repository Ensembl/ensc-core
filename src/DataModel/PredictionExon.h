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
  char *displayLabel;

#define FUNCSTRUCTTYPE PredictionExonFuncs
struct PredictionExonStruct {
  PREDICTIONEXON_DATA
};
#undef FUNCSTRUCTTYPE

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

#define PredictionExon_setContig(predictionExon,c) SeqFeature_setContig((predictionExon),(c))
#define PredictionExon_getContig(predictionExon) SeqFeature_getContig((predictionExon))

#define PredictionExon_free(predictionExon) SeqFeature_free((predictionExon))

PredictionExon *PredictionExon_new();

ECOSTRING PredictionExon_setDisplayLabel(PredictionExon *sf, char *label);
#define PredictionExon_getDisplayLabel(predictionExon) (predictionExon)->displayLabel

void PredictionExon_freeImpl(PredictionExon *sf);

#ifdef __PREDICTIONEXON_MAIN__
 PredictionExonFuncs
   predictionExonFuncs = {
                      PredictionExon_freeImpl,
                      NULL, // shallowCopy
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
