#ifndef __PREDICTIONTRANSCRIPT_H__
#define __PREDICTIONTRANSCRIPT_H__

#include "DataModelTypes.h"

#include "SeqFeature.h"
#include "FeatureSet.h"
#include "StableIdInfo.h"

struct PredictionTranscriptStruct {
  SeqFeature sf;
  FeatureSet fs;
  StableIdInfo si;
  char *type;
};

PredictionTranscript *PredictionTranscript_new(void);

#define PredictionTranscript_setStableId(transcript,sid)  StableIdInfo_setStableId(&((transcript)->si),(sid))
#define PredictionTranscript_getStableId(transcript)  StableIdInfo_getStableId(&((transcript)->si))

#define PredictionTranscript_setStart(transcript,start) SeqFeature_setStart(&((transcript)->sf),start)
#define PredictionTranscript_getStart(transcript) SeqFeature_getStart(&((transcript)->sf))

#define PredictionTranscript_setEnd(transcript,end) SeqFeature_setEnd(&((transcript)->sf),end)
#define PredictionTranscript_getEnd(transcript) SeqFeature_getEnd(&((transcript)->sf))

#define PredictionTranscript_setStrand(transcript,strand) SeqFeature_setStrand(&((transcript)->sf),strand)
#define PredictionTranscript_getStrand(transcript) SeqFeature_getStrand(&((transcript)->sf))

#define PredictionTranscript_addExon(trans,exon) FeatureSet_addFeature(&((trans)->fs),exon)
#define PredictionTranscript_getExonAt(trans,ind) FeatureSet_getFeatureAt(&((trans)->fs),ind)

//#define PredictionTranscript_setExonCount(trans,ec) FeatureSet_setNumFeature(&((trans)->fs))
#define PredictionTranscript_setExonCount(trans,ec)
#define PredictionTranscript_getExonCount(trans) FeatureSet_getNumFeature(&((trans)->fs))

#define PredictionTranscript_removeAllExons(trans) FeatureSet_removeAll(&((trans)->fs))

#define PredictionTranscript_setDbID(transcript,id) SeqFeature_setDbID(&((transcript)->sf),(id))
#define PredictionTranscript_getDbID(transcript) SeqFeature_getDbID(&((transcript)->sf))

#define PredictionTranscript_setAdaptor(transcript,ad) SeqFeature_setAdaptor(&((transcript)->sf),(ad))
#define PredictionTranscript_getAdaptor(transcript) SeqFeature_getAdaptor(&((transcript)->sf))

#define PredictionTranscript_setAnalysis(transcript,an) SeqFeature_setAnalysis(&((transcript)->sf),(an))
#define PredictionTranscript_getAnalysis(transcript) SeqFeature_getAnalysis(&((transcript)->sf))

void PredictionTranscript_flushExons(PredictionTranscript *trans);

void PredictionTranscript_free(PredictionTranscript *trans);

#endif
