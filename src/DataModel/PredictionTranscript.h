#ifndef __PREDICTIONTRANSCRIPT_H__
#define __PREDICTIONTRANSCRIPT_H__

#include "DataModelTypes.h"

#include "AnnotatedSeqFeature.h"
#include "FeatureSet.h"
#include "StableIdInfo.h"

ANNOTATEDSEQFEATUREFUNC_TYPES(PredictionTranscript)

typedef struct PredictionTranscriptFuncsStruct {
  ANNOTATEDSEQFEATUREFUNCS_DATA(PredictionTranscript)
} PredictionTranscriptFuncs;



#define FUNCSTRUCTTYPE PredictionTranscriptFuncs
struct PredictionTranscriptStruct {
  ANNOTATEDSEQFEATURE_DATA
  FeatureSet fs;
  char *type;
};
#undef FUNCSTRUCTTYPE

PredictionTranscript *PredictionTranscript_new(void);

#define PredictionTranscript_setStableId(transcript,sid)  StableIdInfo_setStableId(&((transcript)->si),(sid))
#define PredictionTranscript_getStableId(transcript)  StableIdInfo_getStableId(&((transcript)->si))

#define PredictionTranscript_setStart(transcript,start) AnnotatedSeqFeature_setStart((transcript),start)
#define PredictionTranscript_getStart(transcript) AnnotatedSeqFeature_getStart((transcript))

#define PredictionTranscript_setEnd(transcript,end) AnnotatedSeqFeature_setEnd((transcript),end)
#define PredictionTranscript_getEnd(transcript) AnnotatedSeqFeature_getEnd((transcript))

#define PredictionTranscript_setStrand(transcript,strand) AnnotatedSeqFeature_setStrand((transcript),strand)
#define PredictionTranscript_getStrand(transcript) AnnotatedSeqFeature_getStrand((transcript))

#define PredictionTranscript_addExon(transcript,exon) FeatureSet_addFeature(&((transcript)->fs),exon)
#define PredictionTranscript_getExonAt(transcript,ind) FeatureSet_getFeatureAt(&((transcript)->fs),ind)

//#define PredictionTranscript_setExonCount(transcript,ec) FeatureSet_setNumFeature(&((transcript)->fs))
#define PredictionTranscript_setExonCount(transcript,ec)
#define PredictionTranscript_getExonCount(transcript) FeatureSet_getNumFeature(&((transcript)->fs))

#define PredictionTranscript_removeAllExons(transcript) FeatureSet_removeAll(&((transcript)->fs))

#define PredictionTranscript_setDbID(transcript,id) AnnotatedSeqFeature_setDbID((transcript),(id))
#define PredictionTranscript_getDbID(transcript) AnnotatedSeqFeature_getDbID((transcript))

#define PredictionTranscript_setAdaptor(transcript,ad) AnnotatedSeqFeature_setAdaptor((transcript),(ad))
#define PredictionTranscript_getAdaptor(transcript) AnnotatedSeqFeature_getAdaptor((transcript))

#define PredictionTranscript_setAnalysis(transcript,an) AnnotatedSeqFeature_setAnalysis((transcript),(an))
#define PredictionTranscript_getAnalysis(transcript) AnnotatedSeqFeature_getAnalysis((transcript))

void PredictionTranscript_flushExons(PredictionTranscript *transcript);

void PredictionTranscript_free(PredictionTranscript *transcript);

int PredictionTranscript_getLength(PredictionTranscript *trans);


#ifdef __PREDICTIONTRANSCRIPT_MAIN__
  PredictionTranscriptFuncs 
    predictionTranscriptFuncs = {
 /* Not done yet
                       PredictionTranscript_getStart,
                       PredictionTranscript_setStart,
                       PredictionTranscript_getEnd,
                       PredictionTranscript_setEnd,
 */
                       NULL, // getStart  
                       NULL, // setStart 
                       NULL, // getEnd
                       NULL, // setEnd
                       NULL, // getStrand
                       NULL, // setStrand
                       NULL, // getSeq
                       NULL, // setSeq
                       PredictionTranscript_getLength,
                       NULL, // reverseComplement
                       NULL, // transformToRawContig
                       NULL, // transformToSlice
                       NULL, // transformRawContigToSlice
                       NULL, // transformSliceToRawContig
                       NULL  // transformSliceToSlice
                      };
#else
  extern PredictionTranscriptFuncs predictionTranscriptFuncs;
#endif

#endif
