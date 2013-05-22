#ifndef __PREDICTIONTRANSCRIPT_H__
#define __PREDICTIONTRANSCRIPT_H__

#include "DataModelTypes.h"

#include "EcoString.h"

#include "AnnotatedSeqFeature.h"
#include "FeatureSet.h"
#include "StableIdInfo.h"
#include "Mapper.h"

ANNOTATEDSEQFEATUREFUNC_TYPES(PredictionTranscript)

typedef struct PredictionTranscriptFuncsStruct {
  ANNOTATEDSEQFEATUREFUNCS_DATA(PredictionTranscript)
} PredictionTranscriptFuncs;



#define FUNCSTRUCTTYPE PredictionTranscriptFuncs
struct PredictionTranscriptStruct {
  ANNOTATEDSEQFEATURE_DATA
  Mapper *exonCoordMapper;
  Vector *exons;
  int codingRegionStart;
  int codingRegionEnd;
  char codingRegionStartIsSet;
  char codingRegionEndIsSet;
  char startIsSet;
  char endIsSet;
  Vector *translateableExons;
  char *displayLabel;
  ECOSTRING type;
};
#undef FUNCSTRUCTTYPE

PredictionTranscript *PredictionTranscript_new(void);

ECOSTRING PredictionTranscript_setDisplayLabel(PredictionTranscript *pt, char *label);
#define PredictionTranscript_getDisplayLabel(transcript) (transcript)->displayLabel

#define PredictionTranscript_setStableId(transcript,sid)  StableIdInfo_setStableId(&((transcript)->si),(sid))
#define PredictionTranscript_getStableId(transcript)  StableIdInfo_getStableId(&((transcript)->si))

int PredictionTranscript_setStart(PredictionTranscript *transcript, int start);
int PredictionTranscript_getStart(PredictionTranscript *transcript);

int PredictionTranscript_setEnd(PredictionTranscript *transcript,int end);
int PredictionTranscript_getEnd(PredictionTranscript *transcript);

#define PredictionTranscript_setStrand(transcript,strand) AnnotatedSeqFeature_setStrand((transcript),strand)
#define PredictionTranscript_getStrand(transcript) AnnotatedSeqFeature_getStrand((transcript))

#define PredictionTranscript_setCodingRegionEndIsSet(transcript,flag) (transcript)->codingRegionEndIsSet = (flag)
#define PredictionTranscript_getCodingRegionEndIsSet(transcript) (transcript)->codingRegionEndIsSet

#define PredictionTranscript_setCodingRegionStartIsSet(transcript,flag) (transcript)->codingRegionStartIsSet = (flag)
#define PredictionTranscript_getCodingRegionStartIsSet(transcript) (transcript)->codingRegionStartIsSet

#define PredictionTranscript_setEndIsSet(transcript,flag) (transcript)->endIsSet = (flag)
#define PredictionTranscript_getEndIsSet(transcript) (transcript)->endIsSet

#define PredictionTranscript_setStartIsSet(transcript,flag) (transcript)->startIsSet = (flag)
#define PredictionTranscript_getStartIsSet(transcript) (transcript)->startIsSet

void PredictionTranscript_addExon(PredictionTranscript *transcript,PredictionExon *exon,int *positionP);
#define PredictionTranscript_getExonAt(transcript,ind) Vector_getElementAt((transcript)->exons,ind)
#define PredictionTranscript_getExons(transcript) (transcript)->exons

int PredictionTranscript_setExonCount(PredictionTranscript *trans, int count);
int PredictionTranscript_getExonCount(PredictionTranscript *trans);

#define PredictionTranscript_setSlice(predictionTranscript,slice) SeqFeature_setSlice((predictionTranscript),(slice))
#define PredictionTranscript_getSlice(predictionTranscript) SeqFeature_getSlice((predictionTranscript))

//#define PredictionTranscript_removeAllExons(transcript) FeatureSet_removeAll(&((transcript)->fs))

#define PredictionTranscript_setDbID(transcript,id) AnnotatedSeqFeature_setDbID((transcript),(id))
#define PredictionTranscript_getDbID(transcript) AnnotatedSeqFeature_getDbID((transcript))

#define PredictionTranscript_setAdaptor(transcript,ad) AnnotatedSeqFeature_setAdaptor((transcript),(ad))
#define PredictionTranscript_getAdaptor(transcript) AnnotatedSeqFeature_getAdaptor((transcript))

#define PredictionTranscript_setAnalysis(transcript,an) AnnotatedSeqFeature_setAnalysis((transcript),(an))
#define PredictionTranscript_getAnalysis(transcript) AnnotatedSeqFeature_getAnalysis((transcript))

#define PredictionTranscript_getSeqRegionStart(t) SeqFeature_getSeqRegionStart((t))
#define PredictionTranscript_getSeqRegionEnd(t) SeqFeature_getSeqRegionEnd((t))
#define PredictionTranscript_getSeqRegionStrand(t) SeqFeature_getSeqRegionStrand((t))

void PredictionTranscript_flushExons(PredictionTranscript *transcript);

void PredictionTranscript_free(PredictionTranscript *transcript);

int PredictionTranscript_getLength(PredictionTranscript *trans);

int PredictionTranscript_setCodingRegionStart(PredictionTranscript *trans, int start);
int PredictionTranscript_getCodingRegionStart(PredictionTranscript *trans);

int PredictionTranscript_setCodingRegionEnd(PredictionTranscript *trans, int start);
int PredictionTranscript_getCodingRegionEnd(PredictionTranscript *trans);

Vector *PredictionTranscript_getAllExons(PredictionTranscript *trans, int wishUndefinedExon);
Vector *PredictionTranscript_getAllTranslateableExons(PredictionTranscript *trans);

void PredictionTranscript_sort(PredictionTranscript *trans);

char *PredictionTranscript_translate(PredictionTranscript *trans);

char *PredictionTranscript_getcDNA(PredictionTranscript *trans);

MapperRangeSet *PredictionTranscript_pep2Genomic(PredictionTranscript *trans, int start, int end);

MapperRangeSet *PredictionTranscript_cDNA2Genomic(PredictionTranscript *trans, int start, int end);
MapperRangeSet *PredictionTranscript_genomic2cDNA(PredictionTranscript *trans, int start, int end, int strand, BaseContig *contig);

Mapper *PredictionTranscript_getcDNACoordMapper(PredictionTranscript *trans);

void PredictionTranscript_free(PredictionTranscript *trans);

#ifdef __PREDICTIONTRANSCRIPT_MAIN__
  PredictionTranscriptFuncs 
    predictionTranscriptFuncs = {
                       PredictionTranscript_free,
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
