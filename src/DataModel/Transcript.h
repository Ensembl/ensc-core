#ifndef __TRANSCRIPT_H__
#define __TRANSCRIPT_H__

#include "DataModelTypes.h"

#include "AnnotatedSeqFeature.h"
#include "FeatureSet.h"
#include "Storable.h"
#include "Translation.h"
#include "IDHash.h"
#include "Vector.h"
#include "Mapper.h"

ANNOTATEDSEQFEATUREFUNC_TYPES(Transcript)

typedef struct TranscriptFuncsStruct {
  ANNOTATEDSEQFEATUREFUNCS_DATA(Transcript)
} TranscriptFuncs;

#define FUNCSTRUCTTYPE TranscriptFuncs
struct TranscriptStruct {
  ANNOTATEDSEQFEATURE_DATA
  FeatureSet   fs;
  Vector      *dbLinks;
  Translation *translation;
  IDType       translationId;
  char         codingRegionStartIsSet;
  char         codingRegionEndIsSet;
  char         cDNACodingStartIsSet;
  char         cDNACodingEndIsSet;
  int          codingRegionStart;
  int          codingRegionEnd;
  int          cDNACodingStart;
  int          cDNACodingEnd;
  Mapper      *exonCoordMapper;
  char        *type;
};
#undef FUNCSTRUCTTYPE

#define Transcript_setStableId(transcript,sid)  StableIdInfo_setStableId(&((transcript)->si),(sid))
char *Transcript_getStableId(Transcript *transcript);

char *Transcript_setType(Transcript *transcript,char *type);
#define Transcript_getType(transcript) (transcript)->type

#define Transcript_setVersion(transcript,ver)  StableIdInfo_setVersion(&((transcript)->si),(ver))
int Transcript_getVersion(Transcript *transcript);

#define Transcript_setStart(transcript,start) AnnotatedSeqFeature_setStart((transcript),start)
#define Transcript_getStart(transcript) AnnotatedSeqFeature_getStart((transcript))

#define Transcript_setEnd(transcript,end) AnnotatedSeqFeature_setEnd((transcript),end)
#define Transcript_getEnd(transcript) AnnotatedSeqFeature_getEnd((transcript))

#define Transcript_setStrand(transcript,strand) AnnotatedSeqFeature_setStrand((transcript),strand)
#define Transcript_getStrand(transcript) AnnotatedSeqFeature_getStrand((transcript))

#define Transcript_setCreated(transcript,cd)  StableIdInfo_setCreated(&((transcript)->si),cd)
#define Transcript_getCreated(transcript)  StableIdInfo_getCreated(&((transcript)->si))

#define Transcript_setModified(transcript,mod)  StableIdInfo_setModified(&((transcript)->si),mod)
#define Transcript_getModified(transcript)  StableIdInfo_getModified(&((transcript)->si))

Transcript *Transcript_new(void);

#define Transcript_addExon(transcript,exon) FeatureSet_addFeature(&((transcript)->fs),exon)
#define Transcript_getExonAt(transcript,ind) FeatureSet_getFeatureAt(&((transcript)->fs),ind)

#define Transcript_getExonCount(transcript) FeatureSet_getNumFeature(&((transcript)->fs))
#define Transcript_getExons(transcript) FeatureSet_getFeatures(&((transcript)->fs))

#define Transcript_removeAllExons(transcript) FeatureSet_removeAll(&((transcript)->fs))

#define Transcript_setTranslationId(transcript,tid) (transcript)->translationId = (tid)
#define Transcript_getTranslationId(transcript) (transcript)->translationId

#define Transcript_setTranslation(transcript,tn) (transcript)->translation = (tn)
Translation *Transcript_getTranslation(Transcript *trans);

#define Transcript_setDbID(transcript,id) AnnotatedSeqFeature_setDbID((transcript),(id))
#define Transcript_getDbID(transcript) AnnotatedSeqFeature_getDbID((transcript))

#define Transcript_setAdaptor(transcript,ad) AnnotatedSeqFeature_setAdaptor((transcript),(ad))
#define Transcript_getAdaptor(transcript) AnnotatedSeqFeature_getAdaptor((transcript))

#define Transcript_setCodingRegionStartIsSet(transcript, flag)  (transcript)->codingRegionStartIsSet = (flag)
#define Transcript_getCodingRegionStartIsSet(transcript)  (transcript)->codingRegionStartIsSet

#define Transcript_setCodingRegionEndIsSet(transcript, flag)  (transcript)->codingRegionEndIsSet = (flag)
#define Transcript_getCodingRegionEndIsSet(transcript)  (transcript)->codingRegionEndIsSet

#define Transcript_setcDNACodingStartIsSet(transcript, flag)  (transcript)->cDNACodingStartIsSet = (flag)
#define Transcript_getcDNACodingStartIsSet(transcript)  (transcript)->cDNACodingStartIsSet

#define Transcript_setcDNACodingEndIsSet(transcript, flag)  (transcript)->cDNACodingEndIsSet = (flag)
#define Transcript_getcDNACodingEndIsSet(transcript)  (transcript)->cDNACodingEndIsSet

Exon *Transcript_getStartExon(Transcript *trans);
Exon *Transcript_getEndExon(Transcript *trans);

void Transcript_flushExons(Transcript *transcript);

int Transcript_addDBLink(Transcript *transcript, DBEntry *dbe);

Transcript *Transcript_transform(Transcript *transcript, IDHash *exonTransforms);

void Transcript_sort(Transcript *trans);

int Transcript_getCodingRegionEnd(Transcript *trans);
int Transcript_setCodingRegionEnd(Transcript *trans, int end);

int Transcript_setCodingRegionStart(Transcript *trans, int start);
int Transcript_getCodingRegionStart(Transcript *trans); 

Vector *Transcript_getAllTranslateableExons(Transcript *trans);

int Transcript_setcDNACodingStart(Transcript *trans, int start);
int Transcript_getcDNACodingStart(Transcript *trans);

int Transcript_setcDNACodingEnd(Transcript *trans, int end);
int Transcript_getcDNACodingEnd(Transcript *trans);

MapperRangeSet *Transcript_cDNA2Genomic(Transcript *trans, int start, int end);
MapperRangeSet *Transcript_genomic2cDNA(Transcript *trans, int start, int end, int strand, BaseContig *contig);
MapperRangeSet *Transcript_genomic2Pep(Transcript *trans, int start, int end, int strand, BaseContig *contig);


Mapper *Transcript_getcDNACoordMapper(Transcript *trans);

char *Transcript_translate(Transcript *trans);

void Transcript_free(Transcript *trans);
int Transcript_getLength(Transcript *trans);


#ifdef __TRANSCRIPT_MAIN__
  TranscriptFuncs 
    transcriptFuncs = {
                       Transcript_free,
                       NULL,
                       NULL,
                       NULL,
                       NULL,
                       NULL, // getStrand
                       NULL, // setStrand
                       NULL, // getSeq
                       NULL, // setSeq
                       Transcript_getLength,
                       NULL, // reverseComplement
                       NULL, // transformToRawContig
                       NULL, // transformToSlice
                       NULL, // transformRawContigToSlice
                       NULL, // transformSliceToRawContig
                       NULL  // transformSliceToSlice
                      };
#else
  extern TranscriptFuncs transcriptFuncs;
#endif

#endif
