#ifndef __TRANSCRIPT_H__
#define __TRANSCRIPT_H__

#include "DataModelTypes.h"

#include "AnnotatedSeqFeature.h"
#include "Storable.h"
#include "Translation.h"
#include "IDHash.h"
#include "Vector.h"
#include "Mapper.h"
#include "IntronSupportingEvidence.h"

ANNOTATEDSEQFEATUREFUNC_TYPES(Transcript)

typedef struct TranscriptFuncsStruct {
  ANNOTATEDSEQFEATUREFUNCS_DATA(Transcript)
} TranscriptFuncs;

#define FUNCSTRUCTTYPE TranscriptFuncs
struct TranscriptStruct {
  ANNOTATEDSEQFEATURE_DATA
  ECOSTRING    biotype;
  ECOSTRING    status;
  ECOSTRING    externalDb;
  ECOSTRING    externalStatus;
  char        *externalName;
  char        *description;
  char         editsEnabled;
  char         isCanonical;
  Vector      *exons;
  Vector      *dbLinks;
  Vector      *attributes;
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
  Vector      *iseVector;
  Vector      *supportingEvidence;
};
#undef FUNCSTRUCTTYPE

#define Transcript_isStored(transcript, db) Storable_isStored(&((transcript)->st), (db))

char *Transcript_getSplicedSeq(Transcript *trans);
ECOSTRING Transcript_setBiotype(Transcript *transcript, char *biotype);
#define Transcript_getBiotype(transcript)  (transcript)->biotype

void Transcript_freeAdjustedTranslateableExons(Transcript *transcript, Vector *translateableExons);

Vector *Transcript_getAllSupportingFeatures(Transcript *transcript);
Vector *Transcript_getAllIntronSupportingEvidence(Transcript *transcript);

int Transcript_addIntronSupportingEvidence(Transcript *transcript, IntronSupportingEvidence *ise);

Vector *Transcript_getAllAttributes(Transcript *transcript, char *attribCode);

#define Transcript_setIsCurrent(trans,isC)  StableIdInfo_setIsCurrent(&((trans)->si),(isC))
#define Transcript_getIsCurrent(trans)  StableIdInfo_getIsCurrent(&((trans)->si))

ECOSTRING DNAAlignFeature_setExtraData(DNAAlignFeature *fp, char *extraData);
#define Transcript_setExtraData(trans, ed)  (trans)->extraData = (ed)
#define Transcript_getExtraData(trans)  (trans)->extraData

#define Transcript_setIsCanonical(transcript, flag)  (transcript)->isCanonical = (flag)
int Transcript_getIsCanonical(Transcript *transcript);

#define Transcript_setAnalysis(trans,ana) AnnotatedSeqFeature_setAnalysis((trans),(ana))
#define Transcript_getAnalysis(trans) AnnotatedSeqFeature_getAnalysis((trans))

#define Transcript_setStableId(transcript,sid)  StableIdInfo_setStableId(&((transcript)->si),(sid))
char *Transcript_getStableId(Transcript *transcript);

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

void Transcript_addExon(Transcript *transcript, Exon *exon, int rank);
void Transcript_recalculateCoordinates(Transcript *transcript);
//#define Transcript_addExon(transcript,exon, rank) Vector_addElement((transcript)->exons, (exon))
#define Transcript_getExonAt(transcript,ind) Vector_getElementAt((transcript)->exons, (ind))

#define Transcript_getExonCount(transcript) Vector_getNumElement((transcript)->exons)
#define Transcript_sortExons(transcript, func) Vector_sort((transcript)->exons, (func))

// Note: I've decided to implement this call - not sure about whether thats a good idea, but lets give it a go
#define Transcript_getAllExons(transcript) (transcript)->exons

#define Transcript_removeAllExons(transcript) Vector_removeAll((transcript)->exons)

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

#define Transcript_setSlice(trans,sl) AnnotatedSeqFeature_setSlice((trans),(sl))
#define Transcript_getSlice(trans) AnnotatedSeqFeature_getSlice((trans))

#define Transcript_setScore(trans,s) SeqFeature_setScore((trans),(s))
#define Transcript_getScore(trans) SeqFeature_getScore((trans))

#define Transcript_setEditsEnabled(transcript, flag)  (transcript)->editsEnabled = (flag)
#define Transcript_getEditsEnabled(transcript)  (transcript)->editsEnabled

#define Transcript_setDisplayXref(trans,xref) AnnotatedSeqFeature_setDisplayXref((trans),xref)
#define Transcript_getDisplayXref(trans) AnnotatedSeqFeature_getDisplayXref((trans))

char *Transcript_setDescription(Transcript *t, char *description);
#define Transcript_getDescription(transcript)  (transcript)->description

ECOSTRING Transcript_setStatus(Transcript *t, char *status);
#define Transcript_getStatus(transcript)  (transcript)->status

ECOSTRING Transcript_setExternalDb(Transcript *t, char *externalDb);
ECOSTRING Transcript_setExternalStatus(Transcript *t, char *externalStatus);
char *Transcript_setExternalName(Transcript *t, char *externalName);


#define Transcript_getSeqRegionStart(t) SeqFeature_getSeqRegionStart((t))
#define Transcript_getSeqRegionEnd(t) SeqFeature_getSeqRegionEnd((t))
#define Transcript_getSeqRegionStrand(t) SeqFeature_getSeqRegionStrand((t))
#define Transcript_getSeqRegionName(t) SeqFeature_getSeqRegionName((t))


Exon *Transcript_getStartExon(Transcript *trans);
Exon *Transcript_getEndExon(Transcript *trans);

void Transcript_flushExons(Transcript *transcript);

int Transcript_addDBLink(Transcript *transcript, DBEntry *dbe);

Transcript *Transcript_transform(Transcript *transcript, IDHash *exonTransforms);
Transcript *Transcript_transfer(Transcript *transcript, Slice *slice);

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

void Transcript_swapExons(Transcript *transcript, Exon *oldExon, Exon *newExon);


Mapper *Transcript_getcDNACoordMapper(Transcript *trans);

char *Transcript_translate(Transcript *trans);

void Transcript_free(Transcript *trans);
Transcript *Transcript_shallowCopy(Transcript *trans);
int Transcript_getLength(Transcript *trans);


#ifdef __TRANSCRIPT_MAIN__
  TranscriptFuncs 
    transcriptFuncs = {
                       Transcript_free, // free
                       Transcript_shallowCopy, // shallowCopy
                       NULL, // deepCopy
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
