#ifndef __SEQFEATURE_H__
#define __SEQFEATURE_H__


#include "DataModelTypes.h"
#include "EnsC.h"
#include "Storable.h"
#include "Analysis.h"
#include "BaseContig.h"
#include "Vector.h"

#include "EnsRoot.h"

#define SEQFEATUREFUNC_TYPES(CLASSTYPE) \
typedef int (*CLASSTYPE ## _GetStartFunc)(CLASSTYPE *); \
typedef int (*CLASSTYPE ## _SetStartFunc)(CLASSTYPE *, int start); \
typedef int (*CLASSTYPE ## _GetEndFunc)(CLASSTYPE *); \
typedef int (*CLASSTYPE ## _SetEndFunc)(CLASSTYPE *, int end); \
typedef int (*CLASSTYPE ## _GetStrandFunc)(CLASSTYPE *); \
typedef int (*CLASSTYPE ## _SetStrandFunc)(CLASSTYPE *, int strand); \
typedef Sequence * (*CLASSTYPE ## _GetSeqFunc)(CLASSTYPE *); \
typedef Sequence * (*CLASSTYPE ## _SetSeqFunc)(CLASSTYPE *, Sequence *seq); \
typedef int (*CLASSTYPE ## _GetLengthFunc)(CLASSTYPE *); \
typedef void (*CLASSTYPE ## _ReverseComplementFunc)(CLASSTYPE *); \
typedef Vector * (*CLASSTYPE ## _TransformToRawContigFunc)(CLASSTYPE *sf); \
typedef Vector * (*CLASSTYPE ## _TransformSliceToRawContigFunc)(CLASSTYPE *sf); \
typedef CLASSTYPE * (*CLASSTYPE ## _TransformToSliceFunc)(CLASSTYPE *sf, Slice *slice); \
typedef CLASSTYPE * (*CLASSTYPE ## _TransformRawContigToSliceFunc)(CLASSTYPE *sf, Slice *slice); \
typedef CLASSTYPE * (*CLASSTYPE ## _TransformSliceToSliceFunc)(CLASSTYPE *sf, Slice *slice);

#define SEQFEATUREFUNCS_DATA(CLASSTYPE) \
  CLASSTYPE ## _GetStartFunc getStart; \
  CLASSTYPE ## _SetStartFunc setStart; \
  CLASSTYPE ## _GetEndFunc getEnd; \
  CLASSTYPE ## _SetEndFunc setEnd; \
  CLASSTYPE ## _GetStrandFunc getStrand; \
  CLASSTYPE ## _SetStrandFunc setStrand; \
  CLASSTYPE ## _GetSeqFunc getSeq; \
  CLASSTYPE ## _SetSeqFunc setSeq; \
  CLASSTYPE ## _GetLengthFunc getLength; \
  CLASSTYPE ## _ReverseComplementFunc reverseComplement; \
  CLASSTYPE ## _TransformToRawContigFunc transformToRawContig; \
  CLASSTYPE ## _TransformToSliceFunc transformToSlice; \
  CLASSTYPE ## _TransformRawContigToSliceFunc transformRawContigToSlice; \
  CLASSTYPE ## _TransformSliceToRawContigFunc transformSliceToRawContig; \
  CLASSTYPE ## _TransformSliceToSliceFunc transformSliceToSlice;

SEQFEATUREFUNC_TYPES(SeqFeature)

typedef struct SeqFeatureFuncsStruct {
  SEQFEATUREFUNCS_DATA(SeqFeature)
} SeqFeatureFuncs;

#define SEQFEATURE_DATA \
  ENSROOT_DATA \
  int         start; \
  int         end; \
  signed char phase; \
  signed char endPhase; \
  signed char frame; \
  signed char strand; \
  char *      seqName; \
  Storable    st; \
  Analysis *  analysis; \
  double      score; \
  double      pValue; \
  double      percentId; \
  char        isSplittable; \
  BaseContig *contig;

#define FUNCSTRUCTTYPE SeqFeatureFuncs
struct SeqFeatureStruct {
  SEQFEATURE_DATA
};
#undef FUNCSTRUCTTYPE

SeqFeature *SeqFeature_new(void);

#define SeqFeature_setStart(sf,s) ((sf)->funcs->setStart == NULL ? ((sf)->start = (s)) : \
                                                                   ((sf)->funcs->setStart((sf),(s))))
#define SeqFeature_getStart(sf) ((sf)->funcs->getStart == NULL ? ((sf)->start) : \
                                                                 ((sf)->funcs->getStart((sf))))
#define SeqFeature_setEnd(sf,e) ((sf)->funcs->setEnd == NULL ? ((sf)->end = (e)) : \
                                                                   ((sf)->funcs->setEnd((sf),(e))))
#define SeqFeature_getEnd(sf) ((sf)->funcs->getEnd == NULL ? ((sf)->end) : \
                                                                 ((sf)->funcs->getEnd((sf))))
#define SeqFeature_setStrand(sf,s) ((sf)->funcs->setStrand == NULL ? ((sf)->strand = (s)) : \
                                                                   ((sf)->funcs->setStrand((sf),(s))))
#define SeqFeature_getStrand(sf) ((sf)->funcs->getStrand == NULL ? ((sf)->strand) : \
                                                                 ((sf)->funcs->getStrand((sf))))
#define SeqFeature_getLength(sf) ((sf)->funcs->getLength == NULL ? (SeqFeature_getEnd((sf)) - SeqFeature_getStart((sf)) + 1) : \
                                                                 ((sf)->funcs->getLength((sf))))

#define SeqFeature_setScore(sf,s) (sf)->score = (s)
#define SeqFeature_getScore(sf) (sf)->score

#define SeqFeature_setpValue(sf,s) (sf)->pValue = (s)
#define SeqFeature_getpValue(sf) (sf)->pValue

#define SeqFeature_setPercId(sf,s) (sf)->percentId = (s)
#define SeqFeature_getPercId(sf) (sf)->percentId

#define SeqFeature_setPhase(sf,p) (sf)->phase = (p)
#define SeqFeature_getPhase(sf) (sf)->phase

#define SeqFeature_setEndPhase(sf,ep) (sf)->endPhase = (ep)
#define SeqFeature_getEndPhase(sf) (sf)->endPhase

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

#define SeqFeature_getIsSplittable(sf) (sf)->isSplittable


int SeqFeature_startCompFunc(const void *a, const void *b);
int SeqFeature_reverseStartCompFunc(const void *a, const void *b);


Vector *SeqFeature_transformToRawContig(SeqFeature *sf);
Vector *SeqFeature_transformSliceToRawContig(SeqFeature *sf);

SeqFeature *SeqFeature_transformToSlice(SeqFeature *sf, Slice *slice);
SeqFeature *SeqFeature_transformRawContigToSlice(SeqFeature *sf, Slice *slice);

#ifdef __SEQFEATURE_MAIN__
 SeqFeatureFuncs 
   seqFeatureFuncs = {
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
                      SeqFeature_transformToRawContig,
                      SeqFeature_transformToSlice,
                      SeqFeature_transformRawContigToSlice,
                      SeqFeature_transformSliceToRawContig,
                      NULL // transformSliceToSlice
                     };
#else
 extern SeqFeatureFuncs seqFeatureFuncs;
#endif

#endif
