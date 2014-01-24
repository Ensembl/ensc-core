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

#ifndef __SEQFEATURE_H__
#define __SEQFEATURE_H__


#include "DataModelTypes.h"
#include "EnsC.h"
#include "Storable.h"
#include "Analysis.h"
#include "BaseContig.h"
#include "Vector.h"

#include "EnsRoot.h"

//typedef void (*CLASSTYPE ## _FreeFunc)(CLASSTYPE *); \

#define SEQFEATUREFUNC_TYPES(CLASSTYPE) \
OBJECTFUNC_TYPES(CLASSTYPE) \
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

//  CLASSTYPE ## _FreeFunc free; \

#define SEQFEATUREFUNCS_DATA(CLASSTYPE) \
  OBJECTFUNCS_DATA(CLASSTYPE) \
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
  long         start; \
  long         end; \
  signed char  phase; \
  signed char  endPhase; \
  signed char  frame; \
  signed char  strand; \
  char         isSplittable; \
  ECOSTRING    seqName; \
  Storable     st; \
  Analysis *   analysis; \
  double       pValue; \
  float        score; \
  float        percentId; \
  unsigned int flags; \
  void *       extraData; \
  BaseContig * contig;

#define FUNCSTRUCTTYPE SeqFeatureFuncs
struct SeqFeatureStruct {
  SEQFEATURE_DATA
};
#undef FUNCSTRUCTTYPE

SeqFeature *SeqFeature_new(void);

#define SeqFeature_addFlag(sf, f) (sf)->flags |= (f)
#define SeqFeature_getFlags(sf) (sf)->flags
#define SeqFeature_removeFlag(sf, f) (sf)->flags &= ~((f))

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

#define SeqFeature_setAnalysis(sf,ana) (ana) ? (sf)->analysis = (ana), Object_incRefCount((ana)) : 0//NULL
#define SeqFeature_getAnalysis(sf) (sf)->analysis

char *SeqFeature_setStableId(SeqFeature *sf, char *stableId);
#define SeqFeature_getStableId(sf) (sf)->stableId

#define SeqFeature_setDbID(sf,dbID) Storable_setDbID(&((sf)->st),dbID)
#define SeqFeature_getDbID(sf) Storable_getDbID(&((sf)->st))

#define SeqFeature_setAdaptor(sf,ad) Storable_setAdaptor(&((sf)->st),ad)
#define SeqFeature_getAdaptor(sf) Storable_getAdaptor(&((sf)->st))

#define SeqFeature_setContig(sf,c) (sf)->contig = (BaseContig *)(c)
#define SeqFeature_getContig(sf) (sf)->contig

#define SeqFeature_setSlice(sf,c) (sf)->contig = (BaseContig *)(c)
#define SeqFeature_getSlice(sf) ((Slice *)((sf)->contig))

#define SeqFeature_getIsSplittable(sf) (sf)->isSplittable

SeqFeature *SeqFeature_transfer(SeqFeature *sf, Slice *slice);
Vector *SeqFeature_projectToSlice(SeqFeature *sf, Slice *toSlice);
Vector *SeqFeature_project(SeqFeature *sf, char *csName, char *csVersion);
SeqFeature *SeqFeature_transform(SeqFeature *sf, char *csName, char *csVersion, Slice *toSlice);

int SeqFeature_startCompFunc(const void *a, const void *b);
int SeqFeature_startEndCompFunc(const void *a, const void *b);
int SeqFeature_reverseStartCompFunc(const void *a, const void *b);
int SeqFeature_reverseScoreCompFunc(const void *a, const void *b);
int SeqFeature_startRevEndCompFunc(const void *a, const void *b);


Vector *SeqFeature_transformToRawContigImpl(SeqFeature *sf);
Vector *SeqFeature_transformSliceToRawContigImpl(SeqFeature *sf);

SeqFeature *SeqFeature_transformToSliceImpl(SeqFeature *sf, Slice *slice);
SeqFeature *SeqFeature_transformRawContigToSliceImpl(SeqFeature *sf, Slice *slice);

ECOSTRING SeqFeature_setSeqName(SeqFeature *sf, char *seqName);
ECOSTRING SeqFeature_getSeqName(SeqFeature *sf);

ECOSTRING SeqFeature_getSeqRegionName(SeqFeature *sf);
long SeqFeature_getSeqRegionStart(SeqFeature *sf);
long SeqFeature_getSeqRegionEnd(SeqFeature *sf);
int SeqFeature_getSeqRegionStrand(SeqFeature *sf);

int SeqFeature_overlaps(SeqFeature *sf, SeqFeature *f);


#define SeqFeature_transformToSlice(sf,slice) \
      ((sf)->funcs->transformToSlice == NULL ? \
         (fprintf(stderr,"Error: Null pointer for transformToSlice - bye\n"),  exit(1), (void *)NULL) : \
         ((sf)->funcs->transformToSlice((sf),(slice))))

#define SeqFeature_transformToRawContig(sf) \
      ((sf)->funcs->transformToRawContig == NULL ? \
         (fprintf(stderr,"Error: Null pointer for transformToRawContig - bye\n"),  exit(1), (void *)NULL) : \
         ((sf)->funcs->transformToRawContig((sf))))

#define SeqFeature_transformSliceToRawContig(sf) \
      ((sf)->funcs->transformSliceToRawContig == NULL ? \
         (fprintf(stderr,"Error: Null pointer for transformSliceToRawContig - bye\n"),  exit(1), (Vector *)NULL) : \
         ((sf)->funcs->transformSliceToRawContig((sf))))

#define SeqFeature_transformRawContigToSlice(sf,slice) \
      ((sf)->funcs->transformRawContigToSlice == NULL ? \
         (fprintf(stderr,"Error: Null pointer for transformRawContigToSlice - bye\n"),  exit(1), (void *)NULL) : \
         ((sf)->funcs->transformRawContigToSlice((sf),(slice))))

#define SeqFeature_free(sf) EnsRoot_free((sf))

#define SeqFeature_shallowCopy(sf) EnsRoot_shallowCopy((sf))
#define SeqFeature_deepCopy(sf) EnsRoot_deepCopy((sf))

#ifdef __SEQFEATURE_MAIN__
 SeqFeatureFuncs 
   seqFeatureFuncs = {
                      NULL, // free
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
                      SeqFeature_transformToRawContigImpl,
                      SeqFeature_transformToSliceImpl,
                      SeqFeature_transformRawContigToSliceImpl,
                      SeqFeature_transformSliceToRawContigImpl,
                      NULL // transformSliceToSlice
                     };
#else
 extern SeqFeatureFuncs seqFeatureFuncs;
#endif

#endif
