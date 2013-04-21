#ifndef __SLICE_H__
#define __SLICE_H__

#include "DataModelTypes.h"
#include "Storable.h"
#include "BaseContig.h"
#include "Gene.h"
#include "Vector.h"
#include "EcoString.h"
#include "CoordSystem.h"

BASECONTIGFUNC_TYPES(Slice)

typedef struct SliceFuncsStruct {
  BASECONTIGFUNCS_DATA(Slice)
} SliceFuncs;

#define FUNCSTRUCTTYPE SliceFuncs
struct SliceStruct {
  BASECONTIG_DATA
  int strand;
  char emptyFlag;
  ECOSTRING seqRegionName;
  ECOSTRING assemblyType;
  CoordSystem *coordSystem;
  IDType seqRegionId;
};
#undef FUNCSTRUCTTYPE

#ifdef __SLICE_MAIN__
// NIY: Need to update for new slice
 Slice emptySliceData = {CLASS_SLICE,0,NULL,NULL,NULL,0,{-1,NULL},-1,-1,0,1,NULL,NULL,-1};
 Slice *emptySlice = &emptySliceData;
#else
 extern Slice *emptySlice;
#endif

Slice *Slice_new(char *chr, int start, int end, int strand, char *assemblyType,
                 SliceAdaptor *sa, IDType dbID, int empty);

#define Slice_setDbID(s,dbID) BaseContig_setDbID((s),(dbID))
#define Slice_getDbID(s) BaseContig_getDbID((s))

#define Slice_setAdaptor(s,ad) BaseContig_setAdaptor((s),(ad))
#define Slice_getAdaptor(s) BaseContig_getAdaptor((s))

#define Slice_setEmptyFlag(s,e) (s)->emptyFlag = (e)
#define Slice_getEmptyFlag(s) (s)->emptyFlag

#define Slice_setSeqRegionStart(sl,s) (sl)->start = (s)
#define Slice_getSeqRegionStart(sl) (sl)->start

#define Slice_setSeqRegionEnd(sl,e) (sl)->end = (e)
#define Slice_getSeqRegionEnd(sl) (sl)->end

#define Slice_setSeqRegionId(sl,c) (sl)->seqRegionId = (c)
#define Slice_getSeqRegionId(sl) (sl)->seqRegionId

// SMJS Temporary Chr versions to satisfy linking
#define Slice_setChrStart(sl,s) (sl)->start = (s)
#define Slice_getChrStart(sl) (sl)->start

#define Slice_setChrEnd(sl,e) (sl)->end = (e)
#define Slice_getChrEnd(sl) (sl)->end

#define Slice_setChrId(sl,c) (sl)->seqRegionId = (c)
#define Slice_getChrId(sl) (sl)->seqRegionId

#define Slice_getChrName(sl) (sl)->seqRegionName
// SMJS End temporary Chr versions to satisfy linking

#define Slice_setCoordSystem(sl,cs) (sl)->coordSystem = (cs)
#define Slice_getCoordSystem(sl) (sl)->coordSystem

#define Slice_setStrand(sl,s) (sl)->strand = (s)
#define Slice_getStrand(sl) (sl)->strand

#define Slice_getLength(sl) ((sl)->end - (sl)->start + 1)

ECOSTRING Slice_setAssemblyType(Slice *sl,char *type);
#define Slice_getAssemblyType(sl) (sl)->assemblyType

ECOSTRING Slice_setSeqRegionName(Slice *sl,char *seqRegionName);
#define Slice_getSeqRegionName(sl) (sl)->seqRegionName

ECOSTRING Slice_getName(Slice *sl);
Vector *Slice_getAllGenes(Slice *slice, char *logicName);
Vector *Slice_getAllSimpleFeatures(Slice *slice, char *logicName, double *score);
Vector *Slice_getAllDNAAlignFeatures(Slice *slice, char *logicName, double *score);
Vector *Slice_getAllDNAPepAlignFeatures(Slice *slice, char *logicName, double *score);
Vector *Slice_getAllPredictionTranscripts(Slice *slice, char *logicName);
Vector *Slice_getAllRepeatFeatures(Slice *slice, char *logicName);

Vector *Slice_getAllGenesByType(Slice *slice, char *type);

char *Slice_getSubSeq(Slice *slice, int start, int end, int strand);
char *Slice_getSeq(Slice *slice);


#ifdef __SLICE_MAIN__
  SliceFuncs sliceFuncs = {
                           NULL, // free
                           Slice_getName, 
                           Slice_getSeq,
                           Slice_getSubSeq
                          };
#else
  extern SliceFuncs sliceFuncs;
#endif


#endif
