/*
 * Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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
  ECOSTRING seqRegionName;
  long seqRegionLength;
  IDType seqRegionId;
  CoordSystem *coordSystem;
  long codonTable;

  char isReference;
  char isTopLevel;
  char hasKaryotype;
};
#undef FUNCSTRUCTTYPE

//#ifdef __SLICE_MAIN__
// NIY: Need to update for new slice
// Slice emptySliceData = {CLASS_SLICE,0,NULL,NULL,NULL,0,{-1,NULL},-1,-1,0,1,NULL,NULL,NULL,-1};
// Slice *emptySlice = &emptySliceData;
//#else
// extern Slice *emptySlice;
//#endif

Slice *Slice_new(char *regionName, long start, long end, int strand, long length, CoordSystem *cs, SliceAdaptor *sa);

long Slice_getCodonTableId(Slice *slice);

#define Slice_setStart(sl,s) (sl)->start = (s)
#define Slice_getStart(sl) (sl)->start
#define Slice_setEnd(sl,e) (sl)->end = (e)
#define Slice_getEnd(sl) (sl)->end

#define Slice_setDbID(s,dbID) BaseContig_setDbID((s),(dbID))
#define Slice_getDbID(s) BaseContig_getDbID((s))

#define Slice_setAdaptor(s,ad) BaseContig_setAdaptor((s),(ad))
#define Slice_getAdaptor(s) BaseContig_getAdaptor((s))

#define Slice_setSeqRegionStart(sl,s) (sl)->start = (s)
#define Slice_getSeqRegionStart(sl) (sl)->start

#define Slice_setSeqRegionEnd(sl,e) (sl)->end = (e)
#define Slice_getSeqRegionEnd(sl) (sl)->end

#define Slice_setSeqRegionLength(sl,len) (sl)->seqRegionLength = (len)
#define Slice_getSeqRegionLength(sl) (sl)->seqRegionLength

Slice *Slice_getSeqRegionSlice(Slice *slice);

//#define Slice_setSeqRegionId(sl,c) (sl)->seqRegionId = (c)
//#define Slice_getSeqRegionId(sl) (sl)->seqRegionId

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

//ECOSTRING Slice_setAssemblyType(Slice *sl,char *type);
//#define Slice_getAssemblyType(sl) (sl)->assemblyType

ECOSTRING Slice_setSeqRegionName(Slice *sl,char *seqRegionName);
#define Slice_getSeqRegionName(sl) (sl)->seqRegionName

ECOSTRING Slice_getName(Slice *sl);
Vector *Slice_constrainToRegion(Slice *slice);
Slice *Slice_expand(Slice *slice, long fivePrimeShift, long threePrimeShift, int forceExpand, long *fpRef, long *tpRef);
Vector *Slice_getAllAttributes(Slice *slice, char *attribCode);
Vector *Slice_getAllPredictionTranscripts(Slice *slice, char *logicName, int loadExons, char *dbType);
Vector *Slice_getAllDNAAlignFeatures(Slice *slice, char *logicName, double *score, char *dbType, double *hCoverage);
Vector *Slice_getAllProteinAlignFeatures(Slice *slice, char *logicName, double *score, char *dbType, double *hCoverage);
Vector *Slice_getAllSimilarityFeatures(Slice *slice, char *logicName, double *score);
Vector *Slice_getAllSimpleFeatures(Slice *slice, char *logicName, double *score, char *dbType);
Vector *Slice_getAllRepeatFeatures(Slice *slice, char *logicName, Vector *repeatTypes, char *dbType);
Vector *Slice_getAllGenes(Slice *slice, char *logicName, char *dbType, int loadTranscripts, char *source, char *bioType);
Vector *Slice_getAllGenesByType(Slice *slice, char *type, char *logicName, int loadTranscripts);
Vector *Slice_getAllGenesBySource(Slice *slice, char *source, int loadTranscripts);
Vector *Slice_getAllTranscripts(Slice *slice, int loadExons, char *logicName, char *dbType);
Vector *Slice_getAllExons(Slice *slice, char *dbType);

IDType Slice_getSeqRegionId(Slice *slice);

DBAdaptor *Slice_getSelectedDBAdaptor(Slice *slice, char *dbType);

int Slice_isTopLevel(Slice *slice);




char *Slice_getSubSeq(Slice *slice, int start, int end, int strand);
char *Slice_getSeq(Slice *slice);
Vector *Slice_project(Slice *slice, char *csName, char *csVersion);
Vector *Slice_projectToSlice(Slice *slice, Slice *toSlice);
char *Slice_getCoordSystemName(Slice *slice);
Slice *Slice_expand(Slice *slice, long fivePrimeShift, long threePrimeShift, int forceExpand, long *fpRef, long *tpRef);
Slice *Slice_invert(Slice *slice);



#ifdef __SLICE_MAIN__
  SliceFuncs sliceFuncs = {
                           NULL, // free
                           NULL, // shallowCopy
                           NULL, // deepCopy
                           Slice_getName, 
                           Slice_getSeq,
                           Slice_getSubSeq
                          };
#else
  extern SliceFuncs sliceFuncs;
#endif


#endif
