#ifndef __SLICE_H__
#define __SLICE_H__

#include "Storable.h"
#include "BaseContig.h"
#include "Gene.h"

typedef struct SliceStruct Slice;

struct SliceStruct {
  BASECONTIG_DATA
  int strand;
  Storable st;
  char emptyFlag;
  char *chrName;
  char *assemblyType;
  long chrId;
};

Slice *Slice_new(char *chr, int start, int end, int strand, char *assemblyType,
                 SliceAdaptor *sa, long dbID, int empty);

#define Slice_setDbID(s,dbID) Storable_setDbID(&((s)->st),dbID)
#define Slice_getDbID(s) Storable_getDbID(&((s)->st))

#define Slice_setAdaptor(s,ad) Storable_setAdaptor(&((s)->st),ad)
#define Slice_getAdaptor(s) Storable_getAdaptor(&((s)->st))

#define Slice_setEmptyFlag(s,e) (s)->emptyFlag = (e)
#define Slice_getEmptyFlag(s) (s)->emptyFlag;

#define Slice_setChrStart(sl,s) (sl)->start = (s)
#define Slice_getChrStart(sl) (sl)->start

#define Slice_setChrEnd(sl,e) (sl)->end = (e)
#define Slice_getChrEnd(sl) (sl)->end

#define Slice_setChrId(sl,c) (sl)->chrId = (c)
#define Slice_getChrId(sl) (sl)->chrId

#define Slice_setStrand(sl,s) (sl)->strand = (s)
#define Slice_getStrand(sl) (sl)->strand

char *Slice_setAssemblyType(Slice *sl,char *type);
#define Slice_getAssemblyType(sl) (sl)->assemblyType

char *Slice_setChrName(Slice *sl,char *chrName);
#define Slice_getChrName(sl) (sl)->chrName

char *Slice_getName(Slice *sl, char *retStr);

Gene **Slice_getAllGenes(Slice *slice, char *logicName);


#endif
