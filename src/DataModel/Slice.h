#ifndef __SLICE_H__
#define __SLICE_H__

#include "Storable.h"

typedef struct SliceStruct Slice;

struct SliceStruct {
  Storable st;
  char emptyFlag;
};

Slice *Slice_new(void);

#define Slice_setDbID(s,dbID) Storable_setDbID(&((s)->st),dbID)
#define Slice_getDbID(s) Storable_getDbID(&((s)->st))

#define Slice_setAdaptor(s,ad) Storable_setAdaptor(&((s)->st),ad)
#define Slice_getAdaptor(s) Storable_getAdaptor(&((s)->st))

#define Slice_setEmptyFlag(s,e) (s)->emptyFlag = (e)
#define Slice_getEmptyFlag(s) (s)->emptyFlag;

#endif
