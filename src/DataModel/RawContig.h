#ifndef __RAWCONTIG_H__
#define __RAWCONTIG_H__

#include "Storable.h"
#include "BaseContig.h"

typedef struct RawContigStruct RawContig;

struct RawContigStruct {
  BASECONTIG_DATA
  int length;
  int emblOffset;
  long cloneId;
  char *name;  
  Storable st;
};

RawContig *RawContig_new(void);

#define RawContig_setDbID(rc,id) Storable_setDbID(&((rc)->st),(id))
#define RawContig_getDbID(rc) Storable_getDbID(&((rc)->st))

#define RawContig_setAdaptor(rc,ad) Storable_setAdaptor(&((rc)->st),ad)
#define RawContig_getAdaptor(rc) Storable_getAdaptor(&((rc)->st))

char *RawContig_setName(RawContig *rc, char *name);
char *RawContig_getName(RawContig *rc);

#define RawContig_setCloneID(rc,cid) (rc)->cloneId = (cid)
long RawContig_getCloneID(RawContig *rc);

#define RawContig_setEMBLOffset(rc,eo) (rc)->emblOffset = (eo)
int RawContig_getEMBLOffset(RawContig *rc);

#define RawContig_setLength(rc,l) (rc)->length = (l)
int RawContig_getLength(RawContig *rc);


#endif
