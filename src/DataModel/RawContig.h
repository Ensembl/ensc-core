#ifndef __RAWCONTIG_H__
#define __RAWCONTIG_H__

#include "DataModelTypes.h"
#include "Storable.h"
#include "BaseContig.h"
#include "Vector.h"
#include "SequenceAdaptor.h"

#define FUNCSTRUCTTYPE BaseContigFuncs
struct RawContigStruct {
  BASECONTIG_DATA
  int emblOffset;
  IDType cloneId;
};
#undef FUNCSTRUCTTYPE

RawContig *RawContig_new(void);

#define RawContig_setDbID(rc,id) BaseContig_setDbID((rc),(id))
#define RawContig_getDbID(rc) BaseContig_getDbID((rc))

#define RawContig_setAdaptor(rc,ad) BaseContig_setAdaptor((rc),(ad))
#define RawContig_getAdaptor(rc) BaseContig_getAdaptor((rc))

char *RawContig_setName(RawContig *rc, char *name);
char *RawContig_getName(RawContig *rc);

#define RawContig_setCloneID(rc,cid) (rc)->cloneId = (cid)
long RawContig_getCloneID(RawContig *rc);

#define RawContig_setEMBLOffset(rc,eo) (rc)->emblOffset = (eo)
int RawContig_getEMBLOffset(RawContig *rc);

#define RawContig_setLength(rc,l) (rc)->length = (l)
int RawContig_getLength(RawContig *rc);

Vector *RawContig_getAllSimpleFeatures(RawContig *rc, char *logicName, double *scoreP);
Vector *RawContig_getAllPredictionTranscripts(RawContig *rc, char *logicName);
Vector *RawContig_getAllRepeatFeatures(RawContig *rc, char *logicName);
Vector *RawContig_getAllDNAAlignFeatures(RawContig *rc, char *logicName, double *scoreP);
Vector *RawContig_getAllProteinAlignFeatures(RawContig *rc, char *logicName, double *scoreP);


#endif
