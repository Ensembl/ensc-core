#ifndef __RAWCONTIGADAPTOR_H__
#define __RAWCONTIGADAPTOR_H__

#include "BaseAdaptor.h"
#include "AdaptorTypes.h"
#include "RawContig.h"
#include "IDHash.h"

struct RawContigAdaptorStruct {
  BASEADAPTOR_DATA
  IDHash *rawContigCache;
};

RawContigAdaptor *RawContigAdaptor_new(DBAdaptor *dba);
RawContig *RawContigAdaptor_fetchByDbID(RawContigAdaptor *rca, IDType dbID);
RawContig *RawContigAdaptor_rawContigFromRow(RawContigAdaptor *rca, ResultRow *row);
void RawContigAdaptor_fillRawContigWithRow(RawContigAdaptor *rca, RawContig *rc, ResultRow *row);
void RawContigAdaptor_fetchAttributes(RawContigAdaptor *rca, RawContig *rc);
RawContig *RawContigAdaptor_fetchByName(RawContigAdaptor *rca, char *name);

#endif
