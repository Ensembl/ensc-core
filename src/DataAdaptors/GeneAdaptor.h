#ifndef __GENEADAPTOR_H__
#define __GENEADAPTOR_H__

#include "BaseAdaptor.h"
#include "AdaptorTypes.h"
#include "Gene.h"
#include "Slice.h"

#include "StringHash.h"

struct GeneAdaptorStruct {
  BASEADAPTOR_DATA
  StringHash *sliceGeneCache;
};

GeneAdaptor *GeneAdaptor_new(DBAdaptor *dba);
int GeneAdaptor_listGeneIds(GeneAdaptor *ga, long **geneIds);
Gene *GeneAdaptor_fetchByDbID(GeneAdaptor *ga, long geneId, int chrCoords);
int GeneAdaptor_getStableEntryInfo(GeneAdaptor *ga, Gene *gene);
Gene **GeneAdaptor_fetchAllBySlice(GeneAdaptor *ga, Slice *slice, char *logicName);




#endif
