#ifndef __GENEADAPTOR_H__
#define __GENEADAPTOR_H__

#include "BaseAdaptor.h"
#include "AdaptorTypes.h"
#include "Gene.h"

struct GeneAdaptorStruct {
  BASEADAPTOR_DATA
};

GeneAdaptor *GeneAdaptor_new(DBAdaptor *dba);
int GeneAdaptor_listGeneIds(GeneAdaptor *ga, long **geneIds);
Gene *GeneAdaptor_fetchByDbID(GeneAdaptor *ga, long geneId, int chrCoords);
int GeneAdaptor_getStableEntryInfo(GeneAdaptor *ga, Gene *gene);



#endif
