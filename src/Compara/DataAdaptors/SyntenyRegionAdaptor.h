#ifndef __SYNTENYREGIONADAPTOR_H__
#define __SYNTENYREGIONADAPTOR_H__

#include "BaseComparaAdaptor.h"
#include "ComparaAdaptorTypes.h"
#include "SyntenyRegion.h"
#include "Vector.h"


struct SyntenyRegionAdaptorStruct {
  BASECOMPARAADAPTOR_DATA
};

SyntenyRegionAdaptor *SyntenyRegionAdaptor_new(ComparaDBAdaptor *dba);
SyntenyRegion *SyntenyRegionAdaptor_newRegionFromArray(SyntenyRegionAdaptor *sra, IDType dbID,
                                            IDType cluster, IDType dnaFrag,int start,int end);
SyntenyRegion *SyntenyRegionAdaptor_fetchByDbID(SyntenyRegionAdaptor *sra, IDType dbID);

Vector *SyntenyRegionAdaptor_fetchByClusterId(SyntenyRegionAdaptor *sra, IDType clusterId);
IDType SyntenyRegionAdaptor_store(SyntenyRegionAdaptor *sra, IDType clusterId, SyntenyRegion *region);



#endif
