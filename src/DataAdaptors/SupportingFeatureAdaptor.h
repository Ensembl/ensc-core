#ifndef __SUPPORTINGFEATUREADAPTOR_H__
#define __SUPPORTINGFEATUREADAPTOR_H__

#include "BaseAdaptor.h"
#include "AdaptorTypes.h"
#include "Exon.h"

struct SupportingFeatureAdaptorStruct {
  BASEADAPTOR_DATA
};

SupportingFeatureAdaptor *SupportingFeatureAdaptor_new(DBAdaptor *dba);
Vector *SupportingFeatureAdaptor_fetchAllByExon(SupportingFeatureAdaptor *sfa, Exon *exon);


#endif
