#ifndef __DNAALIGNFEATUREADAPTOR_H__
#define __DNAALIGNFEATUREADAPTOR_H__

#include "BaseFeatureAdaptor.h"
#include "AdaptorTypes.h"

struct DNAAlignFeatureAdaptorStruct {
  BASEFEATUREADAPTOR_DATA
};

DNAAlignFeatureAdaptor *DNAAlignFeatureAdaptor_new(DBAdaptor *dba);

#endif
