#ifndef __PROTEINALIGNFEATUREADAPTOR_H__
#define __PROTEINALIGNFEATUREADAPTOR_H__

#include "BaseFeatureAdaptor.h"
#include "AdaptorTypes.h"

struct ProteinAlignFeatureAdaptorStruct {
  BASEFEATUREADAPTOR_DATA
};

ProteinAlignFeatureAdaptor *ProteinAlignFeatureAdaptor_new(DBAdaptor *dba);

#endif
