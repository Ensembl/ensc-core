#ifndef __METACONTAINER_H__
#define __METACONTAINER_H__

#include "BaseAdaptor.h"
#include "AdaptorTypes.h"
#include "Species.h"


struct MetaContainerStruct {
  BASEADAPTOR_DATA
};

MetaContainer *MetaContainer_new(DBAdaptor *dba);
char *MetaContainer_getDefaultAssembly(MetaContainer *mc);
Species *MetaContainer_getSpecies(MetaContainer *mc);


#endif
