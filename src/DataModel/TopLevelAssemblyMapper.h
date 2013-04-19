#ifndef __TOPLEVELASSEMBLYMAPPER_H__
#define __TOPLEVELASSEMBLYMAPPER_H__

#include "DataModelTypes.h"
#include "EnsC.h"

#include "Vector.h"
#include "IDHash.h"
#include "Mapper.h"
#include "AdaptorTypes.h"
#include "MapperRangeSet.h"

#define MAXCHUNK 2000

struct TopLevelAssemblyMapperStruct {
  AssemblyMapperAdaptor *adaptor;
  Vector *coordSystems;
  CoordSystem *topLevelCs;
  CoordSystem *otherCs;
};

#define TopLevelAssemblyMapper_setAdaptor(tlam, ad) AssemblyMapper_setAdaptor(tlam,ad)
#define TopLevelAssemblyMapper_getAdaptor(tlam) AssemblyMapper_getAdaptor(tlam)


#endif
