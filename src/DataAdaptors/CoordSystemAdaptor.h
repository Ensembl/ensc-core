#ifndef __COORDSYSTEMADAPTOR_H__
#define __COORDSYSTEMADAPTOR_H__

#include "BaseAdaptor.h"
#include "AdaptorTypes.h"
#include "CoordSystem.h"
#include "IDHash.h"
#include "StringHash.h"
#include "Vector.h"

struct CoordSystemAdaptorStruct {
  BASEADAPTOR_DATA
  IDHash *dbIDCache;
  IDHash *rankCache;
  StringHash *nameCache;

  IDHash *isSeqLevelCache;
  IDHash *isDefaultVersionCache;

  StringHash *mappingPaths;

  CoordSystem *topLevel;
};


int CoordSystem_sortByRankFunc(const void *one, const void *two);

void                CoordSystemAdaptor_cacheMappingPaths(CoordSystemAdaptor *csa);
void                CoordSystemAdaptor_cacheSeqRegionMapping(CoordSystemAdaptor *csa);
Vector *            CoordSystemAdaptor_fetchAllByAttrib(CoordSystemAdaptor *csa, char *attrib);
Vector *            CoordSystemAdaptor_fetchAllByName(CoordSystemAdaptor *csa, char *name);
CoordSystem *       CoordSystemAdaptor_fetchByAttrib(CoordSystemAdaptor *csa, char *attrib, char *version);
CoordSystem *       CoordSystemAdaptor_fetchByDbID(CoordSystemAdaptor *csa, IDType dbID);
CoordSystem *       CoordSystemAdaptor_fetchByName(CoordSystemAdaptor *csa, char *name, char *version);
CoordSystem *       CoordSystemAdaptor_fetchByRank(CoordSystemAdaptor *csa, int rank);
Vector *            CoordSystemAdaptor_fetchAll(CoordSystemAdaptor *csa);
CoordSystem *       CoordSystemAdaptor_fetchSeqLevel(CoordSystemAdaptor *csa);
CoordSystem *       CoordSystemAdaptor_fetchTopLevel(CoordSystemAdaptor *csa);
Vector *            CoordSystemAdaptor_getMappingPath(CoordSystemAdaptor *csa, CoordSystem *cs1, CoordSystem *cs2);
CoordSystemAdaptor *CoordSystemAdaptor_new(DBAdaptor *dba);

#endif
