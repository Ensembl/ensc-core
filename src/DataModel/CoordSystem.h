#ifndef __COORDSYSTEM_H__
#define __COORDSYSTEM_H__

#include "DataModelTypes.h"
#include "Storable.h"

#include "EnsRoot.h"

OBJECTFUNC_TYPES(CoordSystem)

typedef struct CoordSystemFuncsStruct {
  OBJECTFUNCS_DATA(CoordSystem)
} CoordSystemFuncs;

#define FUNCSTRUCTTYPE CoordSystemFuncs
struct CoordSystemStruct {
  ENSROOT_DATA
  Storable st;
  char *name;
  char *version;
  int   rank;

  int  isTopLevel;
  int  isSeqLevel;
  int  isDefaultVersion;
};
#undef FUNCSTRUCTTYPE

CoordSystem *CoordSystem_new(char *name, char *version, int rank, IDType dbID, CoordSystemAdaptor *csa, int isDefault, int isSeqLevel, int isTopLevel);

#define CoordSystem_setDbID(cs,dbID) Storable_setDbID(&((cs)->st),dbID)
#define CoordSystem_getDbID(cs) Storable_getDbID(&((cs)->st))

#define CoordSystem_setAdaptor(cs,ad) Storable_setAdaptor(&((cs)->st),ad)
#define CoordSystem_getAdaptor(cs) Storable_getAdaptor(&((cs)->st))

#define CoordSystem_getName(cs) cs->name

#define CoordSystem_getRank(cs) cs->rank

#define CoordSystem_getVersion(cs) cs->version

#define CoordSystem_getIsTopLevel(cs) cs->isTopLevel

#define CoordSystem_getIsSeqLevel(cs) cs->isSeqLevel

#define CoordSystem_getIsDefaultVersion(cs) cs->isDefaultVersion

int CoordSystem_compare(CoordSystem *cs1, CoordSystem *cs2);

void CoordSystem_free(CoordSystem *coordSystem);

#ifdef __COORDSYSTEM_MAIN__
  CoordSystemFuncs
    coordSystemFuncs = {
                    CoordSystem_free,
                   };
#else
  extern CoordSystemFuncs coordSystemFuncs;
#endif


#endif
