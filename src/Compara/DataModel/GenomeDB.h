#ifndef __GENOMEDB_H__
#define __GENOMEDB_H__

#include "ComparaDataModelTypes.h"
#include "EnsRoot.h"

struct GenomeDBStruct {
  ENSROOT_DATA
  Storable st;
  char *name;
  DBAdaptor *dbAdaptor;
  char *assembly;
  char *taxonId;
};

#endif
