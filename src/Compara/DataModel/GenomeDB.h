#ifndef __GENOMEDB_H__
#define __GENOMEDB_H__

#include "ComparaDataModelTypes.h"
#include "EnsRoot.h"
#include "Storable.h"
#include "Vector.h"

OBJECTFUNC_TYPES(GenomeDB)

typedef struct GenomeDBFuncsStruct {
  OBJECTFUNCS_DATA(GenomeDB)
} GenomeDBFuncs;

#define FUNCSTRUCTTYPE GenomeDBFuncs
struct GenomeDBStruct {
  ENSROOT_DATA
  Storable st;
  char *name;
  DBAdaptor *dbAdaptor;
  char *assembly;
  IDType taxonId;
};
#undef FUNCSTRUCTTYPE

#define GenomeDB_setDbID(gd,id) Storable_setDbID(&((gd)->st),(id))
#define GenomeDB_getDbID(gd) Storable_getDbID(&((gd)->st))

#define GenomeDB_setAdaptor(gd,a) Storable_setAdaptor(&((gd)->st),(a))
#define GenomeDB_getAdaptor(gd) Storable_getAdaptor(&((gd)->st))

#define GenomeDB_setDBAdaptor(gd,a) (gd)->dbAdaptor = (a)
#define GenomeDB_getDBAdaptor(gd) (gd)->dbAdaptor

char *GenomeDB_setAssembly(GenomeDB *gd, char *assembly);
#define GenomeDB_getAssembly(gd) (gd)->assembly

char *GenomeDB_setName(GenomeDB *gd, char *name);
#define GenomeDB_getName(gd) (gd)->name

#define GenomeDB_setTaxonId(gd,t) (gd)->taxonId = (t)
#define GenomeDB_getTaxonId(gd) (gd)->taxonId

GenomeDB *GenomeDB_new();

Vector *GenomeDB_linkedGenomesByMethodLinkId(GenomeDB *gdb, IDType methodLinkId);

void GenomeDB_free(GenomeDB *gd);

#ifdef __GENOMEDB_MAIN__
  GenomeDBFuncs
    genomeDBFuncs = {
                     GenomeDB_free
                    };
#else
  extern GenomeDBFuncs enomeDBFuncs;
#endif



#endif
