#define __GENOMEDB_MAIN__
#include "GenomeDB.h"
#undef __GENOMEDB_MAIN__
#include "GenomeDBAdaptor.h"

GenomeDB *GenomeDB_new() {
  GenomeDB *gdb;

  if ((gdb = (GenomeDB *)calloc(1,sizeof(GenomeDB))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for gdb\n");
    return NULL;
  }

  gdb->objectType = CLASS_GENOMEDB;

  gdb->funcs = &genomeDBFuncs;

  Object_incRefCount(gdb);
  return gdb;
}

int GenomeDB_hasConsensus(GenomeDB *gdb, GenomeDB *conGdb, IDType methodLinkId) {
  GenomeDBAdaptor *gda = (GenomeDBAdaptor *)GenomeDB_getAdaptor(gdb);

  // sanity check on the GenomeDB passed in
  if (!conGdb) {
    fprintf(stderr,"Error: No query genome specified\n");
    exit(1);
  }

  // and check that you are not trying to compare the same GenomeDB
  if (conGdb == gdb) {
    fprintf(stderr,"Error: Trying to return consensus/query information from the same db\n");
    exit(1);
  }

  return GenomeDBAdaptor_checkForConsensusDb(gda, gdb, conGdb, methodLinkId);
}

int GenomeDB_hasQuery(GenomeDB *gdb, GenomeDB *queryGdb, IDType methodLinkId) {
  GenomeDBAdaptor *gda = (GenomeDBAdaptor *)GenomeDB_getAdaptor(gdb);

  // sanity check on the GenomeDB passed in
  if (!queryGdb) {
    fprintf(stderr,"Error: No query genome specified\n");
    exit(1);
  }

  // and check that you are not trying to compare the same GenomeDB
  if (queryGdb == gdb) {
    fprintf(stderr,"Error: Trying to return consensus/query information from the same db\n");
    exit(1);
  }

  return GenomeDBAdaptor_checkForQueryDb(gda, gdb, queryGdb, methodLinkId);
}

Vector *GenomeDB_linkedGenomesByMethodLinkId(GenomeDB *gdb, IDType methodLinkId) {
  GenomeDBAdaptor *gda = (GenomeDBAdaptor *)GenomeDB_getAdaptor(gdb);

  return GenomeDBAdaptor_getAllDbLinks(gda, gdb, methodLinkId);
}

char *GenomeDB_setAssembly(GenomeDB *gdb, char *assembly) {
  StrUtil_copyString(&(gdb->assembly), assembly, 0);

  return gdb->assembly;
}

char *GenomeDB_setName(GenomeDB *gdb, char *name) {
  StrUtil_copyString(&(gdb->name), name, 0);

  return gdb->name;
}

void GenomeDB_free(GenomeDB *gdb) {
  Object_decRefCount(gdb);

  if (Object_getRefCount(gdb) > 0) {
    return;
  } else if (Object_getRefCount(gdb) < 0) {
    fprintf(stderr,"Error: Negative reference count for GenomeDB\n"
                   "       Freeing it anyway\n");
  }

  if (gdb->name)     free(gdb->name);
  if (gdb->assembly) free(gdb->assembly);

  free(gdb);
}

