#include "GenomeDB.h"
#include "GenomeDBAdaptor.h"

GenomeDB *GenomeDB_new() {
  GenomeDB *gdb;

  if ((gdb = (GenomeDB *)calloc(1,sizeof(GenomeDB))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for gdb\n");
    return NULL;
  }

  gdb->objectType = CLASS_GENOMEDB;
  return gdb;
}

int GenomeDB_hasConsensus(GenomeDB *gdb, GenomeDB *conGdb, IDType methodLinkId) {
  GenomeDBAdaptor *gda = GenomeDB_getAdaptor(gdb);

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

  return GenomeDBAdaptor_checkForConsensusDB(gda, gdb, consensusGdb, methodLinkId);
}

int GenomeDB_hasQuery(GenomeDB *gdb, GenomeDB *queryGdb, IDType methodLinkId) {
  GenomeDBAdaptor *gda = GenomeDB_getAdaptor(gdb);

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

  return GenomeDBAdaptor_checkForQueryDB(gda, gdb, queryGdb, methodLinkId);
}

Vector *GenomeDB_linkedGenomesByMethodLinkId(GenomeDB *gdb, IDType methodLinkId) {
  GenomeDBAdaptor *gda = GenomeDB_getAdaptor(gdb);

  return GenomeDBAdaptor_getAllDBLinks(gda, gdb, methodLinkId);
}
