#ifndef __BASEADAPTOR_H__
#define __BASEADAPTOR_H__

#include <mysql.h>
#include <string.h>
#include <stdlib.h>

#include "AdaptorTypes.h"
//#include "DBAdaptor.h"
#include "StatementHandle.h"


#ifdef __MAIN_C__
char *Adaptor_TypeStrings[] = {
  "NONE",
  "GENE",
  "TRANSCRIPT",
  "TRANSLATION",
  "EXON",
  "ANALYSIS",
  "CLONE",
  "RAWCONTIG",
  "SLICE",
  "CHROMOSOME",
  "METACONTAINER",
  "SEQUENCE",
  "DNAALIGNFEATURE",
  "ASSEMBLYMAPPER",
  "SIMPLEFEATURE",
  "PROTEINALIGNFEATURE",
  "REPEATFEATURE",
  "REPEATCONSENSUS",
  "PREDICTIONTRANSCRIPT",
  "SUPPORTINGFEATURE",
  "DBENTRY",
  "GENOMEDB",
  "GENOMICALIGN",
  "DNAFRAG",
  "COMPARADNAALIGNFEATURE",
  "SYNTENY",
  "SYNTENYREGION",
  "HOMOLOGY"
};

#else
 extern char *Adaptor_TypeStrings[];
#endif

enum Adaptor_Types {
  NONE,
  GENE_ADAPTOR, 
  TRANSCRIPT_ADAPTOR, 
  TRANSLATION_ADAPTOR, 
  EXON_ADAPTOR,
  ANALYSIS_ADAPTOR,
  CLONE_ADAPTOR,
  RAWCONTIG_ADAPTOR,
  SLICE_ADAPTOR,
  CHROMOSOME_ADAPTOR,
  META_CONTAINER,
  SEQUENCE_ADAPTOR,
  DNAALIGNFEATURE_ADAPTOR,
  ASSEMBLYMAPPER_ADAPTOR,
  SIMPLEFEATURE_ADAPTOR,
  PROTEINALIGNFEATURE_ADAPTOR,
  REPEATFEATURE_ADAPTOR,
  REPEATCONSENSUS_ADAPTOR,
  PREDICTIONTRANSCRIPT_ADAPTOR,
  SUPPORTINGFEATURE_ADAPTOR,
  DBENTRY_ADAPTOR,
  GENOMEDB_ADAPTOR,
  GENOMICALIGN_ADAPTOR,
  DNAFRAG_ADAPTOR,
  COMPARADNAALIGNFEATURE_ADAPTOR,
  SYNTENY_ADAPTOR,
  SYNTENYREGION_ADAPTOR,
  HOMOLOGY_ADAPTOR
};

typedef StatementHandle *(*BaseAdaptor_PrepareFunc)(BaseAdaptor *ba, char *qStr,size_t len);

#define BASEADAPTOR_DEF(DBADAPTOR_TYPE) \
  DBADAPTOR_TYPE *dba; \
  int adaptorType; \
  BaseAdaptor_PrepareFunc prepare;

#define BASEADAPTOR_DATA BASEADAPTOR_DEF(DBAdaptor)

struct BaseAdaptorStruct {
  BASEADAPTOR_DATA
};

void BaseAdaptor_init(BaseAdaptor *ba, DBAdaptor *dba, int adaptorType);
StatementHandle *BaseAdaptor_prepare(BaseAdaptor *ba, char *qStr,size_t len);

#endif
