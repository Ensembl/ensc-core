#ifndef __COMPARADBADAPTOR_H__
#define __COMPARADBADAPTOR_H__

#include "BaseDBAdaptor.h"
#include "DBConnection.h"
#include "ComparaAdaptorTypes.h"
#include "EnsC.h"
#include "StringHash.h"

struct ComparaDBAdaptorStruct {
  BASEDBADAPTOR_DATA
  StringHash *genomes;
};

ComparaDBAdaptor *ComparaDBAdaptor_new(char *host, char *user, char *pass, char *dbname,
                                       unsigned int port, char *confFile);

DBAdaptor *ComparaDBAdaptor_getDBAdaptor(ComparaDBAdaptor *cdba, char *species, char *assembly);
void ComparaDBAdaptor_addDBAdaptor(ComparaDBAdaptor *cdba, DBAdaptor *dba);
DBAdaptor *ComparaDBAdaptor_getDBAdaptor(ComparaDBAdaptor *cdba, char *species, char *assembly);
SyntenyAdaptor *ComparaDBAdaptor_getSyntenyAdaptor(ComparaDBAdaptor *cdba);
GenomeDBAdaptor *ComparaDBAdaptor_getGenomeDBAdaptor(ComparaDBAdaptor *cdba);
DNAFragAdaptor *ComparaDBAdaptor_getDNAFragAdaptor(ComparaDBAdaptor *cdba);
GenomicAlignAdaptor *ComparaDBAdaptor_getGenomicAlignAdaptor(ComparaDBAdaptor *cdba);
HomologyAdaptor *ComparaDBAdaptor_getHomologyAdaptor(ComparaDBAdaptor *cdba);
ComparaDNAAlignFeatureAdaptor *ComparaDBAdaptor_getComparaDNAAlignFeatureAdaptor(ComparaDBAdaptor *cdba);
MetaContainer *ComparaDBAdaptor_getMetaContainer(ComparaDBAdaptor *cdba);




#endif
