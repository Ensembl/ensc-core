#ifndef __DNAFRAGADAPTOR_H__
#define __DNAFRAGADAPTOR_H__

#include "BaseComparaAdaptor.h"
#include "ComparaAdaptorTypes.h"
#include "DNAFrag.h"
#include "Vector.h"

struct DNAFragAdaptorStruct {
  BASECOMPARAADAPTOR_DATA
};

DNAFrag *DNAFragAdaptor_fetchByDbID(DNAFragAdaptor *dfa, IDType dbID);
Vector *DNAFragAdaptor_fetchAllByGenomeDBRegion(DNAFragAdaptor *dfa, GenomeDB *genomeDB,
                      char *dnaFragType, char *name, int *startP, int *endP);
Vector *DNAFragAdaptor_fetchAll(DNAFragAdaptor *dfa);
Vector *DNAFragAdaptor_objectsFromStatementHandle(DNAFragAdaptor *dfa, StatementHandle *sth);
IDType DNAFragAdaptor_store(DNAFragAdaptor *dfa, DNAFrag *dnaFrag);
IDType DNAFragAdaptor_storeIfNeeded(DNAFragAdaptor *dfa, DNAFrag *dnaFrag);


#endif
