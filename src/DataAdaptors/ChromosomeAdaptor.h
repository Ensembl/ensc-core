#ifndef __CHROMOSOMEADAPTOR_H__
#define __CHROMOSOMEADAPTOR_H__

#include "BaseAdaptor.h"
#include "AdaptorTypes.h"
#include "Chromosome.h"
#include "IDHash.h"
#include "StringHash.h"

struct ChromosomeAdaptorStruct {
  BASEADAPTOR_DATA
  IDHash *chrCache;
  StringHash *chrNameCache;
};

ChromosomeAdaptor *ChromosomeAdaptor_new(DBAdaptor *dba);
Chromosome *ChromosomeAdaptor_fetchByDbID(ChromosomeAdaptor *ca, IDType dbID);
Chromosome *ChromosomeAdaptor_fetchByChrName(ChromosomeAdaptor *ca, char *chrName);
Chromosome *ChromosomeAdaptor_chromosomeFromRow(ChromosomeAdaptor *ca, ResultRow  *row);


#endif
