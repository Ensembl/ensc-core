#ifndef __GENOMICALIGNADAPTOR_H__
#define __GENOMICALIGNADAPTOR_H__

#include "BaseComparaAdaptor.h"
#include "ComparaAdaptorTypes.h"
#include "GenomicAlign.h"
#include "Vector.h"

struct GenomicAlignAdaptorStruct {
  BASECOMPARAADAPTOR_DATA
  int maxAlignmentLength;
};

GenomicAlignAdaptor *GenomicAlignAdaptor_new(ComparaDBAdaptor *dba);

Vector *GenomicAlignAdaptor_fetchAllByDNAFragGenomeDB(GenomicAlignAdaptor *gaa,
               DNAFrag *dnaFrag, GenomeDB *targetGenome, int *startP, int *endP,
               char *alignmentType);
Vector *GenomicAlignAdaptor_objectsFromStatementHandle(GenomicAlignAdaptor *gaa, StatementHandle *sth,
                                                       int reverse);
void GenomicAlignAdaptor_addDerivedAlignments(GenomicAlignAdaptor *gaa,
                     Vector *mergedAligns, GenomicAlign *alignA, GenomicAlign *alignB);
IDType GenomicAlignAdaptor_methodLinkIdByAlignmentType(GenomicAlignAdaptor *gaa, char *alignmentType);
char *GenomicAlignAdaptor_alignmentTypeByMethodLinkId(GenomicAlignAdaptor *gaa, IDType methodLinkId);




#endif
