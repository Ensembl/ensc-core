#ifndef __COMPARADNAALIGNFEATUREADAPTOR_H__
#define __COMPARADNAALIGNFEATUREADAPTOR_H__

#include "BaseComparaAdaptor.h"
#include "ComparaAdaptorTypes.h"
#include "DNAAlignFeature.h"
#include "Vector.h"
#include "Cache.h"

struct ComparaDNAAlignFeatureAdaptorStruct {
  BASECOMPARAADAPTOR_DATA
  Cache *regionCache;
};

ComparaDNAAlignFeatureAdaptor *ComparaDNAAlignFeatureAdaptor_new(ComparaDBAdaptor *dba);
Vector *ComparaDNAAlignFeatureAdaptor_fetchAllBySpeciesRegion(ComparaDNAAlignFeatureAdaptor *dafa,
                                                   char *csSpecies, char *csAssembly,
                                                   char *qySpecies, char *qyAssembly,
                                                   char *chrName, int start, int end,
                                                   char *alignmentType);
Vector *ComparaDNAAlignFeatureAdaptor_fetchAllBySlice(ComparaDNAAlignFeatureAdaptor *dafa,
              Slice *slice, char *qySpecies, char *qyAssembly, char *assemblyType);

#endif
