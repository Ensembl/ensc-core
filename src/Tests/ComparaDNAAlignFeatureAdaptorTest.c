#include <stdio.h>

#include "SliceAdaptor.h"
#include "ComparaDBAdaptor.h"
#include "EnsC.h"
#include "ComparaDNAAlignFeatureAdaptor.h"

#include "BaseComparaDBTest.h"

int main(int argc, char *argv[]) {
  ComparaDBAdaptor *cdba;
  ComparaDNAAlignFeatureAdaptor *cdafa;
  Slice *slice = NULL;
  Vector *dnaAligns;
  int i;
  int failed;
  
  initEnsC();

  cdba = Test_initComparaDB();

  slice = Test_getStandardSlice(cdba);

  ok(1, slice!=NULL);

  cdafa = ComparaDBAdaptor_getComparaDNAAlignFeatureAdaptor(cdba);

  ok(2, cdafa!=NULL);

  dnaAligns = ComparaDNAAlignFeatureAdaptor_fetchAllBySlice(cdafa,slice,"mus musculus","NCBIM30","WGA");

  ok(3, dnaAligns!=NULL);
  ok(4, Vector_getNumElement(dnaAligns)!=0);

/*
  for (i=0; i<Vector_getNumElement(dnaAligns); i++) {
    DNAAlignFeature *daf = Vector_getElementAt(dnaAligns,i);
    printf(" %s %d %d and %s %d %d\n", ComparaDNAAlignFeatureRegion_getChrName(daf),
                                       ComparaDNAAlignFeatureRegion_getChrStart(daf),
                                       ComparaDNAAlignFeatureRegion_getChrEnd(daf),
                                       ComparaDNAAlignFeatureRegion_getHitChrName(daf),
                                       ComparaDNAAlignFeatureRegion_getHitChrStart(daf),
                                       ComparaDNAAlignFeatureRegion_getHitChrEnd(daf));
  }
*/

  return 0;
}
