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
  
  initEnsC(argc, argv);

  cdba = Test_initComparaDB();

  slice = Test_getStandardSlice(cdba);

  ok(1, slice!=NULL);

  cdafa = ComparaDBAdaptor_getComparaDNAAlignFeatureAdaptor(cdba);

  ok(2, cdafa!=NULL);

  dnaAligns = ComparaDNAAlignFeatureAdaptor_fetchAllBySlice(cdafa,slice,"mus musculus","NCBIM30","WGA");

  ok(3, dnaAligns!=NULL);
  ok(4, Vector_getNumElement(dnaAligns)!=0);

  for (i=0; i<Vector_getNumElement(dnaAligns); i++) {
    DNAAlignFeature *daf = Vector_getElementAt(dnaAligns,i);
    Slice *slice = (Slice *)DNAAlignFeature_getSlice(daf);
    printf(" %s %ld %ld and %s %d %d\n", DNAAlignFeature_getSeqName(daf),
                                       DNAAlignFeature_getStart(daf),
                                       DNAAlignFeature_getEnd(daf),
                                       DNAAlignFeature_getHitSeqName(daf),
                                       DNAAlignFeature_getHitStart(daf),
                                       DNAAlignFeature_getHitEnd(daf));
    printf(" seq = %s\n",Slice_getSubSeq(slice,
                                         DNAAlignFeature_getStart(daf),
                                         DNAAlignFeature_getEnd(daf), 
                                         DNAAlignFeature_getStrand(daf)));
  }

  return 0;
}
