#include <stdio.h>

#include "SliceAdaptor.h"
#include "DBAdaptor.h"
#include "EnsC.h"
#include "DNAAlignFeature.h"

#include "BaseRODBTest.h"

int main(int argc, char *argv[]) {
  DBAdaptor *dba;
  DNAAlignFeatureAdaptor *dafa;
  Slice *slice;
  Vector *features;
  int i;
  int failed;
  
  initEnsC();

  dba = Test_initROEnsDB();

  slice = Test_getStandardSlice(dba);

  ok(1, slice!=NULL);

  dafa = DBAdaptor_getDNAAlignFeatureAdaptor(dba);

  ok(2, dafa!=NULL);

  features =  Slice_getAllDNAAlignFeatures(slice,"",NULL);

  ok(3, features!=NULL);
  ok(4, Vector_getNumElement(features)!=0);

  failed = 0;
  for (i=0;i<Vector_getNumElement(features) && !failed;i++) {
    DNAAlignFeature *daf = Vector_getElementAt(features,i);
    Vector *ungapped;
    char *oldCigar = DNAAlignFeature_getCigarString(daf);

    printf("in loop hit strand = %d\n",DNAAlignFeature_getHitStrand(daf));
    ungapped = DNAAlignFeature_getUngappedFeatures((BaseAlignFeature *)daf);
    if (!ungapped) failed = 1;
    printf(" cigar = %s num ungapped %d\n",DNAAlignFeature_getCigarString(daf), Vector_getNumElement(ungapped));

    BaseAlignFeature_parseFeatures(daf,ungapped); 
    printf(" cigar now = %s\n",DNAAlignFeature_getCigarString(daf));
    Vector_free(ungapped,NULL);
    if (strcmp(oldCigar,DNAAlignFeature_getCigarString(daf))) {
      printf(" cigars different %s %s\n",oldCigar, DNAAlignFeature_getCigarString(daf));
      failed = 1;
    }
  }
  ok(5, !failed);

  return 0;
}
