#include <stdio.h>

#include "SliceAdaptor.h"
#include "DNAAlignFeatureAdaptor.h"
#include "DBAdaptor.h"
#include "EnsC.h"
#include "DNAAlignFeature.h"

#include "BaseRODBTest.h"
#include "gperftools/tcmalloc.h"

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
  SliceAdaptor *sa = DBAdaptor_getSliceAdaptor(dba);

  ok(2, dafa!=NULL);

  //features =  Slice_getAllDNAAlignFeatures(slice,NULL,NULL, NULL,NULL);

  Slice *slice3 = SliceAdaptor_fetchByRegion(sa,"chromosome","1",2,260000000,1,NULL,0);
  //Slice *slice2 = SliceAdaptor_fetchByRegion(sa,"chromosome","Y",1000000,4000000,1,NULL,0);
  features =  Slice_getAllDNAAlignFeatures(slice3,NULL,NULL, NULL,NULL);

  ok(3, features!=NULL);
  ok(4, Vector_getNumElement(features)!=0);
//  MallocExtension_ReleaseFreeMemory();
//  tc_malloc_stats();
//  exit(1);


  failed = 0;

  for (i=0;i<Vector_getNumElement(features) && !failed;i++) {
    DNAAlignFeature *daf = Vector_getElementAt(features,i);
    Vector *ungapped;
    //fprintf(stderr, "slice start %ld end %ld\n", DNAAlignFeature_getStart(daf), DNAAlignFeature_getEnd(daf));
    char *oldCigar = DNAAlignFeature_getCigarString(daf);

    //printf(" cigar = %s pre ungapped\n",DNAAlignFeature_getCigarString(daf));
    ungapped = DNAAlignFeature_getUngappedFeatures((BaseAlignFeature *)daf);
    if (!ungapped) failed = 1;
    //printf(" cigar = %s num ungapped %d\n",DNAAlignFeature_getCigarString(daf), Vector_getNumElement(ungapped));

    BaseAlignFeature_parseFeatures(daf,ungapped); 
    //printf(" cigar now = %s\n",DNAAlignFeature_getCigarString(daf));
// NIY Check that free func has been set
    Vector_setFreeFunc(ungapped, Object_freeImpl);
    Vector_free(ungapped);
    if (strcmp(oldCigar,DNAAlignFeature_getCigarString(daf))) {
      printf(" cigars different %s %s\n",oldCigar, DNAAlignFeature_getCigarString(daf));
      failed = 1;
    }
  }
  ok(5, !failed);
  fprintf(stderr,"Clearing dafa cache\n");
  DNAAlignFeatureAdaptor_clearCache(dafa);

//  Slice *slice2 = SliceAdaptor_fetchByRegion(sa,"chromosome","Y",1000000,4000000,1,NULL,0);
//  features =  Slice_getAllDNAAlignFeatures(slice2,NULL,NULL, NULL,NULL);
//  Slice *slice3 = SliceAdaptor_fetchByRegion(sa,"chromosome","1",10000000,100000000,1,NULL,0);
//  features =  Slice_getAllDNAAlignFeatures(slice3,NULL,NULL, NULL,NULL);
  features =  Slice_getAllDNAAlignFeatures(slice,NULL,NULL, NULL,NULL);

  fprintf(stderr," \n\n\nTest 6 commented out - NEED TO IMPLEMENT TRANSFORMS ON FEATURES\n\n\n");
/*
  failed = 0;
  for (i=0;i<Vector_getNumElement(features) && !failed;i++) {
    DNAAlignFeature *daf = Vector_getElementAt(features,i);
    int start = DNAAlignFeature_getStart(daf);
    int end   = DNAAlignFeature_getEnd(daf);
    Vector *rdafVector;
    DNAAlignFeature *rdaf;

    // printf("slice start = %d end = %d id " IDFMTSTR "\n",start,end,DNAAlignFeature_getDbID(daf));
    rdafVector = DNAAlignFeature_transformToRawContig(daf);
    if (Vector_getNumElement(rdafVector) > 1) {
      printf("Feature mapped to more than one rawcontig\n");
      failed=1;
    }
    rdaf = Vector_getElementAt(rdafVector,0);

    // printf("rc id " IDFMTSTR " rc start = %d end = %d\n",BaseContig_getDbID(DNAAlignFeature_getContig(rdaf)),
    //       DNAAlignFeature_getStart(rdaf),DNAAlignFeature_getEnd(rdaf));
    daf = DNAAlignFeature_transformToSlice(rdaf, slice);
    if (DNAAlignFeature_getStart(daf) != start ||
        DNAAlignFeature_getEnd(daf) != end) {
      // printf("slice start now = %d end = %d\n",DNAAlignFeature_getStart(daf),DNAAlignFeature_getEnd(daf));
      printf("Remapping to slice produced different coords for " IDFMTSTR "\n", DNAAlignFeature_getDbID(daf));
      failed =1;
    }
  }
  ok(6, !failed);
*/
  return 0;
}
