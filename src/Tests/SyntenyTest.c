#include <stdio.h>

#include "SliceAdaptor.h"
#include "ComparaDBAdaptor.h"
#include "EnsC.h"
#include "SyntenyAdaptor.h"

#include "BaseComparaDBTest.h"

int main(int argc, char *argv[]) {
  ComparaDBAdaptor *cdba;
  SyntenyAdaptor *sa;
  Slice *slice = NULL;
  Vector *synRegions;
  int i;
  int failed;
  
  initEnsC(argc, argv);

  cdba = Test_initComparaDB();

  //slice = Test_getStandardSlice(dba);

  ok(1, slice!=NULL);

  sa = ComparaDBAdaptor_getSyntenyAdaptor(cdba);

  ok(2, sa!=NULL);

  SyntenyAdaptor_setSpecies(sa, "homo sapiens","mus musculus");
  
  synRegions =  SyntenyAdaptor_getSyntenyForChromosome(sa,"1",NULL,NULL);

  ok(3, synRegions!=NULL);
  ok(4, Vector_getNumElement(synRegions)!=0);

  for (i=0; i<Vector_getNumElement(synRegions); i++) {
    SyntenyRegion *sr = Vector_getElementAt(synRegions,i);
    printf(" %s %d %d and %s %d %d\n", SyntenyRegion_getChrName(sr),
                                       SyntenyRegion_getChrStart(sr),
                                       SyntenyRegion_getChrEnd(sr),
                                       SyntenyRegion_getHitChrName(sr),
                                       SyntenyRegion_getHitChrStart(sr),
                                       SyntenyRegion_getHitChrEnd(sr));
  }

  return 0;
}
