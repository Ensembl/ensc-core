#include <stdio.h>

#include "SliceAdaptor.h"
#include "ComparaDBAdaptor.h"
#include "EnsC.h"
#include "HomologyAdaptor.h"

#include "BaseComparaDBTest.h"

int main(int argc, char *argv[]) {
  ComparaDBAdaptor *cdba;
  HomologyAdaptor *ha;
  Slice *slice = NULL;
  Vector *homolList;
  Vector *homols;
  int i;
  int failed;
  
  initEnsC();

  cdba = Test_initComparaDB();

  //slice = Test_getStandardSlice(dba);

  ok(1, slice!=NULL);

  ha = ComparaDBAdaptor_getHomologyAdaptor(cdba);

  ok(2, ha!=NULL);

  homolList =  HomologyAdaptor_listStableIdsFromSpecies(ha,"homo sapiens");

  ok(3, homolList!=NULL);
  ok(4, Vector_getNumElement(homolList)!=0);

  for (i=0;i<Vector_getNumElement(homolList);i++) {
    char *sid = Vector_getElementAt(homolList,i);
    Vector *geneHomols;
    int j;

    printf("sid = %s\n",sid);

    geneHomols = HomologyAdaptor_fetchHomologuesOfGeneInSpecies(ha, "homo sapiens",sid,"mus musculus");

    for (j=0;j<Vector_getNumElement(geneHomols);j++) {
      Homology *hom = Vector_getElementAt(geneHomols,j);
      printf(" homol = %s\n",Homology_getStableId(hom));
    }
    Vector_free(geneHomols,NULL);
  }

  return 0;
}
