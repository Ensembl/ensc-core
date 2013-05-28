#include <stdio.h>

#include "SliceAdaptor.h"
#include "GeneAdaptor.h"
#include "DBAdaptor.h"
#include "EnsC.h"
#include "Gene.h"

#include "BaseRODBTest.h"
#include "BaseRWDBTest.h"
#include "gperftools/tcmalloc.h"

int main(int argc, char *argv[]) {
  DBAdaptor *dba;
  DBAdaptor *writeDba;
  GeneAdaptor *ga;
  Slice *slice;
  Vector *genes;
  int i;
  int failed;
  
  initEnsC(argc, argv);

  dba = Test_initROEnsDB();

  writeDba = Test_initRWEnsDB();

  slice = Test_getStandardSlice(dba);

  ok(1, slice!=NULL);

  ga = DBAdaptor_getGeneAdaptor(writeDba);
  SliceAdaptor *sa = DBAdaptor_getSliceAdaptor(dba);

  ok(2, ga!=NULL);

  //genes =  Slice_getAllGenes(slice,NULL,NULL, NULL,NULL);

//  Slice *slice2 = SliceAdaptor_fetchByRegion(sa,"chromosome","1",1000000,4000000,1,NULL,0);
  genes =  Slice_getAllGenes(slice,NULL,NULL,1,NULL,NULL);

  ok(3, genes!=NULL);
  ok(4, Vector_getNumElement(genes)!=0);

  fprintf(stderr,"Have %d genes to store\n", Vector_getNumElement(genes));
  for (i=0; i<Vector_getNumElement(genes); i++) {
    fprintf(stderr, "storing gene i = %d\n",i);
    Gene *gene = Vector_getElementAt(genes, i);

    GeneAdaptor_store(ga, gene, 1);
  }

  return 0;
}
