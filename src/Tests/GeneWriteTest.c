/*
 * Copyright [1999-2017] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *      http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <stdio.h>

#include "SliceAdaptor.h"
#include "GeneAdaptor.h"
#include "DBAdaptor.h"
#include "EnsC.h"
#include "Gene.h"

#include "BaseRODBTest.h"
#include "BaseRWDBTest.h"
#ifdef HAVE_LIBTCMALLOC
#include "gperftools/tcmalloc.h"
#endif

int main(int argc, char *argv[]) {
  int testResult = 0;
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

  testResult += ok(1, slice!=NULL);

  ga = DBAdaptor_getGeneAdaptor(writeDba);
  SliceAdaptor *sa = DBAdaptor_getSliceAdaptor(dba);

  testResult += ok(2, ga!=NULL);

  //genes =  Slice_getAllGenes(slice,NULL,NULL, NULL,NULL);

//  Slice *slice2 = SliceAdaptor_fetchByRegion(sa,"chromosome","1",1000000,4000000,1,NULL,0);
  genes =  Slice_getAllGenes(slice,NULL,NULL,1,NULL,NULL);

  testResult += ok(3, genes!=NULL);
  testResult += ok(4, Vector_getNumElement(genes)!=0);

  fprintf(stderr,"Have %d genes to store\n", Vector_getNumElement(genes));
  for (i=0; i<Vector_getNumElement(genes); i++) {
    fprintf(stderr, "storing gene i = %d\n",i);
    Gene *gene = Vector_getElementAt(genes, i);

    GeneAdaptor_store(ga, gene, 1);
  }

  return testResult;
}
