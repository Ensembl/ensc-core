/*
 * Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
 * Copyright [2016-2017] EMBL-European Bioinformatics Institute
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
#include "ComparaDBAdaptor.h"
#include "EnsC.h"
#include "GenomicAlignAdaptor.h"

#include "BaseComparaDBTest.h"

int main(int argc, char *argv[]) {
  ComparaDBAdaptor *cdba;
  GenomicAlignAdaptor *gaa;
  Slice *slice = NULL;
  Vector *homolList;
  Vector *homols;
  int i;
  int failed;
  
  initEnsC(argc, argv);

  cdba = Test_initComparaDB();

  slice = Test_getStandardSlice(cdba);

  ok(1, slice!=NULL);

  gaa = ComparaDBAdaptor_getGenomicAlignAdaptor(cdba);

  ok(2, gaa!=NULL);

  homolList =  GenomicAlignAdaptor_listStableIdsFromSpecies(gaa,"homo sapiens");

  ok(3, homolList!=NULL);
  ok(4, Vector_getNumElement(homolList)!=0);

  for (i=0;i<Vector_getNumElement(homolList);i++) {
    char *sid = Vector_getElementAt(homolList,i);
    Vector *geneHomols;
    int j;

    printf("sid = %s\n",sid);

    geneHomols = GenomicAlignAdaptor_fetchHomologuesOfGeneInSpecies(gaa, "homo sapiens",sid,"mus musculus");

    for (j=0;j<Vector_getNumElement(geneHomols);j++) {
      GenomicAlign *hom = Vector_getElementAt(geneHomols,j);
      printf(" homol = %s %s %d %d\n",GenomicAlign_getStableId(hom),
                                      GenomicAlign_getChrName(hom),
                                      GenomicAlign_getChrStart(hom),
                                      GenomicAlign_getChrEnd(hom));
    }
    Vector_free(geneHomols,NULL);
  }

  return 0;
}
