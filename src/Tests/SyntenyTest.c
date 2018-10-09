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
#include "ComparaDBAdaptor.h"
#include "EnsC.h"
#include "SyntenyAdaptor.h"

#include "BaseComparaDBTest.h"

int main(int argc, char *argv[]) {
  int testResult = 0;
  ComparaDBAdaptor *cdba;
  SyntenyAdaptor *sa;
  Vector *synRegions;
  int i;
  int testNum = 0;
  
  initEnsC(argc, argv);

  cdba = Test_initComparaDB();

  sa = ComparaDBAdaptor_getSyntenyAdaptor(cdba);

  testResult += ok(++testNum, sa!=NULL);

  SyntenyAdaptor_setSpecies(sa, "homo sapiens","mus musculus");
  
  synRegions =  SyntenyAdaptor_getSyntenyForChromosome(sa,"1",NULL,NULL);

  testResult += ok(++testNum, synRegions!=NULL);
  testResult += ok(++testNum, Vector_getNumElement(synRegions)!=0);

  for (i=0; i<Vector_getNumElement(synRegions); i++) {
    SyntenyRegion *sr = Vector_getElementAt(synRegions,i);
    printf(" %s %d %d and %s %d %d\n", SyntenyRegion_getChrName(sr),
                                       SyntenyRegion_getChrStart(sr),
                                       SyntenyRegion_getChrEnd(sr),
                                       SyntenyRegion_getHitChrName(sr),
                                       SyntenyRegion_getHitChrStart(sr),
                                       SyntenyRegion_getHitChrEnd(sr));
  }

  return testResult;
}
