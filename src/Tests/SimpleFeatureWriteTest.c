/*
 * Copyright [1999-2024] EMBL-European Bioinformatics Institute
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
#include "SimpleFeatureAdaptor.h"
#include "DBAdaptor.h"
#include "EnsC.h"
#include "SimpleFeature.h"

#include "BaseRODBTest.h"
#include "BaseRWDBTest.h"
#ifdef HAVE_LIBTCMALLOC
#include "gperftools/tcmalloc.h"
#endif

int main(int argc, char *argv[]) {
  int testResult = 0;
  DBAdaptor *dba;
  DBAdaptor *writeDba;
  SimpleFeatureAdaptor *sfa;
  Slice *slice;
  Vector *features;
  int i;
  int failed;
  
  initEnsC(argc, argv);

  dba = Test_initROEnsDB();

  writeDba = Test_initRWEnsDB();

  slice = Test_getStandardSlice(dba);

  testResult += ok(1, slice!=NULL);

  sfa = DBAdaptor_getSimpleFeatureAdaptor(writeDba);
  SliceAdaptor *sa = DBAdaptor_getSliceAdaptor(dba);

  testResult += ok(2, sfa!=NULL);

  //features =  Slice_getAllSimpleFeatures(slice,NULL,NULL, NULL,NULL);

  //Slice *slice3 = SliceAdaptor_fetchByRegion(sa,"chromosome","1",2,260000000,1,NULL,0);
  Slice *slice2 = SliceAdaptor_fetchByRegion(sa,"chromosome","1",1000000,4000000,1,NULL,0);
  features =  Slice_getAllSimpleFeatures(slice2,NULL,NULL, NULL);

  testResult += ok(3, features!=NULL);
  testResult += ok(4, Vector_getNumElement(features)!=0);

  SimpleFeatureAdaptor_store(sfa, features);

  return testResult;
}
