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
#include "DBAdaptor.h"
#include "EnsC.h"
#include "SimpleFeature.h"
#include "BaseFeatureAdaptor.h"

#include "BaseRODBTest.h"

int main(int argc, char *argv[]) {
  DBAdaptor *dba;
  SimpleFeatureAdaptor *sfa;
  Slice *slice;
  Vector *features;
  int i;
  int failed;
  
  initEnsC(argc, argv);

  dba = Test_initROEnsDB();

  slice = Test_getStandardSlice(dba);

  ok(1, slice!=NULL);

  sfa = DBAdaptor_getSimpleFeatureAdaptor(dba);

  ok(2, sfa!=NULL);

  features =  Slice_getAllSimpleFeatures(slice, NULL, NULL, NULL);

  ok(3, features!=NULL);
  ok(4, Vector_getNumElement(features)!=0);

  fprintf(stderr,"Number of features is %d\n", Vector_getNumElement(features));

  failed = 0;
  for (i=0;i<Vector_getNumElement(features) && !failed;i++) {
    SimpleFeature *sf = Vector_getElementAt(features,i);
    long start = SimpleFeature_getStart(sf);
    long end   = SimpleFeature_getEnd(sf);
    Vector *rsfVector;
    SimpleFeature *rsf;

    //printf("slice start = %d end = %d\n",start,end);
    rsf = (SimpleFeature*)SeqFeature_transform((SeqFeature*)sf,"contig",NULL,NULL);

    if (rsf) {
      printf("rc start = %ld end = %ld\n",SimpleFeature_getStart(rsf),SimpleFeature_getEnd(rsf));
    } else {
      printf("no mapped feature\n");
    }

    if (rsf) {
      //sf = SeqFeature_transform(rsf,"chromosome",NULL,slice);
      sf = (SimpleFeature*)SeqFeature_transfer((SeqFeature*)rsf, slice);
      if (SimpleFeature_getStart(sf) != start ||
          SimpleFeature_getEnd(sf) != end) {
        printf("Remapping to slice produced different coords start %ld v %ld   end %ld v %ld\n", 
               SimpleFeature_getStart(sf),start, SimpleFeature_getEnd(sf),end);
        failed =1;
      }
    }
  }
  ok(5, !failed);
  return 0;
}
