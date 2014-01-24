/*
 * Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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
#include "PredictionTranscript.h"

#include "BaseRODBTest.h"

int main(int argc, char *argv[]) {
  DBAdaptor *dba;
  PredictionTranscriptAdaptor *pta;
  Slice *slice;
  Vector *features;
  int i;
  int failed;
  
  initEnsC(argc, argv);

  dba = Test_initROEnsDB();

  slice = Test_getStandardSlice(dba);

  ok(1, slice!=NULL);

  pta = DBAdaptor_getPredictionTranscriptAdaptor(dba);

  ok(2, pta!=NULL);

  features =  Slice_getAllPredictionTranscripts(slice,NULL,1, NULL);

  ok(3, features!=NULL);
  ok(4, Vector_getNumElement(features)!=0);

  failed = 0;
  for (i=0;i<Vector_getNumElement(features) && !failed;i++) {
    PredictionTranscript *pt = Vector_getElementAt(features,i);
    fprintf(stderr, "PT %s\n", PredictionTranscript_getDisplayLabel(pt));
    fprintf(stderr, " from %d to %d strand %d number of exons %d\n",
           PredictionTranscript_getStart(pt),
           PredictionTranscript_getEnd(pt),
           PredictionTranscript_getStrand(pt),
           PredictionTranscript_getExonCount(pt));
    fprintf(stderr, " translation = %s\n",PredictionTranscript_translate(pt));
  }
  ok(5, !failed);
  return 0;
}
