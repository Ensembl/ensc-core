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
  
  initEnsC();

  dba = Test_initROEnsDB();

  slice = Test_getStandardSlice(dba);

  ok(1, slice!=NULL);

  pta = DBAdaptor_getPredictionTranscriptAdaptor(dba);

  ok(2, pta!=NULL);

  features =  Slice_getAllPredictionTranscripts(slice,"");

  ok(3, features!=NULL);
  ok(4, Vector_getNumElement(features)!=0);

  failed = 0;
  for (i=0;i<Vector_getNumElement(features) && !failed;i++) {
    PredictionTranscript *pt = Vector_getElementAt(features,i);
    printf("PT %s\n", PredictionTranscript_getStableId(pt));
    printf(" from %d to %d\n",
           PredictionTranscript_getStart(pt),
           PredictionTranscript_getEnd(pt));
    printf(" translation = %s\n",PredictionTranscript_translate(pt));
  }
  ok(5, !failed);
  return 0;
}
