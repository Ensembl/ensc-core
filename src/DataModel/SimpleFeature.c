#include "SimpleFeature.h"

SimpleFeature *SimpleFeature_new() {
  SimpleFeature *sf;

  if ((sf = (SimpleFeature *)calloc(1,sizeof(SimpleFeature))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for sf\n");
    return NULL;
  }

  return sf;
}
