#include "SimpleFeature.h"
#include "StrUtil.h"

SimpleFeature *SimpleFeature_new() {
  SimpleFeature *sf;

  if ((sf = (SimpleFeature *)calloc(1,sizeof(SimpleFeature))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for sf\n");
    return NULL;
  }

  sf->objectType = CLASS_SIMPLEFEATURE;
 
  return sf;
}

char *SimpleFeature_setDisplayLabel(SimpleFeature *sf, char *label) {
  sf->displayLabel = StrUtil_copyString(&(sf->displayLabel), label, 0);
  return sf->displayLabel;
}
