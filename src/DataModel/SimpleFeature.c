#define __SIMPLEFEATURE_MAIN__
#include "SimpleFeature.h"
#undef __SIMPLEFEATURE__MAIN__
#include "StrUtil.h"

SimpleFeature *SimpleFeature_new() {
  SimpleFeature *sf;

  if ((sf = (SimpleFeature *)calloc(1,sizeof(SimpleFeature))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for sf\n");
    return NULL;
  }

  sf->objectType = CLASS_SIMPLEFEATURE;

  sf->funcs = &simpleFeatureFuncs;
 
  return sf;
}

char *SimpleFeature_setDisplayLabel(SimpleFeature *sf, char *label) {
  sf->displayLabel = StrUtil_copyString(&(sf->displayLabel), label, 0);
  return sf->displayLabel;
}
