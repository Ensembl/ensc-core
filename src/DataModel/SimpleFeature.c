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
  Object_incRefCount(sf);

  sf->funcs = &simpleFeatureFuncs;
 
  return sf;
}

ECOSTRING SimpleFeature_setDisplayLabel(SimpleFeature *sf, char *label) {
  EcoString_copyStr(ecoSTable, &(sf->displayLabel), label, 0);
  return sf->displayLabel;
}

void SimpleFeature_free(SimpleFeature *sf) {
  Object_decRefCount(sf);

  if (Object_getRefCount(sf) > 0) {
    return;
  } else if (Object_getRefCount(sf) < 0) {
    fprintf(stderr,"Error: Negative reference count for SimpleFeature\n"
                   "       Freeing it anyway\n");
  }

  if (sf->displayLabel) EcoString_freeStr(ecoSTable, sf->displayLabel);

  SeqFeature_freePtrs((SeqFeature *)sf);
  free(sf);
}

