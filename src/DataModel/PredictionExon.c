#define __PREDICTIONEXON_MAIN__
#include "PredictionExon.h"
#undef __PREDICTIONEXON__MAIN__
#include "StrUtil.h"

PredictionExon *PredictionExon_new() {
  PredictionExon *pe;

  if ((pe = (PredictionExon *)calloc(1,sizeof(PredictionExon))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for pe\n");
    return NULL;
  }

  pe->objectType = CLASS_PREDICTIONEXON;
  Object_incRefCount(pe);

  pe->funcs = &predictionExonFuncs;
 
  return pe;
}

ECOSTRING PredictionExon_setDisplayLabel(PredictionExon *pe, char *label) {
  EcoString_copyStr(ecoSTable, &(pe->displayLabel), label, 0);
  return pe->displayLabel;
}

void PredictionExon_freeImpl(PredictionExon *pe) {
  Object_decRefCount(pe);

  if (Object_getRefCount(pe) > 0) {
    return;
  } else if (Object_getRefCount(pe) < 0) {
    fprintf(stderr,"Error: Negative reference count for PredictionExon\n"
                   "       Freeing it anyway\n");
  }

  if (pe->displayLabel) EcoString_freeStr(ecoSTable, pe->displayLabel);

  SeqFeature_freePtrs((SeqFeature *)pe);
  free(pe);
}

