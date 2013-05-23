#define __DNAALIGNFEATURE_MAIN__
#include "DNAAlignFeature.h"
#undef __DNAALIGNFEATURE_MAIN__
//#include "ProcUtil.h"

#include <string.h>

DNAAlignFeature *DNAAlignFeature_new() {
  DNAAlignFeature *daf;

  if ((daf = (DNAAlignFeature *)calloc(1,sizeof(DNAAlignFeature))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for dna align feature\n");
    return NULL;
  }

  daf->objectType = CLASS_DNADNAALIGNFEATURE;
  Object_incRefCount(daf);

  daf->funcs = &dnaAlignFeatureFuncs;

  return daf;
}

int DNAAlignFeature_getHitUnit(void) {
  return 1;
}

int DNAAlignFeature_getQueryUnit(void) {
  return 1;
}

ECOSTRING DNAAlignFeature_setExtraData(DNAAlignFeature *daf, char *extraData) {
  EcoString_copyStr(ecoSTable, &(daf->extraData),extraData,0);

  if (daf->extraData == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for extraData\n");
    return NULL;
  }

  return daf->extraData;
}

void DNAAlignFeature_freeImpl(DNAAlignFeature *daf) {
  Object_decRefCount(daf);

  fprintf(stderr,"DAF_freeImpl called\n");
  if (Object_getRefCount(daf) > 0) {
//    ProcUtil_showBacktrace(EnsC_progName);
//    printf("return\n");
    return;
  } else if (Object_getRefCount(daf) < 0) {
    fprintf(stderr,"Error: Negative reference count for DNAAlignFeature\n"
                   "       Freeing it anyway\n");
  }

  BaseAlignFeature_freePtrs((BaseAlignFeature *)daf);
  
  free(daf);
}

DNAAlignFeature *DNAAlignFeature_shallowCopyImpl(DNAAlignFeature *daf) {
  DNAAlignFeature *newDNAAlignFeature = DNAAlignFeature_new();

  memcpy(newDNAAlignFeature,daf,sizeof(DNAAlignFeature));

  return newDNAAlignFeature;
}
