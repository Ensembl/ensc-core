#define __ECOS_MAIN__
#include "EnsC.h"
#undef __ECOS_MAIN__
#include "Stream.h"
#include "StrUtil.h"

void initEnsC(int argc, char **argv) {
  if (!EcoString_initTable(&ecoSTable)) {
    fprintf(stderr, "Failed initialising ecoSTable\n");
    exit(1);
  }

  StrUtil_copyString(&EnsC_progName, argv[0], 0);

  Stream_setDefaults(0);
}


int idTypeCompFunc(const void *one, const void *two) {
  IDType id1 = *((IDType *)one);
  IDType id2 = *((IDType *)two);

  return id1-id2;
}

long *long_new(long val) {
  long *longP;
  if ((longP = (long *)calloc(1,sizeof(long))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for long\n");
    return NULL;
  }
  *longP = val;
  return longP;
}

IDType *IDType_new(IDType val) {
  IDType *IDTypeP;
  if ((IDTypeP = (IDType *)calloc(1,sizeof(IDType))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for IDType\n");
    return NULL;
  }
  *IDTypeP = val;
  return IDTypeP;
}

