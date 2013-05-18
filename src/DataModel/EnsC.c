#define __ECOS_MAIN__
#include "EnsC.h"
#undef __ECOS_MAIN__
#include "Stream.h"

void initEnsC(int argc, char **argv) {
  if (!EcoString_initTable(&ecoSTable)) {
    fprintf(stderr, "Failed initialising ecoSTable\n");
    exit(1);
  }

  StrUtil_copyString(&EnsC_progName, argv[0]);

  Stream_setDefaults(0);
}


int idTypeCompFunc(const void *one, const void *two) {
  IDType id1 = *((IDType *)one);
  IDType id2 = *((IDType *)two);

  return id1-id2;
}

