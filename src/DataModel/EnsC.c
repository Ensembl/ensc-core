#define __ECOS_MAIN__
#include "EnsC.h"
#undef __ECOS_MAIN__
#include "Stream.h"

void initEnsC(void) {
  if (!EcoString_initTable(&ecoSTable)) {
    fprintf(stderr, "Failed initialising ecoSTable\n");
    exit(1);
  }

  Stream_setDefaults(0);
}
