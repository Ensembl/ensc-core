#define __ECOS_MAIN__
#include "EnsC.h"
#undef __ECOS_MAIN__

void initEnsC() {
  if (!EcoString_initTable(&ecoSTable)) {
    fprintf(stderr, "Failed initialising ecoSTable\n");
    exit(1);
  }
}
