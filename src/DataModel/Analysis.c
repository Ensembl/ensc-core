#include "Analysis.h"

Analysis *Analysis_new() {
  Analysis *anal;

  if ((anal = (Analysis *)calloc(1,sizeof(Analysis))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for anal\n");
    return NULL;
  }

  anal->objectType = CLASS_ANALYSIS;

  return anal;
}

#define COMPARE_STRINGS(func,a,b) \
  str1 = func((a)); \
  str2 = func((b)); \
  \
  if (str1 && !str2) { \
    retVal = 1; \
  } else if (!str1 && str2) { \
    return -1; \
  } else if (str1 && str2) { \
    if (strcmp(str1,str2)) return -1; \
  }

#define COMPARE_INTS(func,a,b) \
  i1 = func((a)); \
  i2 = func((b)); \
  \
  if (i1 && !i2) { \
    retVal = 1; \
  } else if (!i1 && i2) { \
    return -1; \
  } else if (i1 && i2) { \
    if (i1 != i2) return -1; \
  }

int Analysis_compare(Analysis *a, Analysis *b) {
  int retVal = 0;
  char *str1,*str2;
  int i1,i2;

  COMPARE_STRINGS(Analysis_getLogicName,a,b);
  COMPARE_STRINGS(Analysis_getProgram,a,b);
  COMPARE_INTS(Analysis_getProgramVersion,a,b);
  COMPARE_STRINGS(Analysis_getProgramFile,a,b);
  COMPARE_STRINGS(Analysis_getDb,a,b);
  COMPARE_INTS(Analysis_getDbVersion,a,b);
  COMPARE_STRINGS(Analysis_getDbFile,a,b);
  COMPARE_STRINGS(Analysis_getGFFSource,a,b);
  COMPARE_STRINGS(Analysis_getGFFFeature,a,b);
  COMPARE_STRINGS(Analysis_getModule,a,b);
  COMPARE_INTS(Analysis_getModuleVersion,a,b);
  COMPARE_STRINGS(Analysis_getParameters,a,b);
  
  return retVal;
}
