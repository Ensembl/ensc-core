#include "DNAFrag.h"
#include "StrUtil.h"

DNAFrag *DNAFrag_new() {
  DNAFrag *df;

  if ((df = (DNAFrag *)calloc(1,sizeof(DNAFrag))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for df\n");
    return NULL;
  }

  df->objectType = CLASS_DNAFRAG;
  return df;
}

char *DNAFrag_setName(DNAFrag *df, char *name) {
  StrUtil_copyString(&(df->name),name,0);

  return df->name;
}
