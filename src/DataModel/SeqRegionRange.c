#include "SeqRegionRange.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

SeqRegionRange *SeqRegionRange_new() {
  SeqRegionRange *srRange;

  if ((srRange = (SeqRegionRange *)calloc(1,sizeof(SeqRegionRange))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for SeqRegionRange\n");
    return NULL;
  }

  return srRange;
}

char *SeqRegionRange_setSeqRegionName(SeqRegionRange *srRange, char *srName) {
  if ((srRange->srName = (char *)malloc(strlen(srName)+1)) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for srName\n");
    exit(1);
  }

  strcpy(sr->srName,srName);

  return sr->srName;
}

void SeqRegionRange_expand(SeqRegionRange *range, int pad) {
  SeqRegionRange_setSeqRegionStart(range,SeqRegionRange_getSeqRegionStart(range) - pad);
  SeqRegionRange_setSeqRegionEnd(range,SeqRegionRange_getSeqRegionEnd(range) + pad);
  
  if (SeqRegionRange_getSeqRegionStart(range) < 1) {
    SeqRegionRange_setSeqRegionStart(range,1); 
  } 
}

void SeqRegionRange_free(SeqRegionRange *range) {
  free(range->srName);
  free(range);
}
