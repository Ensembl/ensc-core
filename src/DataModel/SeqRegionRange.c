/*
 * Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
 * Copyright [2016-2017] EMBL-European Bioinformatics Institute
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *      http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

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

  strcpy(srRange->srName,srName);

  return srRange->srName;
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
