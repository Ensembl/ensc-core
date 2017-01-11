/*
 * Copyright [1999-2017] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

#define __SEQUENCE_MAIN__
#include "Sequence.h"
#undef __SEQUENCE_MAIN__

#include "StrUtil.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

Sequence *Sequence_new() {
  Sequence *seq;

  if ((seq = (Sequence *)calloc(1,sizeof(Sequence))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for seq\n");
    return NULL;
  }

  seq->length = -1;

  seq->objectType = CLASS_SEQUENCE;

  seq->funcs = &sequenceFuncs;

  return seq;
}

ECOSTRING Sequence_setName(Sequence *seq, char *name) {
  if (EcoString_copyStr(ecoSTable,&(seq->name),name,0)) {
    fprintf(stderr,"ERROR: Failed allocating space for seq name\n");
    return NULL;
  }

  return seq->name;
}

char *Sequence_setSeq(Sequence *seq, char *seqstr) {
  if (StrUtil_copyString(&(seq->seq),seqstr,0)) {
    fprintf(stderr,"ERROR: Failed allocating space for seq\n");
    return NULL;
  }

  seq->length = strlen(seq->seq);

  return seq->seq;
}

void Sequence_freePtrs(Sequence *seq) {
  if (seq->name) free(seq->name);
// Is this OK!??
  if (seq->seq) free(seq->seq);

}
