/*
 * Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

#ifndef __SEQREGIONRANGE_H__
#define __SEQREGIONRANGE_H__

#include "DataModelTypes.h"
#include "EnsC.h"

struct SeqRegionRangeStruct {
  char * srName;
  int    srStart;
  int    srEnd;
  IDType srId;
};

SeqRegionRange *SeqRegionRange_new();

#define SeqRegionRange_setSeqRegionStart(srr,s) (srr)->srStart = (s)
#define SeqRegionRange_getSeqRegionStart(srr) (srr)->srStart

#define SeqRegionRange_setSeqRegionEnd(srr,e) (srr)->srEnd = (e)
#define SeqRegionRange_getSeqRegionEnd(srr) (srr)->srEnd

char *setSeqRegionName(SeqRegionRange *srRange, char *srName);
#define SeqRegionRange_getSeqRegionName(srr) (srr)->srName

#define SeqRegionRange_setSeqRegionId(srr,i) (srr)->srId = (i)
#define SeqRegionRange_getSeqRegionId(srr) (srr)->srId

char *SeqRegionRange_setSeqRegionName(SeqRegionRange *srRange, char *srName);

void SeqRegionRange_free(SeqRegionRange *range);

void SeqRegionRange_expand(SeqRegionRange *range, int pad);


#endif
