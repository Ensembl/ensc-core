/*
 * Copyright [1999-2024] EMBL-European Bioinformatics Institute
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

#ifndef __GENOMICRANGE_H__
#define __GENOMICRANGE_H__

#include "DataModelTypes.h"

struct GenomicRangeStruct {
  char *chrName;
  int   chrStart;
  int   chrEnd;
  long  chrId;
};

GenomicRange *GenomicRange_new();

#define GenomicRange_setChrStart(gr,s) (gr)->chrStart = (s)
#define GenomicRange_getChrStart(gr) (gr)->chrStart

#define GenomicRange_setChrEnd(gr,e) (gr)->chrEnd = (e)
#define GenomicRange_getChrEnd(gr) (gr)->chrEnd

char *setChrName(GenomicRange *gr, char *chrName);
#define GenomicRange_getChrName(gr) (gr)->chrName

#define GenomicRange_setChrId(gr,i) (gr)->chrId = (i)
#define GenomicRange_getChrId(gr) (gr)->chrId

char *GenomicRange_setChrName(GenomicRange *gr, char *chrName);

void GenomicRange_free(GenomicRange *range);

void GenomicRange_expand(GenomicRange *range, int pad);


#endif
