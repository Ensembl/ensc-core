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

#ifndef __SEQUENCEADAPTOR_H__
#define __SEQUENCEADAPTOR_H__

#include "BaseAdaptor.h"
#include "AdaptorTypes.h"
#include "RawContig.h"

#include "Cache.h"
#include "IDHash.h"
#include "StringHash.h"
#include "LRUCache.h"

struct SequenceAdaptorStruct {
  BASEADAPTOR_DATA

  LRUCache *seqCache;
//  StringHash *seqCache;
  IDHash *rnaEditsCache;
};

SequenceAdaptor *SequenceAdaptor_new(DBAdaptor *dba);
StatementHandle *SequenceAdaptor_prepare(BaseAdaptor *ba, char *qStr, size_t len);

char *SequenceAdaptor_fetchByRawContigStartEndStrand(SequenceAdaptor *sa,
                                                     //RawContig *rc,
                                                     IDType rcId, 
                                                     int start,
                                                     int end,
                                                     char strand);
char *SequenceAdaptor_fetchByAssemblyLocation(SequenceAdaptor *sa,
          int chrStart, int chrEnd, int strand, char * chrName, char *assemblyType);
char *SequenceAdaptor_fetchBySliceStartEndStrand(SequenceAdaptor *sa,
                                                 Slice *slice, long start, long end,
                                                 int strand);

char *SequenceAdaptor_fetchBySliceStartEndStrandRecursive(SequenceAdaptor *sa,
                                                          Slice *slice, long start, long end,
                                                          int strand, int *recLev);

          

#endif
