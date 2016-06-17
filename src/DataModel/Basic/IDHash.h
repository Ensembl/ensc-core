/*
 * Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

#ifndef __IDHASH_H__
#define __IDHASH_H__

#include "EnsC.h"
#include "Vector.h"

typedef enum IDHashSizesEnum {
  IDHASH_SMALL,
  IDHASH_MEDIUM,
  IDHASH_LARGE
} IDHashSizes;

typedef struct IDKeyValuePairStruct {
  IDType key;
  void *value;
} IDKeyValuePair;

typedef struct IDHashStruct {
  IDKeyValuePair **buckets;
  int  *bucketCounts;
  int   size;
  int   nValue;
} IDHash;


IDHash * IDHash_new(IDHashSizes size);
int      IDHash_add(IDHash *idHash, IDType id, void *val);
int      IDHash_contains(IDHash *idHash, IDType id);
void     IDHash_free(IDHash *idHash, void freeFunc());
int      IDHash_getNumValues(IDHash *idHash);
IDType * IDHash_getKeys(IDHash *idHash);
void *   IDHash_getValue(IDHash *idHash, IDType id);
void **  IDHash_getValues(IDHash *idHash);
int      IDHash_remove(IDHash *idHash, IDType id, void freeFunc());
Vector *IDHash_getValuesVector(IDHash *idHash);

#endif
