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

#ifndef __STRINGHASH_H__
#define __STRINGHASH_H__

#include <string.h>

typedef enum StringHashSizesEnum {
  STRINGHASH_SMALL,
  STRINGHASH_MEDIUM,
  STRINGHASH_LARGE,
  STRINGHASH_HUGE
} StringHashSizes;

typedef struct KeyValuePairStruct {
  char *key;
  int keyLen;
  void *value;
} KeyValuePair;

typedef struct StringHashStruct {
  KeyValuePair **buckets;
  int  *bucketCounts;
  int   size;
  int   nValue;
} StringHash;


void StringHash_printHashStats(StringHash *stringHash, char *hashDesc);
StringHash *StringHash_new(StringHashSizes size);
int     StringHash_add(StringHash *stringHash, char *string, void *val);
int     StringHash_contains(StringHash *stringHash, char *string);
void    StringHash_free(StringHash *stringHash, void freeFunc());
void StringHash_freeNoValFree(StringHash *stringHash);
//int     StringHash_getNumValues(StringHash *stringHash);
#define StringHash_getNumValues(stringHash) stringHash->nValue

void *  StringHash_getValue(StringHash *stringHash, char *string);
void *  StringHash_getValues(StringHash *stringHash);
char ** StringHash_getKeys(StringHash *stringHash);
char ** StringHash_getKeysNoCopy(StringHash *stringHash);
char *  StringHash_getKey(StringHash *stringHash, char *key);
int     StringHash_remove(StringHash *stringHash, char *key, void freeFunc());

unsigned int StringHash_getBucketNum(StringHash *stringHash, char *key, int keyLen);


#endif
