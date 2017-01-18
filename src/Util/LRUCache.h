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

#ifndef __LRUCACHE_H__
#define __LRUCACHE_H__

#include "StringHash.h"
#include <stdio.h>

typedef struct LRUCacheStruct LRUCache;
typedef struct LRUCacheElementStruct LRUCacheElement;

typedef void (*LRUCache_FreeFunc)(void *);

struct LRUCacheStruct {
  int curSize;
  int maxSize;
  StringHash *hash;
  LRUCacheElement *head;
  LRUCacheElement *tail;
};

struct LRUCacheElementStruct {
  void *val;
  char *key;
  int nAccess;
  int size;
  LRUCacheElement *prev;
  LRUCacheElement *next;
  LRUCache_FreeFunc freeFunc; 
};

LRUCache *LRUCache_new(int size);
int LRUCache_contains(LRUCache *cache, char *key);
void LRUCache_empty(LRUCache *cache);
void *LRUCache_get(LRUCache *cache, char *key);
int LRUCache_put(LRUCache *cache, char *key, void *data, LRUCache_FreeFunc freeFunc, int size);
void LRUCache_remove(LRUCache *cache, char *key);
int LRUCache_getSize(LRUCache *cache, char *key);

void LRUCacheElement_free(LRUCacheElement *ce);
LRUCacheElement *LRUCacheElement_new(char *key, void *data, LRUCache_FreeFunc freeFunc, int size);

#endif
