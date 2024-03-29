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

#ifndef __CACHE_H__
#define __CACHE_H__

#include "EcoString.h"
#include <stdio.h>

typedef struct CacheStruct Cache;
typedef struct CacheElementStruct CacheElement;

typedef void (*Cache_FreeFunc)(void *);

struct CacheStruct {
  int size;
  int startPos;
  CacheElement **array;
};

struct CacheElementStruct {
  void *val;
  char *key;
  Cache_FreeFunc freeFunc; 
};

Cache *Cache_new(int size);
void *Cache_findElem(Cache *cache, char *key);
void Cache_empty(Cache *cache);


CacheElement *CacheElement_new(char *key, void *val, Cache_FreeFunc freeFunc);
int Cache_addElement(Cache *cache, char *key, void *data, Cache_FreeFunc freeFunc);

void CacheElement_free(CacheElement *ce);
int Cache_contains(Cache *cache, char *key);

#endif
