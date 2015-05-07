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

#include "Cache.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "StrUtil.h"

Cache *Cache_new(int size) {
  Cache *cache;

  if ((cache = (Cache *)calloc(1,sizeof(Cache))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating cache\n");
    exit(1);
  }

  if ((cache->array = (CacheElement **)calloc(size,sizeof(CacheElement *))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating cache elem array for %d elems\n",size);
    exit(1);
  }

  cache->size = size;

  return cache;
}

int Cache_addElement(Cache *cache, char *key, void *data, Cache_FreeFunc freeFunc) {
  CacheElement *cacheElem;

  if ((cacheElem = CacheElement_new(key,data,freeFunc)) == NULL) {
    fprintf(stderr, "ERROR: Failed adding to cache\n");
    return 0;
  }

  if (cache->array[cache->startPos] != NULL) {
    CacheElement_free(cache->array[cache->startPos]);
  }

  cache->array[cache->startPos] = cacheElem;
  cache->startPos++;

  if (cache->startPos>=cache->size) {
    cache->startPos=0;
  }
  return 1;
}

int Cache_contains(Cache *cache, char *key) {
  int i=0;
  for (i=0;i<cache->size;i++) {
    if (cache->array[i] != NULL &&
        !strcmp(cache->array[i]->key,key)) {
      return 1;
    }
  }

  return 0;
}

void *Cache_findElem(Cache *cache, char *key) {
  int i=0;
  for (i=0;i<cache->size;i++) {
    if (cache->array[i] != NULL &&
        !strcmp(cache->array[i]->key,key)) {
      return cache->array[i]->val;
    }
  }

  return NULL;
}

void Cache_empty(Cache *cache) {
  int i=0;
  for (i=0;i<cache->size;i++) {
    if (cache->array[i] != NULL) {
      CacheElement_free(cache->array[i]);
      cache->array[i] = NULL;
    }
  }

  return;
}

void CacheElement_free(CacheElement *ce) {
  // fprintf(stderr,"CacheElement_free for element with key %s\n", ce->key);
  if (ce->freeFunc) {
    ce->freeFunc(ce->val);
  }
  free(ce->key);
  free(ce);

  return;
}

CacheElement *CacheElement_new(char *key, void *data, Cache_FreeFunc freeFunc) {
  CacheElement *cacheElem;

  if ((cacheElem = (CacheElement *)calloc(1,sizeof(CacheElement))) == NULL) {
    fprintf(stderr, "ERROR: Failed adding to cache\n");
    return NULL;
  }
  StrUtil_copyString(&(cacheElem->key), key, 0);
  cacheElem->val = data;
  cacheElem->freeFunc = freeFunc;

  return cacheElem;
}
