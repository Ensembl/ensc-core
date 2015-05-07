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

#include "LRUCache.h"
#include "StringHash.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "StrUtil.h"

LRUCache *LRUCache_new(int size) {
  LRUCache *cache;

  if ((cache = (LRUCache *)calloc(1,sizeof(LRUCache))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating cache\n");
    exit(1);
  }

  cache->hash = StringHash_new(STRINGHASH_MEDIUM);

  cache->head = cache->tail = NULL;

  cache->maxSize = size;

  cache->curSize = 0;

  return cache;
}

int LRUCache_put(LRUCache *cache, char *key, void *data, LRUCache_FreeFunc freeFunc, int size) {
  LRUCacheElement *cacheElem;

  if ((cacheElem = LRUCacheElement_new(key, data, freeFunc, size)) == NULL) {
    fprintf(stderr, "ERROR: Failed adding to cache\n");
    return 0;
  }

  // Do we already have something with this key - if so, need to remove it
  if (StringHash_contains(cache->hash, key)) {
    LRUCache_remove(cache, key);
  }

  // See if we're over size.
  // If we are remove elements from head until we are sufficiently under to be able to add the new one
  LRUCacheElement *elem = cache->head;
  while (cache->curSize + size > cache->maxSize && elem != NULL) {
    LRUCache_remove(cache, elem->key);

    // Setting this to head may seem odd but its right - we've just removed the head so the head is now
    // the next element in the list - I hate linked lists!!!
    elem = cache->head;
  }

  // If we have an empty list its easy 
  if (cache->head == NULL && cache->tail == NULL) {
    cache->head = cache->tail = cacheElem;
    cacheElem->next = cacheElem->prev = NULL;
  } else {
  // Add the new one to tail
    cache->tail->next = cacheElem;
    cacheElem->prev = cache->tail;
    cache->tail = cacheElem;
    cacheElem->next = NULL;
  }

  // Add the new one to hash
  StringHash_add(cache->hash, key, cacheElem);
  
  // Add to size
  cache->curSize += size;

  return 1;
}

// Contains doesn't count as used
int LRUCache_contains(LRUCache *cache, char *key) {

  // Check hash for key
  return StringHash_contains(cache->hash, key);
}

void *LRUCache_get(LRUCache *cache, char *key) {

  // Check hash for key
  if (StringHash_contains(cache->hash, key)) {
    LRUCacheElement *cacheElem = StringHash_getValue(cache->hash, key);


    // Move entry to tail of list - it is now most recently accessed
    if (cacheElem == cache->tail) {
      // Don't rearrange the linked list if its already the tail
      cacheElem->nAccess++;

    } else { 
      // Need to rearrange
      if (cacheElem == cache->head && cache->head->next != NULL) {
        cache->head = cache->head->next;
      }

      if (cacheElem->prev != NULL) {
        cacheElem->prev->next = cacheElem->next;
      }
      if (cacheElem->next != NULL) {
        cacheElem->next->prev = cacheElem->prev;
      }
      cacheElem->prev = cache->tail;
      cacheElem->next = NULL;

      cache->tail->next = cacheElem;

      // Reset tail to be this entry
      cache->tail = cacheElem;
    }
    cacheElem->nAccess++;

    return cacheElem->val;

  } else {
    // Not in cache
    fprintf(stderr,"Warning: key %s not in LRUCache\n",key);
    return NULL;
  }
}

// Size doesn't count as used
int LRUCache_getSize(LRUCache *cache, char *key) {

  // Check hash for key
  if (StringHash_contains(cache->hash, key)) {
    LRUCacheElement *cacheElem = StringHash_getValue(cache->hash, key);

    return cacheElem->size;

  } else {
    // Not in cache
    fprintf(stderr,"Warning: key %s not in LRUCache\n",key);
    return 0;
  }
}

void LRUCache_empty(LRUCache *cache) {
  // Free the LRUCacheElements using the list
  LRUCacheElement *elem = cache->head;
  while (elem != NULL) {
    LRUCache_remove(cache, elem->key);

    // Setting this to head may seem odd but its right - we've just removed the head so the head is now
    // the next element in the list - I hate linked lists!!!
    elem = cache->head;
  }

  // Free the hash table with no data freeFunc
  StringHash_free(cache->hash, NULL);

  cache->hash = StringHash_new(STRINGHASH_MEDIUM);

  cache->curSize = 0;
  
  return;
}

void LRUCache_remove(LRUCache *cache, char *key) {
  if (StringHash_contains(cache->hash, key)) {
    LRUCacheElement *cacheElem = StringHash_getValue(cache->hash, key);

    // Remove from hash table
    StringHash_remove(cache->hash, key, NULL);
  
    // Remove from linked list
    if (cacheElem == cache->head && cache->head == cache->tail) {
      // Its the only entry in cache - only need to unset head and tail, no list rearranging needed
      cache->head = NULL;
      cache->tail = NULL;
    } else {
      if (cacheElem == cache->head) {
        cache->head = cache->head->next;
      }
      if (cacheElem == cache->tail) {
        cache->tail = cache->tail->prev;
      }
      if (cacheElem->prev != NULL) {
        cacheElem->prev->next = cacheElem->next;
      }
      if (cacheElem->next != NULL) {
        cacheElem->next->prev = cacheElem->prev;
      }
    }
  
    // Reduce size
    cache->curSize -= cacheElem->size;

    // Free element
    LRUCacheElement_free(cacheElem);
    
  } else {
    fprintf(stderr,"Tried to remove nonexistent entry from LRUCache: key %s\n", key);
  }
}

void LRUCacheElement_free(LRUCacheElement *ce) {
  if (ce->freeFunc) {
    ce->freeFunc(ce->val);
  }
  free(ce->key);
  free(ce);

  return;
}

LRUCacheElement *LRUCacheElement_new(char *key, void *data, LRUCache_FreeFunc freeFunc, int size) {
  LRUCacheElement *cacheElem;

  if ((cacheElem = (LRUCacheElement *)calloc(1,sizeof(LRUCacheElement))) == NULL) {
    fprintf(stderr, "ERROR: Failed adding to cache\n");
    return NULL;
  }
  StrUtil_copyString(&(cacheElem->key), key, 0);
  cacheElem->val = data;
  cacheElem->freeFunc = freeFunc;
  cacheElem->size = size;

  return cacheElem;
}
