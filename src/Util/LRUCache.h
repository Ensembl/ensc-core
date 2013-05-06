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

void LRUCacheElement_free(LRUCacheElement *ce);
LRUCacheElement *LRUCacheElement_new(char *key, void *data, LRUCache_FreeFunc freeFunc, int size);

#endif
