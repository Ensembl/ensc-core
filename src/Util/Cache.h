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
