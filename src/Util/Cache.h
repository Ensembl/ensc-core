#ifndef __CACHE_H__
#define __CACHE_H__

#include <stdio.h>

typedef struct CacheStruct Cache;
typedef struct CacheElementStruct CacheElement;

typedef int (*Cache_FreeFunc)();

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
int CacheElement_free(CacheElement *ce);

#endif
