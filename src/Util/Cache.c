#include "Cache.h"

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

int CacheElement_free(CacheElement *ce) {
  if (ce->freeFunc) {
    ce->freeFunc(ce->val);
  }
  free(ce->key);
  free(ce);

  return 1;
}

CacheElement *CacheElement_new(char *key, void *data, int freeFunc()) {
  CacheElement *cacheElem;

  if ((cacheElem = (CacheElement *)calloc(1,sizeof(CacheElement))) == NULL) {
    fprintf(stderr, "ERROR: Failed adding to cache\n");
    return NULL;
  }
  StrUtil_copyString(&(cacheElem->key),key,0);
  cacheElem->val = data;
  cacheElem->freeFunc = freeFunc;

  return cacheElem;
}
