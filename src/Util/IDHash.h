#ifndef __IDHASH_H__
#define __IDHASH_H__

#include "EnsC.h"

typedef enum IDHashSizesEnum {
  IDHASH_SMALL,
  IDHASH_MEDIUM,
  IDHASH_LARGE
} IDHashSizes;

typedef struct IDKeyValuePairStruct {
  int64 key;
  void *value;
} IDKeyValuePair;

typedef struct IDHashStruct {
  IDKeyValuePair **buckets;
  int  *bucketCounts;
  int   size;
  int   nValue;
} IDHash;


IDHash *IDHash_new(IDHashSizes size);
int     IDHash_add(IDHash *idHash, int64 id, void *val);
int     IDHash_contains(IDHash *idHash, int64 id);
void    IDHash_free(IDHash *idHash, int freeFunc());
int     IDHash_getNumValues(IDHash *idHash);
int64 * IDHash_getKeys(IDHash *idHash);
void *  IDHash_getValue(IDHash *idHash, int64 id);
void *  IDHash_getValues(IDHash *idHash);

#endif
