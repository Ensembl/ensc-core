#ifndef __IDHASH_H__
#define __IDHASH_H__

#include "EnsC.h"

typedef enum IDHashSizesEnum {
  IDHASH_SMALL,
  IDHASH_MEDIUM,
  IDHASH_LARGE
} IDHashSizes;

typedef struct IDKeyValuePairStruct {
  IDType key;
  void *value;
} IDKeyValuePair;

typedef struct IDHashStruct {
  IDKeyValuePair **buckets;
  int  *bucketCounts;
  int   size;
  int   nValue;
} IDHash;


IDHash *IDHash_new(IDHashSizes size);
int     IDHash_add(IDHash *idHash, IDType id, void *val);
int     IDHash_contains(IDHash *idHash, IDType id);
void    IDHash_free(IDHash *idHash, int freeFunc());
int     IDHash_getNumValues(IDHash *idHash);
IDType * IDHash_getKeys(IDHash *idHash);
void *  IDHash_getValue(IDHash *idHash, IDType id);
void *  IDHash_getValues(IDHash *idHash);

#endif
