#ifndef __IDHASH_H__
#define __IDHASH_H__

typedef enum IDHashSizesEnum {
  IDHASH_SMALL,
  IDHASH_MEDIUM,
  IDHASH_LARGE
} IDHashSizes;

typedef struct IDKeyValuePairStruct {
  long key;
  void *value;
} IDKeyValuePair;

typedef struct IDHashStruct {
  IDKeyValuePair **buckets;
  int  *bucketCounts;
  int   size;
  int   nValue;
} IDHash;


IDHash *IDHash_new(IDHashSizes size);
int     IDHash_add(IDHash *idHash, long id, void *val);
int     IDHash_contains(IDHash *idHash, long id);
void    IDHash_free(IDHash *idHash, int freeFunc());
int     IDHash_getNumValues(IDHash *idHash);
long *  IDHash_getKeys(IDHash *idHash);
void *  IDHash_getValue(IDHash *idHash, long id);
void *  IDHash_getValues(IDHash *idHash);

#endif
