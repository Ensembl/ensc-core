#ifndef __STRINGHASH_H__
#define __STRINGHASH_H__

#include <string.h>

typedef enum StringHashSizesEnum {
  STRINGHASH_SMALL,
  STRINGHASH_MEDIUM,
  STRINGHASH_LARGE
} StringHashSizes;

typedef struct KeyValuePairStruct {
  char *key;
  void *value;
} KeyValuePair;

typedef struct StringHashStruct {
  KeyValuePair **buckets;
  int  *bucketCounts;
  int   size;
  int   nValue;
} StringHash;


StringHash *StringHash_new(StringHashSizes size);
int     StringHash_add(StringHash *stringHash, char *string, void *val);
int     StringHash_contains(StringHash *stringHash, char *string);
void    StringHash_free(StringHash *stringHash, void freeFunc());
int     StringHash_getNumValues(StringHash *stringHash);
void *  StringHash_getValue(StringHash *stringHash, char *string);
void *  StringHash_getValues(StringHash *stringHash);
char ** StringHash_getKeys(StringHash *stringHash);


#endif
