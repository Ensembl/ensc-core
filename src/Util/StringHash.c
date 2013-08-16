#include <stdio.h>
#include <stdlib.h>

#include "StringHash.h"
#include "StrUtil.h"
#include "Error.h"

#include "EnsC.h"
#include "ProcUtil.h"


unsigned int StringHash_getBucketNum(StringHash *stringHash, char *key, int keyLen) {
  unsigned int hash, i;
  unsigned int hashCode;

  for (hash=keyLen, i=0; i<keyLen; ++i) {
    hash = (hash<<5)^(hash>>27)^key[i];
  }
  hashCode = (hash%stringHash->size);

  return hashCode;
}

StringHash *StringHash_new(StringHashSizes size) {
  StringHash *stringHash;

  if ((stringHash = (StringHash *)calloc(1,sizeof(StringHash))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for stringHash\n");
    return NULL;
  }

  switch (size) {
    case STRINGHASH_SMALL:
      stringHash->size = 257; 
      break;
    case STRINGHASH_LARGE:
      stringHash->size = 104711; 
      break;
    case STRINGHASH_HUGE:
      stringHash->size = 1000039; 
      break;
    case STRINGHASH_MEDIUM:
    default:
      stringHash->size = 32353; 
  }

  if ((stringHash->buckets = (KeyValuePair **)calloc(stringHash->size,sizeof(KeyValuePair *))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for stringHash->buckets\n");
    return NULL;
  }

  if ((stringHash->bucketCounts = (int *)calloc(stringHash->size,sizeof(int))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for stringHash->bucketCounts\n");
    return NULL;
  }

  return stringHash;
}

char **StringHash_getKeys(StringHash *stringHash) {
  int i;
  int j;
  char **keys;
  int keyCnt = 0;

  if (!stringHash->nValue) {
    fprintf(stderr, "Warning: tried to get Keys for empty StringHash\n");
    return NULL;
  }

  if ((keys = (char **)calloc(stringHash->nValue,sizeof(char *))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for keys\n");
    return NULL;
  }

  for (i=0; i<stringHash->size; i++) {
    if (stringHash->bucketCounts[i]) {
      for (j=0; j<stringHash->bucketCounts[i]; j++) {
        StrUtil_copyString(&(keys[keyCnt++]),stringHash->buckets[i][j].key,0);
      }
    }
  }
  if (keyCnt != stringHash->nValue) {
    fprintf(stderr,"ERROR: Internal StringHash error - keyCnt != stringHash->nValue\n");
    //ProcUtil_showBacktrace(EnsC_progName);
  }
  return keys;
}

char **StringHash_getKeysNoCopy(StringHash *stringHash) {
  int i;
  int j;
  char **keys;
  int keyCnt = 0;

  if (!stringHash->nValue) {
    fprintf(stderr, "Warning: tried to get Keys for empty StringHash\n");
    return NULL;
  }

  if ((keys = (char **)calloc(stringHash->nValue,sizeof(char *))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for keys\n");
    return NULL;
  }

  for (i=0; i<stringHash->size; i++) {
    if (stringHash->bucketCounts[i]) {
      for (j=0; j<stringHash->bucketCounts[i]; j++) {
        keys[keyCnt++] = stringHash->buckets[i][j].key;
      }
    }
  }
  if (keyCnt != stringHash->nValue) {
    fprintf(stderr,"ERROR: Internal StringHash error - keyCnt != stringHash->nValue\n");
    //ProcUtil_showBacktrace(EnsC_progName);
  }
  return keys;
}

/*
int StringHash_getNumValues(StringHash *stringHash) {
  return stringHash->nValue;
}
*/

void *StringHash_getValues(StringHash *stringHash) {
  int i;
  int j;
  void **values;
  int valCnt = 0;
  
  if (!stringHash->nValue) {
    return NULL;
  }

  if ((values = (void **)calloc(stringHash->nValue,sizeof(void *))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for values\n");
    return NULL;
  }

  for (i=0; i<stringHash->size; i++) {
    if (stringHash->bucketCounts[i]) {
      for (j=0; j<stringHash->bucketCounts[i]; j++) {
        values[valCnt++] = stringHash->buckets[i][j].value;
      }
    }
  }
  if (valCnt != stringHash->nValue) {
    fprintf(stderr,"ERROR: Internal StringHash error - valCnt != stringHash->nValue\n");
  }
  return values;
}

void StringHash_printHashStats(StringHash *stringHash, char *hashDesc) {
  int nOccupiedBuckets = 0;
  int nMultipleOccupied = 0;
  int maxInBucket = 0;
  int nValues = 0;

  int i;
  for (i=0; i<stringHash->size; i++) {
    if (stringHash->bucketCounts[i]) {
      nOccupiedBuckets++;
      if (stringHash->bucketCounts[i] > 1) {
        nMultipleOccupied++;
      }
      if (stringHash->bucketCounts[i] > maxInBucket) {
        maxInBucket = stringHash->bucketCounts[i];
      }
      nValues += stringHash->bucketCounts[i];
    }
  }

  fprintf(stderr, "\nString Hash statistics for %s hash\n", hashDesc);
  fprintf(stderr, "   Size of hash (number of buckets)    = %d\n", stringHash->size);
  fprintf(stderr, "   Number of values stored             = %d\n", nValues);
  fprintf(stderr, "   Load factor                         = %f\n", (float)nValues/(float)stringHash->size);
  fprintf(stderr, "   Number of buckets used              = %d\n", nOccupiedBuckets);
  fprintf(stderr, "   Number of multiply occupied buckets = %d\n", nMultipleOccupied);
  fprintf(stderr, "   Max num value in a single bucket    = %d\n\n", maxInBucket);
}

void *StringHash_getValue(StringHash *stringHash, char *key) {
  int keyLen = strlen(key);
  int bucketNum = StringHash_getBucketNum(stringHash,key,keyLen);
  int i;

  for (i=0; i<stringHash->bucketCounts[bucketNum]; i++) {
    
    if (keyLen == stringHash->buckets[bucketNum][i].keyLen) {
      //fprintf(stderr,"equal keylen %d %d\n",keyLen, stringHash->buckets[bucketNum][i].keyLen);
      if (!memcmp(key, stringHash->buckets[bucketNum][i].key,keyLen)) {
      //if (!strcmp(key, stringHash->buckets[bucketNum][i].key)) {
        return stringHash->buckets[bucketNum][i].value;
      }
    }
  }

//  fprintf(stderr,"ERROR: Didn't find key %s in StringHash\n",key);
  //fprintf(stderr,"keylen %d\n",keyLen);
  //for (i=0; i<stringHash->bucketCounts[bucketNum]; i++) {
  //  fprintf(stderr,"bucket key %s keylen %d\n",stringHash->buckets[bucketNum][i].key, stringHash->buckets[bucketNum][i].keyLen);
  //}
  return NULL;
}

int StringHash_contains(StringHash *stringHash, char *key) {
  int keyLen = strlen(key);
  int bucketNum = StringHash_getBucketNum(stringHash,key, keyLen);
  int i;

  for (i=0; i<stringHash->bucketCounts[bucketNum]; i++) {
    //if (!strcmp(key,stringHash->buckets[bucketNum][i].key)) {
    if (keyLen == stringHash->buckets[bucketNum][i].keyLen) {
      if (!memcmp(key, stringHash->buckets[bucketNum][i].key,keyLen)) {
        return 1;
      }
    }
  }

  return 0;
}

char *StringHash_getKey(StringHash *stringHash, char *key) {
  int keyLen = strlen(key);
  int bucketNum = StringHash_getBucketNum(stringHash,key, keyLen);
  int i;

  for (i=0; i<stringHash->bucketCounts[bucketNum]; i++) {
    //if (!strcmp(key,stringHash->buckets[bucketNum][i].key)) {
    if (keyLen == stringHash->buckets[bucketNum][i].keyLen) {
      if (!memcmp(key, stringHash->buckets[bucketNum][i].key,keyLen)) {
        return stringHash->buckets[bucketNum][i].key;
      }
    }
  }

  return NULL;
}

int StringHash_add(StringHash *stringHash, char *key, void *val) {
  int keyLen = strlen(key);
  int bucketNum = StringHash_getBucketNum(stringHash,key,keyLen);
  int count = stringHash->bucketCounts[bucketNum];

  if (StringHash_contains(stringHash,key)) {
    int i;

    fprintf(stderr,"WARNING: Duplicate key %s - value will be overwritten\n",key);
    for (i=0; i<stringHash->bucketCounts[bucketNum]; i++) {
      if (!strcmp(key,stringHash->buckets[bucketNum][i].key)) {
        stringHash->buckets[bucketNum][i].value = val;
        return 1;
      }
    }
    
  } else {
    if (!count || !(count%10)) {
      if ((stringHash->buckets[bucketNum] = 
           (KeyValuePair *)realloc(stringHash->buckets[bucketNum],
                                   (count+10) * sizeof(KeyValuePair))) == NULL) {
        fprintf(stderr,"ERROR: Failed allocating space for stringHash bucket\n");
        return 0;
      }
    }
    
    StrUtil_copyString(&(stringHash->buckets[bucketNum][count].key), key, 0);
    if (!stringHash->buckets[bucketNum][count].key) {
      Error_trace("StringHash_add",NULL);
      return 0;
    }
    stringHash->buckets[bucketNum][count].value = val;
    stringHash->buckets[bucketNum][count].keyLen = keyLen;

    stringHash->nValue++;
    stringHash->bucketCounts[bucketNum]++;
  }

  return 1; 
}

void StringHash_free(StringHash *stringHash, void freeFunc()) {
  int i;
  int j;
  
  for (i=0; i<stringHash->size; i++) {
    if (stringHash->bucketCounts[i]) {
      for (j=0; j<stringHash->bucketCounts[i]; j++) {
        free(stringHash->buckets[i][j].key);
        if (freeFunc) {
          freeFunc(stringHash->buckets[i][j].value);
        }
      }
      free(stringHash->buckets[i]);
    }
  }
  
  free(stringHash->buckets);
  free(stringHash->bucketCounts);
  free(stringHash);
}

void StringHash_freeNoValFree(StringHash *stringHash) {
  int i;
  int j;
  
  for (i=0; i<stringHash->size; i++) {
    if (stringHash->bucketCounts[i]) {
      for (j=0; j<stringHash->bucketCounts[i]; j++) {
        free(stringHash->buckets[i][j].key);
      }
      free(stringHash->buckets[i]);
    }
  }
  
  free(stringHash->buckets);
  free(stringHash->bucketCounts);
  free(stringHash);
}

int StringHash_remove(StringHash *stringHash, char *key, void freeFunc()) {
  int keyLen = strlen(key);
  int bucketNum = StringHash_getBucketNum(stringHash, key, keyLen);
  int i;
  int cnt = 0;
  char *toRemoveKey = NULL;
  void *toRemoveVal = NULL;

  for (i=0; i<stringHash->bucketCounts[bucketNum]; i++) {
    if (!strcmp(stringHash->buckets[bucketNum][i].key, key)) {
      if (toRemoveVal) {
        fprintf(stderr,"Error: %s found more than once in hash\n", key);
        exit(1);
      }
      toRemoveKey = stringHash->buckets[bucketNum][i].key;
      toRemoveVal = stringHash->buckets[bucketNum][i].value;
    } else {
      stringHash->buckets[bucketNum][cnt].key   = stringHash->buckets[bucketNum][i].key;
      stringHash->buckets[bucketNum][cnt].value = stringHash->buckets[bucketNum][i].value;
      stringHash->buckets[bucketNum][cnt].keyLen = stringHash->buckets[bucketNum][i].keyLen;
      cnt++;
    }
  }

  if (toRemoveVal) {
    stringHash->bucketCounts[bucketNum]--;

    if (freeFunc) {
      freeFunc(toRemoveVal);
    }
    free(toRemoveKey);
    if (!stringHash->bucketCounts[bucketNum]) {
       free(stringHash->buckets[bucketNum]);
       stringHash->buckets[bucketNum] = NULL;
    }
  }

  stringHash->nValue--;
  
  return 0;
}
