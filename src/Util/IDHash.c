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

#include <stdio.h>
#include "IDHash.h"
#include "Vector.h"

#define IDHash_getBucketNum(idHash,key) ((key)%(idHash)->size);

IDHash *IDHash_new(IDHashSizes size) {
  IDHash *idHash;

  if ((idHash = (IDHash *)calloc(1,sizeof(IDHash))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for idHash\n");
    return NULL;
  }

  switch (size) {
    case IDHASH_SMALL:
      idHash->size = 257; 
      break;
    case IDHASH_LARGE:
      idHash->size = 104711; 
      break;
    case IDHASH_MEDIUM:
    default:
      idHash->size = 32353; 
  }

  if ((idHash->buckets = (IDKeyValuePair **)calloc(idHash->size,sizeof(IDKeyValuePair *))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for idHash->buckets\n");
    return NULL;
  }

  if ((idHash->bucketCounts = (int *)calloc(idHash->size,sizeof(int))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for idHash->buckets\n");
    return NULL;
  }

  return idHash;
}

int IDHash_getNumValues(IDHash *idHash) {
  return idHash->nValue;
}

void **IDHash_getValues(IDHash *idHash) {
  int i;
  int j;
  void **values;
  int valCnt = 0;
  
  if (!idHash->nValue) {
    return NULL;
  }

  if ((values = (void **)calloc(idHash->nValue,sizeof(void *))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for values\n");
    return NULL;
  }

  for (i=0; i<idHash->size; i++) {
    if (idHash->bucketCounts[i]) {
      for (j=0; j<idHash->bucketCounts[i]; j++) {
        values[valCnt++] = idHash->buckets[i][j].value;
      }
    }
  }
  if (valCnt != idHash->nValue) {
    fprintf(stderr,"ERROR: Internal IDHash error - valCnt != idHash->nValue\n");
  }
  return values;
}

Vector *IDHash_getValuesVector(IDHash *idHash) {
  int i;
  int j;
  
  if (!idHash->nValue) {
    return NULL;
  }

  Vector *values = Vector_new();

  for (i=0; i<idHash->size; i++) {
    if (idHash->bucketCounts[i]) {
      for (j=0; j<idHash->bucketCounts[i]; j++) {
        Vector_addElement(values, idHash->buckets[i][j].value);
      }
    }
  }
  if (Vector_getNumElement(values) != idHash->nValue) {
    fprintf(stderr,"ERROR: Internal IDHash error - values vector size != idHash->nValue\n");
  }
  return values;
}

IDType *IDHash_getKeys(IDHash *idHash) {
  int i;
  int j;
  IDType *keys;
  int keyCnt = 0;
  
  if (!idHash->nValue) {
    return NULL;
  }

  if ((keys = (IDType *)calloc(idHash->nValue,sizeof(IDType))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for keys\n");
    return NULL;
  }

  for (i=0; i<idHash->size; i++) {
    if (idHash->bucketCounts[i]) {
      for (j=0; j<idHash->bucketCounts[i]; j++) {
        keys[keyCnt++] = idHash->buckets[i][j].key;
      }
    }
  }
  if (keyCnt != idHash->nValue) {
    fprintf(stderr,"ERROR: Internal IDHash error - keyCnt != idHash->nValue\n");
  }
  return keys;
}

void *IDHash_getValue(IDHash *idHash, IDType id) {
  int bucketNum = IDHash_getBucketNum(idHash,id);
  int i;

  for (i=0; i<idHash->bucketCounts[bucketNum]; i++) {
    if (id == idHash->buckets[bucketNum][i].key) {
      return idHash->buckets[bucketNum][i].value;
    }
  }

//  fprintf(stderr,"ERROR: Didn't find key " IDFMTSTR " in IDHash\n",id);
  return NULL;
}

int IDHash_contains(IDHash *idHash, IDType id) {
  int bucketNum = IDHash_getBucketNum(idHash,id);
  int i;

  for (i=0; i<idHash->bucketCounts[bucketNum]; i++) {
    if (id == idHash->buckets[bucketNum][i].key) {
      return 1;
    }
  }

  return 0;
}

int IDHash_remove(IDHash *idHash, IDType id, void freeFunc()) {
  int bucketNum = IDHash_getBucketNum(idHash,id);
  int i;
  int cnt = 0;
  void *toRemove = NULL;

  for (i=0; i<idHash->bucketCounts[bucketNum]; i++) {
    if (id == idHash->buckets[bucketNum][i].key) {
      if (toRemove) {
        fprintf(stderr,"Error: " IDFMTSTR " found more than once in hash\n",id);
        exit(1);
      }
      toRemove = idHash->buckets[bucketNum][i].value;
    } else {
      idHash->buckets[bucketNum][cnt].key   = idHash->buckets[bucketNum][i].key;
      idHash->buckets[bucketNum][cnt].value = idHash->buckets[bucketNum][i].value;
      cnt++;
    }
  }

  if (toRemove) {
    idHash->bucketCounts[bucketNum]--;
    if (freeFunc) {
      freeFunc(toRemove);
    }
  }
  idHash->nValue--;
  
  return 0;
}

int IDHash_add(IDHash *idHash, IDType id, void *val) {
  int bucketNum = IDHash_getBucketNum(idHash,id);
  int count = idHash->bucketCounts[bucketNum];

  if (IDHash_contains(idHash,id)) {
    int i;

    fprintf(stderr,"WARNING: Duplicate key " IDFMTSTR " - value will be overwritten\n",id);
    for (i=0; i<idHash->bucketCounts[bucketNum]; i++) {
      if (id == idHash->buckets[bucketNum][i].key) {
        idHash->buckets[bucketNum][i].value = val;
        return 1;
      }
    }
    
  } else {
    if (!count || !(count%10)) {
      if ((idHash->buckets[bucketNum] = 
           (IDKeyValuePair *)realloc(idHash->buckets[bucketNum],
                                   (count+10) * sizeof(IDKeyValuePair))) == NULL) {
        fprintf(stderr,"ERROR: Failed allocating space for idHash bucket\n");
        return 0;
      }
    }
    
    idHash->buckets[bucketNum][count].key   = id;
    idHash->buckets[bucketNum][count].value = val;

    idHash->nValue++;
    idHash->bucketCounts[bucketNum]++;
  }

  return 1; 
}

void IDHash_free(IDHash *idHash, void freeFunc()) {
  int i;
  int j;
  
  for (i=0; i<idHash->size; i++) {
    if (idHash->bucketCounts[i]) {
      if (freeFunc) {
        for (j=0; j<idHash->bucketCounts[i]; j++) {
          freeFunc(idHash->buckets[i][j].value);
        }
      }
      free(idHash->buckets[i]);
    }
  }
  
  free(idHash->buckets);
  free(idHash->bucketCounts);
  free(idHash);
}
