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

#include "MapperPairSet.h"
#include <stdio.h>
#include <stdlib.h>

MapperPairSet *MapperPairSet_new() {
  MapperPairSet *mps;

  if ((mps = (MapperPairSet *)calloc(1,sizeof(MapperPairSet))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for mps\n");
    exit(1);
  }

  return mps;
}

void MapperPairSet_addPair(MapperPairSet *mps, MapperPair *pair) {
  if (mps->nPair == 0 || !(mps->nPair%10)) {
    if ((mps->pairs = (MapperPair **)realloc(mps->pairs, 
                                               (mps->nPair+10)*sizeof(MapperPair *))) == NULL) {
      fprintf(stderr,"ERROR: Failed allocating MapperPairSet pairs\n");
      exit(1);
    }
  } 
  mps->pairs[mps->nPair++] = pair;
}

Vector *MapperPairSet_getIds(MapperPairSet *mps, int ind) {
  Vector *idVector = NULL;

  if (mps)
    {
      idVector = Vector_new();

      int i;
      IDType *idP;
      for (i=0; i<mps->nPair; i++) {
        if ((idP = (IDType *)calloc(1,sizeof(IDType))) == NULL) {
          fprintf(stderr,"ERROR: Failed allocating space for id\n");
          exit(1);
        }
        *idP = MapperPair_getUnit(mps->pairs[i], ind)->id;

        Vector_addElement(idVector, idP);
      }
    }

  return idVector;
}

Vector *MapperPairSet_getFromIds(MapperPairSet *mps) {

  return MapperPairSet_getIds(mps, MAPPER_FROM_IND);
}

Vector *MapperPairSet_getToIds(MapperPairSet *mps) {

  return MapperPairSet_getIds(mps, MAPPER_TO_IND);
}

MapperPair *MapperPairSet_removePairAt(MapperPairSet *mps, int ind) {
  MapperPair *removed;

  if (ind < 0) {
    fprintf(stderr,"ERROR: Invalid pair element index %d\n",ind);
    exit(1);
  } else if (ind >= mps->nPair) {
    fprintf(stderr,"ERROR: Invalid pair element index %d\n",ind);
    exit(1);
  }

 removed = mps->pairs[ind];

  int i;
  for (i=ind+1; i<mps->nPair; i++) {
    mps->pairs[i-1] = mps->pairs[i];
  }

  mps->nPair--;

  return removed;
}

static int sortInd;

int MapperPairCompFunc(const void *a, const void *b) {
  MapperPair **mp1 = (MapperPair **)a;
  MapperPair **mp2 = (MapperPair **)b;

  if (MapperPair_getUnit(*mp1,sortInd)->start >
      MapperPair_getUnit(*mp2,sortInd)->start) {
    return 1;
  } else if (MapperPair_getUnit(*mp1,sortInd)->start <
             MapperPair_getUnit(*mp2,sortInd)->start) {
    return -1;
  } else {
    return 0;
  }
}

void MapperPairSet_sort(MapperPairSet *mps, int sInd) {
  sortInd = sInd;
  qsort(mps->pairs,(size_t)mps->nPair,sizeof(MapperPair *),MapperPairCompFunc);
}

void MapperPairSet_free(MapperPairSet *mps, int freePairsFlag) {
  int i;
 
  if (freePairsFlag) {
    for (i=0;i<mps->nPair;i++) {
      MapperPair_free(mps->pairs[i]);
    }
  }
  free(mps->pairs);
  free(mps);
}
