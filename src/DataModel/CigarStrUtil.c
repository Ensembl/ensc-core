/*
 * Copyright [1999-2024] EMBL-European Bioinformatics Institute
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

#include "CigarStrUtil.h"
#include <stdio.h>
#include "StrUtil.h"

char *CigarStrUtil_reverse(char *oldCigarString, int len) {
  char *newCigarString;
  char *oldChP;
  char *newChP;
  char *newPieceStartP;

 if ((newCigarString = (char *)calloc(len+1,sizeof(char))) == NULL) {
    fprintf(stderr, "Error: Failed allocating cigarString\n");
    exit(1);
  }

  oldChP = &(oldCigarString[len-2]);

  newPieceStartP = newChP = newCigarString;
  *newChP = *(oldChP+1);
  newChP++;

  while (oldChP >= oldCigarString-1) {
    if (*oldChP == 'M' || *oldChP == 'D' || *oldChP == 'I' || oldChP == oldCigarString-1) {
      *newChP = '\0';
  
      StrUtil_reverseString(newPieceStartP, newChP-newPieceStartP);
      newPieceStartP = newChP;
    }

    *newChP = *oldChP;
    oldChP--;
    newChP++;
  }


  // NIY does old cigar string need freeing?

  return newCigarString;
}

Vector *CigarStrUtil_getPieces(char *cigarString) {
  Vector *pieces = Vector_new();
  char *chP = cigarString;
  char piece[128];
  char *pieceP = piece;

  Vector_setFreeFunc(pieces,free);

  while (*chP != '\0') {
    if (*chP == 'M' || *chP == 'D' || *chP == 'I') {
      char *tmpP;

      *pieceP = *chP;
      pieceP++;
      *pieceP = '\0';
      tmpP = StrUtil_copyString(&tmpP, piece, 0);
      Vector_addElement(pieces, tmpP);

    } else {
      *pieceP = *chP;
      pieceP++;
    }
    chP++; 
  }
  return pieces;
}
