#include "CigarStrUtil.h"
#include <stdio.h>

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
  *newChP++;

  while (oldChP >= oldCigarString-1) {
    if (*oldChP == 'M' || *oldChP == 'D' || *oldChP == 'I' || oldChP == oldCigarString-1) {
      *newChP = '\0';
  
      StrUtil_reverseString(newPieceStartP, newChP-newPieceStartP);
      newPieceStartP = newChP;
    }

    *newChP = *oldChP;
    *oldChP--;
    *newChP++;
  }
  return newCigarString;
}
