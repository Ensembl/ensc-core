#include "CigarStrUtil.h"
#include "StrUtil.h"
#include <stdlib.h>

#include "BaseTest.h"
#include "EnsC.h"

#define cigar1 "6M3IIMDM" 

int main(int argc, char *argv[]) {
  Vector *pieces;
  char *str;
  char *reverse;

  initEnsC();

  pieces = CigarStrUtil_getPieces(cigar1);

  ok(1, Vector_getNumElement(pieces) == 6); 
  
  str = Vector_getElementAt(pieces,0);

  ok(2, !strcmp(str,"6M"));

  Vector_free(pieces,free);

  StrUtil_copyString(&str,cigar1,0);
  
  ok(3, !strcmp(str,cigar1));
  
  reverse = CigarStrUtil_reverse(str,strlen(str));
  free(str);

  pieces = CigarStrUtil_getPieces(reverse);

  ok(4, Vector_getNumElement(pieces) == 6); 

  str = Vector_getElementAt(pieces,0);

  ok(5, !strcmp(str,"M"));

  return 0;
}
