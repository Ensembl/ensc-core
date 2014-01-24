/*
 * Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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
#include "StrUtil.h"
#include <stdlib.h>

#include "BaseTest.h"
#include "EnsC.h"

#define cigar1 "6M3IIMDM" 

int main(int argc, char *argv[]) {
  Vector *pieces;
  char *str;
  char *reverse;

  initEnsC(argc, argv);

  pieces = CigarStrUtil_getPieces(cigar1);

  ok(1, Vector_getNumElement(pieces) == 6); 
  
  str = Vector_getElementAt(pieces,0);

  ok(2, !strcmp(str,"6M"));

  Vector_free(pieces);

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
