/*
 * Copyright [1999-2017] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

#include "translate.h"

#include "BaseTest.h"

int main(int argc, char *argv[]) {
  int testResult = 0;
  char *frm[6];
  int lengths[6];
  int i;
  char *results[6] = {
    "MMMM",
    "***X",
    "DDD",
    "HHHH",
    "XSSS",
    "III",
  };

  for (i=0;i<6;i++) {
    frm[i]=malloc(200);
  }

  translate("ATGATGATGATG",frm,lengths,1, sizeof("ATGATGATGATG"));
  for (i=0;i<6;i++) {
    testResult += ok(i+1, !strcmp(frm[i], results[i]));
  }

  return testResult;
}
