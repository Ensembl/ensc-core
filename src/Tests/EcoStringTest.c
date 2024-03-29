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

#include <stdio.h>

#include "EcoString.h"

#include "BaseTest.h"

char *TestStr1 = "TestStr";

int main(int argc, char **argv) {
  int testResult = 0;
  ECOSTRTABLE *table;
  ECOSTRING    ecoS1;
  ECOSTRING    ecoS2;
  ECOSTRING    ecoS3;

  if (!EcoString_initTable(&table)) {
    Test_failAndDie("Couldn't create table\n");
  }
  testResult += ok(1, table != NULL);

  EcoString_copyStr(table,&ecoS1, TestStr1, 0);

  testResult += ok(2, ecoS1!=NULL);
  testResult += ok(3, !strcmp(ecoS1,TestStr1));
  
  EcoString_copyStr(table,&ecoS2, TestStr1, 0);
  testResult += ok(4, ecoS1 == ecoS2);
  testResult += ok(5, !EcoString_strcmp(ecoS1,ecoS2));

  EcoString_freeStr(table,ecoS1);
  EcoString_freeStr(table,ecoS2);

  EcoString_copyStr(table,&ecoS3, &(TestStr1[1]), 0);

  EcoString_copyStr(table,&ecoS1, TestStr1, 0);

/* Should have been reallocated */
  testResult += ok(6, EcoString_strcmp(ecoS1,ecoS2));

  testResult += ok(7, !strcmp(ecoS3,&(TestStr1[1])));

  testResult += ok(8, EcoString_strcmp(ecoS1,ecoS3));
 
  return testResult;
}
