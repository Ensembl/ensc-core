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

#include "Vector.h"

#include "BaseTest.h"
#include "EnsC.h"

int main(int argc, char *argv[]) {
  int testResult = 0;
  Vector *v1;
  Vector *v2;
  char *str;

  initEnsC(argc, argv);

  v1 = Vector_new();

  Vector_addElement(v1,"a");
  Vector_addElement(v1,"b");
  Vector_addElement(v1,"c");

  testResult += ok(1,Vector_getNumElement(v1) == 3);

  v2 = Vector_new();

  Vector_addElement(v2,"d");
  Vector_addElement(v2,"e");
  Vector_addElement(v2,"f");
  Vector_addElement(v2,"g");
  
  testResult += ok(2,Vector_getNumElement(v2) == 4);

  Vector_append(v1,v2);

  testResult += ok(3,Vector_getNumElement(v1) == 7);

  Vector_reverse(v1);

  str = Vector_getElementAt(v1,0);

  testResult += ok(4, !strcmp(str,"g"));

  return testResult;
}
