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

#include "StrUtil.h"

#include "BaseTest.h"

char *testStr1 = "the cat sat on the mat";
char *testStr2 = "pog sog grog hog";

int main(int argc, char *argv[]) {
  int testResult = 0;
  int count;
  char token[MAXSTRLEN];
  char *chP;

  chP = testStr1;

  StrUtil_gettok(token,&chP,chP,MAXSTRLEN);
  testResult += ok(1, !strcmp(token,"the"));

  return testResult;
}
