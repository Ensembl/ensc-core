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

#ifndef __BASETEST_H__
#define __BASETEST_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Unusually for me I'm including the complete code so the tests compile easier */

void Test_failAndDie(char *message) {
  fprintf(stderr,message);
  exit(1);
}

int ok(int testNum, int isOk) {
  fprintf(stderr, "Test %d: ",testNum);
  if (isOk) {
    fprintf(stderr, "OK\n");
    return 0;
  } else {
    fprintf(stderr, "Failed\n");
    return 1;
  }
}

#endif 
