/*
 * Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
 * Copyright [2016-2017] EMBL-European Bioinformatics Institute
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

#include <unistd.h>

#include "Stream.h"
#include "FileUtil.h"

#include "BaseTest.h"

char *TestStr1 = "Test string 123";

int main(int argc, char **argv) {
  FILE *fp;
  char str[1024];

  unlink("Hist.log");
  unlink("out.out");

  Stream_setDefaults(1);

  if (!FileUtil_exists("Hist.log")) {
    Test_failAndDie("Failed setting up streams\n");
  }
  ok(1,1);

  Stream_setFile("OutStream","out.out");
  fp = FileUtil_open("out.out","r","main");
  ok(2,fp!=NULL);

  Stream_fprintf(OutStream,TestStr1);
  Stream_flush(OutStream);

  FileUtil_getStrippedLine(str,1024,fp);
  ok(3,!strcmp(str,TestStr1));

  return 0;
}

