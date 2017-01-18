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

#include "Cache.h"

#include "BaseTest.h"

char *testStr1 = "the cat sat on the mat";
char *testStr2 = "pog sog grog hog";
char *testStr3 = "so long, goodbye, farewell";

int main(int argc, char *argv[]) {
  Cache *cache = Cache_new(2);
  char *retStr;
  char *retStr2;

  ok(1, cache!=NULL);

  Cache_addElement(cache,"fred",testStr1,NULL);
  
  retStr = (char *)Cache_findElem(cache,"fred");
  ok(2,retStr!=NULL); 
  ok(3,!strcmp(retStr,testStr1)); 

  
  Cache_addElement(cache,"bloggs",testStr2,NULL);
  Cache_addElement(cache,"micky",testStr3,NULL);

  retStr = Cache_findElem(cache,"fred");
  ok(4,retStr==NULL); 

  retStr = (char *)Cache_findElem(cache,"bloggs");
  retStr2 = (char *)Cache_findElem(cache,"micky");
  ok(5,retStr && retStr2); 
  ok(6,!strcmp(retStr,testStr2) && !strcmp(retStr2,testStr3));

  return 0;
}
