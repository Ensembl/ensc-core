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

#include "LRUCache.h"

#include "BaseTest.h"

char *testStr1 = "the cat sat on the mat";
char *testStr2 = "pog sog grog hog";
char *testStr3 = "so long, goodbye, farewell";

int main(int argc, char *argv[]) {
  LRUCache *cache = LRUCache_new(2);
  char *retStr;
  char *retStr2;

  ok(1, cache!=NULL);

  LRUCache_put(cache,"fred",testStr1,NULL,1);
  
  retStr = (char *)LRUCache_get(cache,"fred");
  ok(2,retStr!=NULL); 
  ok(3,!strcmp(retStr,testStr1)); 

  LRUCache_put(cache,"bloggs",testStr2,NULL,1);
  LRUCache_put(cache,"micky",testStr3,NULL,1);

  retStr = LRUCache_get(cache,"fred");
  ok(4,retStr==NULL); 

  retStr = (char *)LRUCache_get(cache,"bloggs");
  retStr2 = (char *)LRUCache_get(cache,"micky");
  ok(5,retStr && retStr2); 
  ok(6,!strcmp(retStr,testStr2) && !strcmp(retStr2,testStr3));

  // Now access bloggs so its most recently accessed, and then add fred 
  retStr = (char *)LRUCache_get(cache,"bloggs");
  LRUCache_put(cache,"fred",testStr1,NULL,1);
  // Now check that micky is the one that was deleted
  ok(7, ! LRUCache_contains(cache, "micky") && 
          LRUCache_contains(cache, "bloggs") && 
          LRUCache_contains(cache, "fred"));

  LRUCache_empty(cache);
  ok(8, cache->curSize == 0);

  LRUCache_put(cache,"fred",testStr1,NULL,1);
  ok(9, cache->curSize == 1);
  retStr = LRUCache_get(cache,"fred");
  ok(10,!strcmp(retStr,testStr1));

  // One with size > 1
  LRUCache_put(cache,"bloggs",testStr2,NULL,2);
  // Because bloggs is size 2 it should have removed fred
  ok(11, ! LRUCache_contains(cache, "fred") && 
           LRUCache_contains(cache, "bloggs"));


  

  return 0;
}
