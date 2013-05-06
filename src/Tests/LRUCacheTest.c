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

  

  return 0;
}
