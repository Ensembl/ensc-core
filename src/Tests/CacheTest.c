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
