#include "Vector.h"

#include "BaseTest.h"
#include "EnsC.h"

int main(int argc, char *argv[]) {
  Vector *v1;
  Vector *v2;
  char *str;

  initEnsC();

  v1 = Vector_new();

  Vector_addElement(v1,"a");
  Vector_addElement(v1,"b");
  Vector_addElement(v1,"c");

  ok(1,Vector_getNumElement(v1) == 3);

  v2 = Vector_new();

  Vector_addElement(v2,"d");
  Vector_addElement(v2,"e");
  Vector_addElement(v2,"f");
  Vector_addElement(v2,"g");
  
  ok(2,Vector_getNumElement(v2) == 4);

  Vector_append(v1,v2);

  ok(3,Vector_getNumElement(v1) == 7);

  Vector_reverse(v1);

  str = Vector_getElementAt(v1,0);

  ok(4, !strcmp(str,"g"));

  return 0;
}
