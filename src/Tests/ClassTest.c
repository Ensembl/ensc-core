#include <stdio.h>

#include "Class.h"

#include "BaseTest.h"

int main(int argc, char **argv) {
  ok(1, Class_isDescendent(CLASS_OBJECT,CLASS_FEATUREPAIR));
  ok(2,!Class_isDescendent(CLASS_FEATURESET,CLASS_OBJECT));
  ok(3, Class_isDescendent(CLASS_ENSROOT,CLASS_GENE));
  ok(4,!Class_isDescendent(CLASS_STATEMENTHANDLE,CLASS_GENE));

  return 0;
}

