#include "StrUtil.h"

#include "BaseTest.h"

char *testStr1 = "the cat sat on the mat";
char *testStr2 = "pog sog grog hog";

int main(int argc, char *argv[]) {
  int count;
  char token[MAXSTRLEN];
  char *chP;

  chP = testStr1;

  StrUtil_gettok(token,&chP,chP,MAXSTRLEN);
  ok(1, !strcmp(token,"the"));
}
