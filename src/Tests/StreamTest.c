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

