#ifndef __BASETEST_H__
#define __BASETEST_H__

/* Unusually for me I'm including the complete code so the tests compile easier */

void failAndDie(char *message) {
  fprintf(stderr,message);
  exit(1);
}

void ok(int testNum, int isOk) {
  fprintf(stderr, "Test %d: ",testNum);
  if (isOk) {
    fprintf(stderr, "OK\n");
  } else {
    fprintf(stderr, "Failed\n");
  }
}

#endif 
