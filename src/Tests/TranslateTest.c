#include "translate.h"

#include "BaseTest.h"

int main(int argc, char *argv[]) {
  char frm[6][200];
  int lengths[6];

  translate("ATGATGATGTGA",frm,lengths);
  printf("frm 1 = %s\n", frm[0]);
}
