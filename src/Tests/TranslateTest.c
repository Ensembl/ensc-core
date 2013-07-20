#include "translate.h"

#include "BaseTest.h"

int main(int argc, char *argv[]) {
  char *frm[6];
  int lengths[6];
  int i;

  for (i=0;i<6;i++) {
    frm[i]=malloc(200);
  }

  translate("ATGATGATGATG",frm,lengths,1, sizeof("ATGATGATGATG"));
  for (i=0;i<6;i++) {
    printf("frm %d = %s\n", i+1, frm[i]);
  }

  return 0;
}
