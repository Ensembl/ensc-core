#include "SeqUtil.h"

#include "BaseTest.h"

int main(int argc, char *argv[]) {
  char transTab[4][4][4];

  initEnsC();
  SeqUtil_readTransTab("../data/trans0.txt",transTab);
}
