#include <stdio.h>

#include "ComparaDBAdaptor.h"
#include "EnsC.h"

#include "BaseTest.h"

int main(int argc, char *argv[]) {
  ComparaDBAdaptor *cdba;

  initEnsC();

  cdba = ComparaDBAdaptor_new("kaka.sanger.ac.uk","anonymous",NULL,"ensembl_compara_14_1",3306,"Tests/Compara.conf");


  return 0;
}
