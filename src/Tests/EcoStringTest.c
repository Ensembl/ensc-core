#include <stdio.h>

#include "EcoString.h"

#include "BaseTest.h"

char *TestStr1 = "TestStr";

int main(int argc, char **argv) {
  ECOSTRTABLE *table;
  ECOSTRING    ecoS1;
  ECOSTRING    ecoS2;
  ECOSTRING    ecoS3;

  if (!EcoString_initTable(&table)) {
    failAndDie("Couldn't create table\n");
  }
  ok(1, table != NULL);

  EcoString_copyStr(table,&ecoS1, TestStr1, 0);

  ok(2, ecoS1!=NULL);
  ok(3, !strcmp(ecoS1,TestStr1));
  
  EcoString_copyStr(table,&ecoS2, TestStr1, 0);
  ok(4, ecoS1 == ecoS2);
  ok(5, !EcoString_strcmp(ecoS1,ecoS2));

  EcoString_freeStr(table,ecoS1);
  EcoString_freeStr(table,ecoS2);

  EcoString_copyStr(table,&ecoS3, &(TestStr1[1]), 0);

  EcoString_copyStr(table,&ecoS1, TestStr1, 0);

/* Should have been reallocated */
  ok(6, EcoString_strcmp(ecoS1,ecoS2));

  ok(7, !strcmp(ecoS3,&(TestStr1[1])));

  ok(8, EcoString_strcmp(ecoS1,ecoS3));
 
  return 0;
}

