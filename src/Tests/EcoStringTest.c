#include <stdio.h>

#include "EcoString.h"

#include "BaseTest.h"

#define TESTSTR1 "TestStr"

int main(int argc, char **argv) {
  ECOSTRTABLE *table;
  ECOSTRING    ecoS1;
  ECOSTRING    ecoS2;

  if (!EcoString_initTable(&table)) {
    failAndDie("Couldn't create table\n");
  }
  ok(1, table != NULL);

  EcoString_copyStr(table,&ecoS1, TESTSTR1, 0);

  ok(2, ecoS1!=NULL);
  ok(3, !strcmp(ecoS1,TESTSTR1));
  
  EcoString_copyStr(table,&ecoS2, TESTSTR1, 0);
  ok(4, ecoS1 == ecoS2);
  ok(5, !EcoString_strcmp(ecoS1,ecoS2));
   
}

