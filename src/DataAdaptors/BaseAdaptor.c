#define __MAIN_C__
#include "BaseAdaptor.h"
#undef __MAIN_C__

void BaseAdaptor_init(BaseAdaptor *ba, DBAdaptor *dba, int adaptorType) {
  ba->dba = dba;
  ba->adaptorType = adaptorType;
  ba->prepare = BaseAdaptor_prepare;
}

StatementHandle *BaseAdaptor_prepare(BaseAdaptor *ba, char *qStr, size_t len) {
  /*printf("Query = %s len = %d\n",qStr,len);*/
  return DBAdaptor_prepare(ba->dba,qStr,len);
}
