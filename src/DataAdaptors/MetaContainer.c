#include "MetaContainer.h"
#include "StrUtil.h"
#include "BaseAdaptor.h"
#include "MysqlUtil.h"
#include "StatementHandle.h"

MetaContainer *MetaContainer_new(DBAdaptor *dba) {
  MetaContainer *mc;

  if ((mc = (MetaContainer *)calloc(1,sizeof(MetaContainer))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for MetaContainer\n");
    return NULL;
  }
  BaseAdaptor_init((BaseAdaptor *)mc, dba, META_CONTAINER);

  return mc;
}

char *MetaContainer_getDefaultAssembly(MetaContainer *mc) {
  char *qStr = "SELECT meta_value from meta where meta_key = 'assembly.default'";
  StatementHandle *sth = mc->prepare((BaseAdaptor *)mc,qStr,strlen(qStr));
  sth->execute(sth);
  ResultRow *row = sth->fetchRow(sth);
  char *assStr;

  if (row) {
    assStr = row->getStringAt(row,0);
    sth->finish(sth);
    return assStr;
  } else {
    sth->finish(sth);
    return NULL;
  }
}
