#include "MetaContainer.h"
#include "StrUtil.h"
#include "BaseAdaptor.h"
#include "MysqlUtil.h"

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
  MYSQL_RES *results = mc->prepare((BaseAdaptor *)mc,qStr,strlen(qStr));
  MYSQL_ROW row = mysql_fetch_row(results);

  if( row ) {
    mysql_free_result(results);
    return MysqlUtil_getString(row,0);
  } else {
    return NULL;
  }
}
