#ifndef __EXONADAPTOR_H__
#define __EXONADAPTOR_H__

#include "BaseAdaptor.h"
#include "AdaptorTypes.h"
#include "Exon.h"

struct ExonAdaptorStruct {
  BASEADAPTOR_DATA
};

ExonAdaptor *ExonAdaptor_new(DBAdaptor *dba);
Exon *ExonAdaptor_fetchByDbID(ExonAdaptor *ea, long dbID);
Exon *ExonAdaptor_exonFromResults(ExonAdaptor *ea, MYSQL_RES *results, MYSQL_ROW row);
Exon *ExonAdaptor_exonFromRow(ExonAdaptor *ea, MYSQL_ROW row);

#endif
