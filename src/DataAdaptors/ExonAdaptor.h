#ifndef __EXONADAPTOR_H__
#define __EXONADAPTOR_H__

#include "BaseAdaptor.h"
#include "AdaptorTypes.h"
#include "Exon.h"

struct ExonAdaptorStruct {
  BASEADAPTOR_DATA
};

ExonAdaptor *ExonAdaptor_new(DBAdaptor *dba);
Exon *ExonAdaptor_fetchByDbID(ExonAdaptor *ea, IDType dbID);
Exon *ExonAdaptor_exonFromResults(ExonAdaptor *ea, StatementHandle *sth, ResultRow *row);
Exon *ExonAdaptor_exonFromRow(ExonAdaptor *ea, ResultRow *row);

#endif
