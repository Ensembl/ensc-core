#ifndef __ATTRIBUTEADAPTOR_H__
#define __ATTRIBUTEADAPTOR_H__

#include "BaseAdaptor.h"
#include "AdaptorTypes.h"
#include "DataModelTypes.h"

struct AttributeAdaptorStruct {
  BASEADAPTOR_DATA
};

AttributeAdaptor *AttributeAdaptor_new(DBAdaptor *dba);
Vector *AttributeAdaptor_fetchAllByTranscript(AttributeAdaptor *ata, Transcript *transcript, char *code);
Vector *AttributeAdaptor_fetchAllByGene(AttributeAdaptor *ata, Gene *gene, char *code);
Vector *AttributeAdaptor_fetchAllBySlice(AttributeAdaptor *ata, Slice *slice, char *code);
Vector *AttributeAdaptor_fetchAllByTranslation(AttributeAdaptor *ata, Translation *translation, char *code);

Vector *AttributeAdaptor_doFetchAllByTypeAndTableAndID(AttributeAdaptor *ata, char *type, char *table, IDType objectId, char *code);
Vector *AttributeAdaptor_objectsFromStatementHandle(AttributeAdaptor *ata, StatementHandle *sth);


#endif
