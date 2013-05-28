#ifndef __INTRONSUPPORTINGEVIDENCEADAPTOR_H__
#define __INTRONSUPPORTINGEVIDENCEADAPTOR_H__

#include "BaseFeatureAdaptor.h"
#include "BaseAdaptor.h"
#include "AdaptorTypes.h"
#include "IntronSupportingEvidence.h"
#include "Intron.h"

struct IntronSupportingEvidenceAdaptorStruct {
  BASEFEATUREADAPTOR_DATA
};

IntronSupportingEvidenceAdaptor *IntronSupportingEvidenceAdaptor_new(DBAdaptor *dba);
NameTableType *IntronSupportingEvidenceAdaptor_getTables();
char **IntronSupportingEvidenceAdaptor_getColumns();
Vector *IntronSupportingEvidenceAdaptor_listLinkedTranscriptIds(IntronSupportingEvidenceAdaptor *isea, IntronSupportingEvidence *ise);
Vector *IntronSupportingEvidenceAdaptor_fetchAllByTranscript(IntronSupportingEvidenceAdaptor *isea, Transcript *transcript);
IDType *IntronSupportingEvidenceAdaptor_fetchFlankingExonIds(IntronSupportingEvidenceAdaptor *isea, IntronSupportingEvidence *ise, Transcript *transcript, IDType *flanks);
IDType IntronSupportingEvidenceAdaptor_store(IntronSupportingEvidenceAdaptor *isea, IntronSupportingEvidence *ise);
Vector *IntronSupportingEvidenceAdaptor_objectsFromStatementHandle(IntronSupportingEvidenceAdaptor *isea, StatementHandle *sth, AssemblyMapper *assMapper, Slice *destSlice);

void IntronSupportingEvidenceAdaptor_storeTranscriptLinkage(IntronSupportingEvidenceAdaptor *isea, IntronSupportingEvidence *sf, Transcript *transcript, IDType transcriptId); 

#define IntronSupportingEvidenceAdaptor_fetchAllByDbIDList(isea,id,slice)  \
   BaseAdaptor_fetchAllByDbIDList((BaseAdaptor *)(isea), (id), (slice))


#endif
