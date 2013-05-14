#ifndef __TRANSCRIPTADAPTOR_H__
#define __TRANSCRIPTADAPTOR_H__

#include "BaseFeatureAdaptor.h"
#include "BaseFeatureAdaptor.h"
#include "AdaptorTypes.h"
#include "Transcript.h"

struct TranscriptAdaptorStruct {
  BASEFEATUREADAPTOR_DATA
};

TranscriptAdaptor *TranscriptAdaptor_new(DBAdaptor *dba);
IDType             TranscriptAdaptor_store(TranscriptAdaptor *ta, Transcript *transcript, IDType geneId, IDType analysisId);
NameTableType *    TranscriptAdaptor_getTables();
char **            TranscriptAdaptor_getColumns();
NameTableType *    TranscriptAdaptor_leftJoin();
Transcript *       TranscriptAdaptor_fetchByStableId(TranscriptAdaptor *ta, char *stableId);
Vector *           TranscriptAdaptor_fetchAll(TranscriptAdaptor *ta);
Vector *           TranscriptAdaptor_fetchAllVersionsByStableId(TranscriptAdaptor *ta, char *stableId);
Transcript *       TranscriptAdaptor_fetchByTranslationStableId(TranscriptAdaptor *ta, char *translationStableId);
Transcript *       TranscriptAdaptor_fetchByTranslationId(TranscriptAdaptor *ta, IDType translationId);
Vector *           TranscriptAdaptor_fetchAllByGene(TranscriptAdaptor *ta, Gene *gene);
Vector *           TranscriptAdaptor_fetchAllBySlice(TranscriptAdaptor *ta, Slice *slice, int loadExons, char *logicName, char *inputConstraint);
Transcript *       TranscriptAdaptor_fetchByDisplayLabel(TranscriptAdaptor *ta, char *label);
Vector *           TranscriptAdaptor_fetchByExonStableId(TranscriptAdaptor *ta, char *stableId);
Vector *           TranscriptAdaptor_fetchAllByBiotype(TranscriptAdaptor *ta, Vector *biotypes);
void               TranscriptAdaptor_biotypeConstraint(TranscriptAdaptor *ta, Vector *biotypes, char *constraint);
int                TranscriptAdaptor_isTranscriptCanonical(TranscriptAdaptor *ta, Transcript *transcript);
Vector *           TranscriptAdaptor_listDbIDs(TranscriptAdaptor *ta, int ordered);
Vector *           TranscriptAdaptor_listStableIDs(TranscriptAdaptor *ta);
Vector *           TranscriptAdaptor_objectsFromStatementHandle(TranscriptAdaptor *ta, StatementHandle *sth, AssemblyMapper *assMapper, Slice *destSlice);

#define TranscriptAdaptor_genericFetch(ta, constraint, mapper, slice) \
      BaseAdaptor_genericFetch((BaseAdaptor *)(ta), (constraint), (mapper), (slice))

#define TranscriptAdaptor_fetchByDbID(ta,id)  \
   BaseFeatureAdaptor_fetchByDbID((BaseFeatureAdaptor *)(ta), (id))

#define TranscriptAdaptor_fetchAllBySliceConstraint(ta,slice,constraint,logicName)  \
   BaseFeatureAdaptor_fetchAllBySliceConstraint((BaseFeatureAdaptor *)(ta), (slice), (constraint), (logicName))





#endif
