#ifndef __DNAALIGNFEATUREADAPTOR_H__
#define __DNAALIGNFEATUREADAPTOR_H__

#include "BaseFeatureAdaptor.h"
#include "AdaptorTypes.h"

struct DNAAlignFeatureAdaptorStruct {
  BASEFEATUREADAPTOR_DATA
};

DNAAlignFeatureAdaptor *DNAAlignFeatureAdaptor_new(DBAdaptor *dba);
int DNAAlignFeatureAdaptor_store(BaseFeatureAdaptor *bfa, Vector *features);
NameTableType *DNAAlignFeatureAdaptor_getTables(void); 
char **DNAAlignFeatureAdaptor_getColumns(void);
Vector *DNAAlignFeatureAdaptor_objectsFromStatementHandle(BaseFeatureAdaptor *bfa,
                                                       StatementHandle *sth,
                                                       AssemblyMapper *assMapper,
                                                       Slice *slice);
int DNAAlignFeatureAdaptor_store(BaseFeatureAdaptor *bfa, Vector *features);


#define DNAAlignFeatureAdaptor_fetchByDbID(dafa, id) BaseFeatureAdaptor_fetchByDbID((BaseFeatureAdaptor *)(dafa), (id))
#define DNAAlignFeatureAdaptor_fetchAllBySliceAndScore(dafa, slice, score, lname) \
          BaseFeatureAdaptor_fetchAllBySliceAndScore((BaseFeatureAdaptor *)(dafa), (slice), (score), (lname))
#define DNAAlignFeatureAdaptor_fetchAllByRawContigAndScore(dafa, contig, score, lname) \
          BaseFeatureAdaptor_fetchAllByRawContigAndScore((BaseFeatureAdaptor *)(dafa), (contig), (score), (lname))


#endif
