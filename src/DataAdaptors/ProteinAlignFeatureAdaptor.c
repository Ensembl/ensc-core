#include "ProteinAlignFeatureAdaptor.h"


char ***ProteinAlignFeatureAdaptor_tableNames = {{"protein_align_feature","paf"}};

ProteinAlignFeatureAdaptor *ProteinAlignFeatureAdaptor_new(DBAdaptor *dba) {
  ProteinAlignFeatureAdaptor *pafa;

  if ((pafa = (ProteinAlignFeatureAdaptor *)calloc(1,sizeof(ProteinAlignFeatureAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for ProteinAlignFeatureAdaptor\n");
    return NULL;
  }
  BaseFeatureAdaptor_init((BaseFeatureAdaptor *)pafa, dba, PROTEINALIGNFEATURE_ADAPTOR);

  return pafa;
}

int ProteinAlignFeatureAdaptor_store(BaseFeatureAdaptor *baf, Set *features) {
}


char ***ProteinAlignFeatureAdaptor_getTables() {
  return ProteinAlignFeatureAdaptor_tableNames;
}

char *ProteinAlignFeatureAdaptor_getColumns() {
}

Set *ProteinAlignFeatureAdaptor_objectsFromStatementHandle(BaseFeatureAdaptor *baf,
                                                           StatementHandle *sth) {
}
