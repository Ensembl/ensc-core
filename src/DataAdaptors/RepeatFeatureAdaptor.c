#include "RepeatFeatureAdaptor.h"


char ***RepeatFeatureAdaptor_tableNames = {{"repeat_feature","rf"},
                                           {"repeat_consensus","rc"}};

RepeatFeatureAdaptor *RepeatFeatureAdaptor_new(DBAdaptor *dba) {
  RepeatFeatureAdaptor *rfa;

  if ((rfa = (RepeatFeatureAdaptor *)calloc(1,sizeof(RepeatFeatureAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for RepeatFeatureAdaptor\n");
    return NULL;
  }
  BaseFeatureAdaptor_init((BaseFeatureAdaptor *)rfa, dba, REPEATFEATURE_ADAPTOR);

  return rfa;
}

int RepeatFeatureAdaptor_store(BaseFeatureAdaptor *baf, Set *features) {
}


char ***RepeatFeatureAdaptor_getTables() {
  return RepeatFeatureAdaptor_tableNames;
}

char *RepeatFeatureAdaptor_getColumns() {
}

Set *RepeatFeatureAdaptor_objectsFromStatementHandle(BaseFeatureAdaptor *baf,
                                                     StatementHandle *sth) {
}
