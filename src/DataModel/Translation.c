#include "Translation.h"
#include "TranslationAdaptor.h"

Translation *Translation_new() {
  Translation *t;

  if ((t = (Translation *)calloc(1,sizeof(Translation))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for t\n");
    return NULL;
  }

  Translation_setVersion(t,-1);

  t->objectType = CLASS_TRANSLATION;

  return t;
}

char *Translation_getStableId(Translation *translation) {
  TranslationAdaptor *ta = (TranslationAdaptor *)Translation_getAdaptor(translation);

  if (StableIdInfo_getStableId(&(translation->si)) == NULL && ta) {
    TranslationAdaptor_getStableEntryInfo(ta,translation);
  }
  return StableIdInfo_getStableId(&(translation->si));
}

int Translation_getVersion(Translation *translation) {
  TranslationAdaptor *ta = (TranslationAdaptor *)Translation_getAdaptor(translation);

  if (StableIdInfo_getVersion(&(translation->si)) == -1 && ta) {
    TranslationAdaptor_getStableEntryInfo(ta,translation);
  }
  return StableIdInfo_getVersion(&(translation->si));
}

void Translation_transform(Translation *translation, IDHash *exonTransforms) {

  Exon * startExon = Translation_getStartExon(translation);
  Exon * endExon   = Translation_getEndExon(translation);
  IDType startExonRef = (int)startExon;
  IDType endExonRef = (int)endExon;

/* CHECK */
  if (IDHash_contains(exonTransforms,startExonRef)) {
    Translation_setStartExon(translation,IDHash_getValue(exonTransforms,startExonRef));
  } else {
    // do nothing, the start exon wasnt mapped
  }

  if (IDHash_contains(exonTransforms,endExonRef)) {
    Translation_setEndExon(translation,IDHash_getValue(exonTransforms,endExonRef));
  } else {
    // do nothing, the end exon wasnt mapped
  }
}


