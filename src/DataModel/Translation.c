#define __TRANSLATION_MAIN__
#include "Translation.h"
#undef __TRANSLATION_MAIN__
#include "TranslationAdaptor.h"

Translation *Translation_new() {
  Translation *t;

  if ((t = (Translation *)calloc(1,sizeof(Translation))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for t\n");
    return NULL;
  }

  Translation_setVersion(t,-1);

  t->objectType = CLASS_TRANSLATION;

  t->funcs = &translationFuncs;

  Object_incRefCount(t);

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
  IDType startExonRef = (IDType)startExon;
  IDType endExonRef = (IDType)endExon;

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

void Translation_free(Translation *translation) {
  Object_decRefCount(translation);

  if (Object_getRefCount(translation) > 0) {
    return;
  } else if (Object_getRefCount(translation) < 0) {
    fprintf(stderr,"Error: Negative reference count for Translation\n"
                   "       Freeing it anyway\n");
  }

  printf("Translation_free not implemented\n");
}
