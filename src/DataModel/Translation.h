#ifndef __TRANSLATION_H__
#define __TRANSLATION_H__

#include "DataModelTypes.h"
#include "Exon.h"

struct TranslationStruct {
  int start;
  int end;
  Exon *startExon;
  Exon *endExon;
  Storable st; 
  StableIdInfo si;
};

#define Translation_setStart(translation,s) (translation)->start = (s)
#define Translation_getStart(translation) (translation)->start

#define Translation_setEnd(translation,e) (translation)->end = (e)
#define Translation_getEnd(translation) (translation)->end

#define Translation_setStartExon(translation,s) (translation)->startExon = (s)
#define Translation_getStartExon(translation) (translation)->startExon

#define Translation_setEndExon(translation,e) (translation)->endExon = (e)
#define Translation_getEndExon(translation) (translation)->endExon

#define Translation_setDbID(t,dbID) Storable_setDbID(&((t)->st),(dbID))
#define Translation_getDbID(t) Storable_getDbID(&((t)->st))

#define Translation_setAdaptor(t,ad) Storable_setAdaptor(&((t)->st),(ad))
#define Translation_getAdaptor(t) Storable_getAdaptor(&((t)->st))

#define Translation_setStableId(translation,sid)  StableIdInfo_setStableId(&((translation)->si),(sid))
char *Translation_getStableId(Translation *translation);

#define Translation_setVersion(translation,ver)  StableIdInfo_setVersion(&((translation)->si),(ver))
int Translation_getVersion(Translation *translation);

Translation *Translation_new(void);

#endif
