#ifndef __TRANSCRIPT_H__
#define __TRANSCRIPT_H__

#include "DataModelTypes.h"

#include "SeqFeature.h"
#include "FeatureSet.h"
#include "Storable.h"
#include "Translation.h"
#include "IDHash.h"

struct TranscriptStruct {
  SeqFeature sf;
  Translation *translation;
  IDType translationId;
  FeatureSet fs;
  StableIdInfo si;
  char *type;
};

#define Transcript_setStableId(transcript,sid)  StableIdInfo_setStableId(&((transcript)->si),(sid))
char *Transcript_getStableId(Transcript *transcript);

char *Transcript_setType(Transcript *transcript,char *type);
#define Transcript_getType(transcript) (transcript)->type

#define Transcript_setVersion(transcript,ver)  StableIdInfo_setVersion(&((transcript)->si),(ver))
int Transcript_getVersion(Transcript *transcript);

#define Transcript_setStart(transcript,start) SeqFeature_setStart(&((transcript)->sf),start)
#define Transcript_getStart(transcript) SeqFeature_getStart(&((transcript)->sf))

#define Transcript_setEnd(transcript,end) SeqFeature_setEnd(&((transcript)->sf),end)
#define Transcript_getEnd(transcript) SeqFeature_getEnd(&((transcript)->sf))

#define Transcript_setStrand(transcript,strand) SeqFeature_setStrand(&((transcript)->sf),strand)
#define Transcript_getStrand(transcript) SeqFeature_getStrand(&((transcript)->sf))

Transcript *Transcript_new(void);

#define Transcript_addExon(trans,exon) FeatureSet_addFeature(&((trans)->fs),exon)
#define Transcript_getExonAt(trans,ind) FeatureSet_getFeatureAt(&((trans)->fs),ind)

#define Transcript_getExonCount(trans) FeatureSet_getNumFeature(&((trans)->fs))

#define Transcript_removeAllExons(trans) FeatureSet_removeAll(&((trans)->fs))

#define Transcript_setTranslationId(transcript,tid) (transcript)->translationId = (tid)
#define Transcript_getTranslationId(transcript) (transcript)->translationId

#define Transcript_setTranslation(transcript,tn) (transcript)->translation = (tn)
#define Transcript_getTranslation(transcript) (transcript)->translation

#define Transcript_setDbID(transcript,id) SeqFeature_setDbID(&((transcript)->sf),(id))
#define Transcript_getDbID(transcript) SeqFeature_getDbID(&((transcript)->sf))

#define Transcript_setAdaptor(transcript,ad) SeqFeature_setAdaptor(&((transcript)->sf),(ad))
#define Transcript_getAdaptor(transcript) SeqFeature_getAdaptor(&((transcript)->sf))

void Transcript_flushExons(Transcript *trans);


#endif
