#include "Exon.h"
#include "ExonAdaptor.h"

Exon *Exon_new() {
  Exon *exon;

  if ((exon = (Exon *)calloc(1,sizeof(Exon))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for exon\n");
    return NULL;
  }

/* Set to empty values */
  Exon_setModified(exon,0);
  Exon_setCreated(exon,0);
  Exon_setVersion(exon,-1);

  return exon;
}

char *Exon_getStableId(Exon *exon) {
  ExonAdaptor *ea = (ExonAdaptor *)Exon_getAdaptor(exon);

  if (StableIdInfo_getStableId(&(exon->si)) == NULL && ea) {
    ExonAdaptor_getStableEntryInfo(ea,exon);
  }
  return StableIdInfo_getStableId(&(exon->si));
}

time_t Exon_getCreated(Exon *exon) {
  ExonAdaptor *ea = (ExonAdaptor *)Exon_getAdaptor(exon);

  if (StableIdInfo_getCreated(&(exon->si)) == 0 && ea) {
    ExonAdaptor_getStableEntryInfo(ea,exon);
  }
  return StableIdInfo_getCreated(&(exon->si));
}

time_t Exon_getModified(Exon *exon) {
  ExonAdaptor *ea = (ExonAdaptor *)Exon_getAdaptor(exon);

  if (StableIdInfo_getModified(&(exon->si)) == 0 && ea) {
    ExonAdaptor_getStableEntryInfo(ea,exon);
  }
  return StableIdInfo_getModified(&(exon->si));
}

int Exon_getVersion(Exon *exon) {
  ExonAdaptor *ea = (ExonAdaptor *)Exon_getAdaptor(exon);

  if (StableIdInfo_getVersion(&(exon->si)) == -1 && ea) {
    ExonAdaptor_getStableEntryInfo(ea,exon);
  }
  return StableIdInfo_getVersion(&(exon->si));
}

void Exon_sortByStickyRank(Exon *exon) {
  fprintf(stderr,"WARNING: sortByStickRank not implemented\n");
}
