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

int ExonStickyRankCompFunc(const void *a, const void *b) {
  Exon **e1 = (Exon **)a;
  Exon **e2 = (Exon **)b;

  if (Exon_getStickyRank(*e1) > Exon_getStickyRank(*e2)) {
    return -1;
  } else if (Exon_getStickyRank(*e1) < Exon_getStickyRank(*e2)) {
    return 1;
  } else {
    return 0;
  }
}

void Exon_sortByStickyRank(Exon *exon) {
  if (Exon_isSticky(exon)) {
    qsort(Exon_getComponents(exon), Exon_getNumComponentExon(exon), sizeof(void *), 
          ExonStickyRankCompFunc);
  }
  return; 
}

Exon *Exon_transformToSlice(Exon *exon, Slice *slice) {
}
