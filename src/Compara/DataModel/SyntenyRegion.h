#ifndef __SYNTENYREGION_H__
#define __SYNTENYREGION_H__

#include "ComparaDataModelTypes.h"
#include "EnsRoot.h"
#include "Storable.h"

#define FUNCSTRUCTTYPE NoTypeFuncs
struct SyntenyRegionStruct {
  ENSROOT_DATA
  Storable st;
  int start;
  int end;
  IDType clusterId;
  IDType dnaFragId;
  ECOSTRING chrName;
  ECOSTRING hitChrName;
  ECOSTRING seqType;
  ECOSTRING hitSeqType;
  int chrStart;
  int chrEnd;
  int hitChrStart;
  int hitChrEnd;
  int relOri;
};
#undef FUNCSTRUCTTYPE

SyntenyRegion *SyntenyRegion_new();

#define SyntenyRegion_setEnd(sr,e) (sr)->end = (e)
#define SyntenyRegion_getEnd(sr) (sr)->end

#define SyntenyRegion_setStart(sr,s) (sr)->start = (s)
#define SyntenyRegion_getStart(sr) (sr)->start

#define SyntenyRegion_setRelOri(sr,ro) (sr)->relOri = (ro)
#define SyntenyRegion_getRelOri(sr) (sr)->relOri

#define SyntenyRegion_setClusterId(sr,cid) (sr)->clusterId = (cid)
#define SyntenyRegion_getClusterId(sr) (sr)->clusterId

#define SyntenyRegion_setDNAFragId(sr,did) (sr)->dnaFragId = (did)
#define SyntenyRegion_getDNAFragId(sr) (sr)->dnaFragId

#define SyntenyRegion_setDbID(sr,dbid) Storable_setDbID(&((sr)->st),(dbid))
#define SyntenyRegion_getDbID(sr) Storable_getDbID(&((sr)->st))

#define SyntenyRegion_setAdaptor(sr,ad) Storable_setAdaptor((sr),(ad))
#define SyntenyRegion_getAdaptor(sr) Storable_getAdaptor((sr))

#define SyntenyRegion_setChrEnd(sr,e) (sr)->chrEnd = (e)
#define SyntenyRegion_getChrEnd(sr) (sr)->chrEnd

#define SyntenyRegion_setChrStart(sr,s) (sr)->chrStart = (s)
#define SyntenyRegion_getChrStart(sr) (sr)->chrStart

ECOSTRING SyntenyRegion_setChrName(SyntenyRegion *sr, char *chrName);
#define SyntenyRegion_getChrName(sr) (sr)->chrName

#define SyntenyRegion_setHitChrEnd(sr,e) (sr)->hitChrEnd = (e)
#define SyntenyRegion_getHitChrEnd(sr) (sr)->hitChrEnd

#define SyntenyRegion_setHitChrStart(sr,s) (sr)->hitChrStart = (s)
#define SyntenyRegion_getHitChrStart(sr) (sr)->hitChrStart

ECOSTRING SyntenyRegion_setHitChrName(SyntenyRegion *sr, char *hitChrName);
#define SyntenyRegion_getHitChrName(sr) (sr)->hitChrName

ECOSTRING SyntenyRegion_setHitSeqType(SyntenyRegion *sr, char *hitSeqType);

ECOSTRING SyntenyRegion_setSeqType(SyntenyRegion *sr, char *seqType);
#endif
