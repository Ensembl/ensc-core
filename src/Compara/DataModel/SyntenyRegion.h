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
  char *chrName;
  int chrStart;
  int chrEnd;
  char *hitChrName;
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

#define SyntenyRegion_setClusterId(sr,cid) (sr)->clusterId = (cid)
#define SyntenyRegion_getClusterId(sr) (sr)->clusterId

#define SyntenyRegion_setDNAFragId(sr,did) (sr)->dnaFragId = (did)
#define SyntenyRegion_getDNAFragId(sr) (sr)->dnaFragId

#define SyntenyRegion_setDbID(sr,dbid) Storable_setDbID(&((sr)->st),(dbid))
#define SyntenyRegion_getDbID(sr) Storable_getDbID(&((sr)->st))

#define SyntenyRegion_setAdaptor(sr,ad) Storable_setAdaptor((sr),(ad))
#define SyntenyRegion_getAdaptor(sr) Storable_getAdaptor((sr))

#endif
