#ifndef __CLONE_H__
#define __CLONE_H__

#include <time.h>

#include "DataModelTypes.h"
#include "Storable.h"
#include "EnsRoot.h"


#define FUNCSTRUCTTYPE NoTypeFuncs
struct CloneStruct {
  ENSROOT_DATA
  char *emblAcc;
  char *name;
  int  version;
  int  emblVersion;
  char htgPhase;
  time_t created;
  time_t modified;
  Storable st;
};
#undef FUNCSTRUCTTYPE

Clone *Clone_new(void);

#define Clone_setDbID(cl,id) Storable_setDbID(&((cl)->st),(id))
#define Clone_getDbID(cl) Storable_getDbID(&((cl)->st))

#define Clone_setAdaptor(cl,ad) Storable_setAdaptor(&((cl)->st),(ad))
#define Clone_getAdaptor(cl) Storable_getAdaptor(&((cl)->st))

char *Clone_setName(Clone *cl, char *name);
#define Clone_getName(cl) (cl)->name

char *Clone_setEmblAcc(Clone *cl, char *acc);
#define Clone_getEmblAcc(cl) (cl)->emblAcc

#define Clone_setVersion(cl,ver) (cl)->version = (ver)
#define Clone_getVersion(cl) cl->version

#define Clone_setEmblVersion(cl,ev) (cl)->emblVersion = (ev)
#define Clone_getEmblVersion(cl) (cl)->emblVersion

#define Clone_setHTGPhase(cl,htg) (cl)->htgPhase = (htg)
#define Clone_getHTGPhase(cl) (cl)->htgPhase

#define Clone_setCreated(cl,ct) (cl)->created = (ct)
#define Clone_getCreated(cl) (cl)->created

#define Clone_setModified(cl,mod) (cl)->modified = (mod)
#define Clone_getModified(cl) (cl)->modified


#endif
