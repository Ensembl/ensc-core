#ifndef __HOMOLOGY_H__
#define __HOMOLOGY_H__

#include "ComparaDataModelTypes.h"
#include "EnsRoot.h"

#define FUNCSTRUCTTYPE NoTypeFuncs
struct HomologyStruct {
  ENSROOT_DATA
  char *species;
  char *stableId;
  int chrStart;
  int chrEnd;
  char *chrName;
};
#undef FUNCSTRUCTTYPE

Homology *Homology_new();

#define Homology_setChrStart(hom, cs) (hom)->chrStart = (cs)
#define Homology_getChrStart(hom) (hom)->chrStart

#define Homology_setChrEnd(hom, ce) (hom)->chrEnd = (ce)
#define Homology_getChrEnd(hom) (hom)->chrEnd

char *Homology_setChrName(Homology *hom, char *name);
#define Homology_getChrName(hom) (hom)->chrName

char *Homology_setStableId(Homology *hom, char *sid);
#define Homology_getStableId(hom) (hom)->stableId

char *Homology_setSpecies(Homology *hom, char *spec);
#define Homology_getSpecies(hom) (hom)->species


#endif
