#ifndef __HOMOLOGY_H__
#define __HOMOLOGY_H__

#include "ComparaDataModelTypes.h"
#include "EnsC.h"
#include "EnsRoot.h"

OBJECTFUNC_TYPES(Homology)

typedef struct HomologyFuncsStruct {
  OBJECTFUNCS_DATA(Homology)
} HomologyFuncs;

#define FUNCSTRUCTTYPE HomologyFuncs
struct HomologyStruct {
  ENSROOT_DATA
  ECOSTRING species;
  ECOSTRING stableId;
  int chrStart;
  int chrEnd;
  ECOSTRING chrName;
};
#undef FUNCSTRUCTTYPE

Homology *Homology_new();

#define Homology_setChrStart(hom, cs) (hom)->chrStart = (cs)
#define Homology_getChrStart(hom) (hom)->chrStart

#define Homology_setChrEnd(hom, ce) (hom)->chrEnd = (ce)
#define Homology_getChrEnd(hom) (hom)->chrEnd

ECOSTRING Homology_setChromosome(Homology *hom, char *name);
#define Homology_getChromosome(hom) (hom)->chrName

ECOSTRING Homology_setStableId(Homology *hom, char *sid);
#define Homology_getStableId(hom) (hom)->stableId

ECOSTRING Homology_setSpecies(Homology *hom, char *spec);
#define Homology_getSpecies(hom) (hom)->species

void Homology_free(Homology *hom);

#ifdef __HOMOLOGY_MAIN__
  HomologyFuncs
    homologyFuncs = {
                    Homology_free
                   };
#else
  extern HomologyFuncs homologyFuncs;
#endif


#endif
