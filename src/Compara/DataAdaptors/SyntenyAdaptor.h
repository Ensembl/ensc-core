#ifndef __SYNTEYADAPTOR_H__
#define __SYNTEYADAPTOR_H__

#include "BaseComparaAdaptor.h"
#include "ComparaAdaptorTypes.h"
#include "Vector.h"
#include "SyntenyRegion.h"

struct SyntenyAdaptorStruct {
  BASECOMPARAADAPTOR_DATA
  char *speciesMain;
  char *speciesSecondary;
};

SyntenyAdaptor *SyntenyAdaptor_new(ComparaDBAdaptor *dba);

#define SyntenyAdaptor_setSpeciesMain(sa,sm) (sa)->speciesMain = (sm)
#define SyntenyAdaptor_getSpeciesMain(sa) (sa)->speciesMain

#define SyntenyAdaptor_setSpeciesSecondary(sa,ss) (sa)->speciesSecondary = (ss)
#define SyntenyAdaptor_getSpeciesSecondary(sa) (sa)->speciesSecondary

#endif
