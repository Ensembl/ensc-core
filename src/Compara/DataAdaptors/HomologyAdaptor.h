#ifndef __HOMOLOGYADAPTOR_H__
#define __HOMOLOGYADAPTOR_H__

#include "BaseComparaAdaptor.h"
#include "ComparaAdaptorTypes.h"
#include "Homology.h"
#include "Vector.h"

struct HomologyAdaptorStruct {
  BASECOMPARAADAPTOR_DATA
};

HomologyAdaptor *HomologyAdaptor_new(ComparaDBAdaptor *dba);

Vector *HomologyAdaptor_fetchHomologuesBySpeciesRelationshipId(HomologyAdaptor *ha,
                   char *hSpecies, IDType internalId);

IDType HomologyAdaptor_getRelationship(HomologyAdaptor *ha, char *qStr);
int HomologyAdaptor_getRelationships(HomologyAdaptor *ha, char *qStr, IDType **idsP);
Vector *HomologyAdaptor_getHomologues(HomologyAdaptor *ha, char *qStr);



#endif
