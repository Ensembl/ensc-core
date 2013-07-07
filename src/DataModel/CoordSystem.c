#define __COORDSYSTEM_MAIN__
#include "CoordSystem.h"
#include "StrUtil.h"
#undef __COORDSYSTEM_MAIN__

#include <string.h>


/*
=head2 new

  Arg [..]   : List of named arguments:
               -NAME      - The name of the coordinate system
               -VERSION   - (optional) The version of the coordinate system.
                            Note that if the version passed in is undefined,
                            it will be set to the empty string in the
                            resulting CoordSystem object.
               -RANK      - The rank of the coordinate system. The highest
                            level coordinate system should have rank 1, the
                            second highest rank 2 and so on.  An example of
                            a high level coordinate system is 'chromosome' an
                            example of a lower level coordinate system is
                            'clone'.
               -TOP_LEVEL - (optional) Sets whether this is a top-level coord
                            system. Default = 0. This should only be set to
                            true if you are creating an artificial toplevel
                            coordsystem by the name of 'toplevel'
               -SEQUENCE_LEVEL - (optional) Sets whether this is a sequence
                            level coordinate system. Default = 0
               -DEFAULT   - (optional)
                            Whether this is the default version of the 
                            coordinate systems of this name. Default = 0
               -DBID      - (optional) The internal identifier of this
                             coordinate system
               -ADAPTOR   - (optional) The adaptor which provides database
                            interaction for this object
  Example    : $cs = Bio::EnsEMBL::CoordSystem->new(-NAME    => 'chromosome',
                                                    -VERSION => 'NCBI33',
                                                    -RANK    => 1,
                                                    -DBID    => 1,
                                                    -ADAPTOR => adaptor,
                                                    -DEFAULT => 1,
                                                    -SEQUENCE_LEVEL => 0);
  Description: Creates a new CoordSystem object representing a coordinate
               system.
  Returntype : Bio::EnsEMBL::CoordSystem
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/

CoordSystem *CoordSystem_new(char *name, char *version, int rank, IDType dbID, CoordSystemAdaptor *csa, int isDefault, int isSeqLevel, int isTopLevel) {

  if (version == NULL) {
    version = emptyString;
  }

  isTopLevel = isTopLevel  ? 1 : 0;
  isSeqLevel = isSeqLevel ? 1 : 0;
  isDefault  = isDefault  ? 1 : 0;


  if ( isTopLevel == 1 ) {
    if ( rank != 0 ) {
      fprintf(stderr,"RANK argument must be 0 if TOP_LEVEL is 1");
      exit(1);
    }

    if ( name != NULL ) {
      if ( strcmp(name,"toplevel")) {
        fprintf(stderr,"The NAME argument must be 'toplevel' if TOP_LEVEL is 1");
        exit(1);
      }
    } else {
      name="toplevel";
    }

    if ( isSeqLevel == 1 ) {
      fprintf(stderr,"SEQUENCE_LEVEL argument must be 0 if TOP_LEVEL is 1");
      exit(1);
    }
    isDefault = 0;

  } else {
    if ( rank == 0 ) {
      fprintf(stderr,"RANK argument must be non-zero unless TOP_LEVEL is 1");
      exit(1);
    }
    if ( name == NULL ) {
      fprintf(stderr,"The NAME argument is required");
      exit(1);
    } else if ( !strcmp(name,"toplevel")) {
      fprintf( stderr, "Cannot name coord system 'toplevel unless TOP_LEVEL is 1" );
      exit(1);
    }
  }

  if ( rank < 0 ) {
    fprintf(stderr, "The RANK argument must be a non negative integer");
    exit(1);
  }

  CoordSystem *cs;
  if ((cs = (CoordSystem *)calloc(1,sizeof(CoordSystem))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for cs\n");
    return NULL;
  }

  cs->objectType = CLASS_COORDSYSTEM;

  cs->funcs = &coordSystemFuncs;
 
  EcoString_copyStr(ecoSTable,&(cs->version),version,0);
  EcoString_copyStr(ecoSTable,&(cs->name),name,0);
  cs->isTopLevel = isTopLevel;
  cs->isSeqLevel = isSeqLevel;
  cs->isDefaultVersion  = isDefault;
  cs->rank       = rank;

  CoordSystem_setAdaptor(cs, (BaseAdaptor *)csa);
  CoordSystem_setDbID(cs, dbID);

  return cs;
}


char * CoordSystem_getSpecies(CoordSystem *cs) {
  return "CoordSystem_getSpecies not implemented!";
  //return $self->adaptor->db->species;
}



/*
=head2 equals

  Arg [1]    : Bio::EnsEMBL::CoordSystem $cs
               The coord system to compare to for equality.
  Example    : if($coord_sys->equals($other_coord_sys)) { ... }
  Description: Compares 2 coordinate systems and returns true if they are
               equivalent.  The definition of equivalent is sharing the same
               name and version.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/

int CoordSystem_compare(CoordSystem *cs1, CoordSystem *cs2) {

/* Dodgy DAS stuff
  if (!$cs2 || !ref($cs2) || !$cs2->isa('Bio::EnsEMBL::CoordSystem')) {
    if ($cs2->isa('Bio::EnsEMBL::ExternalData::DAS::CoordSystem')) {
      return $cs2->equals($self);
    }
    throw('Argument must be a CoordSystem');
  }
*/

  if (!EcoString_strcmp(CoordSystem_getVersion(cs1), CoordSystem_getVersion(cs2)) && !EcoString_strcmp(CoordSystem_getName(cs1), CoordSystem_getName(cs2))) {
    // Equal, so return 0
    return 0;
  }

  return 1;
}

char *CoordSystem_getNameColonVersion(CoordSystem *cs) {
  if (cs->nameColonVersion == NULL) {
    char tmpStr[1024];
    sprintf(tmpStr,"%s:%s",CoordSystem_getName(cs),CoordSystem_getVersion(cs) ? CoordSystem_getVersion(cs):"");
    StrUtil_copyString(&cs->nameColonVersion, tmpStr, 0);
    cs->lenNameColonVersion = strlen(tmpStr);
  }
 
  return cs->nameColonVersion;
}

int CoordSystem_getLenNameColonVersion(CoordSystem *cs) {
  if (cs->nameColonVersion == NULL) {
    CoordSystem_getNameColonVersion(cs);
  }
  return cs->lenNameColonVersion;
}

void CoordSystem_free(CoordSystem *cs) {
  
  if (cs->name) EcoString_freeStr(ecoSTable, cs->name);
  if (cs->version) EcoString_freeStr(ecoSTable, cs->version);

  free(cs);
}
