/*
 * Copyright [1999-2024] EMBL-European Bioinformatics Institute
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *      http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "CoordSystemAdaptor.h"
#include "BaseAdaptor.h"
#include "MysqlUtil.h"
#include "DBAdaptor.h"
#include "MetaContainer.h"

#include "StrUtil.h"

#include "StatementHandle.h"
#include "ResultRow.h"

CoordSystemAdaptor *CoordSystemAdaptor_new(DBAdaptor *dba) {
  CoordSystemAdaptor *csa;
  char qStr[256];
  StatementHandle *sth;
  ResultRow *row;

  if ((csa = (CoordSystemAdaptor *)calloc(1,sizeof(CoordSystemAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for CoordSystemAdaptor\n");
    return NULL;
  }
  BaseAdaptor_init((BaseAdaptor *)csa, dba, COORDSYSTEM_ADAPTOR);

  //
  // Cache the entire contents of the coord_system table cross-referenced
  // by dbID and name.
  //
  csa->dbIDCache = IDHash_new(IDHASH_SMALL);
  csa->nameCache = StringHash_new(STRINGHASH_SMALL);
  csa->rankCache = IDHash_new(IDHASH_SMALL);

  csa->isSeqLevelCache = IDHash_new(IDHASH_SMALL);
  csa->isDefaultVersionCache = IDHash_new(IDHASH_SMALL);

/* Not implemented
  #cache to store the seq_region_mapping information
  #from internal->external
  $self->{'_internal_seq_region_mapping'} = {};
  #from external->internal
  $self->{'_external_seq_region_mapping'} = {};
*/

  IDType speciesID = 1; // Hardcoded for now

  sprintf(qStr, "SELECT coord_system_id, name, rank, version, attrib "
                    "FROM coord_system "
                    "WHERE species_id = "
                    IDFMTSTR, speciesID); 

  sth = csa->prepare((BaseAdaptor *)csa,qStr,strlen(qStr));
  sth->execute(sth);

  while ((row = sth->fetchRow(sth))) {
    int seqLvl = 0;
    int defaultVer = 0;

    IDType dbID    = row->getLongLongAt(row,0);
    char * name    = row->getStringAt(row,1);
    int    rank    = row->getIntAt(row,2);
    char * version = row->getStringAt(row,3);
    char * attribs = row->getStringAt(row,4);

    if ( attribs && strcmp(attribs,"") ) {
      char **tokens = NULL;
      int ntok;
      int i;

      if (!StrUtil_tokenizeByDelim(&tokens, &ntok, attribs, ",")) {
        fprintf(stderr,"Failed tokenizing attribs string %s\n", attribs);
        free(csa);
        return NULL;
      }
        
      
      for (i=0; i<ntok; i++) {
        char *attrib = tokens[i];

        if ( !strcmp(attrib, "sequence_level") ) {
          seqLvl = 1;
        } else if ( !strcmp(attrib, "default_version")) {
          defaultVer = 1;
        } else {
          fprintf(stderr,"Unknown attrib type %s in CoordSystemAdaptor\n",attrib);
        }
      }

      for (i=0; i<ntok; i++) {
        free(tokens[i]);
      }
      free(tokens);
    }

    CoordSystem *cs = CoordSystem_new( name, version, rank, dbID, csa, seqLvl, defaultVer, 0 /*toplevel*/ );

    IDHash_add(csa->dbIDCache, CoordSystem_getDbID(cs), cs);

    // Moved these down slightly so can use cs as value rather than true/false
    if (seqLvl) {
      IDHash_add(csa->isSeqLevelCache, CoordSystem_getDbID(cs), cs);
    }

    if (defaultVer) {
      IDHash_add(csa->isDefaultVersionCache, CoordSystem_getDbID(cs), cs);
    }

    char lcName[1024];
    strcpy(lcName,name);
    StrUtil_strlwr(lcName);

    Vector *v;
    if (StringHash_contains(csa->nameCache, lcName)) {
      v = StringHash_getValue(csa->nameCache, lcName);
    } else {
      v = Vector_new();
      StringHash_add(csa->nameCache, lcName, v);
    }
    Vector_addElement(v, cs);

    if (IDHash_contains(csa->rankCache, (IDType)CoordSystem_getRank(cs))) {
      fprintf(stderr, "Already have a rank %d coord system\n",CoordSystem_getRank(cs));
    }
    IDHash_add(csa->rankCache, (IDType)CoordSystem_getRank(cs), cs);
  }

  sth->finish(sth);

  CoordSystemAdaptor_cacheMappingPaths(csa);

  CoordSystemAdaptor_cacheSeqRegionMapping(csa);

  return csa;
}

void CoordSystemAdaptor_dumpCachedMappings(CoordSystemAdaptor *csa) {
  char **keys = StringHash_getKeys(csa->mappingPaths);

  int i;
  for (i=0;i<StringHash_getNumValues(csa->mappingPaths);i++) { 
    Vector *path = StringHash_getValue(csa->mappingPaths, keys[i]);
    fprintf(stderr,"Path: %s has %d elements, pointer %p\n",keys[i], Vector_getNumElement(path), path);
    free(keys[i]);
  }
  free(keys);
}

void CoordSystemAdaptor_cacheSeqRegionMapping(CoordSystemAdaptor *csa) {
// Not Implemented Yet

  //
  // This cache will load the information from the seq_region_table, if
  // any, to allow mapping between internal and external seq_region_id.
  //


/*
  # For a given core database, will return the schema_build information.
  my $schema_build = $self->db->_get_schema_build();

  # Prepare the query to get relation for the current database being
  # used.
  my $sql = qq(
  SELECT    s.internal_seq_region_id,
            s.external_seq_region_id
  FROM      seq_region_mapping s,
            mapping_set ms,
            seq_region sr,
            coord_system cs
  WHERE     ms.mapping_set_id = s.mapping_set_id
    AND     ms.internal_schema_build = ?
    AND     s.internal_seq_region_id = sr.seq_region_id
    AND     sr.coord_system_id = cs.coord_system_id
    AND     cs.species_id = ?);

  my $sth = $self->prepare($sql);

  $sth->bind_param( 1, $schema_build,       SQL_VARCHAR );
  $sth->bind_param( 2, $self->species_id(), SQL_INTEGER );

  $sth->execute();

  # Load the cache:
  foreach my $row ( @{ $sth->fetchall_arrayref() } ) {
    # internal->external
    $self->{'_internal_seq_region_mapping'}->{ $row->[0] } = $row->[1];
    # external->internal
    $self->{'_external_seq_region_mapping'}->{ $row->[1] } = $row->[0];
  }

  $sth->finish();
*/

}


void CoordSystemAdaptor_cacheMappingPaths(CoordSystemAdaptor *csa) {
  // Retrieve a list of available mappings from the meta table.  This
  // may eventually be moved a table of its own if this proves too
  // cumbersome.
  StringHash *mappingPaths = StringHash_new(STRINGHASH_SMALL);
  MetaContainer *mc = DBAdaptor_getMetaContainer(csa->dba);
  Vector *mappingStrings = MetaContainer_listValueByKey(mc, "assembly.mapping");

  int i;
  MAP_PATH : for (i=0; i<Vector_getNumElement(mappingStrings); i++) {
    char *mapPath = Vector_getElementAt(mappingStrings, i);

    char **csStrings = NULL;
    int nCsString;

    if (!StrUtil_tokenizeByDelim(&csStrings, &nCsString, mapPath, "|#")) {
      fprintf(stderr, "Failed tokenizing mapPath string %s\n", mapPath);
      return;
    }

    if ( nCsString < 2 ) {
      fprintf(stderr, "Incorrectly formatted assembly.mapping value in meta table: %s\n", mapPath );
      continue;
    }

    Vector *coordSystems = Vector_new();

    int j;
    for (j=0; j<nCsString ; j++) {
      char *csString = csStrings[j];
      char *chP = index(csString,':');
      CoordSystem *cs;

      if (chP == NULL) {
        // csString is just name (no :)
        cs = CoordSystemAdaptor_fetchByName(csa, csString, NULL);
      } else {
        // csString is name and version - replace ':' with \0
        *chP = '\0';
        cs = CoordSystemAdaptor_fetchByName(csa, csString, chP+1);
      }
        

      if (cs == NULL) {
        fprintf(stderr, "Unknown coordinate system specified in meta table assembly.mapping:\n  %s", csString);
  
        goto MAP_PATH;
      }

      Vector_addElement(coordSystems, cs);
    }

    // Free csString tokens
    for (j=0; j<nCsString; j++) { free(csStrings[j]); }
    free(csStrings);

    // If the delimiter is a '#' we want a special case, multiple parts
    // of the same component map to the same assembly part.  As this
    // looks like the "long" mapping, we just make the path a bit longer
    if ( index( mapPath, '#' ) != NULL && Vector_getNumElement(coordSystems) == 2 ) {
      Vector_insertElementAt(coordSystems,1,NULL);
    }

    CoordSystem *cs1 = Vector_getElementAt(coordSystems, 0);
    CoordSystem *cs2 = Vector_getLastElement(coordSystems); //, Vector_getNumElement(coordSystems)-1);

    char key1[1024];
    char key2[1024];
    char keypair[2048];

    sprintf(key1,"%s:%s", CoordSystem_getName(cs1), CoordSystem_getVersion(cs1) ? CoordSystem_getVersion(cs1) : "");
    sprintf(key2,"%s:%s", CoordSystem_getName(cs2), CoordSystem_getVersion(cs2) ? CoordSystem_getVersion(cs2) : "");

    sprintf(keypair,"%s|%s",key1,key2);

    if ( StringHash_contains(mappingPaths,keypair) ) {
      fprintf(stderr,"Meta table specifies multiple mapping paths between coord systems %s and %s.\n"
                     "Choosing shorter path arbitrarily.", key1, key2 );

      Vector *existingMapping = StringHash_getValue(mappingPaths, keypair);

      if ( Vector_getNumElement(existingMapping) < Vector_getNumElement(coordSystems) ) {
        goto MAP_PATH;
      }
    }

    //fprintf(stderr, "Adding to cache: Mapping path %s length = %d\n",keypair, Vector_getNumElement(coordSystems));
    // NOTE: May in fact be a replace - need to think about memory
    StringHash_add(mappingPaths, keypair, coordSystems);
  }

  Vector_setFreeFunc(mappingStrings, free);
  Vector_free(mappingStrings);

  // Create the pseudo coord system 'toplevel' and cache it so that only
  // one of these is created for each database.

  CoordSystem *topLevel = CoordSystem_new("toplevel", "0"/*version*/, 0/*rank*/, 0, csa, 0/*seqlevel*/, 0/*defaultVer*/, 1 /*topLevel*/);
  csa->topLevel = topLevel;

  if (csa->mappingPaths) {
    fprintf(stderr, "Warning: Replacing existing mappingPaths cache in CoordSystemAdaptor\n");
  }

  csa->mappingPaths = mappingPaths;

  return;
}

/*
=head2 fetch_all

  Arg [1]    : none
  Example    : foreach my $cs (@{$csa->fetch_all()}) {
                 print $cs->name(), ' ', $cs->version(), "\n";
               }
  Description: Retrieves every coordinate system defined in the DB.
               These will be returned in ascending order of rank. I.e.
               The highest coordinate system with rank=1 would be first in the
               array.
  Returntype : listref of Bio::EnsEMBL::CoordSystems
  Exceptions : none
  Caller     : general
  Status     : Stable

*/

int CoordSystem_sortByRankFunc(const void *one, const void *two) {
  CoordSystem *cs1 = *((CoordSystem**)one);
  CoordSystem *cs2 = *((CoordSystem**)two);

  return CoordSystem_getRank(cs1) - CoordSystem_getRank(cs2);
}

Vector *CoordSystemAdaptor_fetchAll(CoordSystemAdaptor *csa) {

  Vector *coordSystems = IDHash_getValuesVector(csa->rankCache);

  //order the array by rank in ascending order
  Vector_sort(coordSystems, CoordSystem_sortByRankFunc);

  return coordSystems;
}



/*

  Arg [1]    : int $rank
  Example    : my $cs = $coord_sys_adaptor->fetch_by_rank(1);
  Description: Retrieves a CoordinateSystem via its rank. 0 is a special
               rank reserved for the pseudo coordinate system 'toplevel'.
               undef is returned if no coordinate system of the specified rank
               exists.
  Returntype : Bio::EnsEMBL::CoordSystem
  Exceptions : none
  Caller     : general
  Status     : Stable

*/

CoordSystem *CoordSystemAdaptor_fetchByRank(CoordSystemAdaptor *csa, int rank) {
  if (rank == 0) {
    return CoordSystemAdaptor_fetchTopLevel(csa);
  }

  if (IDHash_contains(csa->rankCache, rank)) {
    return IDHash_getValue(csa->rankCache, rank);
  } else {
    return NULL;
  }
}


/*

  Arg [1]    : string $name
               The name of the coordinate system to retrieve.  Alternatively
               this may be an alias for a real coordinate system.  Valid
               aliases are 'toplevel' and 'seqlevel'.
  Arg [2]    : string $version (optional)
               The version of the coordinate system to retrieve.  If not
               specified the default version will be used.
  Example    : $coord_sys = $csa->fetch_by_name('clone');
               $coord_sys = $csa->fetch_by_name('chromosome', 'NCBI33');
               # toplevel is an pseudo coord system representing the highest
               # coord system in a given region
               # such as the chromosome coordinate system
               $coord_sys = $csa->fetch_by_name('toplevel');
               #seqlevel is an alias for the sequence level coordinate system
               #such as the clone or contig coordinate system
               $coord_sys = $csa->fetch_by_name('seqlevel');
  Description: Retrieves a coordinate system by its name
  Returntype : Bio::EnsEMBL::CoordSystem
  Exceptions : throw if no name argument provided
               warning if no version provided and default does not exist
  Caller     : general
  Status     : Stable

*/

CoordSystem *CoordSystemAdaptor_fetchByName(CoordSystemAdaptor *csa, char *name, char *version) {
  char lcName[1024];
  char lcVersion[1024];

  strcpy(lcName,name);
  StrUtil_strlwr(lcName);

  if (version != NULL) {
    strcpy(lcVersion,version);
    StrUtil_strlwr(lcVersion);
  }

  if (!strcmp(lcName,"seqlevel")) {
    return CoordSystemAdaptor_fetchSeqLevel(csa);
  } else if (!strcmp(lcName,"toplevel")) {
    return CoordSystemAdaptor_fetchTopLevel(csa);
  }


  if (!StringHash_contains(csa->nameCache,lcName)) {
    return NULL;
  }

  Vector *coordSystems = StringHash_getValue(csa->nameCache, lcName);

  int i;
  for (i=0; i<Vector_getNumElement(coordSystems); i++) {
    CoordSystem *cs = Vector_getElementAt(coordSystems, i);
    if (version && CoordSystem_getVersion(cs)) {
      if (!strcasecmp(CoordSystem_getVersion(cs),version)) {
        return cs;
      }
    } else if (IDHash_contains(csa->isDefaultVersionCache, CoordSystem_getDbID(cs))) {
      return cs;
    }
  }

  if (version) {
    //the specific version we were looking for was not found
    return NULL;
  }

  //didn't find a default, just take first one
  CoordSystem *cs = Vector_getElementAt(coordSystems, 0);
  fprintf(stderr, "No default version for coord_system [%s] exists. Using version [%s] arbitrarily", 
          (name ? name : ""), (cs && CoordSystem_getVersion(cs) ? CoordSystem_getVersion(cs) : ""));

  return cs;
}


/*
=head2 fetch_all_by_name

  Arg [1]    : string $name
               The name of the coordinate system to retrieve.  This can be
               the name of an actual coordinate system or an alias for a
               coordinate system.  Valid aliases are 'toplevel' and 'seqlevel'.
  Example    : foreach my $cs (@{$csa->fetch_all_by_name('chromosome')}){
                 print $cs->name(), ' ', $cs->version();
               }
  Description: Retrieves all coordinate systems of a particular name
  Returntype : listref of Bio::EnsEMBL::CoordSystem objects
  Exceptions : throw if no name argument provided
  Caller     : general
  Status     : Stable

*/

Vector *CoordSystemAdaptor_fetchAllByName(CoordSystemAdaptor *csa, char *name) {
  char lcName[1024];

  strcpy(lcName,name);
  StrUtil_strlwr(lcName);

  if (!strcmp(lcName,"seqlevel")) {
    Vector *v = Vector_new();
    Vector_addElement(v, CoordSystemAdaptor_fetchSeqLevel(csa));
    return v;
  } else if (!strcmp(lcName,"toplevel")) {
    Vector *v = Vector_new();
    Vector_addElement(v, CoordSystemAdaptor_fetchTopLevel(csa));
    return v;
  }

  if (StringHash_contains(csa->nameCache,lcName)) {
    Vector *v = (Vector *)StringHash_getValue(csa->nameCache, lcName);
    return v;
  } else {
    return emptyVector;
  }
}


CoordSystem *CoordSystemAdaptor_fetchByDbID(CoordSystemAdaptor *csa, IDType dbID) {

  CoordSystem *cs = IDHash_getValue(csa->dbIDCache, dbID);

  return cs;
}



/*
=head2 fetch_top_level

  Arg [1]    : none
  Example    : $cs = $csa->fetch_top_level();
  Description: Retrieves the toplevel pseudo coordinate system.
  Returntype : Bio::EnsEMBL::CoordSystem object
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/

CoordSystem *CoordSystemAdaptor_fetchTopLevel(CoordSystemAdaptor *csa) {

  return csa->topLevel;
}


/*
=head2 fetch_sequence_level

  Arg [1]    : none
  Example    : ($id, $name, $version) = $csa->fetch_sequence_level();
  Description: Retrieves the coordinate system at which sequence
               is stored at.
  Returntype : Bio::EnsEMBL::CoordSystem
  Exceptions : throw if no sequence_level coord system exists at all
               throw if multiple sequence_level coord systems exists
  Caller     : general
  Status     : Stable

=cut
*/

CoordSystem *CoordSystemAdaptor_fetchSeqLevel(CoordSystemAdaptor *csa) {

  if (IDHash_getNumValues(csa->isSeqLevelCache) == 0) {  
    fprintf(stderr, "No sequence_level coord_system is defined");
    return NULL;
  }


  if (IDHash_getNumValues(csa->isSeqLevelCache) > 1) {  
    fprintf(stderr,"Multiple sequence_level coord_systems are defined. Only one is currently supported");
    return NULL;
  }
  
  void **values = IDHash_getValues(csa->isSeqLevelCache);
  CoordSystem *cs = values[0];
  free(values);

  return cs;
}


Vector *CoordSystemAdaptor_getMappingPath(CoordSystemAdaptor *csa, CoordSystem *cs1, CoordSystem *cs2) {
//  char key1[1024];
//  char key2[1024];
  char *key1;
  char *key2;
  char keypair[2048];
  char revKeypair[2048];
  Vector *path = NULL;

  if (!cs1 || !cs2)
    return emptyVector ;

//  int lenKey1 = sprintf(key1,"%s:%s", CoordSystem_getName(cs1), CoordSystem_getVersion(cs1) ? CoordSystem_getVersion(cs1) : "");
  
  key1 = CoordSystem_getNameColonVersion(cs1);
  int lenKey1 = CoordSystem_getLenNameColonVersion(cs1);

  key2 = CoordSystem_getNameColonVersion(cs2);
  int lenKey2 = CoordSystem_getLenNameColonVersion(cs2);
//  int lenKey2 = sprintf(key2,"%s:%s", CoordSystem_getName(cs2), CoordSystem_getVersion(cs2) ? CoordSystem_getVersion(cs2) : "");

  memcpy(keypair,key1,lenKey1);
  keypair[lenKey1] = '|';
  memcpy(&keypair[lenKey1+1],key2,lenKey2+1);

  //sprintf(keypair,"%s|%s",key1,key2);
  
  if (StringHash_contains(csa->mappingPaths, keypair)) {
    path = StringHash_getValue(csa->mappingPaths, keypair);
    //fprintf(stderr, "first if path %s length = %d\n", keypair, Vector_getNumElement(path));
    if (path) return path; // Not sure if a NULL path should ever be returned from hash but perl was checking for it
  }

  memcpy(revKeypair,key2,lenKey2);
  revKeypair[lenKey2] = '|';
  memcpy(&revKeypair[lenKey2+1],key1,lenKey1+1);

  //sprintf(revKeypair,"%s|%s",key2,key1); // reverse of path

  if (StringHash_contains(csa->mappingPaths, revKeypair)) {
    path = StringHash_getValue(csa->mappingPaths, revKeypair);
    //fprintf(stderr,"second if path %s length = %d\n", revKeypair, Vector_getNumElement(path));
    if (path) return path; // Not sure if a NULL path should ever be returned from hash but perl was checking for it
  }


  if (!path) {
    // No path was explicitly defined, but we might be able to guess a
    // suitable path.  We only guess for missing 2 step paths.

    StringHash *mid1Hash = StringHash_new(STRINGHASH_SMALL);
    StringHash *mid2Hash = StringHash_new(STRINGHASH_SMALL);

    Vector *mappingPathsVec = StringHash_getValues(csa->mappingPaths);

    int i;
    for (i=0; i<Vector_getNumElement(mappingPathsVec); i++) {

      path = Vector_getElementAt(mappingPathsVec, i);

      if (Vector_getNumElement(path) != 2) {
        continue;
      }

      int match = -1;

      if (CoordSystem_compare(Vector_getElementAt(path,0), cs1)) {
        match = 1;
      } else if (CoordSystem_compare(Vector_getElementAt(path,1), cs1)) {
        match = 0;
      }

      if (match > -1) {
        CoordSystem *mid = Vector_getElementAt(path, match);
        char midKey[1024];
        sprintf(midKey,"%s:%s", CoordSystem_getName(mid), CoordSystem_getVersion(mid) ? CoordSystem_getVersion(mid) : "");


        // is the same cs mapped to by other cs?
        if (StringHash_contains(mid2Hash, midKey)) {
          // Make new three way path
          path = Vector_new();
          Vector_addElement(path, cs1); Vector_addElement(path, mid); Vector_addElement(path, cs2);

          // Add to mapping
          StringHash_add(csa->mappingPaths, keypair, path);

          if (!CoordSystem_getVersion(cs1)) key1[strlen(key1)-1]     = '\0';
          if (!CoordSystem_getVersion(cs2)) key2[strlen(key2)-1]     = '\0';
          if (!CoordSystem_getVersion(mid)) midKey[strlen(midKey)-1] = '\0';

          fprintf(stderr, "Using implicit mapping path between '%s' and '%s' coord systems.\n"
                          "An explicit 'assembly.mapping' entry should be added to the meta table.\nExample: '%s|%s|%s'\n", 
                          key1, key2, key1, midKey, key2);

          StringHash_free(mid1Hash, NULL);
          StringHash_free(mid2Hash, NULL);

          return path;
        } else {
          StringHash_add(mid1Hash, midKey, mid);
        }
      }

      match = -1;

      if (CoordSystem_compare(Vector_getElementAt(path,0), cs2)) {
        match = 1;
      } else if (CoordSystem_compare(Vector_getElementAt(path,1), cs2)) {
        match = 0;
      }

      if (match > -1) {
        CoordSystem *mid = Vector_getElementAt(path, match);
        char midKey[1024];
        sprintf(midKey,"%s:%s", CoordSystem_getName(mid), CoordSystem_getVersion(mid) ? CoordSystem_getVersion(mid) : "");


        // is the same cs mapped to by other cs?
        if (StringHash_contains(mid1Hash, midKey)) {
          // Make new three way path
          path = Vector_new();
          Vector_addElement(path, cs2); Vector_addElement(path, mid); Vector_addElement(path, cs1);

          // Add to mapping
          StringHash_add(csa->mappingPaths, revKeypair, path);

          if (!CoordSystem_getVersion(cs1)) key1[strlen(key1)-1]     = '\0';
          if (!CoordSystem_getVersion(cs2)) key2[strlen(key2)-1]     = '\0';
          if (!CoordSystem_getVersion(mid)) midKey[strlen(midKey)-1] = '\0';

          fprintf(stderr, "Using implicit mapping path between '%s' and '%s' coord systems.\n"
                          "An explicit 'assembly.mapping' entry should be added to the meta table.\nExample: '%s|%s|%s'\n", 
                          key1, key2, key1, midKey, key2);

          StringHash_free(mid1Hash, NULL);
          StringHash_free(mid2Hash, NULL);

          return path;
        } else {
          StringHash_add(mid2Hash, midKey, mid);
        }
      }
    }
  }

  if (path) {
    return path;
  } else {
    return emptyVector;
  }
}

/* StoreMappingPath not implemented
=head2 store_mapping_path

  Arg [1]    : Bio::EnsEMBL::CoordSystem $cs1
  Arg [2]    : Bio::EnsEMBL::CoordSystem $cs2
  Arg [3..n] : Bio::EnsEMBL::CoordSystem $cs3..$csN
  Example    : my $pathref = $csa->store_mapping_path($cs1,$cs2);
  Description: Given two or more coordinate systems this will store 
               mapping paths between them in the database. 

               The 'rank' attrib of the CoordSystems is used to
               determine the assembled/component relationships between
               them.

               For example, if $cs1 represents chrs of version
               V1, $cs2 represents contigs, and $cs3 clones then, unless
               they already exist, the following entries will be created 
               in the meta table;
               +------------------+---------------------+ 
               | meta_key         | meta_value          |
               +------------------+---------------------+ 
               | assembly.mapping | chr:V1|clone        |
               | assembly.mapping | clone|contig        |
               | assembly.mapping | chr:V1|clone|contig |
               +------------------+---------------------+


               For a one-step mapping path to be valid there needs to be
               a relationship between the two coordinate systems defined in
               the assembly table.  Two step mapping paths work by building
               on the one-step mapping paths which are already defined.

               The first coordinate system in a one step mapping path must
               be the assembled coordinate system and the second must be
               the component.

  Returntype : reference to a list of lists of new meta_value mapping strings
               created for assembly.mapping
  Exceptions : CoordSystems with no rank/duplicated rank
  Caller     : general
  Status     : Experimental

=cut

sub store_mapping_path{
  my $self = shift;
  my @csystems = @_;

  # Validate and sort the args
  my %seen_ranks;
  @csystems >= 2 or throw('Need two or more CoordSystems');
  my $validate = sub{ 
    ref($_[0]) && $_[0]->isa('Bio::EnsEMBL::CoordSystem') or
        throw('CoordSystem argument expected.');
    my $rank = $_[0]->rank || 
        throw('CoordSystem has no rank: '.$_[0]->name);
    $seen_ranks{$rank} &&
        throw('CoordSystem '.$_[0]->name." shares rank $rank with ".
              $seen_ranks{$rank}->name);
    $seen_ranks{$rank} = $_[0];
  };
  @csystems = sort{$a->rank <=> $b->rank} map{&{$validate}($_)} @csystems;

  # Get a list of all existing assembly.mappings
  #my %mappings = map{$_=>1} @{$meta->list_value_by_key('assembly.mapping')};
  
  # For each pair in the sorted list, store in the DB
  my $meta = $self->db->get_MetaContainer;
  my @retlist;
  for( my $i=1; $i<@csystems; $i++ ){
    for( my $j=0; $j<(@csystems-$i); $j++ ){
      my $mapping = join( "|", 
                          map{join( ':', $_->name, ($_->version||()) )} 
                          @csystems[$j..$j+$i] );
      my $mapping_key = join( "|",
                              map{join( ':', $_->name, ($_->version||'') )}
                              @csystems[$j..$j+$i] );
      # Skip existing
      next if $self->{'_mapping_paths'}->{$mapping_key};
      
      # Update the database
      $meta->store_key_value('assembly.mapping',$mapping);
      push @retlist, $mapping;
    }
  }

  if( @retlist ){
    # Update mapping path cache
    $self->_cache_mapping_paths;
  }

  # Return the mappings that we have just created
  return [@retlist];
}
*/

/*
=head2 fetch_by_attrib

  Arg [1]    : string attrib
  Arg [2]    : (optional) string version
  Example    : $csa->fetch_by_attrib('default_version','NCBIM37');
  Description: Retrieves a CoordSystem object from the database that have the specified
               attrib and version, if no version is specified, returns the default version
  Returntype : Bio::EnsEMBL::CoordSystem object
  Exceptions : throw when attrib not present
  Caller     : general
  Status     : Stable

=cut
*/

CoordSystem *CoordSystemAdaptor_fetchByAttrib(CoordSystemAdaptor *csa, char *attrib, char *version) {
  IDHash *attribCache;

  if ( !strcmp(attrib, "sequence_level") ) {
    attribCache = csa->isSeqLevelCache;
  } else if ( !strcmp(attrib, "default_version")) {
    attribCache = csa->isDefaultVersionCache;
  }

  
  if (!IDHash_getNumValues(attribCache)) {
    fprintf(stderr, "No %s coordinate system defined\n",attrib);
    return NULL;
  }

  int i=0;
  void **values = IDHash_getValues(attribCache);
  for (i=0; i< IDHash_getNumValues(attribCache); i++) {
    CoordSystem *cs = values[i];

    if (version) {
      if (!strcasecmp(version, CoordSystem_getVersion(cs))) {
        free(values);
        return cs;
      }
    } else if (IDHash_contains(csa->isDefaultVersionCache, CoordSystem_getDbID(cs))) {
      free(values);
      return cs;
    }
  }

  //specifically requested attrib system was not found
  if (version) {
    fprintf(stderr,"%s coord_system with version [%s] does not exist\n", attrib, version);
    return NULL;
  }

  CoordSystem *cs = values[0];
  fprintf(stderr, "No default version for %s coord_system exists. Using version [%s] arbitrarily\n",
                  attrib, CoordSystem_getVersion(cs));

  free(values);
  return cs;
}


/*
=head2 fetch_all_by_attrib

  Arg [1]    : string attrib
  Example    : $csa->fetch_all_by_attrib('default_version');
  Description: Retrieves all CoordSystem object from the database that have the specified
               attrib.
  Returntype : reference to a list of Bio::EnsEMBL::CoordSystem objects
  Exceptions : throw when attrib not present
  Caller     : general
  Status     : Stable

=cut
*/

Vector *CoordSystemAdaptor_fetchAllByAttrib(CoordSystemAdaptor *csa, char *attrib) {
  IDHash *attribCache;

  if ( !strcmp(attrib, "sequence_level") ) {
    attribCache = csa->isSeqLevelCache;
  } else if ( !strcmp(attrib, "default_version")) {
    attribCache = csa->isDefaultVersionCache;
  }

  Vector *coordSystems = Vector_newFromArray(IDHash_getValues(attribCache), IDHash_getNumValues(attribCache));

  return coordSystems;
}

/* Store not implemented

=head2 store

  Arg [1]    : Bio::EnsEMBL::CoordSystem
  Example    : $csa->store($coord_system);
  Description: Stores a CoordSystem object in the database.
  Returntype : none
  Exceptions : Warning if CoordSystem is already stored in this database.
  Caller     : none
  Status     : Stable

=cut

sub store {
  my $self = shift;
  my $cs = shift;

  if(!$cs || !ref($cs) || !$cs->isa('Bio::EnsEMBL::CoordSystem')) {
    throw('CoordSystem argument expected.');
  }

  my $db = $self->db();
  my $name = $cs->name();
  my $version = $cs->version();
  my $rank    = $cs->rank();

  my $seqlevel = $cs->is_sequence_level();
  my $default  = $cs->is_default();

  my $toplevel = $cs->is_top_level();

  if($toplevel) {
    throw("The toplevel CoordSystem cannot be stored");
  }

  #
  # Do lots of sanity checking to prevent bad data from being entered
  #

  if($cs->is_stored($db)) {
    warning("CoordSystem $name $version is already in db.\n");
    return;
  }

  if($name eq 'toplevel' || $name eq 'seqlevel' || !$name) {
    throw("[$name] is not a valid name for a CoordSystem.");
  }

  if($seqlevel && keys(%{$self->{'_is_sequence_level'}})) {
    throw("There can only be one sequence level CoordSystem.");
  }

  if(exists $self->{'_name_cache'}->{lc($name)}) {
    my @coord_systems = @{$self->{'_name_cache'}->{lc($name)}};
    foreach my $c (@coord_systems) {
      if(lc($c->version()) eq lc($version)) {
        warning("CoordSystem $name $version is already in db.\n");
        return;
      }
      if($default && $self->{'_is_default_version'}->{$c->dbID()}) {
        throw("There can only be one default version of CoordSystem $name");
      }
    }
  }

  if($rank !~ /^\d+$/) {
    throw("Rank attribute must be a positive integer not [$rank]");
  }
  if($rank == 0) {
    throw("Only toplevel CoordSystem may have rank of 0.");
  }

  if(defined($self->{'_rank_cache'}->{$rank})) {
    throw("CoordSystem with rank [$rank] already exists.");
  }

  my @attrib;

  push @attrib, 'default_version' if($default);
  push @attrib, 'sequence_level' if($seqlevel);

  my $attrib_str = (@attrib) ? join(',', @attrib) : undef;

  #
  # store the coordinate system in the database
  #

  my $sth =
    $db->dbc->prepare(   'INSERT INTO coord_system '
                       . 'SET name = ?, '
                       . 'version = ?, '
                       . 'attrib = ?,'
                       . 'rank = ?,'
                       . 'species_id = ?' );

  $sth->bind_param( 1, $name,               SQL_VARCHAR );
  $sth->bind_param( 2, $version,            SQL_VARCHAR );
  $sth->bind_param( 3, $attrib_str,         SQL_VARCHAR );
  $sth->bind_param( 4, $rank,               SQL_INTEGER );
  $sth->bind_param( 5, $self->species_id(), SQL_INTEGER );

  $sth->execute();
  my $dbID = $sth->{'mysql_insertid'};
  $sth->finish();

  if(!$dbID) {
    throw("Did not get dbID from store of CoordSystem.");
  }

  $cs->dbID($dbID);
  $cs->adaptor($self);

  #
  # update the internal caches that are used for fetching
  #
  $self->{'_is_default_version'}->{$dbID} = 1 if($default);
  $self->{'_is_sequence_level'}->{$dbID} = 1 if($seqlevel);

  $self->{'_name_cache'}->{lc($name)} ||= [];
  push @{$self->{'_name_cache'}->{lc($name)}}, $cs;

  $self->{'_dbID_cache'}->{$dbID} = $cs;
  $self->{'_rank_cache'}->{$rank} = $cs;

  return $cs;
}

*/
