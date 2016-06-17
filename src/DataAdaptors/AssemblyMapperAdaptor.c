/*
 * Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

#include "AssemblyMapperAdaptor.h"
#include "BaseAdaptor.h"
#include "MysqlUtil.h"
#include "AssemblyMapper.h"
#include "ChainedAssemblyMapper.h"
#include "TopLevelAssemblyMapper.h"
#include "Slice.h"
#include "DBAdaptor.h"

#include "StatementHandle.h"
#include "MysqlStatementHandle.h"
#include "CoordSystemAdaptor.h"
#include "SeqRegionCacheEntry.h"
#include "ResultRow.h"
#include "CoordPair.h"
#include "StrUtil.h"

int CHUNKFACTOR = 20;  // 2^20 = approx. 10^6


static char *AMA_FIRST = "first";
static char *AMA_LAST  = "last";

AssemblyMapperAdaptor *AssemblyMapperAdaptor_new(DBAdaptor *dba) {
  AssemblyMapperAdaptor *ama;

  if ((ama = (AssemblyMapperAdaptor *)calloc(1,sizeof(AssemblyMapperAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for AssemblyMapperAdaptor\n");
    return NULL;
  }
  BaseAdaptor_init((BaseAdaptor *)ama, dba, ASSEMBLYMAPPER_ADAPTOR);

  ama->asmMapperCache = StringHash_new(STRINGHASH_SMALL);

  // use a shared cache (for this database) that contains info about
  // seq regions

  ama->srNameCache = DBAdaptor_getSeqRegionNameCache(ama->dba);
  ama->srIdCache   = DBAdaptor_getSeqRegionIdCache(ama->dba);

  ama->multSeqIdCache = NULL;

  return ama;
}


/*
=head2  cache_seq_ids_with_mult_assemblys

  Example    : $self->adaptor->cache_seq_ids_with_mult_assemblys();
  Description: Creates a hash of the component seq region ids that
               map to more than one assembly from the assembly table.
  Retruntype : none
  Exceptions : none
  Caller     : AssemblyMapper, ChainedAssemblyMapper
  Status     : At Risk

=cut
*/

void AssemblyMapperAdaptor_cacheSeqIdsWithMultAssemblies(AssemblyMapperAdaptor *ama) {
  if (ama->multSeqIdCache) {
    return;
  }

  ama->multSeqIdCache = IDHash_new(IDHASH_LARGE);

  IDType speciesID = 1;

  char qStr[1024];
  sprintf(qStr, "SELECT    sra.seq_region_id"
                " FROM      seq_region_attrib sra,"
                  " attrib_type at,"
                  " seq_region sr,"
                  " coord_system cs"
                " WHERE     sra.attrib_type_id = at.attrib_type_id"
                  " AND     code = 'MultAssem'"
                  " AND     sra.seq_region_id = sr.seq_region_id"
                  " AND     sr.coord_system_id = cs.coord_system_id"
                  " AND     cs.species_id = "
                  IDFMTSTR, speciesID);

  StatementHandle *sth = ama->prepare((BaseAdaptor *)ama,qStr,strlen(qStr));
  sth->execute(sth);

  ResultRow *row;
  while ((row = sth->fetchRow(sth))) {
    IDType dbID    = row->getLongLongAt(row,0);

    IDHash_add(ama->multSeqIdCache, dbID, &trueVal);
  }

  sth->finish(sth);
}


char *makeMappingPathKey(Vector *path, char *key) {
  int i;
//  char *key;

//  StrUtil_copyString(&key, "", 0);
  key[0] = '\0';
  

  int pos=0;
  for (i=0; i<Vector_getNumElement(path); i++) {
    CoordSystem *cs = Vector_getElementAt(path, i);
    if (cs == NULL) {
      key[pos++] = '-';
      //key = StrUtil_appendString(key, "-");
    } else {
      char tmp[1024];
      int lenNum = sprintf(tmp,IDFMTSTR,CoordSystem_getDbID(cs));
      //key = StrUtil_appendString(key, tmp );
      memcpy(&key[pos],tmp,lenNum);
      pos+=lenNum;
    }
    if (i < Vector_getNumElement(path)-1) {
      //key = StrUtil_appendString(key, ":");
      key[pos++] = ':';
    }
  }
  key[pos] = '\0';
  //fprintf(stderr, "Made mapping path key %s\n", key);
  
  return key;
}

/*
=head2 fetch_by_CoordSystems

  Arg [1]    : Bio::EnsEMBL::CoordSystem $cs1
               One of the coordinate systems to retrieve the mapper
               between
  Arg [2]    : Bio::EnsEMBL::CoordSystem $cs2
               The other coordinate system to map between
  Description: Retrieves an Assembly mapper for two coordinate
               systems whose relationship is described in the
               assembly table.

               The ordering of the coodinate systems is arbitrary.
               The following two statements are equivalent:
               $mapper = $asma->fetch_by_CoordSystems($cs1,$cs2);
               $mapper = $asma->fetch_by_CoordSystems($cs2,$cs1);
  Returntype : Bio::EnsEMBL::AssemblyMapper
  Exceptions : wrong argument types
  Caller     : general
  Status     : Stable

=cut
*/

AssemblyMapper *AssemblyMapperAdaptor_fetchByCoordSystems(AssemblyMapperAdaptor *ama, CoordSystem *cs1, CoordSystem *cs2) {
  if (cs1 && CoordSystem_getIsTopLevel(cs1)) {
    return (AssemblyMapper *)TopLevelAssemblyMapper_new(ama, cs1, cs2);
  }

  if (cs2 && CoordSystem_getIsTopLevel(cs2)) {
    return (AssemblyMapper *)TopLevelAssemblyMapper_new(ama, cs2, cs1);
  }

  //printf("Srnamecache size at start of fbcs = %d\n",StringHash_getNumValues(ama->srNameCache));

  CoordSystemAdaptor *csa = DBAdaptor_getCoordSystemAdaptor(ama->dba);

  //retrieve the shortest possible mapping path between these systems
  Vector *mappingPath = CoordSystemAdaptor_getMappingPath(csa,cs1,cs2);

  if (!Vector_getNumElement(mappingPath)) {

    // It is perfectly fine not to have a mapping. No warning needed really
    // Just check the return code!!

//#    warning(
//#      "There is no mapping defined between these coord systems:\n" .
//#      $cs1->name() . " " . $cs1->version() . " and " . $cs2->name() . " " .
//#      $cs2->version()
//#    );
    return NULL;
  }

  //char *key = makeMappingPathKey(mappingPath);
  char key[2048];
  makeMappingPathKey(mappingPath, key);

  if (StringHash_contains(ama->asmMapperCache, key)) {
    return StringHash_getValue(ama->asmMapperCache, key);
  }

  switch (Vector_getNumElement(mappingPath)) {
    case 1:
      fprintf(stderr,"Incorrect mapping path defined in meta table. 0 step mapping encountered between:\n"
                     "%s %s and %s %s\n", 
                     CoordSystem_getName(cs1), CoordSystem_getVersion(cs1),
                     CoordSystem_getName(cs2), CoordSystem_getVersion(cs2));
      break;

    case 2:
      {
        // 1 step regular mapping
        //printf("Making a normal AssemblyMapper for path %s\n",key);
        AssemblyMapper *asmMapper = AssemblyMapper_new(ama, mappingPath);
  
        //   If you want multiple pieces on two seqRegions to map to each other
        //   you need to make an assembly.mapping entry that is seperated with a #
        //   instead of an |.
  
        StringHash_add(ama->asmMapperCache, key, asmMapper);
  
        //free(key);
  
        return asmMapper;
      }
      break;
  
    case 3:
      // two step chained mapping
      {
        //printf("Making a ChainedAssemblyMapper for path %s\n",key);
        ChainedAssemblyMapper *casmMapper = ChainedAssemblyMapper_new(ama, mappingPath);
  
        // in multi-step mapping it is possible get requests with the
        // coordinate system ordering reversed since both mappings directions
        // cache on both orderings just in case
        // e.g.   chr <-> contig <-> clone   and   clone <-> contig <-> chr
    
        StringHash_add(ama->asmMapperCache, key, casmMapper);
        //free(key);
  
        // Make a reverse COPY of the mapping path (as its a copy I'm not using Vector_reverse which is inplace reverse 
        Vector *revMappingPath = Vector_new();
        int i;
        for (i=Vector_getNumElement(mappingPath)-1; i>=0; i--) {
          Vector_addElement(revMappingPath, Vector_getElementAt(mappingPath, i));
        }
        //key = makeMappingPathKey(revMappingPath);
        makeMappingPathKey(revMappingPath, key);
        StringHash_add(ama->asmMapperCache, key, casmMapper);
  
        Vector_free(revMappingPath);
        //free(key);
  
        return (AssemblyMapper *)casmMapper;
      }
      break;

    default:
      fprintf(stderr, "Only 1 and 2 step coordinate system mapping is currently\n" 
                      "supported.  Mapping between %s %s and %s %s requires %d steps.\n",
                     CoordSystem_getName(cs1), CoordSystem_getVersion(cs1),
                     CoordSystem_getName(cs2), CoordSystem_getVersion(cs2),
                     Vector_getNumElement(mappingPath));
  }

  return NULL;
}



/*
=head2 register_assembled

  Arg [1]    : Bio::EnsEMBL::AssemblyMapper $asm_mapper
               A valid AssemblyMapper object
  Arg [2]    : integer $asm_seq_region
               The dbID of the seq_region to be registered
  Arg [3]    : int $asm_start
               The start of the region to be registered
  Arg [4]    : int $asm_end
               The end of the region to be registered
  Description: Declares an assembled region to the AssemblyMapper.
               This extracts the relevant data from the assembly
               table and stores it in Mapper internal to the $asm_mapper.
               It therefore must be called before any mapping is
               attempted on that region. Otherwise only gaps will
               be returned.  Note that the AssemblyMapper automatically
               calls this method when the need arises.
  Returntype : none
  Exceptions : throw if the seq_region to be registered does not exist
               or if it associated with multiple assembled pieces (bad data
               in assembly table)
  Caller     : Bio::EnsEMBL::AssemblyMapper
  Status     : Stable

=cut
*/

SeqRegionRange *AssemblyMapperAdaptor_addToRangeVector(Vector *ranges, IDType id, long start, long end, char *name) {
  SeqRegionRange *range = NULL;

  if ((range = (SeqRegionRange *)calloc(1,sizeof(SeqRegionRange))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for SeqRegionRange\n");
  } else {
    SeqRegionRange_setSeqRegionStart(range, start);
    SeqRegionRange_setSeqRegionEnd(range, end);
    SeqRegionRange_setSeqRegionId(range, id);

    if (name) SeqRegionRange_setSeqRegionName(range, name);

    Vector_addElement(ranges, range);
  }

  return range; 
}

void AssemblyMapperAdaptor_registerAssembled(AssemblyMapperAdaptor *ama, AssemblyMapper *asmMapper, IDType asmSeqRegion, long asmStart, long asmEnd) {

  IDType cmpCsId = CoordSystem_getDbID(AssemblyMapper_getComponentCoordSystem(asmMapper));

  //split up the region to be registered into fixed chunks
  //this allows us to keep track of regions that have already been
  //registered and also works under the assumption that if a small region
  //is requested it is likely that other requests will be made in the
  //vicinity (the minimum size registered the chunksize (2^chunkfactor)

  Vector *chunkRegions = Vector_new();
  //determine span of chunks
  //bitwise shift right is fast and easy integer division

  long startChunk, endChunk;

  startChunk = asmStart >> CHUNKFACTOR;
  endChunk   = asmEnd   >> CHUNKFACTOR;

  // inserts have start = end + 1, on boundary condition start_chunk
  // could be less than end chunk
  if(asmStart == asmEnd + 1) {
    long tmp   = startChunk;
    startChunk = endChunk;
    endChunk   = tmp;
  }

  //find regions of continuous unregistered chunks
  long i;
  long beginChunkRegion, endChunkRegion;
  int begun = 0;

  for (i=startChunk; i<=endChunk; i++) {
    if (AssemblyMapper_haveRegisteredAssembled(asmMapper, asmSeqRegion, i)) {
      if (begun) {
        //this is the end of an unregistered region.
        AssemblyMapperAdaptor_addToRangeVector(chunkRegions, 0,  
                                               (beginChunkRegion << CHUNKFACTOR), 
                                               (((endChunkRegion+1) << CHUNKFACTOR)-1), NULL);

        //$begin_chunk_region = $end_chunk_region = undef;
        begun = 0;
      }
    } else {
      if (!begun) beginChunkRegion = i;
      begun = 1;
      endChunkRegion = i+1;
      AssemblyMapper_registerAssembled(asmMapper, asmSeqRegion, i);
    }
  }

  //the last part may have been an unregistered region too
  //if(defined($begin_chunk_region)) 
  if (begun) {
    AssemblyMapperAdaptor_addToRangeVector(chunkRegions, 0,  
                                           (beginChunkRegion << CHUNKFACTOR), 
                                           (((endChunkRegion+1) << CHUNKFACTOR)-1), NULL);
  }

  if (Vector_getNumElement(chunkRegions) == 0) {
    Vector_free(chunkRegions);
    return;
  }

  // keep the Mapper to a reasonable size
  if (AssemblyMapper_getSize(asmMapper) > AssemblyMapper_getMaxPairCount(asmMapper) ) {
    AssemblyMapper_flush(asmMapper);

    // we now have to go and register the entire requested region since we
    // just flushed everything
    chunkRegions->freeElement = SeqRegionRange_free;
    Vector_free(chunkRegions);
    chunkRegions = Vector_new();

    AssemblyMapperAdaptor_addToRangeVector(chunkRegions, 0,  
                                           (startChunk << CHUNKFACTOR), 
                                           (((endChunk+1) << CHUNKFACTOR)-1), NULL);

    long i;
    for (i=startChunk; i <= endChunk; i++ ) {
      AssemblyMapper_registerAssembled(asmMapper, asmSeqRegion, i);
    }
  }

  //  my $asm_seq_region_id =
  //    $self->_seq_region_name_to_id($asm_seq_region,$asm_cs_id);

  // Retrieve the description of how the assembled region is made from
  // component regions for each of the continuous blocks of unregistered,
  // chunked regions


  for (i=0; i<Vector_getNumElement(chunkRegions); i++) {
    SeqRegionRange *region = Vector_getElementAt(chunkRegions, i);

    long regionStart = SeqRegionRange_getSeqRegionStart(region);
    long regionEnd   = SeqRegionRange_getSeqRegionEnd(region);

    char qStr[1024];
    sprintf(qStr, "SELECT"
                    " asm.cmp_start,"
                    " asm.cmp_end,"
                    " asm.cmp_seq_region_id,"
                    " sr.name,"
                    " sr.length,"
                    " asm.ori,"
                    " asm.asm_start,"
                    " asm.asm_end"
                   " FROM"
                   "  assembly asm, seq_region sr"
                   " WHERE asm.asm_seq_region_id = "
                     IDFMTSTR
                    " AND %ld <= asm.asm_end"
                    " AND %ld >= asm.asm_start"
                    " AND asm.cmp_seq_region_id = sr.seq_region_id"
                    " AND sr.coord_system_id = "IDFMTSTR, 
                   asmSeqRegion, regionStart, regionEnd, cmpCsId);
  

    StatementHandle *sth = ama->prepare((BaseAdaptor *)ama,qStr,strlen(qStr));
    sth->execute(sth);


    // 
    // Load the unregistered regions of the mapper
    // 
    ResultRow *row;
    while ((row = sth->fetchRow(sth))) {
      long   cmpStart           = row->getLongAt(row,0);
      long   cmpEnd             = row->getLongAt(row,1);
      IDType cmpSeqRegionId     = row->getLongLongAt(row,2);
      char * cmpSeqRegionName   = row->getStringAt(row,3);
      long   cmpSeqRegionLength = row->getLongAt(row,4);
      int    ori                = row->getIntAt(row,5);
      long   asmStart           = row->getLongAt(row,6);
      long   asmEnd             = row->getLongAt(row,7);

      if (AssemblyMapper_haveRegisteredComponent(asmMapper, cmpSeqRegionId) &&
          !IDHash_contains(ama->multSeqIdCache, cmpSeqRegionId)) {
        continue;
      }
      AssemblyMapper_registerComponent(asmMapper, cmpSeqRegionId);

      Mapper *mapper = AssemblyMapper_getMapper(asmMapper);

      Mapper_addMapCoordinates(mapper,
                 asmSeqRegion, asmStart, asmEnd, ori,
                 cmpSeqRegionId, cmpStart, cmpEnd);

      DBAdaptor_addToSrCaches(ama->dba, cmpSeqRegionId, cmpSeqRegionName, cmpCsId, cmpSeqRegionLength);
    }
    sth->finish(sth);
  }

  chunkRegions->freeElement = SeqRegionRange_free;
  Vector_free(chunkRegions);

}

IDType AssemblyMapperAdaptor_seqRegionNameToId(AssemblyMapperAdaptor *ama, char *srName, IDType csId) {
  if (srName == NULL || csId <= 0) {
    fprintf(stderr, "seq_region_name and coord_system_id args are required\n");
    return 0;
  }

  char key[1024];
  sprintf(key,"%s:"IDFMTSTR, srName, csId);
  if (StringHash_contains(ama->srNameCache, key)) {
    SeqRegionCacheEntry *cacheData = StringHash_getValue(ama->srNameCache, key);
    return cacheData->regionId;
  }

  // Get the seq_region_id via the name.  This would be quicker if we just
  // used internal ids instead but stored but then we lose the ability
  // the transform across databases with different internal ids

  char qStr[1024];
  sprintf(qStr, "SELECT seq_region_id, length "
                  "FROM seq_region "
                 "WHERE name = '%s' AND coord_system_id = "
                 IDFMTSTR, srName, csId);

  StatementHandle *sth = ama->prepare((BaseAdaptor *)ama,qStr,strlen(qStr));
  sth->execute(sth);

  if (sth->numRows(sth) != 1) {
    fprintf(stderr,"Ambiguous or non-existant seq_region [%s] in coord system "IDFMTSTR" (numRowss = "IDFMTSTR")\n", srName,csId,sth->numRows(sth));
    return 0;
  }

  
  ResultRow *row = sth->fetchRow(sth);
  IDType srId   = row->getLongLongAt(row,0);
  long srLength = row->getLongAt(row,1);


  DBAdaptor_addToSrCaches(ama->dba, srId, srName, csId, srLength);

  sth->finish(sth);
  return srId;
}

char *AssemblyMapperAdaptor_seqRegionIdToName(AssemblyMapperAdaptor *ama, IDType srId) {

  if (!srId) {
    fprintf(stderr,"seq_region_id is required");
    return NULL;
  }

  if (IDHash_contains(ama->srIdCache, srId)) {
    SeqRegionCacheEntry *cacheData = IDHash_getValue(ama->srIdCache, srId);
    return cacheData->regionName;
  }

  // Get the seq_region name via the id.  This would be quicker if we just
  // used internal ids instead but stored but then we lose the ability
  // the transform across databases with different internal ids
  char qStr[1024];
  sprintf(qStr, "SELECT name, length, coord_system_id "
                  "FROM seq_region "
                 "WHERE seq_region_id = "
                 IDFMTSTR, srId);

  StatementHandle *sth = ama->prepare((BaseAdaptor *)ama,qStr,strlen(qStr));
  sth->execute(sth);


  if (sth->numRows(sth) != 1) {
    fprintf(stderr, "non-existant seq_region ["IDFMTSTR"]\n", srId);
    return NULL;
  }

  ResultRow *row = sth->fetchRow(sth);
//NIY: Do I need to alloc this?
  char *srName  = row->getStringAt(row,0);
  long srLength = row->getLongAt(row,1);
  IDType csId   = row->getLongLongAt(row,2);


  DBAdaptor_addToSrCaches(ama->dba, srId, srName, csId, srLength);

  sth->finish(sth);
  return srName;
}

/*
=head2 register_component

  Arg [1]    : Bio::EnsEMBL::AssemblyMapper $asm_mapper
               A valid AssemblyMapper object
  Arg [2]    : integer $cmp_seq_region
               The dbID of the seq_region to be registered
  Description: Declares a component region to the AssemblyMapper.
               This extracts the relevant data from the assembly
               table and stores it in Mapper internal to the $asm_mapper.
               It therefore must be called before any mapping is
               attempted on that region. Otherwise only gaps will
               be returned.  Note that the AssemblyMapper automatically
               calls this method when the need arises.
  Returntype : none
  Exceptions : throw if the seq_region to be registered does not exist
               or if it associated with multiple assembled pieces (bad data
               in assembly table)
  Caller     : Bio::EnsEMBL::AssemblyMapper
  Status     : Stable

=cut
*/

void AssemblyMapperAdaptor_registerComponent(AssemblyMapperAdaptor *ama, AssemblyMapper *asmMapper, IDType cmpSeqRegion) {

  IDType asmCsId = CoordSystem_getDbID(AssemblyMapper_getAssembledCoordSystem(asmMapper));

  // do nothing if this region is already registered or special case

  if (AssemblyMapper_haveRegisteredComponent(asmMapper, cmpSeqRegion) &&
      !IDHash_contains(ama->multSeqIdCache, cmpSeqRegion)) {
    return;
  }

//  my $cmp_seq_region_id =
//    $self->_seq_region_name_to_id($cmp_seq_region, $cmp_cs_id);

  // Determine what part of the assembled region this component region makes up

  char qStr[1024];
  sprintf(qStr, "SELECT "
                  " asm.asm_start,"
                  " asm.asm_end,"
                  " asm.asm_seq_region_id,"
                  " sr.name,"
                  " sr.length"
              " FROM"
                  " assembly asm, seq_region sr"
             " WHERE"
                  " asm.cmp_seq_region_id = " IDFMTSTR " AND"
                  " asm.asm_seq_region_id = sr.seq_region_id AND"
                  " sr.coord_system_id = " IDFMTSTR,
                 cmpSeqRegion, asmCsId);

  StatementHandle *sth = ama->prepare((BaseAdaptor *)ama,qStr,strlen(qStr));
  sth->execute(sth);

  if (sth->numRows(sth) == 0) {
    //this component is not used in the assembled part i.e. gap
    // NIY: Chained also??
    AssemblyMapper_registerComponent(asmMapper, cmpSeqRegion);
    return;
  }

  //we do not currently support components mapping to multiple assembled
  // make sure that you've got the correct mapping in the meta-table :
  //   chromosome:EquCab2#contig ( use'#' for multiple mappings )
  //   chromosome:EquCab2|contig ( use '|' delimiter for 1-1 mappings )
  //
  if (sth->numRows(sth) != 1) {
    fprintf(stderr,"Multiple assembled regions for single component region cmp_seq_region_id="IDFMTSTR"\n"
                   "Remember that multiple mappings use the #-operaator in the meta-table (i.e. chromosome:EquCab2#contig\n",
                   cmpSeqRegion);
    return;
  }

  ResultRow *row = sth->fetchRow(sth);
  long asmStart           = row->getLongAt(row,0);
  long asmEnd             = row->getLongAt(row,1);
  IDType asmSeqRegionId   = row->getLongLongAt(row,2);
//NIY: Do I need to alloc this?
  char *asmSeqRegion      = row->getStringAt(row,3);
  long asmSeqRegionLength = row->getLongAt(row,4);

  DBAdaptor_addToSrCaches(ama->dba, asmSeqRegionId, asmSeqRegion, asmCsId, asmSeqRegionLength);


  // Register the corresponding assembled region. This allows a us to
  // register things in assembled chunks which allows us to:
  // (1) Keep track of what assembled regions are registered
  // (2) Use locality of reference (if they want something in same general
  //     region it will already be registered).

  AssemblyMapperAdaptor_registerAssembled(ama, asmMapper, asmSeqRegionId, asmStart, asmEnd);
  sth->finish(sth);
}


/*
=head2 register_chained

  Arg [1]    : Bio::EnsEMBL::ChainedAssemblyMapper $casm_mapper
               The chained assembly mapper to register regions on
  Arg [2]    : string $from ('first' or 'last')
               The direction we are registering from, and the name of the
               internal mapper.
  Arg [3]    : string $seq_region_name
               The name of the seqregion we are registering on
  Arg [4]    : listref $ranges
               A list  of ranges to register (in [$start,$end] tuples).
  Arg [5]    : (optional) $to_slice
               Only register those on this Slice.
  Description: Registers a set of ranges on a chained assembly mapper.
               This function is at the heart of the chained mapping process.
               It retrieves information from the assembly table and
               dynamically constructs the mappings between two coordinate
               systems which are 2 mapping steps apart. It does this by using
               two internal mappers to load up a third mapper which is
               actually used by the ChainedAssemblyMapper to perform the
               mapping.

               This method must be called before any mapping is
               attempted on regions of interest, otherwise only gaps will
               be returned.  Note that the ChainedAssemblyMapper automatically
               calls this method when the need arises.
  Returntype : none
  Exceptions : throw if the seq_region to be registered does not exist
               or if it associated with multiple assembled pieces (bad data
               in assembly table)

               throw if the mapping between the coordinate systems cannot
               be performed in two steps, which means there is an internal
               error in the data in the meta table or in the code that creates
               the mapping paths.
  Caller     : Bio::EnsEMBL::AssemblyMapper
  Status     : Stable

=cut
*/

void AssemblyMapperAdaptor_registerChained(AssemblyMapperAdaptor *ama, ChainedAssemblyMapper *casmMapper, 
                                           char *from, IDType seqRegionId, Vector *ranges, Slice *toSlice ) {

  IDType toSeqRegionId = 0;

  if (toSlice != NULL) {
    if (!CoordSystem_compare(ChainedAssemblyMapper_getFirstCoordSystem(casmMapper), 
                             ChainedAssemblyMapper_getLastCoordSystem(casmMapper))) {
      return AssemblyMapperAdaptor_registerChainedSpecial(ama, casmMapper, from, seqRegionId, ranges, toSlice);
    }

    IDType toSeqRegionId = Slice_getSeqRegionId(toSlice);

    // NIY can we do not defined?? Use 0 for now if (!defined($to_seq_region_id)){
    if (!toSeqRegionId) {
      fprintf(stderr, "Could not get seq_region_id for to_slice %s\n", Slice_getSeqRegionName(toSlice));
      return ;
    }
  }

  char *          startName;
  Mapper *startMidMapper;
  CoordSystem    *startCs;
  RangeRegistry * startRegistry;
  char *          endName;
  Mapper *endMidMapper;
  CoordSystem *   endCs;
  RangeRegistry * endRegistry;

  if (!strcmp(from,"first")) {
    startName       = AMA_FIRST;
    startMidMapper  = ChainedAssemblyMapper_getFirstMiddleMapper(casmMapper);
    startCs         = ChainedAssemblyMapper_getFirstCoordSystem(casmMapper);
    startRegistry   = ChainedAssemblyMapper_getFirstRegistry(casmMapper);
    endMidMapper    = ChainedAssemblyMapper_getLastMiddleMapper(casmMapper);
    endCs           = ChainedAssemblyMapper_getLastCoordSystem(casmMapper);
    endRegistry     = ChainedAssemblyMapper_getLastRegistry(casmMapper);
    endName         = AMA_LAST;
  } else if (!strcmp(from, "last")) {
    startName       = AMA_LAST;
    startMidMapper  = ChainedAssemblyMapper_getLastMiddleMapper(casmMapper);
    startCs         = ChainedAssemblyMapper_getLastCoordSystem(casmMapper);
    startRegistry   = ChainedAssemblyMapper_getLastRegistry(casmMapper);
    endMidMapper    = ChainedAssemblyMapper_getFirstMiddleMapper(casmMapper);
    endCs           = ChainedAssemblyMapper_getFirstCoordSystem(casmMapper);
    endRegistry     = ChainedAssemblyMapper_getFirstRegistry(casmMapper);
    endName         = AMA_FIRST;
  } else {
    fprintf(stderr, "Invalid from argument: [%s], must be 'first' or 'last'",from);
    return;
  }

  Mapper *combinedMapper  = ChainedAssemblyMapper_getFirstLastMapper(casmMapper);
  CoordSystem *midCs      = ChainedAssemblyMapper_getMiddleCoordSystem(casmMapper);
  CoordSystemAdaptor *csa = DBAdaptor_getCoordSystemAdaptor(ama->dba);

  // SMJS Unused my $mid_name        = 'middle'

  // Check for the simple case where the ChainedMapper is short
  if( midCs == NULL) {
    startMidMapper = combinedMapper;
  }

  //#############
  // obtain the first half of the mappings and load them into the start mapper
  //

  //ascertain which is component and which is actually assembled coord system
  Vector *path;

  // check for the simple case, where the ChainedMapper is short
  if (midCs != NULL) {
    path = CoordSystemAdaptor_getMappingPath(csa, startCs, midCs);
  } else {
    path = CoordSystemAdaptor_getMappingPath(csa, startCs, endCs);
  }

  // NIY Free path - NO DON't FREE IT, Its in the cache

  if (Vector_getNumElement(path) != 2 && Vector_getElementAt(path,1) != NULL) {
    //char *pathStr = makeMappingPathKey(path);
    char pathStr[2048];
    makeMappingPathKey(path, pathStr);
    int len  = Vector_getNumElement(path)-1;

    fprintf(stderr, "Unexpected mapping path between start and intermediate coord systems (%s  %s and %s %s).\n"
                    "Expected path length 1, got %d. (path=%s).\n", 
                     CoordSystem_getName(startCs), CoordSystem_getVersion(startCs),
                     CoordSystem_getName(midCs), CoordSystem_getVersion(midCs),
                     len, pathStr);
    return;
  }

  StatementHandle *sth;
  CoordSystem *asmCs;
  CoordSystem *cmpCs;
  asmCs = Vector_getElementAt(path, 0);
  cmpCs = Vector_getLastElement(path);

  //the SQL varies depending on whether we are coming from assembled or
  //component coordinate system

  char asm2CmpBaseStr[1024];
  sprintf(asm2CmpBaseStr, "SELECT "
                          " asm.cmp_start,"
                          " asm.cmp_end,"
                          " asm.cmp_seq_region_id,"
                          " sr.name,"
                          " sr.length,"
                          " asm.ori,"
                          " asm.asm_start,"
                          " asm.asm_end"
                       " FROM"
                          " assembly asm, seq_region sr"
                      " WHERE"
                          " asm.asm_seq_region_id = %"IDFMTSTR" AND"
                          " %%ld <= asm.asm_end AND"
                          " %%ld >= asm.asm_start AND"
                          " asm.cmp_seq_region_id = sr.seq_region_id AND"
                          " sr.coord_system_id = %"IDFMTSTR);


  char cmp2AsmBaseStr[1024];
  sprintf(cmp2AsmBaseStr, "SELECT "
                          " asm.asm_start,"
                          " asm.asm_end,"
                          " asm.asm_seq_region_id,"
                          " sr.name,"
                          " sr.length,"
                          " asm.ori,"
                          " asm.cmp_start,"
                          " asm.cmp_end"
                       " FROM"
                          " assembly asm, seq_region sr"
                      " WHERE"
                          " asm.cmp_seq_region_id = %"IDFMTSTR" AND"
                          " %%ld <= asm.cmp_end AND"
                          " %%ld >= asm.cmp_start AND"
                          " asm.asm_seq_region_id = sr.seq_region_id AND"
                          " sr.coord_system_id = %"IDFMTSTR);

  char cmp2AsmStr[1024];
  char asm2CmpStr[1024];

  if (toSlice != NULL) {
    CoordSystem *toCs = Slice_getCoordSystem(toSlice);
    if (!CoordSystem_compare(asmCs, toCs)) {

      sprintf(asm2CmpStr, "%s", asm2CmpBaseStr);
      sprintf(cmp2AsmStr, "%s AND asm.asm_seq_region_id = "IDFMTSTR, cmp2AsmBaseStr, toSeqRegionId);
    } else if (!CoordSystem_compare(cmpCs, toCs)) {

      sprintf(asm2CmpStr, "%s AND asm.cmp_seq_region_id = "IDFMTSTR, asm2CmpBaseStr, toSeqRegionId);
      sprintf(cmp2AsmStr, "%s", cmp2AsmBaseStr);
    } else {

      sprintf(asm2CmpStr, "%s", asm2CmpBaseStr);
      sprintf(cmp2AsmStr, "%s", cmp2AsmBaseStr);
    }
  } else{

    sprintf(asm2CmpStr, "%s", asm2CmpBaseStr);
    sprintf(cmp2AsmStr, "%s", cmp2AsmBaseStr);
  }

  sth = (!CoordSystem_compare(asmCs, startCs)) ? ama->prepare((BaseAdaptor *)ama, asm2CmpStr, strlen(asm2CmpStr)) :
                                                 ama->prepare((BaseAdaptor *)ama, cmp2AsmStr, strlen(cmp2AsmStr));

  IDType midCsId;

  // check for the simple case where the ChainedMapper is short
  if ( midCs != NULL ) {
    midCsId = CoordSystem_getDbID(midCs);
  } else {
    midCsId = CoordSystem_getDbID(endCs);
  }

  Vector *midRanges   = Vector_new();
  Vector *startRanges = Vector_new();

  //need to perform the query for each unregistered range
  int i;
  for (i=0; i<Vector_getNumElement(ranges); i++) {
    CoordPair *range = Vector_getElementAt(ranges, i);

    long start = CoordPair_getStart(range);
    long end   = CoordPair_getEnd(range);

    sth->execute(sth, seqRegionId, start, end, midCsId);

    //load the start <-> mid mapper with the results and record the mid cs
    //ranges we just added to the mapper

    long midStart, midEnd;
    IDType midSeqRegionId;
    char *midSeqRegion;
    long midLength;
    int ori;
    long startStart, startEnd;

    ResultRow *row;
    while ((row = sth->fetchRow(sth))) {
      midStart       = row->getLongAt(row,0);
      midEnd         = row->getLongAt(row,1);
      midSeqRegionId = row->getLongLongAt(row,2);
      midSeqRegion   = row->getStringAt(row,3);
      midLength      = row->getLongAt(row,4);
      ori            = row->getIntAt(row,5);
      startStart     = row->getLongAt(row,6);
      startEnd       = row->getLongAt(row,7);

      if ( midCs ) {
        Mapper_addMapCoordinates(startMidMapper,
           seqRegionId, startStart, startEnd, ori,
           midSeqRegionId, midStart, midEnd);
      } else {
        if (!strcmp(from,"first")) {
          Mapper_addMapCoordinates(combinedMapper,
             seqRegionId, startStart, startEnd, ori,
             midSeqRegionId, midStart, midEnd);
        } else {
          Mapper_addMapCoordinates(combinedMapper,
             midSeqRegionId, midStart, midEnd, ori,
             seqRegionId, startStart, startEnd);
        }
      }

      //update sr_name cache
      DBAdaptor_addToSrCaches(ama->dba, midSeqRegionId, midSeqRegion,  midCsId, midLength);

      AssemblyMapperAdaptor_addToRangeVector(midRanges, midSeqRegionId, midStart, midEnd, midSeqRegion);  

      AssemblyMapperAdaptor_addToRangeVector(startRanges, seqRegionId, startStart, startEnd, NULL);  

      //the region that we actually register may actually be larger or smaller
      //than the region that we wanted to register.
      //register the intersection of the region so we do not end up doing
      //extra work later

      if (startStart < start || startEnd > end) {
        RangeRegistry_checkAndRegister(startRegistry, seqRegionId, startStart, startEnd, startStart, startEnd, 0);
      }
    }
  }
  sth->finish(sth);

  // in the one step case, we load the mid ranges in the
  // last_registry and we are done
  if (midCs == NULL) {
    int i;
    for (i=0;i<Vector_getNumElement(midRanges); i++) {
      SeqRegionRange *range = Vector_getElementAt(midRanges, i);

      RangeRegistry_checkAndRegister(endRegistry, 
                                     SeqRegionRange_getSeqRegionId(range), 
                                     SeqRegionRange_getSeqRegionStart(range), 
                                     SeqRegionRange_getSeqRegionEnd(range),
                                     SeqRegionRange_getSeqRegionStart(range), 
                                     SeqRegionRange_getSeqRegionEnd(range),
                                     0
                                    );
    }

    // and thats it for the simple case ...
    return;
  }


  //##########
  // now the second half of the mapping
  // perform another query and load the mid <-> end mapper using the mid cs
  // ranges
  // 

  //ascertain which is component and which is actually assembled coord system

  path = CoordSystemAdaptor_getMappingPath(csa, midCs, endCs);
  if (Vector_getNumElement(path) == 2 || ( Vector_getNumElement(path) == 3 && Vector_getElementAt(path, 1) == NULL)) {
    asmCs = Vector_getElementAt(path, 0);
    cmpCs = Vector_getLastElement(path);
  } else {
    //char *pathStr = makeMappingPathKey(path);
    char pathStr[2048];
    makeMappingPathKey(path, pathStr);
    int len  = Vector_getNumElement(path)-1;

    fprintf(stderr, "Unexpected mapping path between intermediate and last coord systems (%s  %s and %s %s).\n"
                    "Expected path length 1, got %d. (path=%s).\n", 
                     CoordSystem_getName(midCs), CoordSystem_getVersion(midCs),
                     CoordSystem_getName(endCs), CoordSystem_getVersion(endCs),
                     len, pathStr);
    return;
  }

  if (toSlice != NULL) {
    CoordSystem *toCs = Slice_getCoordSystem(toSlice);
    if (!CoordSystem_compare(asmCs, toCs)) {

      sprintf(asm2CmpStr, "%s", asm2CmpBaseStr);
      sprintf(cmp2AsmStr, "%s AND asm.asm_seq_region_id = "IDFMTSTR, cmp2AsmBaseStr, toSeqRegionId);
    } else if (!CoordSystem_compare(cmpCs, toCs)) {

      sprintf(asm2CmpStr, "%s AND asm.cmp_seq_region_id = "IDFMTSTR, asm2CmpBaseStr, toSeqRegionId);
      sprintf(cmp2AsmStr, "%s", cmp2AsmBaseStr);
    } else {

      sprintf(asm2CmpStr, "%s", asm2CmpBaseStr);
      sprintf(cmp2AsmStr, "%s", cmp2AsmBaseStr);
    }
  }

  sth = (!CoordSystem_compare(asmCs, midCs)) ? ama->prepare((BaseAdaptor *)ama, asm2CmpStr, strlen(asm2CmpStr)) :
                                               ama->prepare((BaseAdaptor *)ama, cmp2AsmStr, strlen(cmp2AsmStr));

  IDType endCsId = CoordSystem_getDbID(endCs);

  for (i=0; i<Vector_getNumElement(midRanges); i++) {
    SeqRegionRange *midRange = Vector_getElementAt(midRanges, i);

    IDType midSeqRegionId = SeqRegionRange_getSeqRegionId(midRange);
    long start            = SeqRegionRange_getSeqRegionStart(midRange);
    long end              = SeqRegionRange_getSeqRegionEnd(midRange);

    sth->execute(sth, midSeqRegionId, start, end, endCsId);

    //load the end <-> mid mapper with the results and record the mid cs
    //ranges we just added to the mapper

    long endStart, endEnd;
    IDType endSeqRegionId;
    char *endSeqRegion;
    long endLength;
    int ori;
    long midStart, midEnd;

    ResultRow *row;
    while ((row = sth->fetchRow(sth))) {
      endStart       = row->getLongAt(row,0);
      endEnd         = row->getLongAt(row,1);
      endSeqRegionId = row->getLongLongAt(row,2);
      endSeqRegion   = row->getStringAt(row,3);
      endLength      = row->getLongAt(row,4);
      ori            = row->getIntAt(row,5);
      midStart       = row->getLongAt(row,6);
      midEnd         = row->getLongAt(row,7);

      Mapper_addMapCoordinates(endMidMapper,
         endSeqRegionId, endStart, endEnd, ori,
         midSeqRegionId, midStart, midEnd);

      //update sr_name cache
      DBAdaptor_addToSrCaches(ama->dba, endSeqRegionId, endSeqRegion, endCsId, endLength);

      //register this region on the end coord system
      RangeRegistry_checkAndRegister(endRegistry, endSeqRegionId, endStart, endEnd, endStart, endEnd, 0);
    }
  }
  sth->finish(sth);

  //########
  // Now that both halves are loaded
  // Do stepwise mapping using both of the loaded mappers to load
  // the final start <-> end mapper
  //

  AssemblyMapperAdaptor_buildCombinedMapper(ama, startRanges, startMidMapper, endMidMapper,
                                            combinedMapper, startName);
  //all done!
  return;
}


/*
=head2 _register_chained_special

  Arg [1]    : Bio::EnsEMBL::ChainedAssemblyMapper $casm_mapper
               The chained assembly mapper to register regions on
  Arg [2]    : string $from ('first' or 'last')
               The direction we are registering from, and the name of the
               internal mapper.
  Arg [3]    : string $seq_region_name
               The name of the seqregion we are registering on
  Arg [4]    : listref $ranges
               A list  of ranges to register (in [$start,$end] tuples).
  Arg [5]    : (optional) $to_slice
               Only register those on this Slice.
  Description: Registers a set of ranges on a chained assembly mapper.
               This function is at the heart of the chained mapping process.
               It retrieves information from the assembly table and
               dynamically constructs the mappings between two coordinate
               systems which are 2 mapping steps apart. It does this by using
               two internal mappers to load up a third mapper which is
               actually used by the ChainedAssemblyMapper to perform the
               mapping.

               This method must be called before any mapping is
               attempted on regions of interest, otherwise only gaps will
               be returned.  Note that the ChainedAssemblyMapper automatically
               calls this method when the need arises.
  Returntype : none
  Exceptions : throw if the seq_region to be registered does not exist
               or if it associated with multiple assembled pieces (bad data
               in assembly table)

               throw if the mapping between the coordinate systems cannot
               be performed in two steps, which means there is an internal
               error in the data in the meta table or in the code that creates
               the mapping paths.
  Caller     : Bio::EnsEMBL::AssemblyMapper
  Status     : Stable

=cut
*/

void AssemblyMapperAdaptor_registerChainedSpecial(AssemblyMapperAdaptor *ama, ChainedAssemblyMapper *casmMapper, 
                                                  char *from, IDType seqRegionId, Vector *ranges, Slice *toSlice) {
  int found = 0;

  char qStr[1024];
  sprintf(qStr, "SELECT "
                  " asm.cmp_start,"
                  " asm.cmp_end,"
                  " asm.cmp_seq_region_id,"
                  " sr.name,"
                  " sr.length,"
                  " asm.ori,"
                  " asm.asm_start,"
                  " asm.asm_end"
                 " FROM"
                  " assembly asm, seq_region sr"
                " WHERE"
                  " asm.asm_seq_region_id = %"IDFMTSTR" AND"
                  " %%ld <= asm.asm_end AND"
                  " %%ld >= asm.asm_start AND"
                  " asm.cmp_seq_region_id = sr.seq_region_id AND"
                  " sr.coord_system_id = %"IDFMTSTR" AND"
                  " asm.cmp_seq_region_id = %"IDFMTSTR);

  StatementHandle *sth = ama->prepare((BaseAdaptor *)ama,qStr,strlen(qStr));

  char *          startName;
  Mapper *startMidMapper;
  CoordSystem *   startCs;
  RangeRegistry * startRegistry;
  char *          endName;
  Mapper *endMidMapper;
  CoordSystem    *endCs;
  RangeRegistry * endRegistry;

  if (!strcmp(from,"first")) {
    startName       = AMA_FIRST;
    startMidMapper  = ChainedAssemblyMapper_getFirstMiddleMapper(casmMapper);
    startCs         = ChainedAssemblyMapper_getFirstCoordSystem(casmMapper);
    startRegistry   = ChainedAssemblyMapper_getFirstRegistry(casmMapper);
    endMidMapper    = ChainedAssemblyMapper_getLastMiddleMapper(casmMapper);
    endCs           = ChainedAssemblyMapper_getLastCoordSystem(casmMapper);
    endRegistry     = ChainedAssemblyMapper_getLastRegistry(casmMapper);
    endName         = AMA_LAST;
  } else if (!strcmp(from, "last")) {
    startName       = AMA_LAST;
    startMidMapper  = ChainedAssemblyMapper_getLastMiddleMapper(casmMapper);
    startCs         = ChainedAssemblyMapper_getLastCoordSystem(casmMapper);
    startRegistry   = ChainedAssemblyMapper_getLastRegistry(casmMapper);
    endMidMapper    = ChainedAssemblyMapper_getFirstMiddleMapper(casmMapper);
    endCs           = ChainedAssemblyMapper_getFirstCoordSystem(casmMapper);
    endRegistry     = ChainedAssemblyMapper_getFirstRegistry(casmMapper);
    endName         = AMA_FIRST;
  } else {
    fprintf(stderr, "Invalid from argument: [%s], must be 'first' or 'last'",from);
    return;
  }

  Mapper *combinedMapper  = ChainedAssemblyMapper_getFirstLastMapper(casmMapper);
  CoordSystem *midCs      = ChainedAssemblyMapper_getMiddleCoordSystem(casmMapper);
  CoordSystemAdaptor *csa = DBAdaptor_getCoordSystemAdaptor(ama->dba);
  // SMJS Unused my $mid_name        = 'middle';

  // Check for the simple case where the ChainedMapper is short
  if (midCs == NULL) {
    startMidMapper = combinedMapper;
  }

  Vector *path;
  if ( midCs ) {
    path = CoordSystemAdaptor_getMappingPath(csa, startCs, midCs);
  } else {
    path = CoordSystemAdaptor_getMappingPath(csa, startCs, endCs);
  }

  if (midCs == NULL) {
    startMidMapper = combinedMapper;
  }

  if (Vector_getNumElement(path) != 2 && Vector_getElementAt(path,1) != NULL) {
    //char *pathStr = makeMappingPathKey(path);
    char pathStr[2048];
    makeMappingPathKey(path, pathStr);
    int len  = Vector_getNumElement(path)-1;

    fprintf(stderr, "Unexpected mapping path between start and intermediate coord systems (%s  %s and %s %s).\n"
                    "Expected path length 1, got %d. (path=%s).\n", 
                     CoordSystem_getName(startCs), CoordSystem_getVersion(startCs),
                     CoordSystem_getName(midCs), CoordSystem_getVersion(midCs),
                     len, pathStr);
    return;
  }

  CoordSystem *asmCs;
  CoordSystem *cmpCs;
  asmCs = Vector_getElementAt(path, 0);
  cmpCs = Vector_getLastElement(path);

  //Mapper *combinedMapper  = ChainedAssemblyMapper_getFirstLastMapper(casmMapper);
  //CoordSystem *midCs      = ChainedAssemblyMapper_getMiddleCoordSystem(casmMapper);
  //CoordSystemAdaptor *csa = DBAdaptor_getCoordSystemAdaptor(ama->dba);
  // SMJS Unused$mid_name        = 'middle';

  IDType midCsId;

  // Check for the simple case where the ChainedMapper is short
  if (midCs == NULL) {
    startMidMapper = combinedMapper;
  } else {
    midCsId = CoordSystem_getDbID(midCs);
  }

  Vector *midRanges   = Vector_new();
  Vector *startRanges = Vector_new();

  CoordSystem *toCs = Slice_getCoordSystem(toSlice);

  int direction;
  for (direction = 1; direction >= 0; direction--) {
    IDType id1;
    IDType id2;
    if (direction) {
      id1 = seqRegionId;
      id2 = Slice_getSeqRegionId(toSlice);
    } else {
      id1 = Slice_getSeqRegionId(toSlice);
      id2 = seqRegionId;
    }

    int i;
    for (i=0; i<Vector_getNumElement(ranges); i++) {
      CoordPair *range = Vector_getElementAt(ranges, i);

      long start = CoordPair_getStart(range);
      long end   = CoordPair_getEnd(range);

      sth->execute(sth, id1, start, end, CoordSystem_getDbID(toCs), id2);

      long midStart, midEnd;
      IDType midSeqRegionId;
      char *midSeqRegion;
      long midLength;
      int ori;
      long startStart, startEnd;

      ResultRow *row;
      while ((row = sth->fetchRow(sth))) {
        midStart       = row->getLongAt(row,0);
        midEnd         = row->getLongAt(row,1);
        midSeqRegionId = row->getLongLongAt(row,2);
        midSeqRegion   = row->getStringAt(row,3);
        midLength      = row->getLongAt(row,4);
        ori            = row->getIntAt(row,5);
        startStart     = row->getLongAt(row,6);
        startEnd       = row->getLongAt(row,7);

        found = 1;

        if (midCs) {
          Mapper_addMapCoordinates(startMidMapper,
             id1, startStart, startEnd, ori,
             midSeqRegionId, midStart, midEnd);
        } else {
          if (!strcmp(from,"first")) {
            if (direction) {
              Mapper_addMapCoordinates(combinedMapper,
                 id1, startStart, startEnd, ori,
                 midSeqRegionId, midStart, midEnd);
            } else {
              Mapper_addMapCoordinates(combinedMapper,
                 midSeqRegionId, midStart, midEnd, ori,
                 id1, startStart, startEnd);
            }
          } else {
            if (direction) {
              Mapper_addMapCoordinates(combinedMapper,
                 midSeqRegionId, midStart, midEnd, ori,
                 id1, startStart, startEnd);
            } else {
              Mapper_addMapCoordinates(combinedMapper,
                 id1, startStart, startEnd, ori,
                 midSeqRegionId, midStart, midEnd);
            }
          }
        }
        
        //update sr_name cache
        DBAdaptor_addToSrCaches(ama->dba, midSeqRegionId, midSeqRegion, midCsId, midLength);
        
        AssemblyMapperAdaptor_addToRangeVector(midRanges, midSeqRegionId, midStart, midEnd, midSeqRegion);  

        AssemblyMapperAdaptor_addToRangeVector(startRanges, id1, startStart, startEnd, NULL);  
        
        //the region that we actually register may actually be larger or smaller
        //than the region that we wanted to register.
        //register the intersection of the region so we do not end up doing
        //extra work later
        
        if (startStart < start || startEnd > end) {
          RangeRegistry_checkAndRegister(startRegistry, id1, startStart, startEnd, startStart, startEnd, 0);
        }
      }
    }
    sth->finish(sth);

    if (found) {
      if (midCs == NULL ) {
        int j;
        for (j=0;j<Vector_getNumElement(midRanges); j++) {
          SeqRegionRange *range = Vector_getElementAt(midRanges, j);

          RangeRegistry_checkAndRegister(endRegistry, 
                                         SeqRegionRange_getSeqRegionId(range),
                                         SeqRegionRange_getSeqRegionStart(range),
                                         SeqRegionRange_getSeqRegionEnd(range),
                                         SeqRegionRange_getSeqRegionStart(range),
                                         SeqRegionRange_getSeqRegionEnd(range),
                                         0
                                        );
        }
        
        // and thats it for the simple case ...
        Vector_free(midRanges);
        Vector_free(startRanges);
        return;
      }
    }
  }

  Vector_free(midRanges);
  Vector_free(startRanges);
}


/*
=head2 register_all

  Arg [1]    : Bio::EnsEMBL::AssemblyMapper $mapper
  Example    : $mapper = $asm_mapper_adaptor->fetch_by_CoordSystems($cs1,$cs2);

               # make cache large enough to hold all of the mappings
               $mapper->max_pair_count(10e6);
               $asm_mapper_adaptor->register_all($mapper);

               # perform mappings as normal
               $mapper->map($slice->seq_region_name(), $sr_start, $sr_end,
                            $sr_strand, $cs1);
               ...
  Description: This function registers the entire set of mappings between
               two coordinate systems in an assembly mapper.
               This will use a lot of memory but will be much more efficient
               when doing a lot of mapping which is spread over the entire
               genome.
  Returntype : none
  Exceptions : none
  Caller     : specialised prograhsm
  Status     : Stable

=cut
*/

// NIY Not sure if this could be a chained mapper???
void AssemblyMapperAdaptor_registerAll(AssemblyMapperAdaptor *ama, AssemblyMapper *asmMapper) {

  IDType cmpCsId = CoordSystem_getDbID(AssemblyMapper_getComponentCoordSystem(asmMapper));
  IDType asmCsId = CoordSystem_getDbID(AssemblyMapper_getAssembledCoordSystem(asmMapper));

  // retrieve every relevant assembled/component pair from the assembly table

  char qStr[1024];
  sprintf(qStr, "SELECT "
                  " asm.cmp_start,"
                  " asm.cmp_end,"
                  " asm.cmp_seq_region_id,"
                  " cmp_sr.name,"
                  " cmp_sr.length,"
                  " asm.ori,"
                  " asm.asm_start,"
                  " asm.asm_end,"
                  " asm.asm_seq_region_id,"
                  " asm_sr.name,"
                  " asm_sr.length"
              " FROM"
                  " assembly asm, seq_region asm_sr, seq_region cmp_sr"
              " WHERE"
                  " asm.cmp_seq_region_id = cmp_sr.seq_region_id AND"
                  " asm.asm_seq_region_id = asm_sr.seq_region_id AND"
                  " cmp_sr.coord_system_id = "IDFMTSTR" AND"
                  " asm_sr.coord_system_id = "IDFMTSTR, cmpCsId, asmCsId);

  StatementHandle *sth = ama->prepare((BaseAdaptor *)ama,qStr,strlen(qStr));

  sth->execute(sth);

  // load the asmMapper with the assembly information
  IDHash *asmRegistered = IDHash_new(IDHASH_MEDIUM);

  ResultRow *row;
  while ((row = sth->fetchRow(sth))) {
    long cmpStart         = row->getLongAt(row,0);
    long cmpEnd           = row->getLongAt(row,1);
    IDType cmpSeqRegionId = row->getLongLongAt(row,2);
    char *cmpSeqRegion    = row->getStringAt(row,3);
    long cmpLength        = row->getLongAt(row,4);
    int  ori              = row->getIntAt(row,5);
    long asmStart         = row->getLongAt(row,6);
    long asmEnd           = row->getLongAt(row,7);
    IDType asmSeqRegionId = row->getLongLongAt(row,8);
    char *asmSeqRegion    = row->getStringAt(row,9);
    long asmLength        = row->getLongAt(row,10);

    AssemblyMapper_registerComponent(asmMapper, cmpSeqRegionId);

    Mapper *mapper = AssemblyMapper_getMapper(asmMapper);

    Mapper_addMapCoordinates(mapper,
                 asmSeqRegionId, asmStart, asmEnd, ori,
                 cmpSeqRegionId, cmpStart, cmpEnd);

    DBAdaptor_addToSrCaches(ama->dba, cmpSeqRegionId, cmpSeqRegion, cmpCsId, cmpLength);

    // only register each asm seq_region once since it requires some work
    if ( ! IDHash_contains(asmRegistered, asmSeqRegionId)) {
      IDHash_add(asmRegistered, asmSeqRegionId, &trueVal);

      // register all chunks from start of seq region to end
      int endChunk = asmLength >> CHUNKFACTOR;
      int i;
      for (i=0; i<=endChunk; i++) {
        AssemblyMapper_registerAssembled(asmMapper, asmSeqRegionId, i);
      }

      DBAdaptor_addToSrCaches(ama->dba, asmSeqRegionId, asmSeqRegion, asmCsId, asmLength);
    }
  }

  IDHash_free(asmRegistered, NULL);

  sth->finish(sth);

  return;
}


/*

=head2 register_all_chained

  Arg [1]    : Bio::EnsEMBL::ChainedAssemblyMapper $casm_mapper
  Example    : $mapper = $asm_mapper_adaptor->fetch_by_CoordSystems($cs1,$cs2);

               # make the cache large enough to hold all of the mappings
               $mapper->max_pair_count(10e6);
               # load all of the mapping data
               $asm_mapper_adaptor->register_all_chained($mapper);

               # perform mappings as normal
               $mapper->map($slice->seq_region_name(), $sr_start, $sr_end,
                            $sr_strand, $cs1);
               ...
  Description: This function registers the entire set of mappings between
               two coordinate systems in a chained mapper.  This will use a lot
               of memory but will be much more efficient when doing a lot of
               mapping which is spread over the entire genome.
  Returntype : none
  Exceptions : throw if mapper is between coord systems with unexpected
               mapping paths
  Caller     : specialised programs doing a lot of genome-wide mapping
  Status     : Stable

=cut
*/

void AssemblyMapperAdaptor_registerAllChained(AssemblyMapperAdaptor *ama, ChainedAssemblyMapper *casmMapper) {

  CoordSystem *firstCs = ChainedAssemblyMapper_getFirstCoordSystem(casmMapper);
  CoordSystem *midCs   = ChainedAssemblyMapper_getMiddleCoordSystem(casmMapper);
  CoordSystem *lastCs  = ChainedAssemblyMapper_getLastCoordSystem(casmMapper);

  Mapper *startMidMapper = ChainedAssemblyMapper_getFirstMiddleMapper(casmMapper);
  Mapper *endMidMapper   = ChainedAssemblyMapper_getLastMiddleMapper(casmMapper);
  Mapper *combinedMapper = ChainedAssemblyMapper_getFirstLastMapper(casmMapper);

  Vector *ranges = NULL;

  char qStr[1024];
  sprintf(qStr,"SELECT "
                  " asm.cmp_start,"
                  " asm.cmp_end,"
                  " asm.cmp_seq_region_id,"
                  " sr_cmp.name,"
                  " sr_cmp.length,"
                  " asm.ori,"
                  " asm.asm_start,"
                  " asm.asm_end,"
                  " asm.asm_seq_region_id,"
                  " sr_asm.name,"
                  " sr_asm.length"
               " FROM"
                  " assembly asm, seq_region sr_asm, seq_region sr_cmp"
               " WHERE"
                  " sr_asm.seq_region_id = asm.asm_seq_region_id AND"
                  " sr_cmp.seq_region_id = asm.cmp_seq_region_id AND"
                  " sr_asm.coord_system_id = %"IDFMTSTR" AND"
                  " sr_cmp.coord_system_id = %"IDFMTSTR); 


  CoordSystemAdaptor *csa = DBAdaptor_getCoordSystemAdaptor(ama->dba);

  Vector *path;
  if ( midCs == NULL ) {
    path = CoordSystemAdaptor_getMappingPath(csa, firstCs, lastCs);
    if (Vector_getElementAt(path,1) == NULL) {
      Vector_removeElementAt(path, 1);
    }
  } else {
    path = CoordSystemAdaptor_getMappingPath(csa, firstCs, midCs);
    // fix for when we have something like supercontig#contig#chromosome
    if (Vector_getElementAt(path,1) == NULL) {
      Vector_removeElementAt(path, 1);
    }
  }

  if (Vector_getNumElement(path) != 2 && Vector_getElementAt(path,1) != NULL) {
    //char *pathStr = makeMappingPathKey(path);
    char pathStr[2048];
    makeMappingPathKey(path, pathStr);
    int len  = Vector_getNumElement(path)-1;

    fprintf(stderr, "Unexpected mapping path between start and intermediate coord systems (%s  %s and %s %s).\n"
                    "Expected path length 1, got %d. (path=%s).\n", 
                     CoordSystem_getName(firstCs), CoordSystem_getVersion(firstCs),
                     CoordSystem_getName(midCs), CoordSystem_getVersion(midCs),
                     len, pathStr);
    return;
  }

  CoordSystem *asmCs = Vector_getElementAt(path,0);
  CoordSystem *cmpCs = Vector_getElementAt(path,1);

  StatementHandle *sth = ama->prepare((BaseAdaptor *)ama,qStr,strlen(qStr));

  sth->addFlag(sth, MYSQLFLAG_USE_RESULT);

  sth->execute(sth, CoordSystem_getDbID(asmCs), CoordSystem_getDbID(cmpCs));

  long  midStart, midEnd;
  IDType midSeqRegionId;
  char *midSeqRegion;
  long midLength;

  int  ori;

  long startStart, startEnd;
  IDType startSeqRegionId;
  char *startSeqRegion;
  long startLength;


  IDType midCsId;
  IDType startCsId;
  Mapper *mapper;
  RangeRegistry *reg;

  if ( midCs == NULL ) {
    midCsId   = CoordSystem_getDbID(lastCs);
    startCsId = CoordSystem_getDbID(firstCs);
    mapper    = combinedMapper;
  } else {
    midCsId   = CoordSystem_getDbID(midCs);
    startCsId = CoordSystem_getDbID(firstCs);
    mapper    = startMidMapper;
  }

  reg = ChainedAssemblyMapper_getFirstRegistry(casmMapper);

  ResultRow *row;
  while ((row = sth->fetchRow(sth))) {
    if (!CoordSystem_compare(asmCs,firstCs)) {
      midStart           = row->getLongAt(row,0);
      midEnd             = row->getLongAt(row,1);
      midSeqRegionId   = row->getLongLongAt(row,2);
      midSeqRegion      = row->getStringAt(row,3);
      midLength          = row->getLongAt(row,4);
      ori                = row->getIntAt(row,5);
      startStart         = row->getLongAt(row,6);
      startEnd           = row->getLongAt(row,7);
      startSeqRegionId = row->getLongLongAt(row,8);
      startSeqRegion    = row->getStringAt(row,9);
      startLength        = row->getLongAt(row,10);

    } else {
      startStart         = row->getLongAt(row,0);
      startEnd           = row->getLongAt(row,1);
      startSeqRegionId = row->getLongLongAt(row,2);
      startSeqRegion    = row->getStringAt(row,3);
      startLength        = row->getLongAt(row,4);
      ori                = row->getIntAt(row,5);
      midStart           = row->getLongAt(row,6);
      midEnd             = row->getLongAt(row,7);
      midSeqRegionId   = row->getLongLongAt(row,8);
      midSeqRegion      = row->getStringAt(row,9);
      midLength          = row->getLongAt(row,10);
    }

    Mapper_addMapCoordinates(mapper,
       startSeqRegionId, startStart, startEnd, ori,
       midSeqRegionId, midStart, midEnd);

    /* gb10: ranges doesn't seem to get set so surely this can't work?! */
    if (ranges)
      AssemblyMapperAdaptor_addToRangeVector(ranges, startSeqRegionId, startStart, startEnd, NULL);  

    RangeRegistry_checkAndRegister(reg, startSeqRegionId, 1, startLength, 1, startLength, 0 );

    if ( midCs == NULL ) {
      RangeRegistry *lastReg = ChainedAssemblyMapper_getLastRegistry(casmMapper);

      RangeRegistry_checkAndRegister(lastReg, midSeqRegionId, midStart, midEnd, midStart, midEnd, 0);
    }

    DBAdaptor_addToSrCaches(ama->dba, midSeqRegionId, midSeqRegion, midCsId, midLength);

    DBAdaptor_addToSrCaches(ama->dba, startSeqRegionId, startSeqRegion, startCsId, startLength);
  }

  if (midCs == NULL) {
    // thats it for the simple case
    return;
  }


  path = CoordSystemAdaptor_getMappingPath(csa, lastCs, midCs);
  if ( midCs ) {
    if (Vector_getElementAt(path,1) == NULL) {
      Vector_removeElementAt(path, 1);
    }
  }

  if (Vector_getNumElement(path) != 2 && Vector_getElementAt(path,1) != NULL) {
    //char *pathStr = makeMappingPathKey(path);
    char pathStr[2048];
    makeMappingPathKey(path, pathStr);
    int len  = Vector_getNumElement(path)-1;

    fprintf(stderr, "Unexpected mapping path between intermediate and last coord systems (%s  %s and %s %s).\n"
                    "Expected path length 1, got %d. (path=%s).\n", 
                     CoordSystem_getName(midCs), CoordSystem_getVersion(midCs),
                     CoordSystem_getName(lastCs), CoordSystem_getVersion(lastCs),
                     len, pathStr);
    return;
  }

  asmCs = Vector_getElementAt(path,0);
  cmpCs = Vector_getElementAt(path,1);

  sth->execute(sth, CoordSystem_getDbID(asmCs), CoordSystem_getDbID(cmpCs));

  long endStart, endEnd;
  IDType endSeqRegionId;
  char *endSeqRegion;
  long endLength;


  IDType endCsId = CoordSystem_getDbID(lastCs);
  reg = ChainedAssemblyMapper_getLastRegistry(casmMapper);

  while ((row = sth->fetchRow(sth))) {
    if (!CoordSystem_compare(asmCs,midCs)) {
      endStart         = row->getLongAt(row,0);
      endEnd           = row->getLongAt(row,1);
      endSeqRegionId = row->getLongLongAt(row,2);
      endSeqRegion    = row->getStringAt(row,3);
      endLength        = row->getLongAt(row,4);
      ori              = row->getIntAt(row,5);
      midStart         = row->getLongAt(row,6);
      midEnd           = row->getLongAt(row,7);
      midSeqRegionId = row->getLongLongAt(row,8);
      midSeqRegion    = row->getStringAt(row,9);
      midLength        = row->getLongAt(row,10);

    } else {
      midStart         = row->getLongAt(row,0);
      midEnd           = row->getLongAt(row,1);
      midSeqRegionId = row->getLongLongAt(row,2);
      midSeqRegion    = row->getStringAt(row,3);
      midLength        = row->getLongAt(row,4);
      ori              = row->getIntAt(row,5);
      endStart         = row->getLongAt(row,6);
      endEnd           = row->getLongAt(row,7);
      endSeqRegionId = row->getLongLongAt(row,8);
      endSeqRegion    = row->getStringAt(row,9);
      endLength        = row->getLongAt(row,10);
    }

    Mapper_addMapCoordinates(endMidMapper, 
       endSeqRegionId, endStart, endEnd, ori,
       midSeqRegionId, midStart, midEnd);

    RangeRegistry_checkAndRegister(reg, endSeqRegionId, 1, endLength, 1, endLength, 0 );

    DBAdaptor_addToSrCaches(ama->dba, endSeqRegionId, endSeqRegion, endCsId, endLength);
  }

  sth->finish(sth);

  AssemblyMapperAdaptor_buildCombinedMapper(ama,  ranges, startMidMapper, endMidMapper,
                                            combinedMapper, "first" );

  ranges->freeElement = SeqRegionRange_free;
  Vector_free(ranges);
  return;
}



// after both halves of a chained mapper are loaded
// this function maps all ranges in $ranges and loads the
// results into the combined mapper
void AssemblyMapperAdaptor_buildCombinedMapper(AssemblyMapperAdaptor *ama, Vector *ranges, Mapper *startMidMapper, 
                                              Mapper *endMidMapper, Mapper *combinedMapper, char *startName) {
  int ok = 1;
  int i;
  for (i=0; i<Vector_getNumElement(ranges); i++) {
    SeqRegionRange *range = Vector_getElementAt(ranges, i);
    
    IDType seqRegionId = SeqRegionRange_getSeqRegionId(range);
    long   start       = SeqRegionRange_getSeqRegionStart(range);
    long   end         = SeqRegionRange_getSeqRegionEnd(range);

    long sum = 0;

    MapperRangeSet *initialCoords = Mapper_mapCoordinates(startMidMapper, seqRegionId, start, end, 1, startName);

    int j;
    for (j=0; j < MapperRangeSet_getNumRange(initialCoords); j++) {
      MapperRange *icoord = MapperRangeSet_getRangeAt(initialCoords, j);

      //skip gaps
      if (icoord->rangeType == MAPPERRANGE_GAP) {
        sum += MapperRange_getLength(icoord);
        continue;
      }

      if (icoord->rangeType != MAPPERRANGE_COORD) {
        fprintf(stderr,"Belt and Braces check for icoord being a MAPPERRANGE_COORD. Its a %d\n", icoord->rangeType);
        ok = 0;
        break;
      }

      MapperCoordinate *coordICoord = (MapperCoordinate *)icoord;
      //feed the results of the first mapping into the second mapper
      MapperRangeSet *finalCoords =
        Mapper_mapCoordinates(endMidMapper, coordICoord->id, coordICoord->start, coordICoord->end, coordICoord->strand, "middle");


      int k;
      for (k=0; k < MapperRangeSet_getNumRange(finalCoords); k++) {
        MapperRange *fcoord = MapperRangeSet_getRangeAt(finalCoords, k);
        //load up the final mapper

        // Mapper::Coordinate - inheritance means INDEL is also a coordinate
        if (fcoord->rangeType == MAPPERRANGE_COORD || fcoord->rangeType == MAPPERRANGE_INDEL) {
// NIY: Need to check that INDEL is compatible with MapperCoordinate
          MapperCoordinate *coordFCoord = (MapperCoordinate *)fcoord;
          long totalStart = start + sum;
          long totalEnd   = totalStart + MapperRange_getLength(coordFCoord) - 1;
          int  ori = coordFCoord->strand;

          if (!strcmp(startName,"first")) { // add coords in consistent order
            Mapper_addMapCoordinates(combinedMapper,
                             seqRegionId, totalStart, totalEnd, ori,
                             coordFCoord->id, coordFCoord->start, coordFCoord->end);
          } else {
            Mapper_addMapCoordinates(combinedMapper,
                        coordFCoord->id, coordFCoord->start, coordFCoord->end, ori,
                        seqRegionId, totalStart, totalEnd);
          }
        }
        sum += MapperRange_getLength(fcoord);
      }
      MapperRangeSet_free(finalCoords);
    }

    if (ok) {
      MapperRangeSet_free(initialCoords);
    }
  }
  //all done!
}


/*
=head2 seq_regions_to_ids

  Arg [1]    : Bio::EnsEMBL::CoordSystem $coord_system
  Arg [2]    : listref of strings $seq_regions
  Example    : my @ids = @{$asma->seq_regions_to_ids($coord_sys, \@seq_regs)};
  Description: Converts a list of seq_region names to internal identifiers
               using the internal cache that has accumulated while registering
               regions for AssemblyMappers. If any requested regions are
               not  found in the cache an attempt is made to retrieve them
               from the database.
  Returntype : listref of ints
  Exceptions : throw if a non-existant seqregion is provided
  Caller     : general
  Status     : Stable

=cut
*/

Vector *AssemblyMapperAdaptor_seqRegionsToIds(AssemblyMapperAdaptor *ama, CoordSystem *coordSystem, Vector *seqRegions) {

  int ok = 1;
  IDType csId = CoordSystem_getDbID(coordSystem);

  Vector *out = Vector_new();

  int i;
  for (i=0; i< Vector_getNumElement(seqRegions); i++) {
    char *sr = Vector_getElementAt(seqRegions,i);
    char key[1024];

    sprintf(key,"%s:"IDFMTSTR, sr, csId);
    
    if (StringHash_contains(ama->srNameCache, key)) {
      SeqRegionCacheEntry *cacheData = StringHash_getValue(ama->srNameCache, key);
// NIY: Do I need to allocate these?
      Vector_addElement(out, &(cacheData->regionId));
    } else {
      IDType id = AssemblyMapperAdaptor_seqRegionNameToId(ama, sr, csId);
      IDType *idP;

      if ((idP = (IDType *)calloc(1,sizeof(IDType))) == NULL) {
        fprintf(stderr, "ERROR: Failed allocating space for idP\n");
        ok = 0;
        break;
      }
      *idP = id;
      
      Vector_addElement(out, idP);
    }
  }

  if (!ok) {
    Vector_free(out);
    out = NULL;
  }

  return out;
}

IDType AssemblyMapperAdaptor_seqRegionToId(AssemblyMapperAdaptor *ama, CoordSystem *coordSystem, char *seqRegionName) {
  IDType csId = CoordSystem_getDbID(coordSystem);
  char key[1024];

  IDType seqRegionId;

  sprintf(key,"%s:"IDFMTSTR, seqRegionName, csId);
    
  if (StringHash_contains(ama->srNameCache, key)) {
    SeqRegionCacheEntry *cacheData = StringHash_getValue(ama->srNameCache, key);
    seqRegionId = cacheData->regionId;
  } else {
    seqRegionId = AssemblyMapperAdaptor_seqRegionNameToId(ama, seqRegionName, csId);
  }

  return seqRegionId;
}

/*
=head2 seq_ids_to_regions

  Arg [1]    : listref of   seq_region ids
  Example    : my @ids = @{$asma->ids_to_seq_regions(\@seq_ids)};
  Description: Converts a list of seq_region ids to seq region names
               using the internal cache that has accumulated while registering
               regions for AssemblyMappers. If any requested regions are
               not  found in the cache an attempt is made to retrieve them
               from the database.
  Returntype : listref of strings
  Exceptions : throw if a non-existant seq_region_id is provided
  Caller     : general
  Status     : Stable

=cut
*/

Vector *AssemblyMapperAdaptor_seqIdsToRegions(AssemblyMapperAdaptor *ama, Vector *seqRegionIds) {
  Vector *out = Vector_new();

  int i;
  for (i=0; i< Vector_getNumElement(seqRegionIds); i++) {
    IDType *sr = Vector_getElementAt(seqRegionIds,i);

    if (IDHash_contains(ama->srIdCache, *sr)) {
      SeqRegionCacheEntry *cacheData = IDHash_getValue(ama->srIdCache, *sr);
// NIY: Do I need to allocate these?
      Vector_addElement(out, cacheData->regionName);
    } else {
      char *name = AssemblyMapperAdaptor_seqRegionIdToName(ama, *sr);
      Vector_addElement(out, name);
    }
  }

  return out;
}

/*
=head2 delete_cache

 Description: Delete all the caches for the mappings/seq_regions
 Returntype : none
 Exceptions : none
 Caller     : General
 Status     : At risk

=cut
*/


void AssemblyMapperAdaptor_deleteCache(AssemblyMapperAdaptor *ama) {

  fprintf(stderr, "deleteCache is not right yet\n");
  
  StringHash_free(ama->srNameCache, SeqRegionCacheEntry_free);
  IDHash_free(ama->srIdCache, SeqRegionCacheEntry_free);

  Vector *asmMappers = StringHash_getValues(ama->asmMapperCache);
  int i;
  for (i=0; i<Vector_getNumElement(asmMappers); i++) {
    AssemblyMapper *asmMapper = Vector_getElementAt(asmMappers, i);
    // NIY: Need to know what type of assembly mapper it is to call the correct flush (or 'objectize' them)

    AssemblyMapper_flush(asmMapper);

    // Do we also want to free these mappers - not just flush them??
  }

  StringHash_free(ama->asmMapperCache, NULL);

  Vector_free(asmMappers);

  return;
}


/*
=head2 register_region

  Description: DEPRECATED use register_assembled instead

=cut
*/

/* Comment out for now
void AssemblyMapperAdaptor_registerRegion(AssemblyMapperAdaptor *ama, AssemblyMapper *asmMapper, char *type, 
                                          char *chrName, long start, long end){

  fprintf(stderr,"Deprecated: Use register_assembled instead\n");

  AssemblyMapperAdaptor_registerAssembled(ama, asmMapper, chrName, start, end);
 
  return;
}
*/


/*
=head2 register_contig

  Description: DEPRECATED use register_component instead

=cut
*/

/* Comment out for now
void AssemblyMapperAdaptor_registerContig(AssemblyMapperAdaptor *ama, AssemblyMapper *asmMapper, char *type, IDType contigId) {

  fprintf(stderr,"Deprecated: Use register_component instead\n");

   //not sure if the use is passing in a seq_region_name or a
   //seq_region_id...
  AssemblyMapperAdaptor_registerComponent(ama, asmMapper, contigId);
}
*/


/*
=head2 fetch_by_type

  Description: DEPRECATED use fetch_by_CoordSystems instead

=cut
*/

AssemblyMapper *AssemblyMapperAdaptor_fetchByType(AssemblyMapperAdaptor *ama, char *type) {

  fprintf(stderr,"Deprecated: Use fetch_by_CoordSystem instead\n");

  //assume that what the user wanted was a mapper between the sequence coord
  //level and the top coord level

  CoordSystemAdaptor *csa  = DBAdaptor_getCoordSystemAdaptor(ama->dba);

  //CoordSystem *cs1 = CoordSystemAdaptor_fetchTopLevel(csa); // Was , type); but perl now ignores type in this call - not sure if it should
  fprintf(stderr, "Temporary hack in AssemblyMapperAdaptor_fetchByType - use 'chromosome' coord system instead of top level\n");
  CoordSystem *cs1 = CoordSystemAdaptor_fetchByName(csa, "chromosome", NULL); // NIY: Hack for testing
  CoordSystem *cs2 = CoordSystemAdaptor_fetchSeqLevel(csa);

  return AssemblyMapperAdaptor_fetchByCoordSystems(ama, cs1, cs2);
}

