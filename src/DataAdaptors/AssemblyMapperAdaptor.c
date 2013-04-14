#include "AssemblyMapperAdaptor.h"
#include "BaseAdaptor.h"
#include "MysqlUtil.h"
#include "AssemblyMapper.h"
//#include "GenomicRange.h"

#include "StatementHandle.h"
#include "ResultRow.h"

int CHUNKFACTOR = 20;  // 2^20 = approx. 10^6

// May move this
typedef struct seqRegionCacheEntryStruct {
  char * regionName;
  IDType regionId;
  IDType csId;
  long   regionLength;
} SeqRegionCacheEntry;

static char *AMA_FIRST = "first";
static char *AMA_LAST  = "last";

int trueVal = 1;

AssemblyMapperAdaptor *AssemblyMapperAdaptor_new(DBAdaptor *dba) {
  AssemblyMapperAdaptor *ama;

  if ((ama = (AssemblyMapperAdaptor *)calloc(1,sizeof(AssemblyMapperAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for AssemblyMapperAdaptor\n");
    return NULL;
  }
  BaseAdaptor_init((BaseAdaptor *)ama, dba, ASSEMBLYMAPPER_ADAPTOR);

  ama->asmMapperCache = StringHash_new(STRINGHASH_LARGE);

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

void AssemblyMapperAdaptor *cacheSeqIdsWithMultAssemblies(AssemblyMapperAdaptor *ama) {
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

  sth = ama->prepare((BaseAdaptor *)ama,qStr,strlen(qStr));
  sth->execute(sth);

  while (row = sth->fetchRow(sth)) {
    IDType dbID    = row->getLongLongAt(row,0);

    IDHash_add(ama->multSeqIdCache, dbID, trueVal);
  }

  sth->finish(sth);
}


char *makeMappingPathKey(Vector *path) {
  int i;
  char *key;

  StrUtil_copyString(&key, "", 0);
  

  for (i=0; i<Vector_getNumElement(path); i++) {
    CoordSystem *cs = Vector_getElementAt(path, i);
    if (cs == NULL) {
      StrUtil_appendString(&key, "-", 0);
    } else {
      char tmp[1024];
      sprintf(tmp,IDFMTSTR,CoordSystem_getDbID(cs));
      StrUtil_appendString(&key, tmp, 0 );
    }
    if (i < Vector_getNumElement(path)-1) StrUtil_appendString(":");
  }
  
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

AssemblyMapper *AssemblyMapper_fetchByCoordSystems(AssemblyMapperAdaptor *ama, CoordSystem *cs1, CoordSystem *cs2) {
  if (CoordSystem_isTopLevel(cs1)) {
    return TopLevelAssemblyMapper_new(ama, cs1, cs2);
  }

  if (CoordSystem_isTopLevel(cs2)) {
    return TopLevelAssemblyMapper_new(ama, cs2, cs1);
  }

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

  char *key = makeMappingPathKey(mappingPath);

  if (StringHash_contains(ama->asmMapperCache, key)) {
    return StringHash_getValue(ama-asmMapperCache, key);
  }

  switch (Vector_getNumElement(mappingPath)) {
    case 1:
      fprintf(stderr,"Incorrect mapping path defined in meta table. 0 step mapping encountered between:\n"
                     "%s %s and %s %s\n", 
                     CoordSystem_getName(cs1), CoordSystem_getVersion(cs1),
                     CoordSystem_getName(cs2), CoordSystem_getVersion(cs2));
      break;

    case 2:

      // 1 step regular mapping
      AssemblyMapper *asmMapper = AssemblyMapper_new(ama, mappingPath);

      //   If you want multiple pieces on two seqRegions to map to each other
      //   you need to make an assembly.mapping entry that is seperated with a #
      //   instead of an |.

      StringHash_add(ama->asmMapperCache, key, asmMapper);

      free(key);

      return asmMapper;
      break;
  
    case 3:
      // two step chained mapping
      ChainedAssemblyMapper *asmMapper = ChainedAssemblyMapper_new(ama, mappingPath);

      // in multi-step mapping it is possible get requests with the
      // coordinate system ordering reversed since both mappings directions
      // cache on both orderings just in case
      // e.g.   chr <-> contig <-> clone   and   clone <-> contig <-> chr
  
      StringHash_add(ama->asmMapperCache, key, asmMapper);
      free(key);

      Vector *revMappingPath = Vector_reverse(mappingPath);
      key = makeMappingPathKey(revMappingPath);
      StringHash_add(ama->asmMapperCache, key, asmMapper);

      Vector_free(revMappingPath, NULL);
      free(key);

      return (AssemblyMapper *)asmMapper;
      break;

    default:
      fprintf(stderr, "Only 1 and 2 step coordinate system mapping is currently\n" 
                      "supported.  Mapping between %s %s and %s %s requires %d steps.\n",
                     CoordSystem_getName(cs1), CoordSystem_getVersion(cs1),
                     CoordSystem_getName(cs2), CoordSystem_getVersion(cs2)
                     Vector_getNumElement(mappingPath));
      exit(1);
  }
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
  SeqRegionRange *range;

  if ((range = (SeqRegionRange *)calloc(1,sizeof(SeqRegionRange))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for SeqRegionRange\n");
    exit(1);
  }
  SeqRegionRange_setSeqRegionStart(range, start);
  SeqRegionRange_setSeqRegionEnd(range, end);
  SeqRegionRange_setSeqRegionId(range, id);

  if (name) SeqRegionRange_setSeqRegionName(range, name);

  Vector_addElement(range);

  return range; 
}

void AssemblyMapperAdaptor_registerAssembled(AssemblyMapperAdaptor *ama, AssemblyMapper *asm, IDType asmSeqRegion, long asmStart, long asmEnd) {

  IDType asmCsId = CoordSystem_getDbID(AssemblyMapper_getAssembledCoordSystem(asm));
  IDType cmpCsId = CoordSystem_getDbID(AssemblyMapper_getComponentCoordSystem(asm));

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
  if (AssemblyMapper_getSize(asmMapper) > AssemblyMapper_maxPairCount(asmMapper) ) {
    AssemblyMapper_flush(asmMapper);

    // we now have to go and register the entire requested region since we
    // just flushed everything
    Vector_free(chunkRegions, SeqRegionRange_free);
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

    regionStart = SeqRegionRange_getSeqRegionStart(region);
    regionEnd   = SeqRegionRange_getSeqRegionEnd(region);

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
                    " AND %ld <= asm.asm_end AND"
                    " AND %ld >= asm.asm_start AND"
                    " AND asm.cmp_seq_region_id = sr.seq_region_id AND"
                    " AND sr.coord_system_id = "IDFMTSTR, 
                   asmSeqRegion, regionStart, regionEnd, cmpCsId);
  

    sth = ama->prepare((BaseAdaptor *)ama,qStr,strlen(qStr));
    sth->execute(sth);


    // 
    // Load the unregistered regions of the mapper
    // 
    while (row = sth->fetchRow(sth)) {
      long   cmpStart           = row->getLongAt(row,0);
      long   cmpEnd             = row->getLongAt(row,1);
      IDType cmpSeqRegionId     = row->getLongLongAt(row,2);
      char * cmpSeqRegionName   = row->getStringAt(row,3);
      long   cmpSeqRegionLength = row->getLongAt(row,4);
      int    ori                = row->getIntAt(row,5);
      long   asmStart           = row->getLongAt(row,6);
      long   asmEnd             = row->getLongAt(row,7);

      if (AssemblyMapper_haveRegisteredComponent(asmMapper, cmpSeqRegionId) &&
          !IDHash_contains(ama->multiSeqIdCache, cmpSeqRegionId)) {
        continue;
      }
      AssemblyMapper_registerComponent(asmMapper, cmpSeqRegionId);

      Mapper *mapper = AssemblyMapper_getMapper(asmMapper);

      Mapper_addMapCoordinates(mapper,
                 asmSeqRegion, asmStart, asmEnd, ori,
                 cmpSeqRegionId, cmpStart, cmpEnd);

      AssemblyMapperAdaptor_addToSrCaches(cmpSeqRegion, cmpSeqRegionId, cmpCsId, cmpSeqRegionLength);
    }
  }

  Vector_free(chunkRegions, SeqRegionRange_free);

  sth->finish(sth);
}

IDType AssemblyMapperAdaptor_seqRegionNameToId(AssemblyMapperAdaptor *ama, char *srName, IDType csId) {
  if (srName == NULL || csID <= 0) {
    fprintf(stderr, "seq_region_name and coord_system_id args are required\n");
    exit(1);
  }

  char key[1024];
  sprintf(key,"%s:"IDFMTSTR, srName, csID);
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
                 "WHERE name = %s AND coord_system_id = "
                 IDFMTSTR, srName, csId);

  StatementHandle *sth = ama->prepare((BaseAdaptor *)ama,qStr,strlen(qStr));
  sth->execute(sth);

  if (sth->numRows(sth) != 1) {
    fprintf(stderr,"Ambiguous or non-existant seq_region [%s] in coord system "IDFMTSTR, srName,csId);
    exit(1);
  }

  
  ResultRow *row = sth->fetchRow(sth);
  IDType srId   = row->getLongLongAt(row,0);
  long srLength = row->getLongAt(row,1);


  AssemblyMapperAdaptor_addToSrCaches(srName, srId, csId, srLength);

  sth->finish(sth);
  return srId;
}

char *AssemblyMapperAdaptor_seqRegionIdToName(AssemblyMapperAdaptor *ama, IDType srId) {

  if (!srId) {
    fprintf(stderr,"seq_region_id is required");
    exit(1);
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
    exit(1);
  }

  ResultRow *row = sth->fetchRow(sth);
//NIY: Do I need to alloc this?
  char *srName  = row->getStringAt(row,0);
  long srLength = row->getLongAt(row,1);
  IDType csID   = row->getLongLongAt(row,2);


  AssemblyMapperAdaptor_addToSrCaches(srName, srId, csId, srLength);

  sth->finish(sth);
  return srName;
}

void AssemblyMapperAdaptor_addToSrCaches(char *regionName, IDType regionId, IDType csId, long regionLength) {
  char key[1024];
  SeqRegionCacheEntry *cacheData;

  // Allocate the struct
  if ((cacheData = (SeqRegionCacheEntry *)calloc(1,sizeof(SeqRegionCacheEntry))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for SeqRegionCacheEntry\n");
    exit(1);
  }
  
  StrUtil_copyString(cacheData->regionName, regionName, 0);
  cacheData->regionId     = regionId;
  cacheData->csId         = csId;
  cacheData->regionLength = regionLength;
  
  sprintf(key,"%s:"IDFMTSTR, regionName, csId);

  // Do a quick sanity check
  if (StringHash_contains(ama->srNameCache, key)) {
    fprintf(stderr,"Hmm - seq region already in name cache - odd\n");
  }
  StringHash_add(ama->srNameCache, key, cacheData);

  if (IDHash_contains(ama->srIdCache, regionId)) {
    fprintf(stderr,"Hmm - seq region already in id cache - odd\n");
  }
  IDHash_add(ama->srIdCache, regionId, cacheData);

  return;
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

  IDType cmpCsId = CoordSystem_getDbID(AssemblyMapper_getComponentCoordSystem(asm));
  IDType asmCsId = CoordSystem_getDbID(AssemblyMapper_getAssembledCoordSystem(asm));

  // do nothing if this region is already registered or special case

  if (AssemblyMapper_haveRegisteredComponent(asmMapper, cmpSeqRegion) &&
      !IDHash_contains(ama->multiSeqIdCache, cmpSeqRegion)) {
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
              " FROM""
                  " assembly asm, seq_region sr"
             " WHERE""
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
                   "Remember that multiple mappings use the #-operaator in the meta-table (i.e. chromosome:EquCab2\#contig\n",
                   cmpSeqRegion);
    exit(1);
  }

  ResultRow *row = sth->fetchRow(sth);
  long asmStart           = row->getLongAt(row,0);
  long asmEnd             = row->getLongAt(row,1);
  IDType asmSeqRegionId   = row->getLongLongAt(row,2);
//NIY: Do I need to alloc this?
  char *asmSeqRegion      = row->getStringAt(row,3);
  long asmSeqRegionLength = row->getLongAt(row,4);

  AssemblyMapperAdaptor_addToSrCaches(asmSeqRegion, asmSeqRegionId, asmCsId, asmSeqRegionLength);


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

  IDType toSeqRegionId;

  if (toSlice != NULL) {
    if (!CoordSystem_compare(ChainedAssemblyMapper_getFirstCoordSystem(casmMapper), ChainedAssemblyMapper_getLastCoordSystem(casmMapper)) {
      return AssemblyMapperAdaptor_registerChainedSpecial(ama, casmMapper, from, seqRegionId, ranges, toSlice);
    }

    IDType toSeqRegionId = Slice_getSeqRegionId(toSlice);

    // NIY can we do not defined?? Use 0 for now if (!defined($to_seq_region_id)){
    if (!toSeqRegionId)){
      fprintf(stderr, "Could not get seq_region_id for to_slice %s\n", Slice_getSeqRegionName(toSlice));
      exit(1);
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
    exit(1);
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

  // NIY Free path

  if (Vector_getNumElement(path) != 2 && Vector_getElementAt(path,1) != NULL) {
    char *pathStr = makeMappingPathKey(path);
    int len  = Vector_getNumElement(path)-1;

    fprintf(stderr, "Unexpected mapping path between start and intermediate coord systems (%s  %s and %s %s).\n"
                    "Expected path length 1, got %d. (path=%s).\n", 
                     CoordSystem_getName(startCs), CoordSystem_getVersion(startCs),
                     CoordSystem_getName(midCs), CoordSystem_getVersion(midCs),
                     len, pathStr);
    exit(1);
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

  StatementHandle *sth;
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
    SeqRegionRange *range = Vector_getElementAt(ranges, i);

    long start = SeqRegionRange_getSeqRegionStart(range);
    long end   = SeqRegionRange_getSeqRegionEnd(range);

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
    while (row = sth->fetchRow(sth)) {
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
      AssemblyMapperAdaptor_addToSrCaches(midSeqRegion, midSeqRegionId, midCsId, midLength);

      AssemblyMapperAdaptor_addToRangeVector(midRanges, midSeqRegionId, midStart, midEnd, midSeqRegion);  

      AssemblyMapperAdaptor_addToRangeVector(startRanges, seqRegionId, startStart, startEnd, NULL);  

      //the region that we actually register may actually be larger or smaller
      //than the region that we wanted to register.
      //register the intersection of the region so we do not end up doing
      //extra work later

      if (startStart < start || startEnd > end) {
        RangeRegistry_checkAndRegister(startRgistry, seqRegionId, startStart, startEnd);
      }
    }
    sth->finish(sth);
  }

  // in the one step case, we load the mid ranges in the
  // last_registry and we are done
  if (midCs == NULL) {
    int i;
    for (i=0;i<Vector_getNumElement(midRanges); i++) {
      SeqRegionRange *range = Vector_getElementAt(midRanges, i);

      RangeRegistry_checkAndRegister(endRegistry, 
                                     SeqRegionRange_getSeqRegionId(range), 
                                     SeqRegionRange_getSeqRegionStart(range), 
                                     SeqRegionRange_getSeqRegionEnd(range));
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
  Vector *path = CoordSystemAdaptor_getMappingPath(csa, midCs, endCs);
  if (Vector_getNumElement(path) == 2 || ( Vector_getNumElement(path) == 3 && Vector_getElementAt(path, 1) == NULL)) {
    asmCs = Vector_getElementAt(path, 0);
    cmpCs = Vector_getLastElement(path);
  } else {
    char *pathStr = makeMappingPathKey(path);
    int len  = Vector_getNumElement(path)-1;

    fprintf(stderr, "Unexpected mapping path between intermediate and last coord systems (%s  %s and %s %s).\n"
                    "Expected path length 1, got %d. (path=%s).\n", 
                     CoordSystem_getName(midCs), CoordSystem_getVersion(midCs),
                     CoordSystem_getName(endCs), CoordSystem_getVersion(endCs),
                     len, pathStr);
    exit(1);
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
    char *midSeqRegion    = SeqRegionRange_getSeqRegionName(midRange); 
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
    while (row = sth->fetchRow(sth)) {
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
      AssemblyMapperAdaptor_addToSrCaches(endSeqRegion, endSeqRegionId, endCsId, endLength);

      //register this region on the end coord system
      RangeRegistry_checkAndRegister(endRegistry, endSeqRegionId, endStart, endEnd);
    }
    sth->finish(sth);
  }

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
    exit(1);
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
    path = CoordSystemAdaptor_getMappingPath(startCs, midCs);
  } else {
    path = CoordSystemAdaptor_getMappingPath(startCs, endCs);
  }

  if (midCs == NULL) {
    startMidMapper = combinedMapper;
  }

  if (Vector_getNumElement(path) != 2 && Vector_getElementAt(path,1) != NULL) {
    char *pathStr = makeMappingPathKey(path);
    int len  = Vector_getNumElement(path)-1;

    fprintf(stderr, "Unexpected mapping path between start and intermediate coord systems (%s  %s and %s %s).\n"
                    "Expected path length 1, got %d. (path=%s).\n", 
                     CoordSystem_getName(startCs), CoordSystem_getVersion(startCs),
                     CoordSystem_getName(midCs), CoordSystem_getVersion(midCs),
                     len, pathStr);
    exit(1);
  }

  CoordSystem *asmCs;
  CoordSystem *cmpCs;
  asmCs = Vector_getElementAt(path, 0);
  cmpCs = Vector_getLastElement(path);

  Mapper *combinedMapper  = ChainedAssemblyMapper_getFirstLastMapper(casmMapper);
  CoordSystem *midCs      = ChainedAssemblyMapper_getMiddleCoordSystem(casmMapper);
  CoordSystemAdaptor *csa = DBAdaptor_getCoordSystemAdaptor(ama->dba);
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
      SeqRegionRange *range = Vector_getElementAt(ranges, i);

      long start = SeqRegionRange_getSeqRegionStart(range);
      long end   = SeqRegionRange_getSeqRegionEnd(range);

      sth->execute(sth, id1, start, end, CoordSystem_getDbID(toCs), id2);

      long midStart, midEnd;
      IDType midSeqRegionId;
      char *midSeqRegion;
      long midLength;
      int ori;
      long startStart, startEnd;

      ResultRow *row;
      while (row = sth->fetchRow(sth)) {
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
        AssemblyMapperAdaptor_addToSrCaches(midSeqRegion, midSeqRegionId, midCsId, midLength);
        
        AssemblyMapperAdaptor_addToRangeVector(midRanges, midSeqRegionId, midStart, midEnd, midSeqRegion);  

        AssemblyMapperAdaptor_addToRangeVector(startRanges, id1, startStart, startEnd, NULL);  
        
        //the region that we actually register may actually be larger or smaller
        //than the region that we wanted to register.
        //register the intersection of the region so we do not end up doing
        //extra work later
        
        if (startStart < start || startEnd > end) {
          RangeRegistry_checkAndRegister(startRegistry, id1, startStart, startEnd);
        }
      }
      sth->finish(sth);
    }

    if (found) {
      if (midCs == NULL ) {
        int j;
        for (j=0;j<Vector_getNumElement(midRanges); j++) {
          SeqRegionRange *range = Vector_getElementAt(midRanges, j);

          RangeRegistry_checkAndRegister(endRegistry, 
                                         SeqRegionRange_getSeqRegionId(range),
                                         SeqRegionRange_getSeqRegionStart(range),
                                         SeqRegionRange_getSeqRegionEnd(range));
        }
        
        // and thats it for the simple case ...
        Vector_free(midRanges, NULL);
        Vector_free(startRanges, NULL);
        return;
      }
    }
  }

  Vector_free(midRanges, NULL);
  Vector_free(startRanges, NULL);
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
  IDHash *asmRegistered = IDHash_new();

  ResultRow *row;
  while (row = sth->fetchRow(sth)) {
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

    AssemblyMapperAdaptor_addToSrCaches(cmpSeqRegion, cmpSeqRegionId, cmpCsId, cmpLength);

    // only register each asm seq_region once since it requires some work
    if ( ! IDHash_contains(asmRegistered, asmSeqRegionId)) {
      IDHash_add(asmRegistered, asmSeqRegionId, &trueVal);

      // register all chunks from start of seq region to end
      int endChunk = asmLength >> CHUNKFACTOR;
      int i;
      for (i=0; i<=endChunk; i++) {
        AssmeblyMapper_registerAssembled(asmSeqRegionId, i);
      }

      AssemblyMapperAdaptor_addToSrCaches(asmSeqRegion, asmSeqRegionId, asmCsId, asmLength);
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

  Vector *ranges;

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


  CoordSystemAdaptor *csa        = DBAdaptor_getCoordSystemAdaptor(ama->dba);

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
    char *pathStr = makeMappingPathKey(path);
    int len  = Vector_getNumElement(path)-1;

    fprintf(stderr, "Unexpected mapping path between start and intermediate coord systems (%s  %s and %s %s).\n"
                    "Expected path length 1, got %d. (path=%s).\n", 
                     CoordSystem_getName(startCs), CoordSystem_getVersion(startCs),
                     CoordSystem_getName(midCs), CoordSystem_getVersion(midCs),
                     len, pathStr);
    exit(1);
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
  AssemblyMapper *mapper;
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
  while (row = sth->fetchRow(sth)) {
    if (!CoordSystem_compare(asmCs,firstCs)) {
      long midStart           = row->getLongAt(row,0);
      long midEnd             = row->getLongAt(row,1);
      IDType midSeqRegionId   = row->getLongLongAt(row,2);
      char *midSeqRegion      = row->getStringAt(row,3);
      long midLength          = row->getLongAt(row,4);
      int  ori                = row->getIntAt(row,5);
      long startStart         = row->getLongAt(row,6);
      long startEnd           = row->getLongAt(row,7);
      IDType startSeqRegionId = row->getLongLongAt(row,8);
      char *startSeqRegion    = row->getStringAt(row,9);
      long startLength        = row->getLongAt(row,10);

    } else {
      long startStart         = row->getLongAt(row,0);
      long startEnd           = row->getLongAt(row,1);
      IDType startSeqRegionId = row->getLongLongAt(row,2);
      char *startSeqRegion    = row->getStringAt(row,3);
      long startLength        = row->getLongAt(row,4);
      int  ori                = row->getIntAt(row,5);
      long midStart           = row->getLongAt(row,6);
      long midEnd             = row->getLongAt(row,7);
      IDType midSeqRegionId   = row->getLongLongAt(row,8);
      char *midSeqRegion      = row->getStringAt(row,9);
      long midLength          = row->getLongAt(row,10);
    }

    Mapper_addMapCoordinates(mapper,
       startSeqRegionId, startStart, startEnd, ori,
       midSeqRegionId, midStart, midEnd);

    AssemblyMapperAdaptor_addToRangeVector(ranges, startSeqRegionId, startStart, startEnd, NULL);  

    RangeRegistry_checkAndRegister(reg, startSeqRegionId, 1, startLength );

    if ( midCs == NULL ) {
      RangeRegistry *lastReg = ChainedAssemblyMapper_getLastRegistry(casmMapper);

      RangeRegistry_checkAndRegister(lastReg, midSeqRegionId, midStart, midEnd);
    }

    AssemblyMapperAdaptor_addToSrCaches(midSeqRegion, midSeqRegionId, midCsId, midLength);

    AssemblyMapperAdaptor_addToSrCaches(startSeqRegion, startSeqRegionId, startCsId, startLength);
  }

  if (midCs == NULL) {
    // thats it for the simple case
    return;
  }


  Vector *path = CoordSystemAdaptor_getMappingPath(csa, lastCs, midCs);
  if ( midCs ) {
    if (Vector_getElementAt(path,1) == NULL) {
      Vector_removeElementAt(path, 1);
    }
  }

  if (Vector_getNumElement(path) != 2 && Vector_getElementAt(path,1) != NULL) {
    char *pathStr = makeMappingPathKey(path);
    int len  = Vector_getNumElement(path)-1;

    fprintf(stderr, "Unexpected mapping path between intermediate and last coord systems (%s  %s and %s %s).\n"
                    "Expected path length 1, got %d. (path=%s).\n", 
                     CoordSystem_getName(midCs), CoordSystem_getVersion(midCs),
                     CoordSystem_getName(endCs), CoordSystem_getVersion(endCs),
                     len, pathStr);
    exit(1);
  }

  CoordSystem *asmCs = Vector_getElementAt(path,0);
  CoordSystem *cmpCs = Vector_getElementAt(path,1);

  sth->execute(sth, CoordSystem_getDbID(asmCs), CoordSystem_getDbID(cmpCs));

  long endStart, endEnd;
  IDType endSeqRegionId;
  char *endSeqRegion,
  long endLength;


  IDType endCsId = CoordSystem_getDbID(lastCs);
  reg = ChainedAssemblyMapper_getLastRegistry(casmMapper);

  while (row = sth->fetchRow(sth)) {
    if (!CoordSystem_compare(asmCs,midCs)) {
      long endStart         = row->getLongAt(row,0);
      long endEnd           = row->getLongAt(row,1);
      IDType endSeqRegionId = row->getLongLongAt(row,2);
      char *endSeqRegion    = row->getStringAt(row,3);
      long endLength        = row->getLongAt(row,4);
      int  ori              = row->getIntAt(row,5);
      long midStart         = row->getLongAt(row,6);
      long midEnd           = row->getLongAt(row,7);
      IDType midSeqRegionId = row->getLongLongAt(row,8);
      char *midSeqRegion    = row->getStringAt(row,9);
      long midLength        = row->getLongAt(row,10);

    } else {
      long midStart         = row->getLongAt(row,0);
      long midEnd           = row->getLongAt(row,1);
      IDType midSeqRegionId = row->getLongLongAt(row,2);
      char *midSeqRegion    = row->getStringAt(row,3);
      long midLength        = row->getLongAt(row,4);
      int  ori              = row->getIntAt(row,5);
      long endStart         = row->getLongAt(row,6);
      long endEnd           = row->getLongAt(row,7);
      IDType endSeqRegionId = row->getLongLongAt(row,8);
      char *endSeqRegion    = row->getStringAt(row,9);
      long endLength        = row->getLongAt(row,10);
    }

    Mapper_addMapCoordinates(endMidMapper, 
       endSeqRegionId, endStart, endEnd, ori,
       midSeqRegionId, midStart, midEnd);

    RangeRegistry_checkAndRegister(reg, endSeqRegionId, 1, endLength );

    AssemblyMapperAdaptor_addToSrCaches(endSeqRegion, endSeqRegionId, endCsId, endLength);
  }

  sth->finish(sth);

  AssemblyMapperAdaptor_buildCombinedMapper(ama,  ranges, startMidMapper, endMidMapper,
                                            combinedMapper, "first" );

  Vector_free(ranges, SeqRegionRange_free);
  return;
}



// after both halves of a chained mapper are loaded
// this function maps all ranges in $ranges and loads the
// results into the combined mapper
sub AssemblyMapperAdaptor_buildCombinedMapper(AssemblyMapperAdaptor *ama, Vector *ranges, Mapper *startMidMapper, 
                                              Mapper *endMidMapper, Mapper *combinedMapper, char *startName) {
  int i;
  for (i=0; i<Vector_getNumElement(ranges); i++) {
    SeqRegionRange *range = Vector_getElementAt(ranges, i);
    
    IDType seqRegionId = SeqRegionRange_getSeqRegionId(range);
    long   start       = SeqRegionRange_getSeqRegionStart(range);
    long   end         = SeqRegionRange_getSeqRegionEnd(range);

    long sum = 0;

    Vector *initialCoords = Mapper_mapCoordinates(startMidMapper, seqRegionId, start, end, 1, startName);

    foreach my $icoord (@initial_coords) {
      //skip gaps
      if ($icoord->isa('Bio::EnsEMBL::Mapper::Gap')) {
        sum += $icoord->length();
        next;
      }


      //feed the results of the first mapping into the second mapper
      my @final_coords =
        $end_mid_mapper->map_coordinates($icoord->id, $icoord->start,
                                         $icoord->end,
                                         $icoord->strand, "middle");


      foreach my $fcoord (@final_coords) {
        //load up the final mapper

        if ($fcoord->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
          long totalStart = start + sum;
          long totalEnd   = totalStart + $fcoord->length - 1;
          int  ori = $fcoord->strand();

          if (!strcmp(startName,"first")) { // add coords in consistant order
            Mapper_addMapCoordinates(combinedMapper,
                             seqRegionId, totalStart, $totalEnd, ori,
                             $fcoord->id(), $fcoord->start(), $fcoord->end());
          } else {
            Mapper_addMapCoordinates(combinedMapper,
                        $fcoord->id(), $fcoord->start(), $fcoord->end(), ori,
                        seqRegionId, totalStart, totalEnd);
          }

        }
        sum += $fcoord->length();
      }
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

      if ((idP = (IDType *)calloc(1,sizeof(IDType))) == NULL) {
        fprintf(stderr, "ERROR: Failed allocating space for idP\n");
        exit(1);
      }
      *idP = id;
      
      Vector_addElement(out, idP);
    }
  }

  return out;
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

Vector *AssemblyMapperAdaptor_seqIdsToRegions(AssemblyMapperAdaptor *ama, CoordSystem *coordSystem, Vector *seqRegionIds) {
  Vector *out = Vector_new();

  int i;
  for (i=0; i< Vector_getNumElement(seqRegionIds); i++) {
    IDType *sr = Vector_getElementAt(seqRegionIds,i);

    if (IDHash_contains(ama->srIDCache, *sr)) {
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

  fprintf(stderr,"AMA_deletCache not implemented yet\n");
  exit(1);
  
/*

  %{$self->{'sr_name_cache'}}     = ();
  %{$self->{'sr_id_cache'}}       = ();

  foreach my $key (keys %{$self->{'_asm_mapper_cache'}}){
    $self->{'_asm_mapper_cache'}->{$key}->flush();
  }
  %{$self->{'_asm_mapper_cache'}} = ();
*/

  return;
}


/*
=head2 register_region

  Description: DEPRECATED use register_assembled instead

=cut
*/

void AssemblyMapperAdaptor_registerRegion(AssemblyMapperAdaptor *ama, AssemblyMapper *asmMapper, char *type??, 
                                          char *chrName, long start, long end){

  fprintf(stderr,"Deprecated: Use register_assembled instead\n");

  AssemblyMapperAdaptor_registerAssembled(ama, asmMapper, chrName, start, end);
 
  return;
}


/*
=head2 register_contig

  Description: DEPRECATED use register_component instead

=cut
*/

void AssemblyMapperAdaptor_registerContig(AssemblyMapperAdaptor *ama, AssemblyMapper *asmMapper, char *type, IDType contigId) {

  fprintf(stderr,"Deprecated: Use register_component instead\n");

   //not sure if the use is passing in a seq_region_name or a
   //seq_region_id...
  AssemblyMapperAdaptor_registerComponent(ama, asmMapper, contigId);
}


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

  CoordSystem *cs1 = CoordSystemAdaptor_fetchTopLevel(csa, type);
  CoordSystem *cs2 = CoordSystemAdaptor_fetchSequenceLevel(csa);

  return AssemblyMapperAdaptor_fetchByCoordSystems(ama, cs1, cs2);
}

