#include "MetaCoordContainer.h"
#include "StrUtil.h"
#include "BaseAdaptor.h"
#include "MysqlUtil.h"
#include "StatementHandle.h"
#include "CoordSystemAdaptor.h"

MetaCoordContainer *MetaCoordContainer_new(DBAdaptor *dba) {
  MetaCoordContainer *mcc;

  if ((mcc = (MetaCoordContainer *)calloc(1,sizeof(MetaCoordContainer))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for MetaCoordContainer\n");
    return NULL;
  }
  BaseAdaptor_init((BaseAdaptor *)mcc, dba, META_CONTAINER);

  mcc->featureCache = StringHash_new(STRINGHASH_SMALL);
  mcc->maxLenCache = IDHash_new(IDHASH_SMALL);
  // 
  // Retrieve the list of the coordinate systems that features are stored
  // in and cache them.
  //

// Uses dnadb
  CoordSystemAdaptor *csa = DBAdaptor_getCoordSystemAdaptor(mcc->dba->dnadb);

  Vector *coordSystems = CoordSystemAdaptor_fetchAll(csa);

  StatementHandle *sth;
  char tmpStr[1024];
  char qStr[1024];

  sprintf(qStr,
          "SELECT mc.table_name, mc.coord_system_id, mc.max_length "
                 "FROM meta_coord mc "
                 "WHERE mc.coord_system_id in (");

  int i;
  for (i=0; i<Vector_getNumElement(coordSystems); i++) {
    CoordSystem *cs = Vector_getElementAt(coordSystems, i);
    if (i!=0) {
      strcat(qStr,", ");
    }
    sprintf(tmpStr,IDFMTSTR,CoordSystem_getDbID(cs));
    strcat(qStr, tmpStr);
  }
  strcat(qStr,")");

  sth = mcc->prepare((BaseAdaptor *)mcc,qStr,strlen(qStr)); 
  sth->execute(sth);

  ResultRow *row;
  while ((row = sth->fetchRow(sth))) {
    char *tableName = row->getStringAt(row, 0);
    IDType csId     = row->getLongLongAt(row, 1);
    long maxLen     = row->getLongAt(row, 2);
    
    StrUtil_strlwr(tableName);

    if (!StringHash_contains(mcc->featureCache, tableName)) {
      StringHash_add(mcc->featureCache, tableName, Vector_new());
    }

    Vector *vec = StringHash_getValue(mcc->featureCache, tableName);
    
    IDType *csIdP;
    if ((csIdP = calloc(1,sizeof(IDType))) == NULL) {
      fprintf(stderr,"Failed allocating space for id\n");
      exit(1);
    }
    *csIdP = csId;
    Vector_addElement(vec, csIdP);

    if (!IDHash_contains(mcc->maxLenCache, csId)) {
      IDHash_add(mcc->maxLenCache, csId, StringHash_new(STRINGHASH_SMALL));
    }
    StringHash *hash = IDHash_getValue(mcc->maxLenCache, csId);

    long *maxLenP;
    if ((maxLenP = calloc(1,sizeof(long))) == NULL) {
      fprintf(stderr,"Failed allocating space for maxLen\n");
      exit(1);
    }
    *maxLenP = maxLen;
    StringHash_add(hash, tableName, maxLenP);
  }
  sth->finish(sth);

  return mcc;
}


/*
=head2 fetch_all_CoordSystems_by_feature_type

  Arg [1]    : string $table - the name of the table to retrieve coord systems
               for.  E.g. 'gene', 'exon', 'dna_align_feature'
  Example    : @css = @{$mcc->fetch_all_CoordSystems_by_feature_type('gene')};
  Description: This retrieves the list of coordinate systems that features
               in a particular table are stored.  It is used internally by
               the API to perform queries to these tables and to ensure that
               features are only stored in appropriate coordinate systems.
  Returntype : listref of Bio::EnsEMBL::CoordSystem objects
  Exceptions : throw if name argument not provided
  Caller     : BaseFeatureAdaptor
  Status     : Stable

=cut
*/
Vector *MetaCoordContainer_fetchAllCoordSystemsByFeatureType(MetaCoordContainer *mcc, char *origTable) {
  char table[1024];

  if (origTable == NULL || !origTable[0]) {
    fprintf(stderr, "Name argument is required\n");
    exit(1);
  }

  strcpy(table, origTable);
  StrUtil_strlwr(table); //case insensitive matching
 
  if (!StringHash_contains(mcc->featureCache, table)) {
    return Vector_new();
  }

  Vector *csIds = StringHash_getValue(mcc->featureCache, table);
  Vector *coordSystems = Vector_new();

  CoordSystemAdaptor *csa = DBAdaptor_getCoordSystemAdaptor(mcc->dba);

  int i;
  for (i=0; i<Vector_getNumElement(csIds); i++) {
    IDType csId = *(IDType *)(Vector_getElementAt(csIds, i));
    CoordSystem *cs = CoordSystemAdaptor_fetchByDbID(csa, csId);

    if (cs == NULL) {
      fprintf(stderr,"meta_coord table refers to non-existant coord_system with id "IDFMTSTR, csId);
    }

    Vector_addElement(coordSystems, cs);
  }

  return coordSystems;
}

/*
=head2 fetch_max_length_by_CoordSystem_feature_type

  Arg [1]    : Bio::EnsEMBL::CoordSystem $cs
  Arg [2]    : string $table
  Example    : $max_len = $mcc->fetch_max_length_by_CoordSystem_feature_type($cs,'gene');
  Description: Returns the maximum length of features of a given type in
               a given coordinate system.
  Returntype : int or undef
  Exceptions : throw on incorrect argument
  Caller     : BaseFeatureAdaptor
  Status     : Stable

=cut
*/
long MetaCoordContainer_fetchMaxLengthByCoordSystemFeatureType(MetaCoordContainer *mcc, CoordSystem *cs, char *table) {
  char lcTable[1024];

  if (cs == NULL) {
    fprintf(stderr, "CoordSystem argument is NULL - bye!\n");
    exit(1);
  }

  if (table == NULL) {
    fprintf(stderr, "Table name argument is NULL - bye!\n");
    exit(1);
  }
  strcpy(lcTable, table);
  StrUtil_strlwr(lcTable);

  if (!IDHash_contains(mcc->maxLenCache, CoordSystem_getDbID(cs))) {
    fprintf(stderr,"Coord system with id "IDFMTSTR" not found in maxLenCache - returning -1\n", CoordSystem_getDbID(cs));
    return -1;
  }
 
  StringHash *hash = IDHash_getValue(mcc->maxLenCache, CoordSystem_getDbID(cs));

  if (!StringHash_contains(hash, lcTable)) {
    fprintf(stderr,"Table %s in coord system with id "IDFMTSTR" not found in maxLenCache - returning -1\n", table, CoordSystem_getDbID(cs));
    return -1;
  }

  return *((long *)StringHash_getValue(hash, lcTable));
}


/*

=head2 add_feature_type

  Arg [1]    : Bio::EnsEMBL::CoordSystem $cs
               The coordinate system to associate with a feature table
  Arg [2]    : string $table - the name of the table in which features of
               a given coordinate system will be stored in
  Arg [3]    : int $length
               This length is used to update the max_length in the database
               and the internal cache. 
  Example    : $csa->add_feature_table($chr_coord_system, 'gene');
  Description: This function tells the coordinate system adaptor that
               features from a specified table will be stored in a certain
               coordinate system.  If this information is not already stored
               in the database it will be added.
  Returntype : none
  Exceptions : none
  Caller     : BaseFeatureAdaptor
  Status     : Stable

=cut

*/
void MetaCoordContainer_addFeatureType(MetaCoordContainer *mcc, CoordSystem *cs, char *table, long length) {

  if (cs == NULL) {
    fprintf(stderr, "CoordSystem argument is required.\n");
    exit(1);
  }

  if (table == NULL) {
    fprintf(stderr, "Table argument is required.\n");
    exit(1);
  }

  // I really don't like having to do this lc here - it should be done before calling this routine
  // This routine is called for every feature stored so should be efficient
  char lcTable[1024];
  strcpy(lcTable, table);
  StrUtil_strlwr(lcTable);
  
  // Slightly different logic so I don't need to track if I've created a temporary Vector which would need freeing
  // If featureCache doesn't exist then add entry to cache and fetch
  if (!StringHash_contains(mcc->featureCache, lcTable)) {
    StringHash_add(mcc->featureCache, lcTable, Vector_new());
  }
  Vector *csIds = StringHash_getValue(mcc->featureCache, lcTable);

  //my $cs_ids = $self->{'_feature_cache'}->{$table} || [];

  int i;
  int exists = 0;
  for (i=0; i<Vector_getNumElement(csIds)  && !exists; i++) {
    IDType csId = *(IDType *)(Vector_getElementAt(csIds, i));
    if (csId == CoordSystem_getDbID(cs)) {
      exists = 1;
    }
  }

  if (!IDHash_contains(mcc->maxLenCache, CoordSystem_getDbID(cs))) {
    IDHash_add(mcc->maxLenCache, CoordSystem_getDbID(cs), StringHash_new(STRINGHASH_SMALL));
  }
  StringHash *tabMaxLenHash = IDHash_getValue(mcc->maxLenCache, CoordSystem_getDbID(cs));
  int inHash = StringHash_contains(tabMaxLenHash, lcTable);
  char qStr[1024];
 
  if (exists) {
   
    if ( ! inHash || *((long *)StringHash_getValue(tabMaxLenHash, lcTable)) < length ) {

      sprintf(qStr, "UPDATE meta_coord "
                    "SET max_length = %ld "
                    "WHERE coord_system_id = " IDFMTSTR " " 
                    "AND table_name = '%s' "
                    "AND (max_length < %ld OR max_length is null)", 
                    length, CoordSystem_getDbID(cs), lcTable, length);

      StatementHandle *sth = mcc->prepare((BaseAdaptor *)mcc,qStr,strlen(qStr)); 
      //$sth->execute( $cs->dbID(), $table );
      sth->execute(sth);

      long *maxLenP;
      if ((maxLenP = calloc(1,sizeof(long))) == NULL) {
        fprintf(stderr,"Failed allocating space for maxLen\n");
        exit(1);
      }
      *maxLenP = length;

      if (inHash) {
        StringHash_remove(tabMaxLenHash, lcTable, free);
      }
      StringHash_add(tabMaxLenHash, lcTable, maxLenP);

      //$self->{'_max_len_cache'}->{$cs->dbID()}->{$table} = $length;

      sth->finish(sth);
    }
    return;
  }

  // store the new tablename -> coord system relationship in the db
  // ignore failures b/c during the pipeline multiple processes may try
  // to update this table and only the first will be successful
  sprintf(qStr, "INSERT IGNORE INTO meta_coord "
                "SET coord_system_id = "IDFMTSTR", table_name = '%s', max_length = %ld",
			  CoordSystem_getDbID(cs), lcTable, length);

  StatementHandle *sth = mcc->prepare((BaseAdaptor *)mcc,qStr,strlen(qStr)); 

  sth->execute(sth);

  //$sth->execute($cs->dbID, $table, $length );

  //update the internal cache
  // First the featureCache
  IDType *csIdP;
  if ((csIdP = calloc(1,sizeof(IDType))) == NULL) {
    fprintf(stderr,"Failed allocating space for id\n");
    exit(1);
  }
  *csIdP = CoordSystem_getDbID(cs);
  Vector_addElement(csIds, csIdP);


  // Next the maxLenCache
  long *maxLenP;
  if ((maxLenP = calloc(1,sizeof(long))) == NULL) {
    fprintf(stderr,"Failed allocating space for maxLen\n");
    exit(1);
  }
  *maxLenP = length;
// Note - would be good to have a 'replace' function but I need to refactor IDHash slightly to do that - add a setFreeFunc method to it, so can free the old hash value
// For now just remove and then add
  if (inHash) {
    StringHash_remove(tabMaxLenHash, lcTable, free);
  }
  StringHash_add(tabMaxLenHash, lcTable, maxLenP);

  //$self->{'_max_len_cache'}->{$cs->dbID()}->{$table} = $length;

  sth->finish(sth);

  return;
}

