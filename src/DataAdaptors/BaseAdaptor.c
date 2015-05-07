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

#define __MAIN_C__
#include "BaseAdaptor.h"
#undef __MAIN_C__

#include "DBAdaptor.h"

void BaseAdaptor_init(BaseAdaptor *ba, DBAdaptor *dba, int adaptorType) {
  ba->dba = dba;
  ba->adaptorType = adaptorType;
  ba->prepare = BaseAdaptor_prepare;

  ba->objectsFromStatementHandle = BaseAdaptor_objectsFromStatementHandle;
  ba->getTables                  = BaseAdaptor_getTables;
  ba->getColumns                 = BaseAdaptor_getColumns;
  ba->finalClause                = BaseAdaptor_finalClause;
  ba->leftJoin                   = BaseAdaptor_leftJoin;
  ba->defaultWhereClause         = BaseAdaptor_defaultWhereClause;
  ba->store                      = BaseAdaptor_store;

  BaseAdaptor_setSpeciesId(ba, DBAdaptor_getSpeciesId(dba));
}

StatementHandle *BaseAdaptor_prepare(BaseAdaptor *ba, char *qStr, size_t len) {
  /*printf("Query = %s len = %d\n",qStr,len);*/
  return DBAdaptor_prepare(ba->dba,qStr,len);
}


/*
# list primary keys for a particular table
# args are table name and primary key field
# if primary key field is not supplied, tablename_id is assumed
# returns listref of IDs
*/
// For ordered, the default should be 0 (if you just need to fill out the args)
// Note ONLY stable_id can be char, all other pk's must be IDType (see code)
Vector *BaseAdaptor_listDbIDs(BaseAdaptor *ba, char *table, char *pk, int ordered) {
  char colName[1024];

  if (pk == NULL) {
    sprintf(colName, "%s_id", table);
  } else {
    strcpy(colName, pk);
  }

  char qStr[1024];
  sprintf(qStr,"SELECT `%s` FROM `%s`", colName, table );

  int joinWithCs = 0;

  if ( BaseAdaptor_isMultiSpecies(BaseAdaptor *ba)
      // For now just the multi species because I don't have adaptors in the Class hierarchy
      // && $self->isa('Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor')
      // && !$self->isa('Bio::EnsEMBL::DBSQL::UnmappedObjectAdaptor') 
     ) {
    char tmpStr[1024];
    
    sprintf(tmpStr, "JOIN seq_region USING (seq_region_id) "
                    "JOIN coord_system cs USING (coord_system_id) "
                    "WHERE cs.species_id = "IDFMTSTR, BaseAdaptor_getSpeciesId(ba));

    sprintf(qStr, "%s %s", qStr, tmpStr);
  }

  if (ordered) {
    sprintf(qStr, "%s ORDER BY seq_region_id, seq_region_start", qStr);
  }

  StatementHandle *sth = ba->prepare(ba,qStr,strlen(qStr));

  sth->execute(sth);

  Vector *out = Vector_new();

  if (strcmp(pk, "stable_id")) {
    ResultRow *row;
    while ((row = sth->fetchRow(sth))) {
      char *stableId = row->getStringCopyAt(row, 0);
  
      Vector_addElement(out, stableId);
    }
  } else  {
    IDType *idP;
    ResultRow *row;
    while ((row = sth->fetchRow(sth))) {
      IDType id = row->getLongLongAt(row, 0);
  
      if ((idP = calloc(1,sizeof(IDType))) == NULL) {
        fprintf(stderr, "Failed allocating space for a id\n");      
        exit(1);
      } 
  
      *idP = id;
      Vector_addElement(out, idP);
    }
  }

  return out;
}


/*
# _straight_join

#   Arg [1]    : (optional) boolean $new_val
#   Example    : $self->_straight_join(1);
#                $self->generic_fetch($constraint);
#                $self->_straight_join(0);
#   Description: PROTECTED Getter/Setter that turns on/off the use of 
#                a straight join in queries.
#   Returntype : boolean
#   Exceptions : none
#   Caller     : general
*/

/*
int _straight_join {
  my $self = shift;
  if(@_) {
    $self->{'_straight_join'} = shift;
  }

  return $self->{'_straight_join'};
}
*/


/*
=head2 bind_param_generic_fetch

 Arg [1]   : (optional)  scalar $param
              This is the parameter to bind
 Arg [2]   : (optional) int $sql_type
              Type of the parameter (from DBI (:sql_types))
 Example   :  $adaptor->bind_param_generic_fetch($stable_id,SQL_VARCHAR);
              $adaptor->generic_fetch();
 Description:  When using parameters for the query, will call the bind_param to avoid
               some security issues. If there are no arguments, will return the bind_parameters
 ReturnType : listref
 Exceptions:  if called with one argument

=cut

sub bind_param_generic_fetch{
    my $self = shift;
    my $param = shift;
    my $sql_type = shift;

    if (defined $param && !defined $sql_type){
	throw("Need to specify sql_type for parameter $param\n");
    }
    elsif (defined $param && defined $sql_type){
	#check when there is a SQL_INTEGER type that the parameter is really a number
	if ($sql_type eq SQL_INTEGER){
	    throw "Trying to assign a non numerical parameter to an integer value in the database" if ($param !~ /^\d+$/);
	}
	#both paramters have been entered, push it to the bind_param array
	push @{$self->{'_bind_param_generic_fetch'}},[$param,$sql_type];
    }
    elsif (!defined $param && !defined $sql_type){
	#when there are no arguments, return the array
	return $self->{'_bind_param_generic_fetch'};
    }
	
}

# Used to reset the params without circumventing scope
sub _bind_param_generic_fetch {
  my ($self, $_bind_param_generic_fetch) = @_;
  $self->{'_bind_param_generic_fetch'} = $_bind_param_generic_fetch if $_bind_param_generic_fetch;
  return $self->{_bind_param_generic_fetch};
}

*/


/*
=head2 generic_fetch

  Arg [1]    : (optional) string $constraint
               An SQL query constraint (i.e. part of the WHERE clause)
  Arg [2]    : (optional) Bio::EnsEMBL::AssemblyMapper $mapper
               A mapper object used to remap features
               as they are retrieved from the database
  Arg [3]    : (optional) Bio::EnsEMBL::Slice $slice
               A slice that features should be remapped to
  Example    : $fts = $a->generic_fetch('contig_id in (1234, 1235)');
  Description: Performs a database fetch and returns feature objects in
               contig coordinates.
  Returntype : listref of Bio::EnsEMBL::SeqFeature in contig coordinates
  Exceptions : Thrown if there is an issue with querying the data
  Caller     : BaseFeatureAdaptor, ProxyDnaAlignFeatureAdaptor::generic_fetch
  Status     : Stable

=cut
*/
Vector *BaseAdaptor_genericFetch(BaseAdaptor *ba, char *constraint, AssemblyMapper *mapper, Slice *slice) {
  char qStr[655500];
//  char *qStr = calloc(5655500, sizeof(char));
  qStr[0] = '\0';

  BaseAdaptor_generateSql(ba, constraint, NULL, qStr);

  StatementHandle *sth = ba->prepare((BaseAdaptor *)ba,qStr,strlen(qStr));

  sth->execute(sth);

  Vector *res = ba->objectsFromStatementHandle(ba, sth, mapper, slice);
  sth->finish(sth);

  //free(qStr);
  return res;
}

/*
=head2 generic_count

  Arg [1]    : (optional) string $constraint
               An SQL query constraint (i.e. part of the WHERE clause)
  Example    : $number_feats = $a->generic_count('contig_id in (1234, 1235)');
  Description: Performs a database fetch and returns a count of those features
               found. This is analagous to C<generic_fetch()>
  Returntype : Integer count of the elements.
  Exceptions : Thrown if there is an issue with querying the data

=cut
*/

char *countCols[] = {"count(*)", NULL};

int BaseAdaptor_genericCount(BaseAdaptor *ba, char *constraint) {
  char qStr[655500];
  qStr[0] = '\0';

  BaseAdaptor_generateSql(ba, constraint, countCols, qStr);

  StatementHandle *sth = ba->prepare(ba,qStr,strlen(qStr));
  sth->execute(sth);

  if (sth->numRows(sth) != 1) {
    fprintf(stderr, "genericCount didn't return a row - bye!\n");
    exit(1);
  }
  ResultRow *row = sth->fetchRow(sth);
  int count = row->getLongAt(row, 0);

  return count;
}

void BaseAdaptor_generateSql(BaseAdaptor *ba, char *constraint, char **inputColumns, char *sql) {
  NameTableType *tables = ba->getTables();
  int nTable = 0;
  char tmpStr[1024];
  char extraDefaultWhere[1024];
  int i;
  NameTableType *tabs;
  int needCsTab = 0;
  int needSrTab = 0;

  extraDefaultWhere[0] = '\0';

  // Make a copy of tables so it can be resized if necessary
/* This is way too painful - I'm going to use a couple of flags to indicate that table names need to be
   added to the tableNamesStr
  while ((*tables)[i][0] != NULL) {
    nTable++;
  }

  if ((tabs = calloc((nTable+1), 2*sizeof(char *))) == NULL) {
    fprintf(stderr, "");
    exit(1);
  }
  for (i=0;i<nTable;i++) {
    StrUtil_copyString((*tabs)[i][0], (*tables)[i][0], 0);
    StrUtil_copyString((*tabs)[i][1], (*tables)[i][1], 0);
  }
*/
  
  char srAlias[128];
  char csAlias[128];
  srAlias[0] = '\0';
  csAlias[0] = '\0';

  // Hack for feature types that needs to be restricted to species_id (in
  // coord_system).
  if (BaseAdaptor_isMultiSpecies(BaseAdaptor *ba)
      // For now just the multi species because I don't have adaptors in the Class hierarchy
      // && $self->isa('Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor')
      // && !$self->isa('Bio::EnsEMBL::DBSQL::UnmappedObjectAdaptor') 
     ) {

    // We do a check to see if there is already seq_region
    // and coord_system defined to ensure we get the right
    // alias.  We then do the extra query irrespectively of
    // what has already been specified by the user.

    while ((*tables)[i][0] != NULL) {
      char **t = (*tables)[i];
      if (!strcmp(t[NAME], "seq_region")) {
        strcpy(srAlias, t[SYN]);
      } else if (!strcmp(t[NAME], "coord_system")) {
        strcpy(csAlias, t[SYN]);
      }
    }
    if (!srAlias[0]) {
      strcpy(srAlias, "sr");
      needSrTab = 1;
    }
    if (!csAlias[0]) {
      strcpy(csAlias, "cs");
      needCsTab = 1;
    }

/* Done with needCsTab and needSrTab above
    if ( !exists( $thash{seq_region} ) ) {
      push( @tabs, [ 'seq_region', $sr_alias ] );
    }
    if ( !exists( $thash{coord_system} ) ) {
      push( @tabs, [ 'coord_system', $cs_alias ] );
    }
*/
    sprintf(extraDefaultWhere,"%s.seq_region_id = %s.seq_region_id "
                         "AND %s.coord_system_id = %s.coord_system_id "
                         "AND %s.species_id = "IDFMTSTR,
                         (*tables)[0][1], srAlias, srAlias,
                         csAlias, csAlias, BaseAdaptor_getSpeciesId(ba));

  }

  if (inputColumns == NULL) {
    inputColumns = ba->getColumns();
  }

  char columnsStr[4128];
  columnsStr[0] = '\0';
  i = 0;
  while (inputColumns[i] != NULL) {
    if (columnsStr[0]) {
      strcat(columnsStr, ", ");
    }
    strcat(columnsStr, inputColumns[i]);
    i++;
  }

  // 
  // Construct a left join statement if one was defined, and remove the
  // left-joined table from the table list
  // 
  NameTableType *ljList = ba->leftJoin();

  char leftJoinPrefix[128];
  char leftJoin[4128];
  char tableNamesStr[4128];
  leftJoin[0] = '\0';
  leftJoinPrefix[0] = '\0';
  tableNamesStr[0] = '\0';

  if (ljList != NULL) {
    StringHash *ljHash = StringHash_new(STRINGHASH_SMALL);
    
    i = 0;
    while ((*ljList)[i][0]) {
      StringHash_add(ljHash, (*ljList)[i][0], (*ljList)[i][1]);
      i++;
    }

    i = 0;
    while ((*tables)[i][0]) {
      char **t = (*tables)[i];

      char tAliasKey[1024];
      sprintf(tAliasKey,"%s %s", t[NAME],t[SYN]);
      
      if (StringHash_contains(ljHash, t[NAME]) || StringHash_contains(ljHash, tAliasKey)) {
        char *condition;
        if (StringHash_contains(ljHash, t[NAME])) {
          condition = StringHash_getValue(ljHash, t[NAME]);
        } else {
          condition = StringHash_getValue(ljHash, tAliasKey);
        }
        sprintf(tmpStr,"\n LEFT JOIN %s %s ON %s ) ", t[NAME], t[SYN], condition);
        strcat(leftJoin, tmpStr);
        strcat(leftJoinPrefix,"(");
      } else {
        if (tableNamesStr[0]) {
          sprintf(tmpStr,", %s %s",t[NAME],t[SYN]);
        } else {
          sprintf(tmpStr,"%s %s",t[NAME],t[SYN]);
        }
        strcat(tableNamesStr,tmpStr);
      }
      i++;
    }
    StringHash_free(ljHash, NULL);
  } else {
    i=0;
    while ((*tables)[i][0]) {
      char **t = (*tables)[i];
      if (tableNamesStr[0]) {
        sprintf(tmpStr,", %s %s",t[NAME],t[SYN]);
      } else {
        sprintf(tmpStr,"%s %s",t[NAME],t[SYN]);
      }
      strcat(tableNamesStr,tmpStr);
      i++;
    }
  }

  if (needCsTab) {
    if (tableNamesStr[0]) {
      strcat(tableNamesStr, ",");
    }
    sprintf(tableNamesStr, "%s coord_system %s", tableNamesStr, csAlias);
  }
  if (needSrTab) {
    if (tableNamesStr[0]) {
      strcat(tableNamesStr, ",");
    }
    sprintf(tableNamesStr, "%s seq_region %s", tableNamesStr, srAlias);
  }

  char straightJoin[128];
  straightJoin[0] = '\0';

  // Just a variable in C not a function
  // Only actually set for repeats as far as I can see
  if(ba->straightJoinFlag) {
    strcpy(straightJoin, "STRAIGHT_JOIN");
  }

  sprintf(sql, "SELECT %s %s\n"
               " FROM %s (%s) %s",
               straightJoin, columnsStr, leftJoinPrefix, tableNamesStr, leftJoin);

  char defaultWhereStr[4128];
  strcpy(defaultWhereStr, ba->defaultWhereClause());

  if (extraDefaultWhere[0]) {
    if (defaultWhereStr[0]) {
      sprintf(tmpStr, "\n AND %s", extraDefaultWhere);
      strcat(defaultWhereStr, tmpStr);
    } else {
      strcpy(defaultWhereStr, extraDefaultWhere);
    }
  }

  // append a where clause if it was defined
  if (constraint) {
    sprintf(sql, "%s\n WHERE %s", sql, constraint);
    if (defaultWhereStr[0]) {
      sprintf(sql, "%s AND \n    %s", sql, defaultWhereStr);
    }
  } else if (defaultWhereStr[0]) {
    sprintf(sql, "%s\n WHERE %s", sql, defaultWhereStr);
  }

  //append additional clauses which may have been defined
  if ((ba->finalClause())[0]) strcat(sql, ba->finalClause());
  
  // FOR DEBUG:
  //fprintf(stderr, "SQL:\n%s\n", sql);
  
  return;
}


/*
=head2 fetch_by_dbID

  Arg [1]    : int $id
               The unique database identifier for the feature to be obtained
  Example    : $feat = $adaptor->fetch_by_dbID(1234));
               $feat = $feat->transform('contig');
  Description: Returns the feature created from the database defined by the
               the id $id.  The feature will be returned in its native
               coordinate system.  That is, the coordinate system in which it
               is stored in the database.  In order to convert it to a
               particular coordinate system use the transfer() or transform()
               method.  If the feature is not found in the database then
               undef is returned instead
  Returntype : Bio::EnsEMBL::Feature or undef
  Exceptions : thrown if $id arg is not provided
               does not exist
  Caller     : general
  Status     : Stable

=cut
*/
SeqFeature *BaseAdaptor_fetchByDbID(BaseAdaptor *ba, IDType id) {

  if (BaseAdaptor_hasNoIdCache(ba)) {
    SeqFeature *sf = BaseAdaptor_uncachedFetchByDbID(ba, id);
    return sf;
  }
  fprintf(stderr,"Cached fetching not implemented yet - didn't seem to be used much in perl so I didn't bother - ask Steve\n");
  exit(1);

  // Please compiler
  return NULL;
  //return IDLRUCache_get(BaseAdaptor_getIdCache(ba), id);
}

/*
# The actual implementation moved sideways to allow for uncached access
# otherwise we'd constantly loop
*/

SeqFeature *BaseAdaptor_uncachedFetchByDbID(BaseAdaptor *ba, IDType id) {
  char constraint[1024];

  //construct a constraint like 't1.table1_id = 123'
  NameTableType *tables = ba->getTables();
  char **t = (*tables)[0];
  sprintf(constraint, "%s.%s_id = "IDFMTSTR, t[SYN], t[NAME], id); 

  //Should only be one
  Vector *vec = BaseAdaptor_genericFetch(ba, constraint, NULL, NULL);

  if (Vector_getNumElement(vec) > 1) {
    fprintf(stderr, "Got more than one feature back in fetch ID call - bye!\n");
    exit(1);
  }

  SeqFeature *feat = NULL;
  if (Vector_getNumElement(vec) == 1) {
    feat = Vector_getElementAt(vec, 0);
    Object_incRefCount(feat);
  }
  
// NIY May want to set a free func???
  Vector_free(vec);

  return feat;
}

/*
=head2 fetch_all_by_dbID_list

  Arg [1]    : listref of integers $id_list
               The unique database identifiers for the features to
               be obtained.
  Arg [2]    : optional - Bio::EnsEMBL::Slice to map features onto.
  Example    : @feats = @{$adaptor->fetch_all_by_dbID_list([1234, 2131, 982]))};
  Description: Returns the features created from the database
               defined by the the IDs in contained in the provided
               ID list $id_list.  The features will be returned
               in their native coordinate system.  That is, the
               coordinate system in which they are stored in the
               database.  In order to convert the features to a
               particular coordinate system use the transfer() or
               transform() method.  If none of the features are
               found in the database a reference to an empty list is
               returned.
  Returntype : listref of Bio::EnsEMBL::Features
  Exceptions : thrown if $id arg is not provided
               does not exist
  Caller     : general
  Status     : Stable

=cut
*/


Vector *BaseAdaptor_fetchAllByDbIDList(BaseAdaptor *ba, Vector *idList, Slice *slice) {
  if (BaseAdaptor_hasNoIdCache(ba)) {
    return BaseAdaptor_uncachedFetchAllByDbIDList(ba, idList, slice);
  }
  fprintf(stderr,"Cached fetching not implemented yet - didn't seem to be used much in perl so I didn't bother - ask Steve\n");
  exit(1);

  // Please compiler
  return NULL;
  //return IDCache_getByList(BaseAdaptor_getIDCache(ba), idList, slice);
}

/*
# The actual implmenetation moved sideways to allow for uncached access
# otherwise we'd constantly loop
*/
// Note I didn't implement the stable id fetching uggliness here. I'll probably make a separate method for that
// if necessary
Vector *BaseAdaptor_uncachedFetchAllByDbIDList(BaseAdaptor *ba, Vector *idList, Slice *slice) {
  if ( idList == NULL) {
    fprintf(stderr, "id_list list reference argument is required - bye!");
    exit(1);
  }
  char constraintPref[1024];
  

  if (!Vector_getNumElement(idList)) {
    return Vector_new();
  }

  NameTableType *tables = ba->getTables();
  char **t = (*tables)[0];

  sprintf(constraintPref, "%s.%s_id ", t[SYN], t[NAME] ); 

  // Ensure that we do not exceed MySQL's max_allowed_packet (defaults to
  // 1 MB) splitting large queries into smaller queries of at most 256 KB.
  // Assuming a (generous) average dbID string
  // length of 16, this means 16384 dbIDs in each query.
  int maxSize = 16384;

  // Uniquify the list
  IDHash *idListHash = IDHash_new(IDHASH_MEDIUM);

  int i;
  for (i=0; i<Vector_getNumElement(idList); i++) {
    IDType id = *(IDType *)(Vector_getElementAt(idList, i));
    if (!IDHash_contains(idListHash, id)) {
      IDHash_add(idListHash, id, &trueVal);
    }
  }

  IDType *uniqueIds = IDHash_getKeys(idListHash);
  int nUniqueId = IDHash_getNumValues(idListHash);

  IDHash_free(idListHash, NULL);

  Vector *out = Vector_new();

  int endPoint;
  int lenNum;
  for (i=0; i<nUniqueId; i+=maxSize) {
    char constraint[655500];
    strcpy(constraint, constraintPref);
  
    // Special case for one remaining Id
    if (i == nUniqueId-1) {
      sprintf(constraint, "%s = "IDFMTSTR, constraint, uniqueIds[i]);
    } else {
      char tmpStr[1024];
      int endPoint = sprintf(constraint, "%s IN (", constraint);
      int j;
      for (j=0; j<maxSize && j+i<nUniqueId; j++) {
        if (j!=0) {
          constraint[endPoint++] = ',';
          constraint[endPoint++] = ' ';
        }
        lenNum = sprintf(tmpStr, IDFMTSTR, uniqueIds[i+j]);
        memcpy(&(constraint[endPoint]), tmpStr, lenNum);
        endPoint+=lenNum;
      }
      constraint[endPoint++] = ')';
      constraint[endPoint] = '\0';
    }

    Vector *resChunk = BaseAdaptor_genericFetch(ba, constraint, NULL, slice);

    Vector_append(out, resChunk);

    Vector_free(resChunk);
  }
  free(uniqueIds);

  return out;
}

/*
# might not be a good idea, but for convenience
# shouldnt be called on the BIG tables though
*/

Vector *BaseAdaptor_fetchAll(BaseAdaptor *ba) {
  return BaseAdaptor_genericFetch(ba, NULL, NULL, NULL);
}

/*
=head2 last_insert_id

  Arg [1]     : (optional) $field the name of the field the inserted ID was pushed 
                into
  Arg [2]     : (optional) HashRef used to pass extra attributes through to the 
                DBD driver
  Arg [3]     : (optional) $table the name of the table to use if the adaptor
                does not implement C<_tables()>
  Description : Delegating method which uses DBI to extract the last inserted 
                identifier. If using MySQL we just call the DBI method 
                L<DBI::last_insert_id()> since MySQL ignores any extra
                arguments. See L<DBI> for more information about this 
                delegated method. 
  Example     : my $id = $self->last_insert_id('my_id'); my $other_id = $self->last_insert_id();
  Returntype  : Scalar or undef
  
=cut
*/

/* Note I'm not using this method in C - I'm only supporting mysql, and there's a method in 
   MysqlStatementHandle to get the insert id
sub last_insert_id {
  my ($self, $field, $attributes, $table) = @_;
  my $dbc = $self->dbc();
  my $dbh = $dbc->db_handle();
  my @args;
  if($dbc->driver() eq 'mysql') {
    @args = (undef,undef,undef,undef);
  }
  else {
    if(!$table) {
      ($table) = $self->_tables();
    }
    @args = (undef, $dbc->dbname(), $table->[0], $field);
  }
  $attributes ||= {};
  return $dbh->last_insert_id(@args, $attributes);
}
*/

// This ID caching not implemented so remove all the complex stuff around this
int BaseAdaptor_hasNoIdCache(BaseAdaptor *ba)  {
  return 1;
}

/*
=head2 _id_cache

  Description : Used to return an instance of a support BaseCache module
                which can be used to speed up object access. The method
                also respects the DBAdaptor's no_cache() flag and will
                return undef in those situations
  Example     : my $cache = $self->_id_cache();
  Returntype  : Bio:EnsEMBL::DBSQL::Support::BaseCache
  
=cut
*/

/* This ID caching barely used in perl as far as I can see so not implemented
sub _id_cache {
  my ($self) = @_;
  return if $self->db()->no_cache() && !$self->ignore_cache_override;
  if(! exists $self->{_id_cache}) {
    $self->{_id_cache} = $self->_build_id_cache();
  }
  return $self->{_id_cache};
}
*/


/*
=head2 _no_id_cache

  Description : Flags if the ID based caching is active or not. This could be
                due to the adaptor not wanting to cache or because of
                a global no_cache() flag on the DBAdaptor instance
  Returntype  : Boolean
  
=cut
*/

/* This ID caching barely used in perl as far as I can see so not implemented
sub _no_id_cache {
  my ($self) = @_;
  return 1 if ! $self->_id_cache();
  return 0;
}
*/

/*
=head2 ignore_cache_override

    Description : Method to interfere with no_cache directive from Registry on
                  a per adaptor basis. This method should be called after new()
                  in order to trigger the _build_id_cache at first query.                  
    Example     : $adaptor->ignore_cache_override(1);              
    Returntype  : Boolean
=cut
*/

/* This ID caching barely used in perl as far as I can see so not implemented - maybe its done somewhere evil like the registry config, but who cares  about that
sub ignore_cache_override {
    my $self = shift;
    $self->{'_override'} = shift if(@_);
    unless (defined($self->{'_override'})) {return}
    return $self->{'_override'}; 
}
*/

/*
#_tables
#
#  Args       : none
#  Example    : $tablename = $self->_table_name()
#  Description: ABSTRACT PROTECTED
#               Subclasses are responsible for implementing this
#               method.  It should list of [tablename, alias] pairs.
#               Additionally the primary table (with the dbID,
#               analysis_id, and score) should be the first table in
#               the list. e.g:
#               ( ['repeat_feature',   'rf'],
#                 ['repeat_consensus', 'rc']);
#               used to obtain features.  
#  Returntype : list of [tablename, alias] pairs
#  Exceptions : thrown if not implemented by subclass
#  Caller     : BaseFeatureAdaptor::generic_fetch
#
*/

NameTableType *BaseAdaptor_getTables(void) {
  fprintf(stderr,"ERROR: Abstract method getTables not defined by implementing subclass\n");
  exit(1);

// Please the compiler
  return NULL;
}

/*
#_columns
#
#  Args       : none
#  Example    : $tablename = $self->_columns()
#  Description: ABSTRACT PROTECTED
#               Subclasses are responsible for implementing this
#               method.  It should return a list of columns to be
#               used for feature creation.
#  Returntype : list of strings
#  Exceptions : thrown if not implemented by subclass
#  Caller     : BaseFeatureAdaptor::generic_fetch
#
*/

char **BaseAdaptor_getColumns(void) {
  fprintf(stderr,"ERROR: Abstract method getColumns not defined by implementing subclass\n");
  exit(1);

// Please the compiler
  return NULL;
}

/*
# _default_where_clause
#
#  Arg [1]    : none
#  Example    : none
#  Description: May be overridden to provide an additional where
#               constraint to the SQL query which is generated to
#               fetch feature records.  This constraint is always
#               appended to the end of the generated where clause
#  Returntype : string
#  Exceptions : none
#  Caller     : generic_fetch
#
*/

char *BaseAdaptor_defaultWhereClause(void) {
  return "";
}


/*
# _left_join

#  Arg [1]    : none
#  Example    : none
#  Description: Can be overridden by a subclass to specify any left
#               joins which should occur.  The table name specigfied
#               in the join must still be present in the return
#               values of.
#  Returntype : a {'tablename' => 'join condition'} pair
#  Exceptions : none
#  Caller     : general
#
*/

// Abuse of NameTableType (basically a pair of strings)
NameTableType *BaseAdaptor_leftJoin(void) {
  return NULL;
}


/*
#_final_clause

#  Arg [1]    : none
#  Example    : none
#  Description: May be overriden to provide an additional clause
#               to the end of the SQL query used to fetch feature
#               records.  This is useful to add a required ORDER BY
#               clause to the query for example.
#  Returntype : string
#  Exceptions : none
#  Caller     : generic_fetch
*/
char *BaseAdaptor_finalClause(void) {
  return "";
}

/*
#_objs_from_sth

#  Arg [1]    : DBI::row_hashref $hashref containing key-value pairs 
#               for each of the columns specified by the _columns method
#  Example    : my @feats = $self->_obj_from_hashref
#  Description: ABSTRACT PROTECTED
#               The subclass is responsible for implementing this
#               method.  It should take in a DBI row hash reference
#               and return a list of created features in contig
#               coordinates.
#  Returntype : list of Bio::EnsEMBL::*Features in contig coordinates
#  Exceptions : thrown if not implemented by subclass
#  Caller     : BaseFeatureAdaptor::generic_fetch
*/
Vector *BaseAdaptor_objectsFromStatementHandle(BaseAdaptor *bfa, StatementHandle *sth,
                                               AssemblyMapper *mapper, Slice *slice) {
  fprintf(stderr,"ERROR: Abstract method objectsFromStatementHandle not defined by implementing subclass\n");
  exit(1);

// Please the compiler
  return NULL;
}

/*
=head2 store

  Arg [1]    : list of Bio::EnsEMBL::SeqFeature
  Example    : $adaptor->store(@feats);
  Description: ABSTRACT  Subclasses are responsible for implementing this
               method.  It should take a list of features and store them in
               the database.
  Returntype : none
  Exceptions : thrown method is not implemented by subclass
  Caller     : general
  Status     : At Risk
             : throws if called.

=cut
*/

int BaseAdaptor_store(BaseAdaptor *ba, Vector *sfs) {

  fprintf(stderr, "Abstract method store not defined by implementing subclass\n");
  exit(1);
}

/*
#_build_id_cache

#  Example    : my $id_cache = $self->_build_id_cache
#  Description: ABSTRACT PROTECTED
#               The subclass is responsible for returning an instance
#               of the Bio::EnsEMBL::DBSQL::Support::BaseCache
#               which can be used to speed up ID based fetch operations
#  Returntype : Instance of Bio::EnsEMBL::DBSQL::Support::BaseCache
#  Exceptions : Could be thrown by the implementing sub-class 
#  Caller     : BaseAdaptor::_id_cache
*/
/* This ID caching barely used in perl as far as I can see so not implemented - maybe its done somewhere evil like the registry config, but who cares  about that
sub _build_id_cache {
  return;
}
*/

/*
sub dump_data {
  my $self = shift;
  my $data = shift;

  my $dumper = Data::Dumper->new([$data]);
  $dumper->Indent(0);
  $dumper->Terse(1);
   my $dump = $dumper->Dump();
# $dump =~ s/'/\\'/g; 
 # $dump =~ s/^\$VAR1 = //;
  return $dump;
}

sub get_dumped_data {
    my $self = shift;
    my $data = shift;

    $data =~ s/\n|\r|\f|\\//g;
    return eval ($data);
}


*/

