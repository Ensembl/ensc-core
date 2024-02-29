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

#define __STORABLE_MAIN__
#include "Storable.h"
#undef __STORABLE_MAIN__
#include "DBAdaptor.h"
#include "DBConnection.h"
#include "BaseAdaptor.h"

// NOTE: Currently I've done Storable as a "has a" member of the various storable object classes - that may have been an error, but for now I'm going to live with it

/*
=head2 is_stored

  Arg [1]    : Bio::EnsEMBL::DBSQL::DBConnection
             : or Bio::EnsEMBL::DBSQL::DBAdaptor
  Example    : do_something if($object->is_stored($db));
  Description: Returns true if this object is stored in the provided database.
               This works under the assumption that if the adaptor and dbID are
               set and the database of the adaptor shares the port, dbname and
               hostname with the provided database, this object is stored in
               that database.
  Returntype : 1 or 0
  Exceptions : throw if dbID is set but adaptor is not
               throw if adaptor is set but dbID is not
               throw if incorrect argument is passed
  Caller     : store methods
  Status     : Stable

=cut
*/
//New
static int messageOnlyOnce = 1;
int Storable_isStored(Storable *storable, DBAdaptor *db) {

/* Perl seems to allow either DBAdaptor or DBConnection, not sure why.
   I'm going to restrict to DBAdaptor for now
  if($db and $db->isa('Bio::EnsEMBL::DBSQL::DBAdaptor')) {
    $db = $db->dbc();
  }
  if(!$db || !ref($db) || !$db->isa('Bio::EnsEMBL::DBSQL::DBConnection')) {
    throw('db argument must be a Bio::EnsEMBL::DBSQL::DBConnection');
  }
*/

  BaseAdaptor *adaptor = Storable_getAdaptor(storable);
  IDType dbID = Storable_getDbID(storable);

  if (dbID && (adaptor == NULL)) {
    if (messageOnlyOnce) {
      fprintf(stderr, "Storable object has a dbID but not an adaptor.\n" 
                      "Storable objects must have neither OR both.\n");
      messageOnlyOnce = 0;
    }
    return 0;
  }

  if ((adaptor != NULL) && ! dbID) {
    if (messageOnlyOnce) {
      fprintf(stderr, "Storable object has an adaptor but not a dbID.\n"
                       "Storable objects must have neither OR both.\n");
      messageOnlyOnce = 0;
    }
    return 0;
  }

  if (adaptor == NULL && dbID == 0) {
    return 0;
  }

  DBConnection *curDbc = adaptor->dba->dbc;
  DBConnection *dbc    = db->dbc;

  // 
  // Databases are the same if they share the same port, host and dbname
  //
  if (   DBConnection_getPort(dbc) == DBConnection_getPort(curDbc) &&
       ! EcoString_strcmp(DBConnection_getHost(dbc), DBConnection_getHost(curDbc)) &&
       ! EcoString_strcmp(DBConnection_getDbName(dbc), DBConnection_getDbName(curDbc))) {
    return 1;
  }

  return 0;
}


