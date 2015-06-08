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

#include "SliceAdaptor.h"

#include "DBAdaptor.h"
#include "BaseAdaptor.h"
#include "CoordSystemAdaptor.h"
#include "MysqlUtil.h"
#include "SeqRegionCacheEntry.h"
#include "Vector.h"
#include "StrUtil.h"
#include "Slice.h"
#include "CoordSystem.h"
#include "ProjectionSegment.h"
#include "StringHash.h"
#include "AssemblyMapperAdaptor.h"
#include "ExonAdaptor.h"
#include "GeneAdaptor.h"
#include "TranscriptAdaptor.h"


/*
=head1 DESCRIPTION

This module is responsible for fetching Slices representing genomic
regions from a database.  A Details on how slices can be used are in the
Bio::EnsEMBL::Slice module.

=head1 METHODS

=cut
*/

// Private struct type for exception cache
typedef struct ExceptionCacheDataStruct {
  IDType seqRegionId;
  long seqRegionStart;
  long seqRegionEnd;
  char *excType;
  IDType excSeqRegionId;
  long excSeqRegionStart;
  long excSeqRegionEnd;
} ExceptionCacheData;

ExceptionCacheData *ExceptionCacheData_new(IDType seqRegionId, long seqRegionStart, long seqRegionEnd, char *excType,
                                           IDType excSeqRegionId, long excSeqRegionStart, long excSeqRegionEnd) {
  ExceptionCacheData *ecd;
  if ((ecd = (ExceptionCacheData *)calloc(1,sizeof(ExceptionCacheData))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for ExceptionCacheData\n");
    return NULL;
  }

  ecd->seqRegionId       = seqRegionId;
  ecd->seqRegionStart    = seqRegionStart;
  ecd->seqRegionEnd      = seqRegionEnd;
  StrUtil_copyString(&ecd->excType, excType, 0);
  ecd->excSeqRegionId    = excSeqRegionId;
  ecd->excSeqRegionStart = excSeqRegionStart;
  ecd->excSeqRegionEnd   = excSeqRegionEnd;
  
  return ecd;
}

void ExceptionCacheData_free(ExceptionCacheData *ecd) {
  free(ecd->excType);
  free(ecd);
}


SliceAdaptor *SliceAdaptor_new(DBAdaptor *dba) {
  SliceAdaptor *sa;

  if ((sa = (SliceAdaptor *)calloc(1,sizeof(SliceAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for SliceAdaptor\n");
    return NULL;
  }
  BaseAdaptor_init((BaseAdaptor *)sa, dba, SLICE_ADAPTOR);

  sa->prepare = SliceAdaptor_prepare;

  // use a shared cache (for this database) that contains info about
  // seq regions
  sa->srNameCache = DBAdaptor_getSeqRegionNameCache(sa->dba);
  sa->srIdCache   = DBAdaptor_getSeqRegionIdCache(sa->dba);

// Not implementing LRGs
//  $self->{'lrg_region_test'} = undef;
//  my $meta_container = $self->db->get_MetaContainer();
//  my @values = $meta_container->list_value_by_key("LRG");
//  if(scalar(@values) and $values[0]->[0]){
//    $self->{'lrg_region_test'} = $values[0]->[0];
//  }

  return sa;
}

IDType SliceAdaptor_getSpeciesID(SliceAdaptor *sa) {
  return 1;
}

/*
=head2 fetch_by_region

  Arg [1]    : string $coord_system_name (optional)
               The name of the coordinate system of the slice to be created
               This may be a name of an actual coordinate system or an alias
               to a coordinate system.  Valid aliases are 'seqlevel' or
               'toplevel'.
  Arg [2]    : string $seq_region_name
               The name of the sequence region that the slice will be
               created on.
  Arg [3]    : int $start (optional, default = 1)
               The start of the slice on the sequence region
  Arg [4]    : int $end (optional, default = seq_region length)
               The end of the slice on the sequence region
  Arg [5]    : int $strand (optional, default = 1)
               The orientation of the slice on the sequence region
  Arg [6]    : string $version (optional, default = default version)
               The version of the coordinate system to use (e.g. NCBI33)
  Arg [7]    : boolean $no_fuzz (optional, default = undef (false))
               If true (non-zero), do not use "fuzzy matching" (see below).
  Example    : $slice = $slice_adaptor->fetch_by_region('chromosome', 'X');
               $slice = $slice_adaptor->fetch_by_region('clone', 'AC008066.4');
  Description: Retrieves a slice on the requested region.  At a minimum the
               name the name of the seq_region to fetch must be provided.

               If no coordinate system name is provided than a slice on the
               highest ranked coordinate system with a matching
               seq_region_name will be returned.  If a version but no
               coordinate system name is provided, the same behaviour will
               apply, but only coordinate systems of the appropriate version
               are considered.  The same applies if the 'toplevel' coordinate
               system is specified, however in this case the version is
               ignored.  The coordinate system should always be specified if
               it is known, since this is unambiguous and faster.

               Some fuzzy matching is performed if no exact match for
               the provided name is found.  This allows clones to be
               fetched even when their version is not known.  For
               example fetch_by_region('clone', 'AC008066') will
               retrieve the sequence_region with name 'AC008066.4'.

               The fuzzy matching can be turned off by setting the
               $no_fuzz argument to a true value.

               If the requested seq_region is not found in the database undef
               is returned.

  Returntype : Bio::EnsEMBL::Slice or undef
  Exceptions : throw if no seq_region_name is provided
               throw if invalid coord_system_name is provided
               throw if start > end is provided
  Caller     : general
  Status     : Stable

=cut
*/


//
// ARNE: This subroutine needs simplification!! 
//
// STEVE: I guess no one's listening!!!!!
//

Slice *SliceAdaptor_fetchByRegion(SliceAdaptor *sa, char *coordSystemName, char *inputSeqRegionName, long start, long end, int strand, char *version, int noFuzz) {

  if (start == POS_UNDEF)     start  = 1;
  if (strand == STRAND_UNDEF) strand = 1;

  if ( !inputSeqRegionName)  {
    fprintf(stderr,"seq_region_name argument is required in fetchByRegion\n");
    exit(1);
  }

  // Note make a copy of the seq region name to make memory management easier
  char *seqRegionName = StrUtil_copyString(&seqRegionName, inputSeqRegionName, 0);

  CoordSystem *cs = NULL;
  CoordSystemAdaptor *csa = DBAdaptor_getCoordSystemAdaptor(sa->dba);

  if (coordSystemName != NULL ) {
    cs = CoordSystemAdaptor_fetchByName(csa, coordSystemName, version );

    // REMOVE THESE THREE LINES WHEN STICKLEBACK DB IS FIXED!
    // Anne/ap5 (2007-10-09):
    // The problem was that the stickleback genebuild called the
    // chromosomes 'groups', which meant they weren't being picked out by
    // the karyotype drawing code.  Apparently they are usually called
    // 'groups' in the stickleback community, even though they really are
    // chromosomes!

    if (cs == NULL && !strcmp(coordSystemName, "chromosome") ) {
      cs = CoordSystemAdaptor_fetchByName(csa, "group", version );
    }

    if (cs == NULL ) {
      fprintf(stderr, "Unknown coordinate system:\nname='%s' version='%s'\n",
                      coordSystemName, version==NULL ? "" : version );
    } else {
      // fetching by toplevel is same as fetching w/o name or version
      if ( CoordSystem_getIsTopLevel(cs) ) {
        cs      = NULL;
        version = NULL;
      }
    }
  } 

  char constraint[1024];
  char sql[2048];
  char tmpStr[1024];
  char key[1024];

  // Empty key
  key[0] = '\0';

  if (cs) {
    sprintf(sql,"SELECT sr.name, sr.seq_region_id, sr.length, "IDFMTSTR" "
                         "FROM seq_region sr ",
                         CoordSystem_getDbID(cs));

    sprintf(constraint, "AND sr.coord_system_id = "IDFMTSTR, CoordSystem_getDbID(cs));

    sprintf(key,"%s:"IDFMTSTR, seqRegionName, CoordSystem_getDbID(cs));
  } else {

    sprintf(sql, "SELECT sr.name, sr.seq_region_id, sr.length, cs.coord_system_id "
                   "FROM seq_region sr, coord_system cs ");

    sprintf(constraint, "AND sr.coord_system_id = cs.coord_system_id "
                       "AND cs.species_id = "IDFMTSTR, SliceAdaptor_getSpeciesID(sa));

    if (version != NULL) {
      sprintf(tmpStr," AND cs.version = '%s'",version);
      strcat(constraint,tmpStr);
    }

    strcat(constraint, " ORDER BY cs.rank ASC");
  }

  // check the cache so we only go to the db if necessary
  long length;

  SeqRegionCacheEntry *cacheData = NULL;

  if (strlen(key) && StringHash_contains(sa->srNameCache, key)) {
    cacheData = StringHash_getValue(sa->srNameCache, key);
  }

  if ( cacheData != NULL ) {
    length = cacheData->regionLength;
  } else {
    char qStr[4196];    

    sprintf(qStr,"%s WHERE sr.name = '%s' %s", sql, seqRegionName, constraint);

    StatementHandle *sth = sa->prepare((BaseAdaptor *)sa,qStr,strlen(qStr));

    sth->execute(sth);

    if ( sth->numRows(sth) == 0 ) {
      sth->finish(sth);

      // try synonyms
      sprintf(qStr,"select s.name, cs.name, cs.version "
                   "from seq_region s join seq_region_synonym ss "
                   "using (seq_region_id) join coord_system cs using (coord_system_id) "
                   "where ss.synonym = '%s' and cs.species_id ="IDFMTSTR, seqRegionName, SliceAdaptor_getSpeciesID(sa));
      
      StatementHandle *synSqlSth = sa->prepare((BaseAdaptor *)sa,qStr,strlen(qStr));

      synSqlSth->execute(synSqlSth);

      ResultRow *row;
      if ((row = synSqlSth->fetchRow(synSqlSth))) {
        char *newName        = row->getStringAt(row, 0);
        char *newCoordSystem = row->getStringAt(row, 1);
        char *newVersion     = row->getStringAt(row, 2);

        Slice *retSlice;
        
        // Need to change logic slightly to ease memory management (don't want to do synSqlSth->finish before this 
        // conditional because it need the row results which would have been freed if I'd done the finish)
        if (cs == NULL || !strcmp(CoordSystem_getName(cs), newCoordSystem)) {
          retSlice = SliceAdaptor_fetchByRegion(sa, newCoordSystem, newName, start, end, strand, newVersion, noFuzz);
        } else if (strcmp(CoordSystem_getName(cs), newCoordSystem)) {
          fprintf(stderr,"Searched for a known feature on coordinate system: "IDFMTSTR" but found it on: %s\n"
                         "No result returned, consider searching without coordinate system or use toplevel.\n", 
                         CoordSystem_getDbID(cs), newCoordSystem);
          retSlice = NULL;
        }
        synSqlSth->finish(synSqlSth);

// NIY:: Any other memory tidying

        return retSlice;
      }
      synSqlSth->finish(synSqlSth);

      if (noFuzz) {
        return NULL;
      }

      // Do fuzzy matching, assuming that we are just missing a version
      // on the end of the seq_region name.

      // NOTE: Started on this statement (%s.%%%% added from bind_param below it)
      // Note I need 4 '%' characters to end up with one because it goes through two sprintfs
      // which each need an escaped % character
      fprintf(stderr,"Before fuzzy\n");
      sprintf(qStr, "%s WHERE sr.name LIKE '%s.%%%%' %s", sql, seqRegionName, constraint);

      sth = sa->prepare((BaseAdaptor *)sa,qStr,strlen(qStr));

      sth->execute(sth);

      int prefixLen = strlen(seqRegionName) + 1;
      long highVer  = 0L;
      CoordSystem *highCs = cs;

      // Find the fuzzy-matched seq_region with the highest postfix
      // (which ought to be a version).

      int i = 0;
      int hadVer = 0;

      while ((row = sth->fetchRow(sth))) {
        char *tmpName  = row->getStringAt(row,0);
        IDType id      = row->getLongLongAt(row,1);
        long tmpLength = row->getLongAt(row,2);
        IDType csId    = row->getLongLongAt(row,3);

        CoordSystem *tmpCs = cs != NULL ? cs : CoordSystemAdaptor_fetchByDbID(csa, csId);

        // cache values for future reference
        DBAdaptor_addToSrCaches(sa->dba, id, tmpName, csId, tmpLength);

        // skip versions which are non-numeric and apparently not
        // versions
        long tmpVer;

        if (!StrUtil_isLongInteger(&tmpVer, &(tmpName[prefixLen]))) {
          continue;
        }

        fprintf(stderr,"Here tmpVer = %ld\n",tmpVer);
        // take version with highest num, if two versions match take one
        // with highest ranked coord system (lowest num)
        if ( !hadVer ||
             tmpVer > highVer ||
             (tmpVer == highVer && CoordSystem_getRank(tmpCs) < CoordSystem_getRank(highCs))) {
// NIY Do we need to string copy??
          free(seqRegionName);
          seqRegionName = StrUtil_copyString(&seqRegionName, tmpName, 0);
          length        = tmpLength;
          highVer       = tmpVer;
          highCs        = tmpCs;
          hadVer        = 1;
        }
        i++;
      }

      // warn if fuzzy matching found more than one result
      if (i>1) {
        fprintf(stderr, "Fuzzy matching of seq_region_name returned more than one result.\n"
                        "You might want to check whether the returned seq_region\n"
                        "(%s:%s) is the one you intended to fetch.\n",
                        CoordSystem_getName(highCs), seqRegionName);
      }

      cs = highCs;

      fprintf(stderr,"After fuzzy hadVer = %d\n",hadVer);
      sth->finish(sth);

      // return if we did not find any appropriate match:

      if (!hadVer) return NULL;

    } else {
      ResultRow *row = sth->fetchRow(sth);

      free(seqRegionName);
      seqRegionName = StrUtil_copyString(&seqRegionName, row->getStringAt(row,0),0); // Note not local to block
      IDType id     = row->getLongLongAt(row,1);
      length        = row->getLongAt(row,2); // Note not local to block
      IDType csId   = row->getLongLongAt(row,3);

      // cache to speed up for future queries
      DBAdaptor_addToSrCaches(sa->dba, id, seqRegionName, csId, length);

      cs = CoordSystemAdaptor_fetchByDbID(csa, csId );
      sth->finish(sth);
    }
  }

  // end set to POS_UNDEF equivalent to not defined - may come back to haunt me!!
  if ( end == POS_UNDEF) {
    end = length;
  }
  //if (!defined($end) ) { $end = $length }

  //If this was given then check if we've got a circular seq region otherwise
  //let it fall through to the normal Slice method
/* No circular stuff
  if ( $end + 1 < $start ) {
    my $cs_id = $cs->dbID();
    my $seq_region_id = $self->{'sr_name_cache'}->{"$seq_region_name:$cs_id"}->[0];
    if($self->is_circular($seq_region_id)) {
      my $new_sl =
        Bio::EnsEMBL::CircularSlice->new(
                                     -COORD_SYSTEM    => $cs,
                                     -SEQ_REGION_NAME => $seq_region_name,
                                     -SEQ_REGION_LENGTH => $length,
                                     -START             => $start,
                                     -END               => $end,
                                     -STRAND            => 1,
                                     -ADAPTOR           => $self );
  
      return $new_sl;
    }
  }
*/

/* No LRG Stuff
  if ( defined( $self->{'lrg_region_test'} )
       and substr( $cs->name, 0, 3 ) eq $self->{'lrg_region_test'} )
*/
  if (0) {
/*
    return
      Bio::EnsEMBL::LRGSlice->new( -COORD_SYSTEM    => $cs,
                                   -SEQ_REGION_NAME => $seq_region_name,
                                   -SEQ_REGION_LENGTH => $length,
                                   -START             => $start,
                                   -END               => $end,
                                   -STRAND            => $strand,
                                   -ADAPTOR           => $self );
*/
  } else {
    Slice *slice = Slice_new(seqRegionName, start, end, strand, length, cs, sa);

    free(seqRegionName);

    return slice;
  }
}

/*
=head2 fetch_by_toplevel_location

  Arg [1]     : string $location
                Ensembl formatted location. Can be a format like 
                C<name:start-end>, C<name:start..end>, C<name:start:end>, 
                C<name:start>, C<name>. We can also support strand 
                specification as a +/- or 1/-1. 
                
                Location names must be separated by a C<:>. All others can be
                separated by C<..>, C<:> or C<->.
  Arg[2]      : boolean $no_warnings
                Suppress warnings from this method
  Arg[3]      : boolean $no_fuzz
                Stop fuzzy matching of sequence regions from occuring
  Example     : my $slice = $sa->fetch_by_toplevel_location('X:1-10000')
                my $slice = $sa->fetch_by_toplevel_location('X:1-10000:-1')
  Description : Converts an Ensembl location/region into the sequence region
                name, start and end and passes them onto C<fetch_by_region()>. 
                The code assumes that the required slice is on the top level
                coordinate system. The code assumes that location formatting
                is not perfect and will perform basic cleanup before parsing.
  Returntype  : Bio::EnsEMBL::Slice
  Exceptions  : If $location is false otherwise see C<fetch_by_location()>
                or C<fetch_by_region()>
  Caller      : General
  Status      : Beta

=cut
*/

Slice *SliceAdaptor_fetchByTopLevelLocation(SliceAdaptor *sa, char *location, int noWarnings, int noFuzz) {

  return SliceAdaptor_fetchByLocation(sa, location, "toplevel", NULL, noWarnings, noFuzz);
}

/*
=head2 fetch_by_location

  Arg [1]     : string $location
                Ensembl formatted location. Can be a format like 
                C<name:start-end>, C<name:start..end>, C<name:start:end>, 
                C<name:start>, C<name>. We can also support strand 
                specification as a +/- or 1/-1. 
                
                Location names must be separated by a C<:>. All others can be
                separated by C<..>, C<:> or C<->.
  Arg[2]      : String $coord_system_name
                The coordinate system to retrieve
  Arg[3]      : String $coord_system_version
                Optional parameter. Version of the coordinate system to fetch
  Arg[4]      : boolean $no_warnings
                Suppress warnings from this method
  Arg[5]      : boolean $no_fuzz
                Stop fuzzy matching of sequence regions from occuring
  Example     : my $slice = $sa->fetch_by_toplevel_location('X:1-10000')
                my $slice = $sa->fetch_by_toplevel_location('X:1-10000:-1')
  Description : Converts an Ensembl location/region into the sequence region
                name, start and end and passes them onto C<fetch_by_region()>. 
                The code assumes that the required slice is on the top level
                coordinate system. The code assumes that location formatting
                is not perfect and will perform basic cleanup before parsing.
  Returntype  : Bio::EnsEMBL::Slice
  Exceptions  : If $location or coordinate system is false otherwise 
                see C<fetch_by_region()>
  Caller      : General
  Status      : Beta

=cut
*/

Slice *SliceAdaptor_fetchByLocation(SliceAdaptor *sa, char *location, char *coordSystemName, char *coordSystemVersion, int noWarnings, int noFuzz) {
  char *seqRegionName;
  long start;
  long end;
  int strand;
  
  if (coordSystemName == NULL) {
    fprintf(stderr,"No coordinate system name specified\n");
    exit(1);
  }
  
  SliceAdaptor_parseLocationToValues(sa, location, noWarnings, 0 /*no errors */, &seqRegionName, &start, &end, &strand);

  // Perl just did a return which seemed a bit odd
  if (seqRegionName == NULL) {
    fprintf(stderr,"No seqRegionName specified - bye\n");
    exit(1);
    //return;
  }
    
  if (start!=POS_UNDEF && end!=POS_UNDEF && start > end) {
    fprintf(stderr,"Cannot request a slice whose start is greater than its end. Start: %ld. End: %ld\n", start, end);
    exit(1);
  }
  
  Slice *slice = SliceAdaptor_fetchByRegion(sa, coordSystemName, seqRegionName, start, end, strand, coordSystemVersion, noFuzz);

  if (slice == NULL) return NULL;
  
  long srl  = Slice_getSeqRegionLength(slice);
  char *name = Slice_getSeqRegionName(slice);

  if (start!=POS_UNDEF && start > srl) {
    fprintf(stderr, "Cannot request a slice whose start (%ld) is greater than %ld for %s.\n", start, srl, name);
    exit(1);
  }
  if (end!=POS_UNDEF && end > srl) {
    if (!noWarnings) {
      fprintf(stderr,"Requested end (%ld) is greater than %ld for %s. Resetting to %ld\n", end, srl, name, srl);
    }
    Slice_setEnd(slice, srl);
  }
  
  return slice;
}

/*
=head2 parse_location_to_values

  Arg [1]     : string $location
                Ensembl formatted location. Can be a format like 
                C<name:start-end>, C<name:start..end>, C<name:start:end>, 
                C<name:start>, C<name>. We can also support strand 
                specification as a +/- or 1/-1. 
                
                Location names must be separated by a C<:>. All others can be
                separated by C<..>, C<:> or C<->.
  Arg[2]      : boolean $no_warnings
                Suppress warnings from this method
  Arg[3]      : boolean $no_errors
                Supress errors being thrown from this method
  Example                        : my ($name, $start, $end, $strand) = $sa->parse_location_to_values('X:1..100:1);
  Description	: Takes in an Ensembl location String and returns the parsed
                values
  Returntype 	: List. Contains name, start, end and strand 

=cut
*/
/* Formats I'm going to allow
    1:1

    1:1-10
    1:1-10:1
    1:1-10:-1

    1:1:10
    1:1:10:1
    1:1:10:-1
    1:1:10:+
    1:1:10:-

    1:1..10
    1:1..10:1
    1:1..10:-1

Numbers can have commas eg 1,000 or underscores 1_000
   eg 1:1..1,000 is OK

There can be spaces anywhere but they will be ignored - they have no meaning
   1 : 1 .. 10 is fine but completely equivalent to 1:1-10
   1 : 1 .. 1 0  is also fine BUT is also equivalent to 1:1-10

I'm NOT currently allowing a ':' separated range which is just one position and a strand, because that's difficult to distinguish from a 
start end range
*/

void SliceAdaptor_parseLocationToValues(SliceAdaptor *sa, char *location, int noWarnings, int noErrors, 
                                        char **retSeqRegionName, long *retStart, long *retEnd, int *retStrand ) {
  char *tmpLoc;
  char *seqRegionName;
  long start = POS_UNDEF;
  long end   = POS_UNDEF;
  int strand = 1;
  int i;

  
  if (location == NULL) {
    fprintf(stderr,"You must specify a location\n");
    exit(1);
  }   
  
  StrUtil_copyString(&tmpLoc, location, 0);

  StrUtil_rmspace(tmpLoc);


  char **tokens = NULL;
  int nTok;

  StrUtil_tokenizeByDelim(&tokens, &nTok, tmpLoc, ":");

  if (nTok > 4 || nTok < 1) {
    fprintf(stderr, "Unable to parse location string %s\n", location);
    exit(1);
  }

  seqRegionName = tokens[0];

  long num;
  if (nTok >= 2) {
    // Is it a .. or - type range

    char *dotChP = NULL;
    char *dashChP = NULL;

    // Allow for a negative start position eg 1:-10..10:1
    int startInd = 0;
    if (tokens[1][0] == '-') startInd++;

    if ((dotChP = strstr(&(tokens[1][startInd]),"..")) || (dashChP = strstr(&(tokens[1][startInd]),"-"))) {

      // Can't be a ':' separated range so 4 tokens is an error
      if (nTok == 4) {
        fprintf(stderr, "Unable to parse location string %s (4 tok)\n", location);
        exit(1);
      }

      // can't use dots and dashs together as a range separator!
      if (dotChP && dashChP) {
        fprintf(stderr, "Unable to parse location string %s (dot and dash)\n", location);
        exit(1);
      }
        
      // should be two numbers around the ".." or "-"
      // Make them in to elements in tokens so can parse in same way for all
      char *chP;
      int sepLen;
      if (dotChP) { 
        sepLen = 2; 
        chP = dotChP;
      } else {
        sepLen = 1; 
        chP = dashChP;
      }

      *chP = '\0'; // make tokens[1] end at start of delim;
      char *newTok;
      StrUtil_copyString(&newTok, chP+sepLen, 0);

      if ((tokens = realloc(tokens, (nTok+1)*sizeof(char *))) == NULL) {
        fprintf(stderr, "Failed allocating space for tokens\n");
        exit(1);
      }
      for (i=nTok-1; i>1; i--) {
        tokens[i+1] = tokens[i];
      }
      tokens[2] = newTok;
      nTok++;
    }

    for (i=1; i<nTok; i++) {
      StrUtil_rmChar(tokens[i],',');
      StrUtil_rmChar(tokens[i],'_');
    }

    // if tokens[1] isn't an integer error
    if ( ! StrUtil_isLongInteger(&num, tokens[1]))  {
      fprintf(stderr, "Unable to parse location string %s (tok1 not long)\n", location);
      exit(1);
    } else {
      start = num;
    }

    if (nTok == 2) {
      end = start;
    // If nTok > 2) then may be an end position or maybe just a strand
    } else {
      int lastInd = nTok-1;
      int haveSetStrand = 0;
      if (strlen(tokens[lastInd]) == 1 && (tokens[lastInd][0] == '+' || tokens[lastInd][0] == '-')) {
        haveSetStrand = 1;
        strand = (tokens[lastInd][0] == '+' ? 1 : -1);
      // is it an integer
      } else if (StrUtil_isLongInteger(&num, tokens[lastInd])) {
        // if its -1 or 1 its a strand
        if (num ==1 || num == -1) {
          end  = start;
          strand = num;
          haveSetStrand = 1;
        // else its the end pos
        } else {
          end = num;
        }
      } else {
        fprintf(stderr, "Unable to parse location string %s (last tok)\n", location);
        exit(1);
      }

      if (nTok == 4) {
        if (haveSetStrand) {
          //end is tokens[2]
          if (StrUtil_isLongInteger(&num, tokens[2])) {
            end = num;
          } else {
            fprintf(stderr, "Unable to parse location string %s (nTok==4 if)\n", location);
            exit(1);
          }
        // else its an error because we have four tokens but the last one isn't a strand
        } else {
          fprintf(stderr, "Unable to parse location string %s (4 tok, no strand)\n", location);
          exit(1);
        }
      }
    }
  }
  printf("Parsed %s to range %s %ld %ld %d\n", location, seqRegionName, start, end, strand);

  // Make the new range object before freeing temporary data
  StrUtil_copyString(retSeqRegionName, seqRegionName, 0);
  *retStart  = start;
  *retEnd    = end;
  *retStrand = strand;
  
  // Free up temporary data
  free(tmpLoc);
  for (i=0;i<nTok;i++) {
    free(tokens[i]);
  }
  free(tokens);
  
  return;

/*
  //cleanup any nomenclature like 1_000 or 1 000 or 1,000
  my $number_seps_regex = qr/\s+|,|_/;
  my $separator_regex = qr/(?:-|[.]{2}|\:)?/;
  my $number_regex = qr/[0-9,_ E]+/xms;
  my $strand_regex = qr/[+-1]|-1/xms;
  
  my $regex = qr/^((?:\w|\.|_|-)+) \s* :? \s* ($number_regex)? $separator_regex ($number_regex)? $separator_regex ($strand_regex)? $/xms;
  my ($seq_region_name, $start, $end, $strand);
  if(($seq_region_name, $start, $end, $strand) = $location =~ $regex) {
    
    if(defined $strand) {
      if(!looks_like_number($strand)) {
        $strand = ($strand eq '+') ? 1 : -1;
      }
    }
    
    if(defined $start) {
      $start =~ s/$number_seps_regex//g; 
      if($start < 1) {
        warning "Start was less than 1 (${start}) which is not allowed. Resetting to 1"  if ! $no_warnings;
        $start = 1;
      }
    }
    if(defined $end) {
      $end =~ s/$number_seps_regex//g;
      if($end < 1) {
        throw "Cannot request negative or 0 end indexes through this interface. Given $end but expected something greater than 0" unless $no_errors;
      }
    }
    
    if(defined $start && defined $end && $start > $end) {
      throw "Cannot request a slice whose start is greater than its end. Start: $start. End: $end" unless $no_errors;
    }
  }
  

*/

  
//  return ($seq_region_name, $start, $end, $strand);
}

/*
=head2 fetch_by_region_unique

  Arg [1]    : string $coord_system_name (optional)
               The name of the coordinate system of the slice to be created
               This may be a name of an actual coordinate system or an alias
               to a coordinate system.  Valid aliases are 'seqlevel' or
               'toplevel'.
  Arg [2]    : string $seq_region_name
               The name of the sequence region that the slice will be
               created on.
  Arg [3]    : int $start (optional, default = 1)
               The start of the slice on the sequence region
  Arg [4]    : int $end (optional, default = seq_region length)
               The end of the slice on the sequence region
  Arg [5]    : int $strand (optional, default = 1)
               The orientation of the slice on the sequence region
  Arg [6]    : string $version (optional, default = default version)
               The version of the coordinate system to use (e.g. NCBI33)
  Arg [7]    : boolean $no_fuzz (optional, default = undef (false))
               If true (non-zero), do not use "fuzzy matching" (see below).
  Example    : $slice = $slice_adaptor->fetch_by_region_unique('chromosome', 'HSCHR6_MHC_COX');
  Description: Retrieves a slice on the requested region but returns only the unique
               parts of the slice.  At a minimum the
               name the name of the seq_region to fetch must be provided.

               If no coordinate system name is provided than a slice on the
               highest ranked coordinate system with a matching
               seq_region_name will be returned.  If a version but no
               coordinate system name is provided, the same behaviour will
               apply, but only coordinate systems of the appropriate version
               are considered.  The same applies if the 'toplevel' coordinate
               system is specified, however in this case the version is
               ignored.  The coordinate system should always be specified if
               it is known, since this is unambiguous and faster.

               Some fuzzy matching is performed if no exact match for
               the provided name is found.  This allows clones to be
               fetched even when their version is not known.  For
               example fetch_by_region('clone', 'AC008066') will
               retrieve the sequence_region with name 'AC008066.4'.

               The fuzzy matching can be turned off by setting the
               $no_fuzz argument to a true value.

               If the requested seq_region is not found in the database undef
               is returned.

  Returntype : listref Bio::EnsEMBL::Slice
  Exceptions : throw if no seq_region_name is provided
               throw if invalid coord_system_name is provided
               throw if start > end is provided
  Caller     : general
  Status     : Stable

=cut
*/

Vector *SliceAdaptor_fetchByRegionUnique(SliceAdaptor *sa, char *coordSystemName, char *seqRegionName, long start, long end, int strand, char *version, int noFuzz) {

  Vector *out = Vector_new();

  Slice *slice = SliceAdaptor_fetchByRegion(sa, coordSystemName, seqRegionName, start, end, strand, version, noFuzz);


  if ( sa->asmExcCache == NULL ) {
    SliceAdaptor_buildExceptionCache(sa);
  }

  IDType seqRegionId = SliceAdaptor_getSeqRegionId(sa, slice);

  if (seqRegionId && IDHash_contains(sa->asmExcCache, seqRegionId)) {
    // Dereference symlinked assembly regions.  Take out any regions
    // which are symlinked because these are duplicates.

    Vector *projection = SliceAdaptor_fetchNormalizedSliceProjection(sa, slice, 0);
    Vector_setFreeFunc(projection, ProjectionSegment_free);

    int i;
    for (i=0; i<Vector_getNumElement(projection); i++) {
      ProjectionSegment *segment = Vector_getElementAt(projection, i);
      Slice *toSlice = ProjectionSegment_getToSlice(segment);

      if ( !strcmp(Slice_getSeqRegionName(toSlice),  Slice_getSeqRegionName(slice)) && 
           !CoordSystem_compare( Slice_getCoordSystem(toSlice), Slice_getCoordSystem(slice)) ) {

        Vector_addElement(out, toSlice);
      }
    }
    Vector_free(projection);
  } else if (seqRegionId) {
    Vector_addElement(out, slice);
  } else {
    fprintf(stderr, "Error getting sequence region ID for slice [%s %s %d %d]", 
            (coordSystemName ? coordSystemName : ""), (seqRegionName ? seqRegionName : ""), start, end);
  }

  return out;
}

/*
=head2 fetch_by_name

  Arg [1]    : string $name
  Example    : $name  = 'chromosome:NCBI34:X:1000000:2000000:1';
               $slice = $slice_adaptor->fetch_by_name($name);
               $slice2 = $slice_adaptor->fetch_by_name($slice3->name());
  Description: Fetches a slice using a slice name (i.e. the value returned by
               the Slice::name method).  This is useful if you wish to 
               store a unique identifier for a slice in a file or database or
               pass a slice over a network.
               Slice::name allows you to serialise/marshall a slice and this
               method allows you to deserialise/unmarshal it.

               Returns undef if no seq_region with the provided name exists in
               the database.

  Returntype : Bio::EnsEMBL::Slice or undef
  Exceptions : throw if incorrent arg provided
  Caller     : Pipeline
  Status     : Stable

=cut
*/

Slice *SliceAdaptor_fetchByName(SliceAdaptor *sa, char *name) {

  if (name == NULL) {
    fprintf(stderr,"name argument is required\n");
    exit(1);
  }

  int nTok;
  char **tokens;
  
  StrUtil_tokenizeByDelim(&tokens, &nTok, name, ":");

  if (nTok < 3 || nTok > 6) {
    fprintf(stderr,"Malformed slice name [%s].  Format is coord_system:version:name:start:end:strand\n", name);
    exit(1);
  }

  // Rearrange arguments to suit fetch_by_region
  // Example    : $name  = 'chromosome:NCBI34:X:1000000:2000000:1';
  char *coordSystemName = tokens[0];
  char *version         = tokens[1][0] != '\0' ? tokens[1] : NULL;
  char *seqRegionName   = tokens[2][0] != '\0' ? tokens[2] : NULL;

  long start = POS_UNDEF;
  long end   = POS_UNDEF;
  int strand = STRAND_UNDEF;

  if (nTok > 3 && tokens[3][0] != '\0') {
    if (!StrUtil_isLongInteger(&start, tokens[3])) {
      fprintf(stderr,"Malformed slice name [%s].  Format is coord_system:version:name:start:end:strand\n", name);
      exit(1);
    }
  }
  if (nTok > 4 && tokens[4][0] != '\0') {
    if (!StrUtil_isLongInteger(&end, tokens[4])) {
      fprintf(stderr,"Malformed slice name [%s].  Format is coord_system:version:name:start:end:strand\n", name);
      exit(1);
    }
  }
  if (nTok > 5 && tokens[5][0] != '\0') {
    long tmp;
    if (!StrUtil_isLongInteger(&tmp, tokens[5])) {
      fprintf(stderr,"Malformed slice name [%s].  Format is coord_system:version:name:start:end:strand\n", name);
      exit(1);
    }
    strand = (int)tmp;
  }

  Slice *slice = SliceAdaptor_fetchByRegion(sa, coordSystemName, seqRegionName, start, end, strand, version, 0 /* No fuzz */);

  // Free up temporary data
  int i;
  for (i=0;i<nTok;i++) {
    free(tokens[i]);
  }
  free(tokens);
  
  return slice;
}



/*
=head2 fetch_by_seq_region_id

  Arg [1]    : string $seq_region_id
               The internal identifier of the seq_region to create this slice
               on
  Arg [2]    : optional start
  Arg [3]    : optional end
  Arg [4]    : optional strand
  Example    : $slice = $slice_adaptor->fetch_by_seq_region_id(34413);
  Description: Creates a slice object of an entire seq_region using the
               seq_region internal identifier to resolve the seq_region.
               Returns undef if no such slice exists.
  Returntype : Bio::EnsEMBL::Slice or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/

Slice *SliceAdaptor_fetchBySeqRegionId(SliceAdaptor *sa, IDType seqRegionId, long start, long end, int strand) {
  char *name;
  long length;
  CoordSystem *cs;
  IDType csId;

  CoordSystemAdaptor *csa = DBAdaptor_getCoordSystemAdaptor(sa->dba);

  if ( !IDHash_contains(sa->srIdCache, seqRegionId)) {
    char qStr[1024];
    sprintf(qStr, "SELECT sr.name, sr.coord_system_id, sr.length "
                  "FROM seq_region sr "
                  "WHERE sr.seq_region_id = "IDFMTSTR, seqRegionId );

    StatementHandle *sth = sa->prepare((BaseAdaptor *)sa,qStr,strlen(qStr));
    sth->execute(sth);

    if ( sth->numRows(sth) == 0 ) {
      sth->finish(sth);
      return NULL;
    }

    ResultRow *row = sth->fetchRow(sth);
// NOTE: To make easier to handle memory management of name have changed logic compared to perl to always fetch from cache after added to cache
    name = row->getStringAt(row,0);
    csId = row->getLongLongAt(row,1);
    length = row->getLongAt(row,2);

    cs = CoordSystemAdaptor_fetchByDbID(csa, csId);

    //cache results to speed up repeated queries
    DBAdaptor_addToSrCaches(sa->dba, seqRegionId, name, csId, length);
    sth->finish(sth);
  }

  SeqRegionCacheEntry *cacheData = IDHash_getValue(sa->srIdCache, seqRegionId);

  name   = cacheData->regionName;
  length = cacheData->regionLength;
  csId   = cacheData->csId;

  cs = CoordSystemAdaptor_fetchByDbID(csa, csId);

  start = (start == POS_UNDEF) ? 1 : start;
  end = (end == POS_UNDEF) ?     length : end;
  strand = (strand == STRAND_UNDEF) ? 1 : strand;

  Slice *slice = Slice_new(name, start, end, strand, length, cs, sa);

  return slice;
}



/*
=head2 get_seq_region_id

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               The slice to fetch a seq_region_id for
  Example    : $srid = $slice_adaptor->get_seq_region_id($slice);
  Description: Retrieves the seq_region id (in this database) given a slice
               Seq region ids are not stored on the slices themselves
               because they are intended to be somewhat database independent
               and seq_region_ids vary accross databases.
  Returntype : int
  Exceptions : throw if the seq_region of the slice is not in the db
               throw if incorrect arg provided
  Caller     : BaseFeatureAdaptor
  Status     : Stable

=cut
*/

IDType SliceAdaptor_getSeqRegionId(SliceAdaptor *sa, Slice *slice) {

  if (slice==NULL) { // || !ref($slice) || !($slice->isa('Bio::EnsEMBL::Slice') or $slice->isa('Bio::EnsEMBL::LRGSlice'))) {
    fprintf(stderr, "Slice argument is required\n");
    exit(1);
  }
  
  char *seqRegionName = Slice_getSeqRegionName(slice);

  if (!seqRegionName)
    return 0;

  char key[1024];
//  sprintf(key, "%s:"IDFMTSTR, seqRegionName, CoordSystem_getDbID(Slice_getCoordSystem(slice)));

  char *endP = stpcpy(key,seqRegionName);
  *endP++ = ':';
  stpcpy(endP, CoordSystem_getDbIDStr(Slice_getCoordSystem(slice)));


  //if (StringHash_contains(sa->srNameCache, key)) {
    SeqRegionCacheEntry *cacheData = StringHash_getValue(sa->srNameCache, key);
    if (cacheData) return cacheData->regionId;
  //}

  IDType csId = CoordSystem_getDbID(Slice_getCoordSystem(slice));

  char qStr[1024];
  sprintf(qStr,"SELECT seq_region_id, length " 
               "FROM seq_region " 
               "WHERE name = '%s' AND coord_system_id = "IDFMTSTR, seqRegionName, csId);

  StatementHandle *sth = sa->prepare((BaseAdaptor *)sa,qStr,strlen(qStr));

  //force seq_region_name cast to string so mysql cannot treat as int
  sth->execute(sth);

  if (sth->numRows(sth) != 1) {
    fprintf(stderr, "Non existent or ambiguous seq_region:\n  coord_system=["IDFMTSTR"],\n   name=[%s]\n", csId, seqRegionName);
    exit(1);
  }

  ResultRow *row = sth->fetchRow(sth);

  IDType seqRegionId = row->getLongLongAt(row,0);
  long length        = row->getLongAt(row,1);

  sth->finish(sth);

  //cache information for future requests
  DBAdaptor_addToSrCaches(sa->dba, seqRegionId, seqRegionName, csId, length);

  return seqRegionId;
}


/*
=head2 fetch_all

  Arg [1]    : string $coord_system_name
               The name of the coordinate system to retrieve slices of.
               This may be a name of an acutal coordinate system or an alias
               to a coordinate system.  Valid aliases are 'seqlevel' or
               'toplevel'.
  Arg [2]    : string $coord_system_version (optional)
               The version of the coordinate system to retrieve slices of
  Arg [3]    : bool $include_non_reference (optional)
               If this argument is not provided then only reference slices
               will be returned. If set, both reference and non refeference
               slices will be rerurned.
  Arg [4]    : int $include_duplicates (optional)
               If set duplicate regions will be returned.
               
               NOTE: if you do not use this option and you have a PAR
               (pseudo-autosomal region) at the beginning of your seq_region
               then your slice will not start at position 1, so coordinates
               retrieved from this slice might not be what you expected.

  Arg[5]     : bool $include_lrg (optional)  (default 0)
               If set lrg regions will be returned aswell.


  Example    : @chromos = @{$slice_adaptor->fetch_all('chromosome','NCBI33')};
               @contigs = @{$slice_adaptor->fetch_all('contig')};

               # get even non-reference regions
               @slices = @{$slice_adaptor->fetch_all('toplevel',undef,1)};

               # include duplicate regions (such as pseudo autosomal regions)
               @slices = @{$slice_adaptor->fetch_all('toplevel', undef,0,1)};

  Description: Retrieves slices of all seq_regions for a given coordinate
               system.  This is analagous to the methods fetch_all which were
               formerly on the ChromosomeAdaptor, RawContigAdaptor and
               CloneAdaptor classes.  Slices fetched span the entire
               seq_regions and are on the forward strand.
               If the coordinate system with the provided name and version
               does not exist an empty list is returned.
               If the coordinate system name provided is 'toplevel', all
               non-redundant toplevel slices are returned (note that any
               coord_system_version argument is ignored in that case).

               Retrieved slices can be broken into smaller slices using the
               Bio::EnsEMBL::Utils::Slice module.

  Returntype : listref of Bio::EnsEMBL::Slices
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/

Vector *SliceAdaptor_fetchAll(SliceAdaptor *sa, char *csName, char *csVersion, int flags) {

  Vector *out = Vector_new();
  // 
  // verify existence of requested coord system and get its id
  // 
  CoordSystemAdaptor *csa = DBAdaptor_getCoordSystemAdaptor(sa->dba);

  CoordSystem *origCs = CoordSystemAdaptor_fetchByName(csa, csName, csVersion);

  
  if (origCs == NULL) return out;

  IDHash *badVals = IDHash_new(IDHASH_MEDIUM);

  // 
  // Get a hash of non reference seq regions
  // 
  if ( ! (flags & SA_INCLUDE_NON_REFERENCE) ) {
    char qStr[1024];
    sprintf(qStr, "SELECT sr.seq_region_id "
                  "FROM seq_region sr, seq_region_attrib sra, "
                       "attrib_type at, coord_system cs "
                  "WHERE at.code = 'non_ref' "
                  "AND sra.seq_region_id = sr.seq_region_id "
                  "AND at.attrib_type_id = sra.attrib_type_id "
                  "AND sr.coord_system_id = cs.coord_system_id "
                  "AND cs.species_id = "IDFMTSTR, SliceAdaptor_getSpeciesID(sa) );

    StatementHandle *sth = sa->prepare((BaseAdaptor *)sa,qStr,strlen(qStr));
    sth->execute(sth);

    ResultRow *row;
    while ((row = sth->fetchRow(sth))) {
      IDType seqRegionId = row->getLongLongAt(row,0);
      
      IDHash_add(badVals, seqRegionId, &trueVal);
    }

    sth->finish(sth);
  }

  //
  // if we do not want lrg's then add them to the bad list;
  //
  // LRGs are always bad in my book, so just get rid of them
  if (1) { 
  //if ( !$include_lrg ) {
    char qStr[1024];
    sprintf(qStr, "SELECT sr.seq_region_id "
                  "FROM seq_region sr, seq_region_attrib sra, "
                       "attrib_type at, coord_system cs "
                  "WHERE at.code = 'LRG' "
                  "AND sra.seq_region_id = sr.seq_region_id "
                  "AND at.attrib_type_id = sra.attrib_type_id "
                  "AND sr.coord_system_id = cs.coord_system_id "
                  "AND cs.species_id = "IDFMTSTR, SliceAdaptor_getSpeciesID(sa) );

    StatementHandle *sth = sa->prepare((BaseAdaptor *)sa,qStr,strlen(qStr));
    sth->execute(sth);

    ResultRow *row;
    while ((row = sth->fetchRow(sth))) {
      IDType seqRegionId = row->getLongLongAt(row,0);
      
      if (!IDHash_contains(badVals, seqRegionId)) {
        IDHash_add(badVals, seqRegionId, &trueVal);
      }
    }

    sth->finish(sth);
  }

  //
  // Retrieve the seq_regions from the database
  //

  StatementHandle *sth;
  char qStr[1024];

  if ( CoordSystem_getIsTopLevel(origCs)) {
    sprintf(qStr, "SELECT sr.seq_region_id, sr.name, "
                         "sr.length, sr.coord_system_id "
                  "FROM seq_region sr, seq_region_attrib sra, "
                       "attrib_type at, coord_system cs "
                  "WHERE at.code = 'toplevel' "
                  "AND sra.seq_region_id = sr.seq_region_id "
                  "AND at.attrib_type_id = sra.attrib_type_id "
                  "AND sr.coord_system_id = cs.coord_system_id "
                  "AND cs.species_id = "IDFMTSTR, SliceAdaptor_getSpeciesID(sa) );

    sth = sa->prepare((BaseAdaptor *)sa,qStr,strlen(qStr));
    sth->execute(sth);
  } else {
    sprintf(qStr,"SELECT sr.seq_region_id, sr.name, "
                        "sr.length, sr.coord_system_id "
                 "FROM seq_region sr "
                 "WHERE sr.coord_system_id = "IDFMTSTR, CoordSystem_getDbID(origCs) );

    sth = sa->prepare((BaseAdaptor *)sa,qStr,strlen(qStr));
    sth->execute(sth);
  }

  int cacheCount = 0;

  ResultRow *row;
  while ((row = sth->fetchRow(sth))) {
    IDType seqRegionId = row->getLongLongAt(row,0);
    char * name        = row->getStringAt(row,1);
    long   length      = row->getLongAt(row,2);
    IDType csId        = row->getLongLongAt(row,3);

// NIY: Should not use defined here, should use exists (I'll use 'contains')
    //if (!defined($bad_vals{$seq_region_id})){
    if ( ! IDHash_contains(badVals, seqRegionId)) {
      CoordSystem *cs = CoordSystemAdaptor_fetchByDbID(csa, csId);

      if (cs == NULL) {
        fprintf(stderr,"seq_region %s references non-existent coord_system "IDFMTSTR".\n", name, csId);
        exit(1);
      }

      //cache values for future reference, but stop adding to the cache once we
      //we know we have filled it up

      // Ughhh - this is the only place we check for this max. Is it worth it????

      // Just cache everything - its C, how big can it be
      //if (cacheCount < $Bio::EnsEMBL::Utils::SeqRegionCache::SEQ_REGION_CACHE_SIZE) {
        DBAdaptor_addToSrCaches(sa->dba, seqRegionId, name, csId, length);

        cacheCount++;
      //}


      Slice *slice = Slice_new(name, 1, length, 1, length, cs, sa);

      if (! (flags & SA_INCLUDE_DUPLICATES)) {
        // test if this slice *could* have a duplicate (exception) region
        if (sa->asmExcCache == NULL) SliceAdaptor_buildExceptionCache(sa);

        if (IDHash_contains(sa->asmExcCache, seqRegionId)) {
          // Dereference symlinked assembly regions.  Take out
          // any regions which are symlinked because these are duplicates

          Vector *projection = SliceAdaptor_fetchNormalizedSliceProjection(sa, slice, 0);
          Vector_setFreeFunc(projection, ProjectionSegment_free);

          int i;
          for (i=0; i<Vector_getNumElement(projection); i++) {
            ProjectionSegment *segment = Vector_getElementAt(projection, i);

            Slice *toSlice = ProjectionSegment_getToSlice(segment);

            if (!strcmp(Slice_getSeqRegionName(toSlice), Slice_getSeqRegionName(slice)) &&
                !CoordSystem_compare(Slice_getCoordSystem(toSlice), Slice_getCoordSystem(slice))) {
              Vector_addElement(out, toSlice);
            }
          }
          Vector_free(projection);
        } else {
          // no duplicate regions
          Vector_addElement(out, slice);
        }
      } else {
        // we want duplicates anyway so do not do any checks
        Vector_addElement(out, slice);
      }
    }
  }

  IDHash_free(badVals,NULL);

  return out;
}


/*
=head2 fetch_all_karyotype
  Example    : my $top = $slice_adptor->fetch_all_karyotype()
  Description: returns the list of all slices which are part of the karyotype
  Returntype : listref of Bio::EnsEMBL::Slices
  Caller     : general
  Status     : At Risk

=cut
*/
Vector *SliceAdaptor_fetchAllKaryotype(SliceAdaptor *sa) {
  CoordSystemAdaptor *csa = DBAdaptor_getCoordSystemAdaptor(sa->dba);

  char qStr[1024];
  sprintf(qStr, "SELECT sr.seq_region_id, sr.name, "
                      "sr.length, sr.coord_system_id "
                 "FROM seq_region sr, seq_region_attrib sra, "
                      "attrib_type at, coord_system cs "
                 "WHERE at.code = 'karyotype_rank' "
                 "AND at.attrib_type_id = sra.attrib_type_id "
                 "AND sra.seq_region_id = sr.seq_region_id "
                 "AND sr.coord_system_id = cs.coord_system_id "
                 "AND cs.species_id = "IDFMTSTR, SliceAdaptor_getSpeciesID(sa) );

  StatementHandle *sth = sa->prepare((BaseAdaptor *)sa,qStr,strlen(qStr));

  sth->execute(sth);

  Vector *out = Vector_new();

  ResultRow *row;
  while ((row = sth->fetchRow(sth))) {
    char * name        = row->getStringAt(row,1);
    long   length      = row->getLongAt(row,2);
    IDType csId        = row->getLongLongAt(row,3);

    CoordSystem *cs = CoordSystemAdaptor_fetchByDbID(csa, csId);

    Slice *slice = Slice_new(name, 1, length, 1, length, cs, sa);

    Vector_addElement(out, slice);
  }

  return out;
}

/*
=head2 is_toplevel
  Arg        : int seq_region_id 
  Example    : my $top = $slice_adptor->is_toplevel($seq_region_id)
  Description: Returns 1 if slice is a toplevel slice else 0
  Returntype : int
  Caller     : Slice method is_toplevel
  Status     : At Risk

=cut
*/

int SliceAdaptor_isTopLevel(SliceAdaptor *sa, IDType id) {
  char qStr[1024];

  sprintf(qStr, "SELECT at.code from seq_region_attrib sra, attrib_type at "
                "WHERE sra.seq_region_id = "IDFMTSTR" "
                "AND at.attrib_type_id = sra.attrib_type_id "
                "AND at.code = 'toplevel'", id);

  StatementHandle *sth = sa->prepare((BaseAdaptor *)sa,qStr,strlen(qStr));

  sth->execute(sth);

  int isTop = 0;
  if (sth->numRows(sth) > 0) {
    isTop = 1;
  }

  sth->finish(sth);
  return isTop;
}

/*
=head2 has_karyotype
  Arg        : int seq_region_id 
  Example    : my $karyotype = $slice_adptor->has_karyotype($seq_region_id)
  Description: Returns 1 if slice is a part of a karyotype else 0
  Returntype : int
  Caller     : Slice method has_karyotype
  Status     : At Risk

=cut
*/

int SliceAdaptor_hasKaryotype(SliceAdaptor *sa, IDType id) {
  char qStr[1024];

  sprintf(qStr,"SELECT at.code from seq_region_attrib sra, attrib_type at "
               "WHERE sra.seq_region_id = "IDFMTSTR" "
               "AND at.attrib_type_id = sra.attrib_type_id "
               "AND at.code = 'karyotype_rank'", id);

  StatementHandle *sth = sa->prepare((BaseAdaptor *)sa,qStr,strlen(qStr));
  sth->execute(sth);

  int hasKary = 0;
  if (sth->numRows(sth) > 0) {
    hasKary = 1;
  }

  sth->finish(sth);
  return hasKary;
}

/*
=head2 get_karyotype_rank
  Arg        : int seq_region_id 
  Example    : my $rank = $slice_adptor->get_karyotype_rank($seq_region_id)
  Description: Returns the rank of a slice if it is part of the karyotype else 0
  Returntype : int
  Caller     : Slice method get_karyotype_rank
  Status     : At Risk

=cut
*/

int SliceAdaptor_getKaryotypeRank(SliceAdaptor *sa, IDType id) {
  char qStr[1024];

  sprintf(qStr, "SELECT sra.value from seq_region_attrib sra, attrib_type at "
                "WHERE sra.seq_region_id = "IDFMTSTR" "
                "AND at.attrib_type_id = sra.attrib_type_id "
                "AND at.code = 'karyotype_rank'", id);

  StatementHandle *sth = sa->prepare((BaseAdaptor *)sa,qStr,strlen(qStr));
  sth->execute(sth);

// NIY: perl has these lines which shouldn't be there
//  my $code;
//  $sth->bind_columns( \$code );

  if (sth->numRows(sth) == 0) {
    fprintf(stderr, "No karyotype rank for id "IDFMTSTR"\n", id);
    exit(1);
  }

  ResultRow *row = sth->fetchRow(sth);
  int rank = row->getIntAt(row,0);

  sth->finish(sth);

  return rank;
}

/*
=head2 is_reference
  Arg        : int seq_region_id 
  Example    : my $reference = $slice_adptor->is_reference($seq_region_id)
  Description: Returns 1 if slice is a reference slice else 0
  Returntype : int
  Caller     : Slice method is_reference
  Status     : At Risk

=cut
*/

int SliceAdaptor_isReference(SliceAdaptor *sa, IDType id) {
  char qStr[1024];
  sprintf(qStr,"SELECT at.code from seq_region_attrib sra, attrib_type at "
               "WHERE sra.seq_region_id = "IDFMTSTR" "
               "AND at.attrib_type_id = sra.attrib_type_id "
               "AND at.code = 'non_ref'", id);

  StatementHandle *sth = sa->prepare((BaseAdaptor *)sa,qStr,strlen(qStr));
  sth->execute(sth);

  int isRef = 1;
  if (sth->numRows(sth) > 0) {
    isRef = 0;
  }

  sth->finish(sth);
  return isRef;
}

/* Not implementing circular stuff
=head2 is_circular

  Arg[1]      : int seq_region_id
  Example     : my $circular = $slice_adptor->is_circular($seq_region_id);
  Description : Indicates if the sequence region was circular or not
  Returntype  : Boolean
  
=cut

sub is_circular {
  my ($self, $id) = @_;
  
  if (! defined $self->{is_circular}) {
    $self->_build_circular_slice_cache();
  }
  
  return 0 if $self->{is_circular} == 0;
  return (exists $self->{circular_sr_id_cache}->{$id}) ? 1 : 0;
}
*/

/*
=head2 fetch_by_band

 Title   : fetch_by_band
 Usage   :
 Function: Does not work please use fetch_by_chr_band
 Example :
 Returns : Bio::EnsEMBL::Slice
 Args    : the band name
 Status     : AT RISK

=cut
*/

/* Comment above claims does not work, so don't implement
sub fetch_by_band {
  my ($self,$band) = @_;

  my $sth = $self->dbc->prepare
        ("select s.name,max(k.seq_region_id)-min(k.seq_region_id, min(k.seq_region_start), max(k.seq_region_id) " .
         "from karyotype as k " .
         "where k.band like ? and k.seq_region_id = s.seq_region_id");

  $sth->bind_param(1,"$band%",SQL_VARCHAR);
  $sth->execute();
  my ( $seq_region_name, $discrepancy, $seq_region_start, $seq_region_end) = $sth->fetchrow_array;

  if($seq_region_name && $discrepancy>0) {
    throw("Band maps to multiple seq_regions");
  } else {
    return $self->fetch_by_region('toplevel',$seq_region_name,$seq_region_start,$seq_region_end);
  }
  throw("Band not recognised in database");
}
*/

/*
=head2 fetch_by_chr_band

 Title   : fetch_by_chr_band
 Usage   :
 Function: create a Slice representing a series of bands
 Example :
 Returns : Bio::EnsEMBL::Slice
 Args    : the band name
 Status     : Stable

=cut
*/
Slice *SliceAdaptor_fetchByChrBand(SliceAdaptor *sa, char *chr, char *band) {

  Slice *chrSlice = SliceAdaptor_fetchByRegion(sa, "toplevel", chr, POS_UNDEF, POS_UNDEF, STRAND_UNDEF, NULL, 0);
  IDType seqRegionId = SliceAdaptor_getSeqRegionId(sa, chrSlice);

  if (seqRegionId) {
    char qStr[1024];
    sprintf(qStr,"SELECT MIN(k.seq_region_start), "
                        "MAX(k.seq_region_end) "
                        "FROM karyotype k "
                        "WHERE k.seq_region_id = "IDFMTSTR
                        "AND k.band LIKE '%s%%%%'", seqRegionId, band );

  StatementHandle *sth = sa->prepare((BaseAdaptor *)sa,qStr,strlen(qStr));
  sth->execute(sth);

    // How can we have undefined ??? When no row???
    if (sth->numRows(sth) == 1) {
      ResultRow *row = sth->fetchRow(sth);
      long sliceStart = row->getLongAt(row,0);
      long sliceEnd   = row->getLongAt(row,1);

      // Need to free chrSlice
      Slice_free(chrSlice);
      sth->finish(sth);

      return SliceAdaptor_fetchByRegion(sa, "toplevel", chr, sliceStart, sliceEnd, STRAND_UNDEF, NULL, 0);
    }

    sth->finish(sth);
  } else {
    fprintf(stderr, "Error getting sequence region ID for slice");
  }

  Slice_free(chrSlice);

  fprintf(stderr,"Band not recognised in database\n");
  exit(1);
} 


/*
=head2 fetch_by_exon_stable_id

  Arg [1]    : string $exonid
               The stable id of the exon around which the slice is 
               desired
  Arg [2]    : (optional) int $size
               The length of the flanking regions the slice should encompass 
               on either side of the exon (0 by default)
  Example    : $slc = $sa->fetch_by_exon_stable_id('ENSE00000302930',10);
  Description: Creates a slice around the region of the specified exon. 
               If a context size is given, the slice is extended by that 
               number of basepairs on either side of the exon.
               
               The slice will be created in the exon's native coordinate system
               and in the forward orientation.
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : Thrown if the exon is not in the database.
  Caller     : general
  Status     : Stable

=cut
*/
Slice *SliceAdaptor_fetchByExonStableId(SliceAdaptor *sa, char *exonId, int size, int isPercent) {
  if (exonId == NULL) {
    fprintf(stderr,"Exon argument is required.");
    exit(1);
  }

  ExonAdaptor *ea = DBAdaptor_getExonAdaptor(sa->dba);
  Exon *exon = ExonAdaptor_fetchByStableId(ea, exonId);

  if (exon == NULL) {
    fprintf(stderr,"Exon [%s] does not exist in DB.", exonId);
    exit(1);
  }

  return SliceAdaptor_fetchByFeature(sa, (SeqFeature *)exon, size, isPercent);
}

/*
=head2 fetch_by_transcript_stable_id

  Arg [1]    : string $transcriptid
               The stable id of the transcript around which the slice is 
               desired
  Arg [2]    : (optional) int $size
               The length of the flanking regions the slice should encompass 
               on either side of the transcript (0 by default)
  Example    : $slc = $sa->fetch_by_transcript_stable_id('ENST00000302930',10);
  Description: Creates a slice around the region of the specified transcript. 
               If a context size is given, the slice is extended by that 
               number of basepairs on either side of the
               transcript.
               
               The slice will be created in the transcript's native coordinate
               system and in the forward orientation.
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : Thrown if the transcript is not in the database.
  Caller     : general
  Status     : Stable

=cut
*/
Slice *SliceAdaptor_fetchByTranscriptStableId(SliceAdaptor *sa, char *transcriptId, int size, int isPercent) {
  if (transcriptId == NULL) {
    fprintf(stderr,"Transcript argument is required.");
    exit(1);
  }

  TranscriptAdaptor *ta = DBAdaptor_getTranscriptAdaptor(sa->dba);
  Transcript *transcript = TranscriptAdaptor_fetchByStableId(ta, transcriptId);

  if (transcript == NULL) {
    fprintf(stderr,"Transcript [%s] does not exist in DB.", transcriptId);
    exit(1);
  }

  return SliceAdaptor_fetchByFeature(sa, (SeqFeature *)transcript, size, isPercent);
}

/*
=head2 fetch_by_transcript_id

  Arg [1]    : int $transcriptid
               The unique database identifier of the transcript around which 
               the slice is desired
  Arg [2]    : (optional) int $size
               The length of the flanking regions the slice should encompass 
               on either side of the transcript (0 by default)
  Example    : $slc = $sa->fetch_by_transcript_id(24, 1000);
  Description: Creates a slice around the region of the specified transcript.
               If a context size is given, the slice is extended by that
               number of basepairs on either side of the
               transcript.
               
               The slice will be created in the transcript's native coordinate
               system and in the forward orientation.
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : throw on incorrect args
               throw if transcript is not in database
  Caller     : general
  Status     : Stable

=cut
*/

Slice *SliceAdaptor_fetchByTranscriptId(SliceAdaptor *sa, IDType transcriptId, int size, int isPercent) {
  TranscriptAdaptor *ta = DBAdaptor_getTranscriptAdaptor(sa->dba);
  Transcript *transcript = (Transcript *)TranscriptAdaptor_fetchByDbID(ta, transcriptId);

  if (transcript == NULL) {
    fprintf(stderr,"Transcript ["IDFMTSTR"] does not exist in DB.", transcriptId);
    exit(1);
  }

  return SliceAdaptor_fetchByFeature(sa, (SeqFeature *)transcript, size, isPercent);
}

/*

=head2 fetch_by_gene_stable_id

  Arg [1]    : string $geneid
               The stable id of the gene around which the slice is
               desired
  Arg [2]    : (optional) int $size
               The length of the flanking regions the slice should encompass
               on either side of the gene (0 by default)
  Example    : $slc = $sa->fetch_by_gene_stable_id('ENSG00000012123',10);
  Description: Creates a slice around the region of the specified gene.
               If a context size is given, the slice is extended by that
               number of basepairs on either side of the gene.
               
               The slice will be created in the gene's native coordinate system
               and in the forward orientation.
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : throw on incorrect args
               throw if transcript does not exist
  Caller     : general
  Status     : Stable

=cut
*/

Slice *SliceAdaptor_fetchByGeneStableId(SliceAdaptor *sa, char *geneId, int size, int isPercent) {
  if (geneId == NULL) {
    fprintf(stderr,"Gene argument is required.");
    exit(1);
  }

  GeneAdaptor *ga = DBAdaptor_getGeneAdaptor(sa->dba);
  Gene *gene = GeneAdaptor_fetchByStableId(ga, geneId);

  if (gene == NULL) {
    fprintf(stderr,"Gene [%s] does not exist in DB.", geneId);
    exit(1);
  }

  return SliceAdaptor_fetchByFeature(sa, (SeqFeature *)gene, size, isPercent);
}

/*

=head2 fetch_by_Feature

  Arg [1]    : Bio::EnsEMBL::Feature $feat
               The feature to fetch the slice around
  Arg [2]    : int size (optional)
               The desired number of flanking basepairs around the feature.
               The size may also be provided as a percentage of the feature 
               size such as 200% or 80.5%.
  Example    : $slice = $slice_adaptor->fetch_by_Feature($feat, 100);
  Description: Retrieves a slice around a specific feature.  All this really
               does is return a resized version of the slice that the feature
               is already on. Note that slices returned from this method
               are always on the forward strand of the seq_region regardless of
               the strandedness of the feature passed in.
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : throw if the feature does not have an attached slice
               throw if feature argument is not provided
  Caller     : fetch_by_gene_stable_id, fetch_by_transcript_stable_id,
               fetch_by_gene_id, fetch_by_transcript_id
  Status     : Stable

=cut
*/

// Note: In C can't do percent magic - size is a number
// In perl if size not specified its 0, so use that when you call this function and don't want anything special for size
Slice *SliceAdaptor_fetchByFeature(SliceAdaptor *sa, SeqFeature *feature, int size, int isPercent) {

  Slice *slice = (Slice *)SeqFeature_getSlice(feature);
  if (slice == NULL) { // NIY: Do we need to do this check in C? || !($slice->isa('Bio::EnsEMBL::Slice') or $slice->isa('Bio::EnsEMBL::LRGSlice') )) {
    fprintf(stderr,"Feature must be attached to a valid slice.\n");
    exit(1);
  }

  long fStart = SeqFeature_getStart(feature);
  long fEnd   = SeqFeature_getEnd(feature);

/* NIY: Not sure if this will every happen - hopefully not
  if(!defined($fstart) || !defined($fend)) {
    throw('Feature must have defined start and end.');
  }
*/

  //convert the feature slice coordinates to seq_region coordinates
  long sliceStart  = Slice_getStart(slice);
  long sliceEnd    = Slice_getEnd(slice);
  int  sliceStrand = Slice_getStrand(slice);

  if (sliceStart != 1 || sliceStrand != 1) {
    if (sliceStrand == 1) {
      fStart = fStart + sliceStart - 1;
      fEnd   = fEnd   + sliceStart - 1;
    } else {
      long tmpStart = fStart;
      fStart = sliceEnd - fEnd     + 1;
      fEnd   = sliceEnd - tmpStart + 1;
    }
  }

  // Size may be stored as a %age of the length of the feature
  // Size = 100% gives no context
  // Size = 200% gives context - 50% the size of the feature either side of 
  // feature

  if (isPercent) {
    size = (int)( ((double)size-100.0)/200.0 * (double)(fEnd-fStart+1) );
  }

  //return a new slice covering the region of the feature

  // NIY: I called it newSlice instead of 'S'
  //      What is raw_feature_strand below??

  Slice *newSlice = Slice_new(Slice_getSeqRegionName(slice), fStart - size, fEnd + size, 1, 
                              Slice_getSeqRegionLength(slice), Slice_getCoordSystem(slice), sa);

// Hhhhh! What's this??? Doesn't seem to be used so its gone
  //$S->{'_raw_feature_strand'}  = $feature->strand * $slice_strand if $feature->can('strand');

  return newSlice;
}

/*
=head2 fetch_by_misc_feature_attribute

  Arg [1]    : string $attribute_type
               The code of the attribute type
  Arg [2]    : (optional) string $attribute_value
               The value of the attribute to fetch by
  Arg [3]    : (optional) int $size
               The amount of flanking region around the misc feature desired.
  Example    : $slice = $sa->fetch_by_misc_feature_attribute('superctg',
                                                             'NT_030871');
               $slice = $sa->fetch_by_misc_feature_attribute('synonym',
                                                             'AL00012311',
                                                             $flanking);
  Description: Fetches a slice around a MiscFeature with a particular
               attribute type and value. If no value is specified then
               the feature with the particular attribute is used.
               If no size is specified then 0 is used.
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : Throw if no feature with the specified attribute type and value
               exists in the database
               Warning if multiple features with the specified attribute type
               and value exist in the database.
  Caller     : webcode
  Status     : Stable

=cut
*/

/* Need to implement MiscFeatureAdaptor before enabling this one
Slice *SliceAdaptor_fetchByMiscFeatureAttribute(SliceAdaptor *sa, char *attribTypeCode, char *attribValue, int size, int isPercent)  {

  MiscFeatureAdaptor *mfa = DBAdaptor_getMiscFeatureAdaptor(sa->dba);

  Vector *feats = MiscFeatureAdaptor_fetchAllByAttributeTypeValue(mfa, attribTypeCode, attribValue);

  if (Vector_getNumElement(feats) == 0) {
    fprintf(stderr,"MiscFeature with $attrib_type_code=$attrib_value does not exist in DB.\n");
    exit(1);
  }

  if (Vector_getNumElement(feats) > 1) {
    fprintf(stderr,"Warning: MiscFeature with %s=%s is ambiguous - using first one found.\n",
            attribTypeCode, attribValue);
  }

  SeqFeature *feat = Vector_getElementAt(feats, 0);
  
  Slice *slice = SliceAdaptor_fetchByFeature(sa, feat, size, isPercent);

//NIY: Freeing the actual feature in the vector
  Vector_free(feats);

  return slice;
}
*/

/*
=head2 fetch_normalized_slice_projection

  Arg [1]    : Bio::EnsEMBL::Slice $slice
  Arg [2]    : boolean $filter_projections 
               Optionally filter the projections to remove anything 
               which is the same sequence region as the given slice
  Example    :  ( optional )
  Description: gives back a project style result. The returned slices 
               represent the areas to which there are symlinks for the 
               given slice. start, end show which area on given slice is 
               symlinked
  Returntype : [[start,end,$slice][]]
  Exceptions : none
  Caller     : BaseFeatureAdaptor
  Status     : Stable

=cut
*/
typedef struct AssemblyExceptionUnitStruct {
  long seqRegionStart;
  long seqRegionEnd;
  IDType excSeqRegionId;
  long excSeqRegionStart;
  long excSeqRegionEnd;
} AssemblyExceptionUnit;

AssemblyExceptionUnit *AssemblyExceptionUnit_new(long seqRegionStart, long seqRegionEnd, IDType excSeqRegionId,
                                                 long excSeqRegionStart, long excSeqRegionEnd) {
  AssemblyExceptionUnit *aeu;
  if ((aeu = (AssemblyExceptionUnit *)calloc(1,sizeof(AssemblyExceptionUnit))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for AssemblyExceptionUnit\n");
    return NULL;
  }

  aeu->seqRegionStart    = seqRegionStart;
  aeu->seqRegionEnd      = seqRegionEnd;
  aeu->excSeqRegionId    = excSeqRegionId;
  aeu->excSeqRegionStart = excSeqRegionStart;
  aeu->excSeqRegionEnd   = excSeqRegionEnd;
  
  return aeu;
}

int AssemblyExceptionUnit_endCompFunc(const void *one, const void *two) {
  AssemblyExceptionUnit *aeu1 = *((AssemblyExceptionUnit**)one);
  AssemblyExceptionUnit *aeu2 = *((AssemblyExceptionUnit**)two);

  return aeu1->seqRegionEnd - aeu2->seqRegionEnd;
}

void AssemblyExceptionUnit_free(AssemblyExceptionUnit *aeu) {
  free(aeu);
}


Vector *SliceAdaptor_fetchNormalizedSliceProjection(SliceAdaptor *sa, Slice *slice, int filterProjections) {
  IDType sliceSeqRegionId = SliceAdaptor_getSeqRegionId(sa, slice);

  if (sa->asmExcCache == NULL) SliceAdaptor_buildExceptionCache(sa);

  Vector *result = NULL;

  if (sliceSeqRegionId && IDHash_contains(sa->asmExcCache, sliceSeqRegionId)) {
    result = IDHash_getValue(sa->asmExcCache, sliceSeqRegionId);
  }

  Vector *haps = Vector_new();
  Vector *pars = Vector_new();

  if (result != NULL) {
    int i;
    for (i=0; i<Vector_getNumElement(result); i++) {
      ExceptionCacheData *ecd = Vector_getElementAt(result,i);
  
      // need overlapping PAR and all HAPs if any
      if (!strcmp(ecd->excType,"PAR")) {
        if( ecd->seqRegionStart <= Slice_getEnd(slice) && 
            ecd->seqRegionEnd >= Slice_getStart(slice) ) {
          Vector_addElement(pars, AssemblyExceptionUnit_new(ecd->seqRegionStart, ecd->seqRegionEnd, ecd->excSeqRegionId,
                                                            ecd->excSeqRegionStart, ecd->excSeqRegionEnd));
        }
      } else {
        Vector_addElement(haps, AssemblyExceptionUnit_new(ecd->seqRegionStart, ecd->seqRegionEnd, ecd->excSeqRegionId,
                                                          ecd->excSeqRegionStart, ecd->excSeqRegionEnd));
      }
    }
  }

  Vector *out = Vector_new();

  if (!Vector_getNumElement(pars) && !Vector_getNumElement(haps)) {
    //just return this slice, there were no haps or pars
    ProjectionSegment *segment = ProjectionSegment_new(1, Slice_getLength(slice), slice);
    Vector_addElement(out, segment);
    Vector_free(haps);
    Vector_free(pars);
    return out;
  }

  Vector *syms = Vector_new();

  if ( Vector_getNumElement(haps) > 0 ) {
    Vector_sort(haps, AssemblyExceptionUnit_endCompFunc);

    AssemblyExceptionUnit **sortHaps = (AssemblyExceptionUnit **)Vector_toArray(haps);
    int nSortHap = Vector_getNumElement(haps);

    int count =0;
    long chrStart = 1;
    long hapStart = 1;
    int last = 0;

    if (sliceSeqRegionId) {
      Slice *seqRegSlice = SliceAdaptor_fetchBySeqRegionId(sa, sliceSeqRegionId, POS_UNDEF, POS_UNDEF, STRAND_UNDEF);
      Slice *excSlice    = SliceAdaptor_fetchBySeqRegionId(sa, sortHaps[0]->excSeqRegionId, POS_UNDEF, POS_UNDEF, STRAND_UNDEF);

      long len1 = Slice_getLength(seqRegSlice);
      long len2 = Slice_getLength(excSlice);

      while (count <= nSortHap && !last) { // Note goes one past end of sortHaps array 
        long chrEnd;
        long hapEnd;

        // Can we really have undefined's here???
        //if (defined($sort_haps[$count]) and defined($sort_haps[$count][0]) )
        if (count < nSortHap) {
          hapEnd = sortHaps[count]->seqRegionStart-1;
          chrEnd = sortHaps[count]->excSeqRegionStart-1;
          
        } else {
          last = 1;
          hapEnd = len1;
          chrEnd = len2;
          long diff = (hapEnd-hapStart)-(chrEnd-chrStart);
          if (diff > 0){
            AssemblyExceptionUnit *aeu = AssemblyExceptionUnit_new(hapStart, hapEnd, sortHaps[0]->excSeqRegionId, 
                                                                   chrStart, chrEnd + diff);  
            Vector_addElement(syms, aeu);
          } else if(diff < 0) {
            AssemblyExceptionUnit *aeu = AssemblyExceptionUnit_new(hapStart, hapEnd - diff, sortHaps[0]->excSeqRegionId, 
                                                                   chrStart, chrEnd);  
            Vector_addElement(syms, aeu);
          } else {
            AssemblyExceptionUnit *aeu = AssemblyExceptionUnit_new(hapStart, hapEnd, sortHaps[0]->excSeqRegionId, 
                                                                   chrStart, chrEnd);  
            Vector_addElement(syms, aeu);
          }        
          continue;
        }

        // NIY: Check if possible for hapEnd not to be set (perl was if ($hap_end and ...)
        if (hapEnd && hapStart < len1){ // if hap at start or end of chromosome
          AssemblyExceptionUnit *aeu = AssemblyExceptionUnit_new(hapStart, hapEnd, sortHaps[0]->excSeqRegionId, 
                                                                 chrStart, chrEnd);  
          Vector_addElement(syms, aeu);
        }
        chrStart = chrEnd + (sortHaps[count]->excSeqRegionEnd - sortHaps[count]->excSeqRegionStart) + 2;
        hapStart = hapEnd + (sortHaps[count]->seqRegionEnd - sortHaps[count]->seqRegionStart) + 2;
        count++;
      }
    }
  }


  // for now haps and pars should not be both there, but in theory we 
  // could handle it here by cleverly merging the pars into the existing syms,
  // for now just:
  Vector_append(syms, pars);

  Mapper *mapper = Mapper_new( "sym", "org", NULL, NULL );
  int i;
  for (i=0;i<Vector_getNumElement(syms);i++) {
    AssemblyExceptionUnit *sym = Vector_getElementAt(syms, i);
    Mapper_addMapCoordinates(mapper, 
                             sliceSeqRegionId, sym->seqRegionStart, sym->seqRegionEnd, 1, 
                             sym->excSeqRegionId, sym->excSeqRegionStart, sym->excSeqRegionEnd);
  }

  MapperRangeSet *linked = Mapper_mapCoordinates(mapper, sliceSeqRegionId, Slice_getStart(slice), Slice_getEnd(slice),
                                                 Slice_getStrand(slice), "sym" );

  // gaps are regions where there is no mapping to another region
  long relStart = 1;

  // if there was only one coord and it is a gap, we know it is just the
  // same slice with no overlapping symlinks
  // Note: For C rearranged logic slightly so can tidy up mapper etc in one place
  if (MapperRangeSet_getNumRange(linked) == 1 && MapperRangeSet_getRangeAt(linked, 0)->rangeType == MAPPERRANGE_GAP) {

    ProjectionSegment *segment = ProjectionSegment_new(1, Slice_getLength(slice), slice);
    Vector_addElement(out, segment);

  } else {
    for (i=0; i<MapperRangeSet_getNumRange(linked); i++) {
      MapperRange *coord = MapperRangeSet_getRangeAt(linked, i);
  
      if ( coord->rangeType == MAPPERRANGE_GAP) {
  
        Slice *excSlice = Slice_new(Slice_getSeqRegionName(slice), coord->start, coord->end, Slice_getStrand(slice), 
                                    Slice_getSeqRegionLength(slice), Slice_getCoordSystem(slice), sa);
  
        ProjectionSegment *segment = ProjectionSegment_new(relStart, MapperRange_getLength(coord)+relStart-1, excSlice);
        Vector_addElement(out, segment);
      } else {
        MapperCoordinate *mc = (MapperCoordinate *)coord;
  
        Slice *excSlice = SliceAdaptor_fetchBySeqRegionId(sa,  mc->id, POS_UNDEF, POS_UNDEF, STRAND_UNDEF);
  
        Slice *exc2Slice = Slice_new(Slice_getSeqRegionName(excSlice), mc->start, mc->end, mc->strand, 
                                     Slice_getSeqRegionLength(excSlice), Slice_getCoordSystem(excSlice), sa);
          
        ProjectionSegment *segment = ProjectionSegment_new(relStart, MapperRange_getLength(coord)+relStart-1, exc2Slice);
        Vector_addElement(out, segment);

        // NIY: May need to free excSlice
        Slice_free(excSlice);
      }
      relStart += MapperRange_getLength(coord);
    }
    
    if (filterProjections) {
      // Was return but do we want to just in place filter out so can do clean up before return
      SliceAdaptor_filterSliceProjections(sa, slice, out);
    }
  }

  // NIY: Tidy up before return
  MapperRangeSet_free(linked);
  Vector_free(pars);
  Vector_free(haps);
  Vector_free(syms);
  Mapper_free(mapper);

  return out;
}

/*
=head2 _filter_Slice_projections

    Arg [1]     : Bio::EnsEMBL::Slice The slice the projections were made from
    Arg [2]     : Array The projections which were fetched from the previous slice
    Description : Removes any projections which occur within the same sequence 
                  region as the given Slice object
    Returntype  : ArrayRef Bio::EnsEMBL::ProjectionSegment; Returns an array
                  of projected segments
=cut
*/

Vector *SliceAdaptor_filterSliceProjections(SliceAdaptor *sa, Slice *slice, Vector *projections) {
  if ( ! Vector_getNumElement(projections) ) {
    fprintf(stderr, "Was not given any projections to filter. Database may have incorrect assembly_exception information loaded\n");
    exit(1);
  }
  
  // Want to get features on the FULL original slice as well as any
  // symlinked slices.
  
  // Filter out partial slices from projection that are on same
  // seq_region as original slice.

  IDType srId = Slice_getSeqRegionId(slice);
  int i;

  for (i=0; i<Vector_getNumElement(projections); i++) {
    ProjectionSegment *segment = Vector_getElementAt(projections,i);
  
    if (Slice_getSeqRegionId(ProjectionSegment_getToSlice(segment)) == srId) {
      Vector_removeElementAt(projections, i);
      i--;
    }
  }
      
  ProjectionSegment *segment = ProjectionSegment_new(1,Slice_getLength(slice), slice);

  Vector_addElement(projections, segment);

  return projections;
}


/* Not implementing store functions yet
=head2 store

  Arg [1]    : Bio::EnsEMBL::Slice $slice
  Arg [2]    : (optional) $seqref reference to a string
               The sequence associated with the slice to be stored.
  Example    : $slice = Bio::EnsEMBL::Slice->new(...);
               $seq_region_id = $slice_adaptor->store($slice, \$sequence);
  Description: This stores a slice as a sequence region in the database
               and returns the seq region id. The passed in slice must
               start at 1, and must have a valid seq_region name and coordinate
               system. The attached coordinate system must already be stored in
               the database.  The sequence region is assumed to start at 1 and
               to have a length equalling the length of the slice.  The end of
               the slice must equal the seq_region_length.
               If the slice coordinate system is the sequence level coordinate
               system then the seqref argument must also be passed.  If the
               slice coordinate system is NOT a sequence level coordinate
               system then the sequence argument cannot be passed.
  Returntype : int 
  Exceptions : throw if slice has no coord system.
               throw if slice coord system is not already stored.
               throw if slice coord system is seqlevel and no sequence is 
                     provided.
               throw if slice coord system is not seqlevel and sequence is
                     provided.
               throw if slice does not start at 1
               throw if sequence is provided and the sequence length does not
                     match the slice length.
               throw if the SQL insert fails (e.g. on duplicate seq region)
               throw if slice argument is not passed
               throw if the slice end is not equal to seq_region_length
  Caller     : database loading scripts
  Status     : Stable

=cut

sub store {
  my $self = shift;
  my $slice = shift;
  my $seqref = shift;

  #
  # Get all of the sanity checks out of the way before storing anything
  #

  if(!ref($slice) || !($slice->isa('Bio::EnsEMBL::Slice') or $slice->isa('Bio::EnsEMBL::LRGSlice'))) {
    throw('Slice argument is required');
  }

  my $cs = $slice->coord_system();
  throw("Slice must have attached CoordSystem.") if(!$cs);

  my $db = $self->db();
  if(!$cs->is_stored($db)) {
    throw("Slice CoordSystem must already be stored in DB.") 
  }

  if($slice->start != 1 || $slice->strand != 1) {
    throw("Slice must have start==1 and strand==1.");
  }

  if($slice->end() != $slice->seq_region_length()) {
    throw("Slice must have end==seq_region_length");
  }

  my $sr_len = $slice->length();
  my $sr_name  = $slice->seq_region_name();

  if(!$sr_name) {
    throw("Slice must have valid seq region name.");
  }

  if($cs->is_sequence_level()) {
    if(!$seqref) {
      throw("Must provide sequence for sequence level coord system.");
    }
    if(ref($seqref) ne 'SCALAR') {
      throw("Sequence must be a scalar reference.");
    }
    my $seq_len = length($$seqref);

    if($seq_len != $sr_len) {
      throw("Sequence length ($seq_len) must match slice length ($sr_len).");
    }
  } else {
    if($seqref) {
      throw("Cannot provide sequence for non-sequence level seq regions.");
    }
  }

  #store the seq_region

  my $sth = $db->dbc->prepare("INSERT INTO seq_region " .
                         "SET    name = ?, " .
                         "       length = ?, " .
                         "       coord_system_id = ?" );

  $sth->bind_param(1,$sr_name,SQL_VARCHAR);
  $sth->bind_param(2,$sr_len,SQL_INTEGER);
  $sth->bind_param(3,$cs->dbID,SQL_INTEGER);

  $sth->execute();

  my $seq_region_id = $sth->{'mysql_insertid'};

  if(!$seq_region_id) {
    throw("Database seq_region insertion failed.");
  }

  if($cs->is_sequence_level()) {
    #store sequence if it was provided
    my $seq_adaptor = $db->get_SequenceAdaptor();
    $seq_adaptor->store($seq_region_id, $$seqref);
  }

  #synonyms
  if(defined($slice->{'synonym'})){
    foreach my $syn (@{$slice->{'synonym'}} ){
      $syn->seq_region_id($seq_region_id); # set the seq_region_id
      $syn->adaptor->store($syn);
    }
  }
  
  
  $slice->adaptor($self);

  return $seq_region_id;
}


=head2 store_assembly

  Arg [1]    : Bio::EnsEMBL::Slice $asm_slice
  Arg [2]    : Bio::EnsEMBL::Slice $cmp_slice
  Example    : $asm = $slice_adaptor->store_assembly( $slice1, $slice2 );
  Description: Creates an entry in the analysis table based on the 
               coordinates of the two slices supplied. Returns a string 
               representation of the assembly that gets created.
  Returntype : string
  Exceptions : throw if either slice has no coord system (cs).
               throw unless the cs rank of the asm_slice is lower than the 
               cmp_slice.
               throw if there is no mapping path between coord systems
               throw if the lengths of each slice are not equal
               throw if there are existing mappings between either slice
               and the oposite cs
  Caller     : database loading scripts
  Status     : Experimental

=cut

sub store_assembly{
  my $self = shift;
  my $asm_slice = shift;
  my $cmp_slice = shift;

  #
  # Get all of the sanity checks out of the way before storing anything
  #

  if(!ref($asm_slice) || !($asm_slice->isa('Bio::EnsEMBL::Slice') or $asm_slice->isa('Bio::EnsEMBL::LRGSlice'))) {
    throw('Assembled Slice argument is required');
  }
  if(!ref($cmp_slice) || !($cmp_slice->isa('Bio::EnsEMBL::Slice') or $cmp_slice->isa('Bio::EnsEMBL::LRGSlice')) ) {
    throw('Assembled Slice argument is required');
  }

  my $asm_cs = $asm_slice->coord_system();
  throw("Slice must have attached CoordSystem.") if(!$asm_cs);
  my $cmp_cs = $cmp_slice->coord_system();
  throw("Slice must have attached CoordSystem.") if(!$cmp_cs);

  unless( $asm_cs->rank < $cmp_cs->rank ){
    throw("Assembled Slice CoordSystem->rank must be lower than ".
          "the component Slice Coord_system" );
  }

  my @path =
    @{ $asm_cs->adaptor()->get_mapping_path( $asm_cs, $cmp_cs ) };

  if ( !@path ) {
    throw("No mapping path defined between "
        . $asm_cs->name() . " and "
        . $cmp_cs->name() );
  }

  if( $asm_slice->length != $cmp_slice->length ){
    throw("The lengths of the assembled and component slices are not equal" );
  }

  # For now we disallow any existing mappings between the asm slice and cmp
  # CoordSystem and vice-versa. 
  # Some cases of multiple mappings may be allowable by the API, but their 
  # logic needs to be coded below.

  my $asm_proj = $asm_slice->project( $cmp_cs->name, $cmp_cs->version );
  if( @$asm_proj ){
    throw("Regions of the assembled slice are already assembled ".
          "into the component CoordSystem" ); 
  }
  my $cmp_proj = $cmp_slice->project( $asm_cs->name, $asm_cs->version );
  if( @$cmp_proj ){
    throw("Regions of the component slice are already assembled ".
          "into the assembled CoordSystem" ); 
  }

  #
  # Checks complete. Store the data
  #
  my $sth = $self->db->dbc->prepare
      ("INSERT INTO assembly " .
       "SET     asm_seq_region_id = ?, " .
       "        cmp_seq_region_id = ?, " .
       "        asm_start = ?, " .
       "        asm_end   = ?, " .
       "        cmp_start = ?, " .
       "        cmp_end   = ?, " .
       "        ori       = ?" );

  my $asm_seq_region_id = $self->get_seq_region_id( $asm_slice );
  my $cmp_seq_region_id = $self->get_seq_region_id( $cmp_slice );
  my $ori = $asm_slice->strand * $cmp_slice->strand;

  $sth->bind_param(1,$asm_seq_region_id,SQL_INTEGER);
  $sth->bind_param(2,$cmp_seq_region_id,SQL_INTEGER);
  $sth->bind_param(3,$asm_slice->start,SQL_INTEGER);
  $sth->bind_param(4,$asm_slice->end,SQL_INTEGER);
  $sth->bind_param(5,$cmp_slice->start,SQL_INTEGER);
  $sth->bind_param(6,$cmp_slice->end,SQL_INTEGER);
  $sth->bind_param(7,$ori,SQL_INTEGER);

  $sth->execute();

  #use Data::Dumper qw( Dumper );
  #warn Dumper( $self->db->{seq_region_cache} );
  #$self->db->{seq_region_cache} = undef;
  #$self->_cache_seq_regions();

  my $ama = $self->db->get_AssemblyMapperAdaptor();
  $ama->delete_cache();


  return $asm_slice->name . "<>" . $cmp_slice->name;

}
*/

/*
=head2 prepare

  Arg [1]    : String $sql
  Example    :  ( optional )
  Description: overrides the default adaptor prepare method.
               All slice sql will usually use the dna_db.
  Returntype : DBD::sth 
  Exceptions : none
  Caller     : internal, convenience method
  Status     : Stable

=cut
*/

StatementHandle *SliceAdaptor_prepare(BaseAdaptor *ba, char *qStr, size_t len) {
  //printf("Query = %s len = %ld\n",qStr,len);
  return DBAdaptor_prepare(ba->dba->dnadb,qStr,len);
}


void SliceAdaptor_buildExceptionCache(SliceAdaptor *sa) {

  // build up a cache of the entire assembly exception table
  // it should be small anyway
  char qStr[1024];
  sprintf(qStr, "SELECT ae.seq_region_id, ae.seq_region_start, "
                    "ae.seq_region_end, ae.exc_type, ae.exc_seq_region_id, "
                    "ae.exc_seq_region_start, ae.exc_seq_region_end "
                    "FROM assembly_exception ae, "
                    "seq_region sr, coord_system cs "
                    "WHERE sr.seq_region_id = ae.seq_region_id "
                    "AND sr.coord_system_id = cs.coord_system_id "
                    "AND cs.species_id = "IDFMTSTR, SliceAdaptor_getSpeciesID(sa));

  StatementHandle *sth = sa->prepare((BaseAdaptor *)sa,qStr,strlen(qStr));

  sth->execute(sth);

  sa->asmExcCache = IDHash_new(IDHASH_SMALL);

  ResultRow *row;
  while ((row = sth->fetchRow(sth))) {
    IDType seqRegionId       = row->getLongLongAt(row,0);
    long   seqRegionStart    = row->getLongAt(row,1);
    long   seqRegionEnd      = row->getLongAt(row,2);
    char * excType           = row->getStringAt(row,3);
    IDType excSeqRegionId    = row->getLongLongAt(row,4);
    long   excSeqRegionStart = row->getLongAt(row,5);
    long   excSeqRegionEnd   = row->getLongAt(row,6);

    ExceptionCacheData *ecd = ExceptionCacheData_new(seqRegionId, seqRegionStart, seqRegionEnd, excType, excSeqRegionId, excSeqRegionStart, excSeqRegionEnd);

    if (! IDHash_contains(sa->asmExcCache, seqRegionId)) {
      IDHash_add(sa->asmExcCache, seqRegionId, Vector_new()); 
    }
    Vector *vec = IDHash_getValue(sa->asmExcCache, seqRegionId);
    Vector_addElement(vec, ecd);
  }

  sth->finish(sth);
}

/*
=head2 cache_toplevel_seq_mappings

  Args       : none
  Example    : $slice_adaptor->cache_toplevel_seq_mappings();
  Description: caches all the assembly mappings needed for genes
  Returntype : None
  Exceptions : None
  Caller     : general
  Status     : At Risk
             : New experimental code

=cut
*/

void SliceAdaptor_cacheTopLevelSeqMappings(SliceAdaptor *sa) {

  // Get the sequence level to map to

  char qStr[1024];
  sprintf(qStr,"SELECT name "
               "FROM  coord_system "
               "WHERE attrib like '%%%%sequence_level%%%%' "
               "AND   species_id = "IDFMTSTR, SliceAdaptor_getSpeciesID(sa));

  StatementHandle *sth = sa->prepare((BaseAdaptor *)sa,qStr,strlen(qStr));
  sth->execute(sth);

  ResultRow *row = sth->fetchRow(sth);
// think I need to make a copy of this
  char *sequenceLevel = StrUtil_copyString(&sequenceLevel, row->getStringAt(row,0), 0);

  sth->finish(sth);

  CoordSystemAdaptor *   csa = DBAdaptor_getCoordSystemAdaptor(sa->dba);
  AssemblyMapperAdaptor *ama = DBAdaptor_getAssemblyMapperAdaptor(sa->dba);

  CoordSystem *cs1 = CoordSystemAdaptor_fetchByName(csa, sequenceLevel, NULL);

  //get level to map to

  sprintf(qStr, "SELECT DISTINCT(cs.name) "
                "FROM  seq_region sr, "
                "      seq_region_attrib sra, "
                "      attrib_type at, "
                "      coord_system cs "
                "WHERE sra.seq_region_id = sr.seq_region_id "
                "AND   sra.attrib_type_id = at.attrib_type_id "
                "AND   at.code = 'toplevel' "
                "AND   cs.coord_system_id = sr.coord_system_id "
                "AND   cs.species_id = "IDFMTSTR, SliceAdaptor_getSpeciesID(sa));

  sth = sa->prepare((BaseAdaptor *)sa,qStr,strlen(qStr));
  sth->execute(sth);

// csn in perl is now csName
  while ((row = sth->fetchRow(sth))) {
    char *csName = row->getStringAt(row,0);

    if (strcmp(csName,sequenceLevel)) { // if not sequence level
      CoordSystem *cs2 = CoordSystemAdaptor_fetchByName(csa, sequenceLevel, NULL);
      AssemblyMapper *am = AssemblyMapperAdaptor_fetchByCoordSystems(ama, cs1, cs2);
      AssemblyMapper_registerAll(am);
    }
  }

  free(sequenceLevel);
}


/* No circular stuff
sub _build_circular_slice_cache {
  my $self = shift;

  # build up a cache of circular sequence region ids
  my $sth =
            $self->prepare( "SELECT sra.seq_region_id FROM seq_region_attrib sra "
                          . "INNER JOIN attrib_type at ON sra.attrib_type_id = at.attrib_type_id "
                        . "INNER JOIN seq_region sr ON sra.seq_region_id = sr.seq_region_id "
                        . "INNER JOIN coord_system cs ON sr.coord_system_id = cs.coord_system_id "
                        . "WHERE code = 'circular_seq' and cs.species_id = ?");

  $sth->bind_param( 1, $self->species_id(), SQL_INTEGER );
  $sth->execute();

  my $id;
  my %hash;
  if ( ($id) = $sth->fetchrow_array() ) {
          $self->{'circular_sr_id_cache'} = \%hash;
        $self->{'is_circular'} = 1;
        $hash{ $id } = $id;
         while ( ($id) = $sth->fetchrow_array() ) {
                    $hash{ $id } = $id;
          }
  } else {
        $self->{'is_circular'} = 0;
  }
  $sth->finish();
} ## end _build_circular_slice_cache

*/

/* Not implementing for now
#####################################
# sub DEPRECATED METHODs
#####################################

=head2 fetch_by_mapfrag

 Function: DEPRECATED use fetch_by_misc_feature_attribute('synonym',$mapfrag)

=cut

sub fetch_by_mapfrag{
   my ($self,$mymapfrag,$flag,$size) = @_;
   deprecate('Use fetch_by_misc_feature_attribute instead');
   $flag ||= 'fixed-width'; # alt.. 'context'
   $size ||= $flag eq 'fixed-width' ? 100000 : 0;
   return $self->fetch_by_misc_feature_attribute('synonym',$mymapfrag,$size);
}



=head2 fetch_by_chr_start_end

  Description: DEPRECATED use fetch_by_region instead

=cut

sub fetch_by_chr_start_end {
  my ($self,$chr,$start,$end) = @_;
  deprecate('Use fetch_by_region() instead');

  #assume that by chromosome the user actually meant top-level coord
  #system since this is the old behaviour of this deprecated method
  my $csa = $self->db->get_CoordSystemAdaptor();
  my ($cs) = @{$csa->fetch_all()}; # get the highest coord system

  return $self->fetch_by_region($cs->name,$chr,$start,$end,1,$cs->version);
}



=head2 fetch_by_contig_name

  Description: Deprecated. Use fetch_by_region(), Slice::project(),
               Slice::expand() instead

=cut

sub fetch_by_contig_name {
  my ($self, $name, $size) = @_;

  deprecate('Use fetch_by_region(), Slice::project() and Slice::expand().');

  #previously wanted chromosomal slice on a given contig.  Assume this means
  #a top-level slice on a given seq_region in the seq_level coord system
  my $csa = $self->db()->get_CoordSystemAdaptor();
  my $seq_level = $csa->fetch_sequence_level();

  my $seq_lvl_slice = $self->fetch_by_region($seq_level->name(), $name);

  if(!$seq_lvl_slice) {
    return undef;
  }

  my @projection = @{$seq_lvl_slice->project('toplevel')};

  if(@projection != 1) {
    warning("$name is mapped to multiple toplevel locations.");
  }

  return $projection[0]->[2]->expand($size, $size);
}


=head2 fetch_by_clone_accession

  Description: DEPRECATED.  Use fetch_by_region, Slice::project, Slice::expand
               instead.

=cut

sub fetch_by_clone_accession{
  my ($self,$name,$size) = @_;

  deprecate('Use fetch_by_region(), Slice::project() and Slice::expand().');

  my $csa = $self->db()->get_CoordSystemAdaptor();
  my $clone_cs = $csa->fetch_by_name('clone');

  if(!$clone_cs) {
    warning('Clone coordinate system does not exist for this species');
    return undef;
  }

  #this unfortunately needs a version on the end to work
  if(! ($name =~ /\./)) {
    my $sth =
      $self->prepare(  "SELECT sr.name "
                     . "FROM   seq_region sr, coord_system cs "
                     . "WHERE  cs.name = 'clone' "
                     . "AND    cs.coord_system_id = sr.coord_system_id "
                     . "AND    sr.name LIKE '$name.%'"
                     . "AND    cs.species_id = ?" );

    $sth->bind_param( 1, $self->species_id(), SQL_INTEGER );
    $sth->execute();

    if(!$sth->rows()) {
      $sth->finish();
      throw("Clone $name not found in database");
    }

    ($name) = $sth->fetchrow_array();

    $sth->finish();
  }

  my $clone = $self->fetch_by_region($clone_cs->name(), $name);
  return undef if(!$clone);

  my @projection = @{$clone->project('toplevel')};

  if(@projection != 1) {
    warning("$name is mapped to multiple locations.");
  }

  return $projection[0]->[2]->expand($size, $size);
}


=head2 fetch_by_supercontig_name

  Description: DEPRECATED. Use fetch_by_region(), Slice::project() and
               Slice::expand() instead

=cut

sub fetch_by_supercontig_name {
  my ($self,$name, $size) = @_;

  deprecate('Use fetch_by_region(), Slice::project() and Slice::expand().');

  my $csa = $self->db()->get_CoordSystemAdaptor();
  my $sc_level = $csa->fetch_by_name('supercontig');

  if(!$sc_level) {
    warning('No supercontig coordinate system exists for this species.');
    return undef;
  }

  my $sc_slice = $self->fetch_by_region($sc_level->name(),$name);

  return undef if(!$sc_slice);

  my @projection = @{$sc_slice->project('toplevel')};

  if(@projection > 1) {
    warning("$name is mapped to multiple locations in toplevel");
  }

  return $projection[0]->[2]->expand($size, $size);
}




=head2 list_overlapping_supercontigs

  Description: DEPRECATED use Slice::project instead

=cut

sub list_overlapping_supercontigs {
   my ($self,$slice) = @_;

   deprecate('Use Slice::project() instead.');

   my $csa = $self->db()->get_CoordSystemAdaptor();
   my $sc_level = $csa->fetch_by_name('supercontig');

   if(!$sc_level) {
     warning('No supercontig coordinate system exists for this species.');
     return undef;
   }

   my @out;
   foreach my $seg (@{$slice->project($sc_level->name(), $sc_level->version)}){
     push @out, $seg->[2]->seq_region_name();
   }

   return \@out;
}


*/

/*
=head2 fetch_by_chr_name

  Description: DEPRECATED. Use fetch by region instead

=cut
*/

Slice *SliceAdaptor_fetchByChrName(SliceAdaptor *sa, char *chrName) {
  fprintf(stderr,"deprecated: Use fetch_by_region() instead.\n");

  CoordSystemAdaptor *csa = DBAdaptor_getCoordSystemAdaptor(sa->dba);

  CoordSystem *topCs = Vector_getElementAt(CoordSystemAdaptor_fetchAll(csa), 0);

  return SliceAdaptor_fetchByRegion(sa, CoordSystem_getName(topCs), chrName, POS_UNDEF, POS_UNDEF, STRAND_UNDEF,
                                    CoordSystem_getVersion(topCs), 0);
}

