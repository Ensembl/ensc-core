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

#include "Mapper.h"
#include "MapperCoordinate.h"
#include "IndelCoordinate.h"
#include "StrUtil.h"
#include <stdio.h>
#include <string.h>

/*
=head2 new

  Arg [1]    : string $from
               The name of the 'from' coordinate system
  Arg [2]    : string $to
               The name of the 'to' coordinate system
  Arg [3]    : (optional) Bio::EnsEMBL::CoordSystem $from_cs
               The 'from' coordinate system
  Arg [4]    : (optional) Bio::EnsEMBL::CoordSystem $to_cs
  Example    : my $mapper = Bio::EnsEMBL::Mapper->new('FROM', 'TO');
  Description: Constructor.  Creates a new Bio::EnsEMBL::Mapper object.
  Returntype : Bio::EnsEMBL::Mapper
  Exceptions : none
  Caller     : general

=cut
*/


Mapper *Mapper_new(char *from, char *to, CoordSystem *fromCs, CoordSystem *toCs) {
  Mapper *m;

  if ((m = (Mapper *)calloc(1,sizeof(Mapper))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for mapper\n");
    exit(1);
  }

  if ( to == NULL || from == NULL ) {
    fprintf(stderr, "Must supply 'to' and 'from' tags");
    exit(1);
  }

  Mapper_setFrom(m, from);
  Mapper_setTo(m, to);

  Mapper_setFromCoordSystem(m, fromCs);
  Mapper_setToCoordSystem(m, toCs);

  Mapper_setPairCount(m, 0);

  Mapper_setPairHash(m, MAPPER_FROM_IND, IDHash_new(IDHASH_MEDIUM));
  Mapper_setPairHash(m, MAPPER_TO_IND,   IDHash_new(IDHASH_MEDIUM));

  return m;
}

char *Mapper_setFrom(Mapper *m, char *from) {
  return StrUtil_copyString(&(m->from), from, 0);
  //return m->from;
}

char *Mapper_setTo(Mapper *m, char *to) {
  return StrUtil_copyString(&(m->to), to, 0);
}

/*
=head2 flush

  Args       : none
  Example    : none
  Description: removes all cached information out of this mapper
  Returntype : none
  Exceptions : none
  Caller     : AssemblyMapper, ChainedAssemblyMapper

=cut
*/

void Mapper_flush(Mapper *m)  {
// NIY: What should free func be for these
  IDHash_free(Mapper_getPairHash(m, MAPPER_FROM_IND), NULL);
  IDHash_free(Mapper_getPairHash(m, MAPPER_TO_IND), NULL);

  Mapper_setPairHash(m, MAPPER_FROM_IND, IDHash_new(IDHASH_MEDIUM));
  Mapper_setPairHash(m, MAPPER_TO_IND,   IDHash_new(IDHASH_MEDIUM));

  Mapper_setPairCount(m, 0);
}

/*
=head2 map_coordinates

    Arg  1      string $id
                id of 'source' sequence
    Arg  2      int $start
                start coordinate of 'source' sequence
    Arg  3      int $end
                end coordinate of 'source' sequence
    Arg  4      int $strand
                raw contig orientation (+/- 1)
    Arg  5      string $type
                nature of transform - gives the type of
                coordinates to be transformed *from*
    Function    generic map method
    Returntype  array of Bio::EnsEMBL::Mapper::Coordinate
                and/or   Bio::EnsEMBL::Mapper::Gap
    Exceptions  none
    Caller      Bio::EnsEMBL::Mapper

=cut
*/

MapperRangeSet *Mapper_mapCoordinates(Mapper *m, IDType id, long start, long end, int strand, char *type) {


  // special case for handling inserts:
  if ( start == end+1 ) {
    return Mapper_mapInsert(m, id, start, end, strand, type, 0 /*fastmap flag */);
  } else if (start > end+1) {
    fprintf(stderr,"ERROR: Start is greater than end for id " IDFMTSTR ", start %ld, end %ld\n",id,start,end);
    exit(1);
  }


  IDHash *hash;
  CoordSystem *cs;
  if( Mapper_getIsSorted(m) == 0 ) {
    Mapper_sort(m);
  }

  int from, to;
  if (!Mapper_compareType(type,Mapper_getTo(m))) {
    from = MAPPER_TO_IND;
    to   = MAPPER_FROM_IND;
    cs   = Mapper_getFromCoordSystem(m);
  } else if (!Mapper_compareType(type,Mapper_getFrom(m))) {
    from = MAPPER_FROM_IND;
    to   = MAPPER_TO_IND;
    cs   = Mapper_getToCoordSystem(m);
  } else {
    fprintf(stderr, "Invalid type [%s] in mapper (not from [%s] or to [%s])\n", type, Mapper_getFrom(m), Mapper_getTo(m));
    exit(1);
  }

  hash = Mapper_getPairHash(m, from);

  if (!hash) {
    fprintf(stderr,"ERROR: Type %s is neither to or from coordinate system\n",type);
    exit(1);
  }


  MapperRangeSet *results = MapperRangeSet_new();


// Was upcasing the id - its a number in C, I haven't found a case yet where its a string
  if (!IDHash_contains(hash, id)) {
    // one big gap!
    MapperRange *gap = (MapperRange *)MapperGap_new(start,end,0); // Perl didn't set rank so use 0
    MapperRangeSet_addRange(results,gap);
    return results;
  }


  MapperPairSet *pairs = IDHash_getValue(hash,id); //my $lr = $hash->{ uc($id) };


  MapperPair *lastUsedPair = NULL;

  int startIdx, endIdx, midIdx;
  MapperPair *pair;
  MapperUnit *selfCoord;


  startIdx = 0;
  endIdx   = MapperPairSet_getNumPair(pairs)-1;

  // binary search the relevant pairs
  // helps if the list is big
  while ( ( endIdx - startIdx ) > 1 ) {
    midIdx = ( startIdx + endIdx ) >> 1;
    pair   = MapperPairSet_getPairAt(pairs, midIdx);

    selfCoord = MapperPair_getUnit(pair, from);

    if ( selfCoord->end < start ) {
      startIdx = midIdx;
    } else {
      endIdx = midIdx;
    }
  }

  int rank       = 0;
  long origStart = start;
  IDType lastTargetCoord;
  int lastTargetCoordIsSet = 0;

  int i;
  for (i=startIdx; i<MapperPairSet_getNumPair(pairs); i++) {
    MapperPair *pair = MapperPairSet_getPairAt(pairs,i);
    MapperUnit *selfCoord   = MapperPair_getUnit(pair, from);
    MapperUnit *targetCoord = MapperPair_getUnit(pair, to);

    //
    // But not the case for haplotypes!! need to test for this case???
    // so removing this till a better solution is found
    // 
    //
    //     if($self_coord->{'start'} < $start){
    //       $start = $orig_start;
    //       $rank++;
    //     }

    if ( lastTargetCoordIsSet && targetCoord->id != lastTargetCoord ) {
      if ( selfCoord->start < start ) {    // i.e. the same bit is being mapped to another assembled bit
        start = origStart;
      }
    } else {
      lastTargetCoord = targetCoord->id;
      lastTargetCoordIsSet = 1;
    }

    // if we haven't even reached the start, move on
    if (selfCoord->end < origStart) {
      continue;
    }

    // if we have over run, break
    if (selfCoord->start > end) {
      break;
    }

// Check is start not origStart
    if (start < selfCoord->start) {
      // gap detected
      MapperRange *gap = (MapperRange *)MapperGap_new(start, selfCoord->start-1, rank);
      MapperRangeSet_addRange(results,gap);
      start = gap->end+1;
    }

    long targetStart, targetEnd;

    MapperRange *res;

    if ( MapperPair_isIndel(pair) ) {
      // When next pair is an IndelPair and not a Coordinate, create the
      // new mapping Coordinate, the IndelCoordinate.
      targetStart = targetCoord->start;
      targetEnd   = targetCoord->end;

      // create a Gap object
      MapperGap *gap = MapperGap_new(start,
                                     selfCoord->end < end ? selfCoord->end : end,
                                     0); // Perl didn't set rank - don't know if need to

      // create the Coordinate object
      MapperCoordinate *coord = MapperCoordinate_new(targetCoord->id,
                                                     targetStart,
                                                     targetEnd,
                                                     pair->ori * strand, 
                                                     cs, 
                                                     0); // Perl didn't set rank - don't know if need to

      //and finally, the IndelCoordinate object with
      res = (MapperRange *)IndelCoordinate_new(gap, coord);
    } else {
      // start is somewhere inside the region
      if (pair->ori == 1) {
        targetStart = targetCoord->start + (start - selfCoord->start);
      } else {
        targetEnd = targetCoord->end - (start - selfCoord->start);
      }

      // Either we are enveloping this map or not.  If yes, then end
      // point (self perspective) is determined solely by target.  If
      // not we need to adjust.
      if (end > selfCoord->end) {
        // enveloped
        if( pair->ori == 1 ) {
          targetEnd = targetCoord->end;
        } else {
          targetStart = targetCoord->start;
        }
      } else {
        // need to adjust end
        if (pair->ori == 1) {
          targetEnd = targetCoord->start + (end - selfCoord->start);
        } else {
          targetStart = targetCoord->end - (end - selfCoord->start);
        }
      }

      res = (MapperRange *)MapperCoordinate_new(targetCoord->id,
                                                targetStart,
                                                targetEnd,
                                                pair->ori * strand, 
                                                cs, 
                                                rank);
    } // end else [ if ( exists $pair->{'indel'...})]

    MapperRangeSet_addRange(results, res);

    lastUsedPair = pair;
    start = selfCoord->end+1;
  }

  if (lastUsedPair == NULL) {
    MapperRange *gap = (MapperRange *)MapperGap_new(start,end, 0); // Perl doesn't set rank, so use 0
    MapperRangeSet_addRange(results,gap);

  } else if (MapperPair_getUnit(lastUsedPair, from)->end < end) {
    // gap at the end
    MapperRange *gap = (MapperRange *)MapperGap_new(
                           MapperPair_getUnit(lastUsedPair,from)->end + 1,
                           end, 0 ); // Perl didn't set rank so use 0
    MapperRangeSet_addRange(results,gap);
  }

  if (strand == -1) {
    MapperRangeSet_reverse(results);
  }

  return results;
}


/*
=head2 map_insert

  Arg [1]    : string $id
  Arg [2]    : int $start - start coord. Since this is an insert should always
               be one greater than end.
  Arg [3]    : int $end - end coord. Since this is an insert should always
               be one less than start.
  Arg [4]    : int $strand (0, 1, -1)
  Arg [5]    : string $type - the coordinate system name the coords are from.
  Arg [6]    : boolean $fastmap - if specified, this is being called from
               the fastmap call. The mapping done is not any faster for
               inserts, but the return value is different.
  Example    : 
  Description: This is in internal function which handles the special mapping
               case for inserts (start = end +1).  This function will be called
               automatically by the map function so there is no reason to
               call it directly.
  Returntype : list of Bio::EnsEMBL::Mapper::Coordinate and/or Gap objects
  Exceptions : none
  Caller     : map_coordinates()

=cut
*/

MapperRangeSet *Mapper_mapInsert(Mapper *m, IDType id, long start, long end, int strand, char *type, int fastmap) {

  // swap start/end and map the resultant 2bp coordinate
  long tmp;
  tmp = start;
  start = end;
  end = tmp;
  
  MapperRangeSet *coords = Mapper_mapCoordinates(m, id, start, end, strand, type);

  MapperRangeSet *retSet = MapperRangeSet_new(); 

  if (MapperRangeSet_getNumRange(coords) == 1) {
// Note assuming its a mapper coordinate - not sure if this is always true
    MapperCoordinate *c = (MapperCoordinate *)MapperRangeSet_getRangeAt(coords,0);

    // swap start and end to convert back into insert
    MapperCoordinate *newC = MapperCoordinate_new(c->id, c->end, c->start, c->strand, c->coordSystem, 0); // Perl didn't set rank so use 0
    MapperRangeSet_addRange(retSet, (MapperRange *)newC);

  } else {
    if (MapperRangeSet_getNumRange(coords) != 2) {
      fprintf(stderr, "Unexpected: Got %d expected 2.\n", MapperRangeSet_getNumRange(coords));
      exit(1);
    }

    // adjust coordinates, remove gaps
    MapperRange *c1, *c2;

    if (strand == -1) {
      c1 = MapperRangeSet_getRangeAt(coords,1);
      c2 = MapperRangeSet_getRangeAt(coords,0);
    } else {
      c1 = MapperRangeSet_getRangeAt(coords,0);
      c2 = MapperRangeSet_getRangeAt(coords,1);
    }

    //@coords = ();

    // Was ref($c1) eq 'Bio::EnsEMBL::Mapper::Coordinate' so think this WON'T include MAPPERRANGE_INDEL

    MapperCoordinate *newC1 = NULL;
    if (c1->rangeType == MAPPERRANGE_COORD) {
      MapperCoordinate *mcC1 = (MapperCoordinate *)c1;
      // insert is after first coord
      if (mcC1->strand * strand == -1) {
        mcC1->end--;
      } else {
        mcC1->start++;
      }
      newC1 = MapperCoordinate_new(mcC1->id, mcC1->start, mcC1->end, mcC1->strand, mcC1->coordSystem, 0); // Perl didn't set rank, so use 0
    }

    // (see above for note on this condition if (ref($c2) eq 'Bio::EnsEMBL::Mapper::Coordinate') 
    MapperCoordinate *newC2 = NULL;
    if (c2->rangeType == MAPPERRANGE_COORD) {
      MapperCoordinate *mcC2 = (MapperCoordinate *)c2;
      // insert is before second coord
      if(mcC2->strand * strand == -1) {
        mcC2->start++;
      } else {
        mcC2->end--;
      }
      newC2 = MapperCoordinate_new(mcC2->id, mcC2->start, mcC2->end, mcC2->strand, mcC2->coordSystem, 0); // Perl didn't set rank, so use 0
    }
   
    if (strand == -1) { // Add in 2, 1 order
      if (newC2) MapperRangeSet_addRange(retSet, (MapperRange *)newC2);
      if (newC1) MapperRangeSet_addRange(retSet, (MapperRange *)newC1);
    } else {  // Add in 1, 2 order
      if (newC1) MapperRangeSet_addRange(retSet, (MapperRange *)newC1);
      if (newC2) MapperRangeSet_addRange(retSet, (MapperRange *)newC2);
    }

  }

  if (fastmap) {
    if (MapperRangeSet_getNumRange(coords) != 1) {
      MapperRangeSet_free(retSet);
      retSet = NULL;
    } else {
      MapperRange *c = MapperRangeSet_getRangeAt(coords, 0);

      // Type check - not sure what we're supposed to have here so belt and braces seem appropriate
      if (c->rangeType != MAPPERRANGE_COORD) {
        fprintf(stderr, "Expected a Coordinate range, got a %d range (look up in MapperRange.h)\n", c->rangeType);
        exit(1);
      }
      
      MapperCoordinate *mcC = (MapperCoordinate *)c;

      MapperCoordinate *newC = MapperCoordinate_new(mcC->id, mcC->start, mcC->end, mcC->strand, mcC->coordSystem, 0); // Perl didn't set rank, so use 0
      MapperRangeSet_addRange(retSet, (MapperRange *)newC);

/*
      return ($c->{'id'}, $c->{'start'}, $c->{'end'},
              $c->{'strand'}, $c->{'coord_system'});
*/
    }
  }

  MapperRangeSet_free(coords);

  return retSet;
}

/*
=head2 fastmap

    Arg  1      string $id
                id of 'source' sequence
    Arg  2      int $start
                start coordinate of 'source' sequence
    Arg  3      int $end
                end coordinate of 'source' sequence
    Arg  4      int $strand
                raw contig orientation (+/- 1)
    Arg  5      int $type
                nature of transform - gives the type of
                coordinates to be transformed *from*
    Function    inferior map method. Will only do ungapped unsplit mapping.
                Will return id, start, end strand in a list.
    Returntype  list of results
    Exceptions  none
    Caller      Bio::EnsEMBL::AssemblyMapper

=cut
*/
// NIY: May need some reworking to handle mapInsert because I'd changed the way it returns data
// Change back to returning MapperRangeSet
MapperRangeSet *Mapper_fastMap(Mapper *m, IDType id, long start, long end, int strand, char *type) {
  MapperPairSet *pairs;
  int i;
  IDHash *hash;
  int from, to;
  CoordSystem *cs;

  if(end+1 == start) {
    return Mapper_mapInsert(m, id, start, end, strand, type, 1);
  }

  if(!Mapper_compareType(type, Mapper_getTo(m))) {
    from = MAPPER_TO_IND;
    to   = MAPPER_FROM_IND;
    cs   = Mapper_getFromCoordSystem(m);

  } else {
    from = MAPPER_FROM_IND;
    to   = MAPPER_TO_IND;
    cs   = Mapper_getToCoordSystem(m);
  }

  hash = Mapper_getPairHash(m, from);

  if (!hash) {
    fprintf(stderr,"ERROR: Type %s is neither to or from coordinate system\n",type);
    exit(1);
  }

  if (Mapper_getIsSorted(m) == 0) {
    Mapper_sort(m);
  }

  if (!IDHash_contains(hash, id)) {
    fprintf(stderr,"ERROR: Fastmap expects to be able to find an id. It couldnt for " IDFMTSTR "\n",id);
    exit(1);
  }

  pairs = IDHash_getValue(hash,id);

  MapperRangeSet *retSet = MapperRangeSet_new();

  for (i=0;i<MapperPairSet_getNumPair(pairs);i++) {
    MapperPair *pair = MapperPairSet_getPairAt(pairs,i);
    MapperUnit *selfCoord   = MapperPair_getUnit(pair, from);
    MapperUnit *targetCoord = MapperPair_getUnit(pair, to);

    // only super easy mapping is done
    if(start < selfCoord->start ||
       end   > selfCoord->end ) {
      continue;
    }

    if (pair->ori == 1) {
      MapperCoordinate *retRange = MapperCoordinate_new(targetCoord->id, 
                                                        targetCoord->start + start - selfCoord->start, 
                                                        targetCoord->start + end   - selfCoord->start,
                                                        strand,
                                                        cs, 
                                                        0); // Perl didn't set rank, so use 0 

/*
      retRange->id     = targetCoord->id;
      retRange->start  = targetCoord->start + start - selfCoord->start;
      retRange->end    = targetCoord->start + end   - selfCoord->start;
      retRange->strand = strand;
      retRange->coordSystem = cs;
*/

      MapperRangeSet_addRange(retSet, (MapperRange *)retRange);
      break;
    } else {
      MapperCoordinate *retRange = MapperCoordinate_new(targetCoord->id, 
                                                        targetCoord->end - (end - selfCoord->start),
                                                        targetCoord->end - (start - selfCoord->start),
                                                        -strand,
                                                        cs, 
                                                        0); // Perl didn't set rank, so use 0 

/*
      retRange->id     = targetCoord->id;
      retRange->start  = targetCoord->end - (end - selfCoord->start);
      retRange->end    = targetCoord->end - (start - selfCoord->start);
      retRange->strand = -strand;
      retRange->coordSystem = cs;
*/

      MapperRangeSet_addRange(retSet, (MapperRange *)retRange);
      break;
    }
  }

  // NIY: Here we return empty set, in mapInsert it returns NULL for empty fastmap - need to work out which is right
  return retSet;
}


/*
=head2 add_map_coordinates

    Arg  1      int $id
                id of 'source' sequence
    Arg  2      int $start
                start coordinate of 'source' sequence
    Arg  3      int $end
                end coordinate of 'source' sequence
    Arg  4      int $strand
                relative orientation of source and target (+/- 1)
    Arg  5      int $id
                id of 'target' sequence
    Arg  6      int $start
                start coordinate of 'target' sequence
    Arg  7      int $end
                end coordinate of 'target' sequence
    Function    Stores details of mapping between
                'source' and 'target' regions.
    Returntype  none
    Exceptions  none
    Caller      Bio::EnsEMBL::Mapper

=cut
*/
void Mapper_addMapCoordinates(Mapper *m, IDType contigId, int contigStart, int contigEnd,
                              int contigOri, IDType chrId, int chrStart, int chrEnd) {
  MapperPair *pair;
  MapperUnit *from;
  MapperUnit *to;
  IDHash     *fromHash;
  IDHash     *toHash;
  MapperPairSet *mps;

  if ((contigEnd - contigStart) != (chrEnd - chrStart)) {
    fprintf(stderr,"ERROR: Cannot deal with mis-lengthed mappings so far\n");
    fprintf(stderr,"Contig %d to %d and chromosome %d to %d\n",contigStart,contigEnd,
            chrStart,chrEnd);
    exit(1);
  }

  pair = MapperPair_new();

  from = MapperUnit_new();

  from->start = contigStart;
  from->end = contigEnd;
  from->id = contigId;

  to = MapperUnit_new();

  to->start = chrStart;
  to->end = chrEnd;
  to->id = chrId;

  MapperPair_setUnit(pair,MAPPER_TO_IND,to);
  MapperPair_setUnit(pair,MAPPER_FROM_IND,from);

  pair->ori = contigOri;

  // place into hash on both ids
  fromHash = Mapper_getPairHash(m, MAPPER_FROM_IND);
  toHash   = Mapper_getPairHash(m, MAPPER_TO_IND);

  if (!IDHash_contains(toHash,chrId)) {
    IDHash_add(toHash,chrId,MapperPairSet_new());
  }
  mps = (MapperPairSet *)IDHash_getValue(toHash,chrId);
  MapperPairSet_addPair(mps,pair);

  if (!IDHash_contains(fromHash,contigId)) {
    IDHash_add(fromHash,contigId,MapperPairSet_new());
  }
  mps = (MapperPairSet *)IDHash_getValue(fromHash,contigId);
  MapperPairSet_addPair(mps,pair);

  Mapper_incPairCount(m);
  Mapper_setIsSorted(m,0);
}


/*
=head2 add_indel_coordinates

    Arg  1      int $id
                id of 'source' sequence
    Arg  2      int $start
                start coordinate of 'source' sequence
    Arg  3      int $end
                end coordinate of 'source' sequence
    Arg  4      int $strand
                relative orientation of source and target (+/- 1)
    Arg  5      int $id
                id of 'targe' sequence
    Arg  6      int $start
                start coordinate of 'targe' sequence
    Arg  7      int $end
                end coordinate of 'targe' sequence
    Function    stores details of mapping between two regions:
                'source' and 'target'. Returns 1 if the pair was added, 0 if it
                was already in. Used when adding an indel
    Returntype  int 0,1
    Exceptions  none
    Caller      Bio::EnsEMBL::Mapper

=cut
*/

// This is almost identical to Mapper_addCoordinates (just the isIndel line as far as I can see!) - I should refactor this
int Mapper_addIndelCoordinates( Mapper *m, IDType contigId, long contigStart, long contigEnd, int contigOri, IDType chrId, long chrStart, long chrEnd) {
  MapperPair *pair;
  MapperUnit *from;
  MapperUnit *to;
  IDHash     *fromHash;
  IDHash     *toHash;
  MapperPairSet *mps;


  //we need to create the IndelPair object to add to both lists, to and from
  pair = MapperPair_new();

  from = MapperUnit_new();

  from->start = contigStart;
  from->end = contigEnd;
  from->id = contigId;

  to = MapperUnit_new();

  to->start = chrStart;
  to->end = chrEnd;
  to->id = chrId;

  MapperPair_setUnit(pair,MAPPER_TO_IND,to);
  MapperPair_setUnit(pair,MAPPER_FROM_IND,from);


  pair->ori = contigOri;
  pair->isIndel = 1;

  // place into hash on both ids
  fromHash = Mapper_getPairHash(m, MAPPER_FROM_IND);
  toHash   = Mapper_getPairHash(m, MAPPER_TO_IND);

  if (!IDHash_contains(toHash,chrId)) {
    IDHash_add(toHash,chrId,MapperPairSet_new());
  }
  mps = (MapperPairSet *)IDHash_getValue(toHash,chrId);
  MapperPairSet_addPair(mps,pair);

  if (!IDHash_contains(fromHash,contigId)) {
    IDHash_add(fromHash,contigId,MapperPairSet_new());
  }
  mps = (MapperPairSet *)IDHash_getValue(fromHash,contigId);
  MapperPairSet_addPair(mps,pair);

  Mapper_incPairCount(m);
  Mapper_setIsSorted(m,0);

  return 1;
}


/*
=head2 map_indel

  Arg [1]    : string $id
  Arg [2]    : int $start - start coord. Since this is an indel should always
               be one greater than end.
  Arg [3]    : int $end - end coord. Since this is an indel should always
               be one less than start.
  Arg [4]    : int $strand (0, 1, -1)
  Arg [5]    : string $type - the coordinate system name the coords are from.
  Example    : @coords = $mapper->map_indel();
  Description: This is in internal function which handles the special mapping
               case for indels (start = end +1). It will be used to map from
               a coordinate system with a gap to another that contains an
               insertion. It will be mainly used by the Variation API.
  Returntype : Bio::EnsEMBL::Mapper::Unit objects
  Exceptions : none
  Caller     : general

=cut
*/

/*  Only used by alignstrainslice as far as I can see, so don't bother implementing for now.
sub Mapper_mapIndel {
  my ( $self, $id, $start, $end, $strand, $type ) = @_;

  // swap start/end and map the resultant 2bp coordinate
  ( $start, $end ) = ( $end, $start );

  if ( !$self->{'_is_sorted'} ) { $self->_sort() }

  my $hash = $self->{"_pair_$type"};

  my ( $from, $to, $cs );

  if ( $type eq $self->{'to'} ) {
    $from = 'to';
    $to   = 'from';
    $cs   = $self->{'from_cs'};
  } else {
    $from = 'from';
    $to   = 'to';
    $cs   = $self->{'to_cs'};
  }

  unless ( defined $hash ) {
    throw("Type $type is neither to or from coordinate systems");
  }
  my @indel_coordinates;

  my ( $start_idx, $end_idx, $mid_idx, $pair, $self_coord );
  my $lr = $hash->{ uc($id) };

  $start_idx = 0;
  $end_idx   = $#{$lr};

  // binary search the relevant pairs
  // helps if the list is big
  while ( ( $end_idx - $start_idx ) > 1 ) {
    $mid_idx    = ( $start_idx + $end_idx ) >> 1;
    $pair       = $lr->[$mid_idx];
    $self_coord = $pair->{$from};
    if ( $self_coord->{'end'} <= $start ) {
      $start_idx = $mid_idx;
    } else {
      $end_idx = $mid_idx;
    }
  }

  for ( my $i = $start_idx; $i <= $#{$lr}; $i++ ) {
    $pair = $lr->[$i];
    my $self_coord   = $pair->{$from};
    my $target_coord = $pair->{$to};

    if ( exists $pair->{'indel'} ) {
      #need to return unit coordinate
      my $to =
        Bio::EnsEMBL::Mapper::Unit->new( $target_coord->{'id'},
                                         $target_coord->{'start'},
                                         $target_coord->{'end'}, );
      push @indel_coordinates, $to;
      last;
    }
  }

  return @indel_coordinates;
}
*/


/*
=head2 add_Mapper

    Arg  1      Bio::EnsEMBL::Mapper $mapper2
    Example     $mapper->add_Mapper($mapper2)
    Function    add all the map coordinates from $mapper to this mapper.
                This object will contain mapping pairs from both the old
                object and $mapper2.
    Returntype  int 0,1
    Exceptions  throw if 'to' and 'from' from both Bio::EnsEMBL::Mappers
                are incompatible
    Caller      $mapper->methodname()

=cut
*/

// Changed type to void because only ever returned 1 (or exited)
/* Only used by alignslice so don't bother implementing for now
void Mapper_addMapper(Mapper *m, Mapper *mToAdd) {

  my $mapper_to = $mapper->{'to'};
  my $mapper_from = $mapper->{'from'};

  if ($mapper_to ne $self->{'to'} or $mapper_from ne $self->{'from'}) {
    fprintf(stderr,"Trying to add an incompatible Mapper\n");
    exit(1);
  }

  int countA = 0;
  foreach my $seq_name (keys %{$mapper->{"_pair_$mapper_to"}}) {

    push(@{$self->{"_pair_$mapper_to"}->{$seq_name}}, @{$mapper->{"_pair_$mapper_to"}->{$seq_name}});

    $count_a += scalar(@{$mapper->{"_pair_$mapper_to"}->{$seq_name}});

  }

  int countB = 0;
  foreach my $seq_name (keys %{$mapper->{"_pair_$mapper_from"}}) {

    push(@{$self->{"_pair_$mapper_from"}->{$seq_name}}, @{$mapper->{"_pair_$mapper_from"}->{$seq_name}});

    $count_b += scalar(@{$mapper->{"_pair_$mapper_from"}->{$seq_name}});

  }

  if (countA == countB) {
    $self->{'pair_count'} += countA;
  } else {
    fprintf(stderr, "Trying to add a funny Mapper\n");
    exit(1);
  }

  Mapper_setIsSorted(m, 0);
  return 1;
}
*/


/*
=head2 list_pairs

    Arg  1      int $id
                id of 'source' sequence
    Arg  2      int $start
                start coordinate of 'source' sequence
    Arg  3      int $end
                end coordinate of 'source' sequence
    Arg  4      string $type
                nature of transform - gives the type of
                coordinates to be transformed *from*
    Function    list all pairs of mappings in a region
    Returntype  list of Bio::EnsEMBL::Mapper::Pair
    Exceptions  none
    Caller      Bio::EnsEMBL::Mapper

=cut
*/
MapperPairSet *Mapper_listPairs(Mapper *m, IDType id, long start, long end, char *type) {
  MapperPairSet *pairs;
  IDHash *hash;
  int from, to;
  MapperPairSet *retSet;
  int i;

  if (start > end) {
    fprintf(stderr,"ERROR: Start is greater than end for id " IDFMTSTR ", start %ld, end %ld\n",id,start,end);
  }

  if( Mapper_getIsSorted(m) == 0 ) {
    Mapper_sort(m);
  }

  if (!Mapper_compareType(type, Mapper_getTo(m))) {
    from = MAPPER_TO_IND;
    to   = MAPPER_FROM_IND;
  } else {
    from = MAPPER_FROM_IND;
    to   = MAPPER_TO_IND;
  }

  hash = Mapper_getPairHash(m, from);

  if (!hash) {
    fprintf(stderr,"ERROR: Type %s is neither to or from coordinate system\n",type);
    exit(1);
  }


  if (!IDHash_contains(hash, id)) {
    return NULL;
  }

  pairs = IDHash_getValue(hash,id);


  retSet = MapperPairSet_new();
  //Mapper_dump(m,NULL);
  //fprintf(stderr, "listPairs with %d %d %d\n",id,start,end);

  if (start == -1 && end == -1) {
    for (i=0;i<MapperPairSet_getNumPair(pairs);i++) {
      MapperPairSet_addPair(retSet,MapperPairSet_getPairAt(pairs,i));
    }
  } else {

    for (i=0;i<MapperPairSet_getNumPair(pairs);i++) {
      MapperPair *pair = MapperPairSet_getPairAt(pairs,i);
      MapperUnit *fromCoord   = MapperPair_getUnit(pair, from);
      //fprintf(stderr," unit %d %d\n",fromCoord->start,fromCoord->end);

      if( fromCoord->end < start ) {
        continue;
      }
      if( fromCoord->start > end ) {
        break;
      }
      MapperPairSet_addPair(retSet,pair);
    }
  }
  return retSet;
}

/*
#  _dump
#
#     Arg  1      *FileHandle $fh
#     Function    convenience dump function
#                 possibly useful for debugging
#     Returntype  none
#     Exceptions  none
#     Caller      internal
#
*/

void Mapper_dump(Mapper *m, FILE *fp) {
  IDHash *fromHash;
  IDType *keys;
  int   nKey;
  int   i;

  if (fp==NULL) {
    fp = stderr;
  }

  fromHash = Mapper_getPairHash(m, MAPPER_TO_IND);

  keys = IDHash_getKeys(fromHash);
  nKey = IDHash_getNumValues(fromHash);

  for (i=0;i<nKey;i++) {
    IDType id = keys[i];
    MapperPairSet *set = IDHash_getValue(fromHash,id);
    int j;

    fprintf(fp, "From Hash " IDFMTSTR " with %d pairs\n",id, MapperPairSet_getNumPair(set));

    for (j=0; j<MapperPairSet_getNumPair(set); j++) {
      MapperPair *pair = MapperPairSet_getPairAt(set,j);
      MapperUnit *fromCoord = MapperPair_getUnit(pair, MAPPER_FROM_IND);
      MapperUnit *toCoord   = MapperPair_getUnit(pair, MAPPER_TO_IND);

      fprintf(fp, "    %ld %ld:%ld %ld " IDFMTSTR "\n",fromCoord->start,fromCoord->end,
              toCoord->start,toCoord->end,toCoord->id);
    }
  }

  free(keys);
}


/*
# _sort
#
#    Function    sort function so that all
#                mappings are sorted by
#                chromosome start
#    Returntype  none
#    Exceptions  none
#    Caller      internal
#
*/

void Mapper_sort(Mapper *m) {

  IDHash *fromHash = Mapper_getPairHash(m, MAPPER_FROM_IND);
  IDHash *toHash   = Mapper_getPairHash(m, MAPPER_TO_IND);
  MapperPairSet **sets;
  int nSet;
  int i;

  sets = (MapperPairSet **)IDHash_getValues(fromHash);
  nSet = IDHash_getNumValues(fromHash);
  for (i=0;i<nSet;i++) {
    MapperPairSet_sort(sets[i],MAPPER_FROM_IND);
  }
  free(sets);

  sets = (MapperPairSet **)IDHash_getValues(toHash);
  nSet = IDHash_getNumValues(toHash);
  for (i=0;i<nSet;i++) {
    MapperPairSet_sort(sets[i],MAPPER_TO_IND);
  }
  free(sets);

// MergePairs yet to do
  Mapper_mergePairs(m);

  Mapper_setIsSorted(m, 1);

}

// this function merges pairs that are adjacent into one
// This function is a pain in the arse to implement in C
void Mapper_mergePairs(Mapper *m) {

  int to   = MAPPER_TO_IND;
  int from = MAPPER_FROM_IND;

  Mapper_setPairCount(m, 0);

  IDHash *toPairHash   = Mapper_getPairHash(m, MAPPER_TO_IND);
  IDHash *fromPairHash = Mapper_getPairHash(m, MAPPER_FROM_IND);
 
  MapperPairSet **toPairValues = (MapperPairSet **)IDHash_getValues(toPairHash);

  
  int pairInd;
  for (pairInd = 0; pairInd<IDHash_getNumValues(toPairHash); pairInd++) {
    MapperPairSet *pairs = toPairValues[pairInd];

    int i = 0;
    int next = 1;
    int length = MapperPairSet_getNumPair(pairs)-1; //$#{$lr};

    while (next <= length) {
      MapperPair *currentPair = MapperPairSet_getPairAt(pairs, i);
      MapperPair *nextPair    = MapperPairSet_getPairAt(pairs, next);
      MapperPair *delPair     = NULL;

      if (MapperPair_isIndel(currentPair) || MapperPair_isIndel(nextPair)) {
        //necessary to modify the merge function to not merge indels
        next++;
        i++;

      } else {
        // duplicate filter
        if ( MapperPair_getUnit(currentPair,to)->start == MapperPair_getUnit(nextPair,to)->start  &&
             MapperPair_getUnit(currentPair,from)->id == MapperPair_getUnit(nextPair,from)->id ) {

          delPair = nextPair;

        } else if (( MapperPair_getUnit(currentPair,from)->id == MapperPair_getUnit(nextPair,from)->id ) &&
                   ( nextPair->ori == currentPair->ori ) &&
                   ( MapperPair_getUnit(nextPair,to)->start-1 == MapperPair_getUnit(currentPair,to)->end )) {

          if ( currentPair->ori == 1 ) {

            // check forward strand merge
            if ( MapperPair_getUnit(nextPair,from)->start-1 == MapperPair_getUnit(currentPair,from)->end) {
              // normal merge with previous element
              MapperPair_getUnit(currentPair,to)->end = MapperPair_getUnit(nextPair,to)->end;
              MapperPair_getUnit(currentPair,from)->end = MapperPair_getUnit(nextPair,from)->end;
              delPair = nextPair;
            }
          } else {

            // check backward strand merge
            if ( MapperPair_getUnit(nextPair,from)->end+1  == MapperPair_getUnit(currentPair,from)->start ) {

              // yes its a merge
              MapperPair_getUnit(currentPair,to)->end = MapperPair_getUnit(nextPair,to)->end;
              MapperPair_getUnit(currentPair,from)->start = MapperPair_getUnit(nextPair,from)->start;
              delPair = nextPair;
            }
          }
        }

        if (delPair != NULL) { // Have a pair to delete
          // Remove from the to pair set
          MapperPairSet_removePairAt(pairs, next); //splice( @$lr, $next, 1 );

          MapperPairSet *fromPairs = IDHash_getValue(fromPairHash, MapperPair_getUnit(delPair, from)->id); //$self->{"_pair_$map_from"}->{uc($del_pair->{'from'}->{'id'})}; 

          int j;
          for (j=0; j < MapperPairSet_getNumPair(fromPairs); j++) {
            MapperPair *fromPair = MapperPairSet_getPairAt(fromPairs, j);
            if ( fromPair == delPair) { // Is this really going to be an equality ??? //$lr_from->[$j] == $del_pair )
              MapperPairSet_removePairAt(fromPairs, j); //splice( @$lr_from, $j, 1 );
              break;
            }
          }

          // NIY: Do we need to free delPair???

          length--;
          if ( length < next ) break;
        } else {
          next++;
          i++;
        }
      }
    }

    Mapper_addToPairCount(m, MapperPairSet_getNumPair(pairs)); //    $self->{'pair_count'} += scalar( @$lr );
  }
}

void Mapper_free(Mapper *mapper) {
  Mapper_flush(mapper);

  // NIY: Probably more cleanup needed

  free(mapper);
}
