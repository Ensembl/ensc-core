#include "Mapper.h"
#include "MapperCoordinate.h"
#include <stdio.h>

Mapper *Mapper_new(CoordSystem from, CoordSystem to) {
  Mapper *m;

  if ((m = (Mapper *)calloc(1,sizeof(Mapper))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for mapper\n");
    exit(1);
  }

  Mapper_setFrom(m, from);
  Mapper_setTo(m, to);
  
  Mapper_setPairHash(m, MAPPER_FROM_IND, IDHash_new(IDHASH_MEDIUM));
  Mapper_setPairHash(m, MAPPER_TO_IND,   IDHash_new(IDHASH_MEDIUM));
  
  return m;
}

/*
=head2 map_coordinates

    Arg  1      int $id
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
    Function    generic map method
    Returntype  array of Bio::EnsEMBL::Mapper::Coordinate
                and/or   Bio::EnsEMBL::Mapper::Gap
    Exceptions  none
    Caller      Bio::EnsEMBL::Mapper

=cut
*/

MapperRangeSet *Mapper_mapCoordinates(Mapper *m, long id, int start, int end, int strand, CoordSystem type) {
  IDHash *hash;
  MapperPair *lastUsedPair;
  MapperRangeSet *results;
  MapperPairSet *pairs;
  int i;
  int from, to;

  if(type == Mapper_getTo(m)) {
    from = MAPPER_TO_IND;
    to   = MAPPER_FROM_IND;
  } else {
    from = MAPPER_FROM_IND;
    to   = MAPPER_TO_IND;
  }

  hash = Mapper_getPairHash(m, from);

  if (!hash) {
    fprintf(stderr,"ERROR: Type %d is neither to or from coordinate system\n",type);
    exit(1);
  }

  if( Mapper_getIsSorted(m) == 0 ) {
    Mapper_sort(m);
  }

  results = MapperRangeSet_new();

  if (!IDHash_contains(hash, id)) {
    // one big gap!
    MapperRange *gap = (MapperRange *)MapperGap_new(start,end);
    MapperRangeSet_addRange(results,gap);
    return results;
  }

  pairs = IDHash_getValue(hash,id);
  
  for (i=0;i<MapperPairSet_getNumPair(pairs);i++) {
    MapperPair *pair = MapperPairSet_getPairAt(pairs,i);
    MapperUnit *selfCoord   = MapperPair_getUnit(pair, from);
    MapperUnit *targetCoord = MapperPair_getUnit(pair, to);
    MapperRange *range;
    int targetStart,targetEnd;

    // if we haven't even reached the start, move on
    if (selfCoord->end < start) {
      continue;
    }

    // if we have over run, break
    if (selfCoord->start > end) {
      break;
    }

    if (start < selfCoord->start) {
      // gap detected
      MapperRange *gap = (MapperRange *)MapperGap_new(start, selfCoord->start-1);
      MapperRangeSet_addRange(results,gap);	
      start = gap->end+1;
    }


    // start is somewhere inside the region
    if (pair->ori == 1) {
      targetStart = 
	   targetCoord->start + (start - selfCoord->start);
    } else {
      targetEnd = 
	   targetCoord->end - (start - selfCoord->start);
    }

    // either we are enveloping this map or not. If yes, then end
    // point (self perspective) is determined solely by target. If not
    // we need to adjust

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
        targetEnd = 
               targetCoord->start + (end - selfCoord->start);
      } else {
        targetStart = 
               targetCoord->end - (end - selfCoord->start);
      }
    }

    range = (MapperRange *)MapperCoordinate_new(targetCoord->id,
                                                targetStart,
                                                targetEnd,
                                                pair->ori * strand);
    MapperRangeSet_addRange(results, range);

    lastUsedPair = pair;
    start = selfCoord->end+1;
  }

  if (lastUsedPair == NULL) {

    MapperRange *gap = (MapperRange *)MapperGap_new(start,end);
    MapperRangeSet_addRange(results,gap);

  } else if (MapperPair_getUnit(lastUsedPair, from)->end < end ) {
    // gap at the end
    MapperRange *gap = (MapperRange *)MapperGap_new(
			   MapperPair_getUnit(lastUsedPair,from)->end + 1,
			   end);
    MapperRangeSet_addRange(results,gap);
  }

  if (strand == -1) {
    MapperRangeSet_reverse(results);
  }

  return results;
}

/*
=head2 fastmap

    Arg  1      int $id
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

int Mapper_fastMap(Mapper *m, long id, int start, int end, int strand, CoordSystem type, MapperCoordinate *retRange) {
  MapperRangeSet *results;
  MapperPairSet *pairs;
  int i;
  IDHash *hash;
  int from, to;

  if(type == Mapper_getTo(m)) {
    from = MAPPER_TO_IND;
    to   = MAPPER_FROM_IND;
  } else {
    from = MAPPER_FROM_IND;
    to   = MAPPER_TO_IND;
  }

  hash = Mapper_getPairHash(m, from);

  if (!hash) {
    fprintf(stderr,"ERROR: Type %d is neither to or from coordinate system\n",type);
    exit(1);
  }

  if (Mapper_getIsSorted(m) == 0) {
    Mapper_sort(m);
  }

  if (!IDHash_contains(hash, id)) {
    fprintf(stderr,"ERROR: Fastmap expects to be able to find an id. It couldnt for %d\n",id);
    exit(1);
  }

  pairs = IDHash_getValue(hash,id);
  
  for (i=0;i<MapperPairSet_getNumPair(pairs);i++) {
    MapperPair *pair = MapperPairSet_getPairAt(pairs,i);
    MapperUnit *selfCoord   = MapperPair_getUnit(pair, from);
    MapperUnit *targetCoord = MapperPair_getUnit(pair, to);
    int targetStart,targetEnd;

    // only super easy mapping is done 
    if(start < selfCoord->start ||
       end   > selfCoord->end ) {
      continue;
    }
     
    if (pair->ori == 1) {
      retRange->id     = targetCoord->id;
      retRange->start  = targetCoord->start + start - selfCoord->start;
      retRange->end    = targetCoord->start + end   - selfCoord->start;
      retRange->strand = strand;
      return 1;
    } else {
      retRange->id     = targetCoord->id;
      retRange->start  = targetCoord->end - (end - selfCoord->start);
      retRange->end    = targetCoord->end - (start - selfCoord->start);
      retRange->strand = -strand;
      return 1;
    }
  }

  return 0;
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
                id of 'targe' sequence
    Arg  6      int $start
                start coordinate of 'targe' sequence
    Arg  7      int $end
                end coordinate of 'targe' sequence
    Function    stores details of mapping between two regions:
                'source' and 'target'
    Returntype  none
    Exceptions  none
    Caller      Bio::EnsEMBL::Mapper

=cut
*/

void Mapper_addMapCoordinates(Mapper *m, long contigId, int contigStart, int contigEnd,
                              int contigOri, long chrId, int chrStart, int chrEnd) {
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

  Mapper_setIsSorted(m,0);
}

/*
=head2 list_pairs

    Arg  1      int $id
                id of 'source' sequence
    Arg  2      int $start
                start coordinate of 'source' sequence
    Arg  3      int $end
                end coordinate of 'source' sequence
    Arg  4      int $type
                nature of transform - gives the type of
                coordinates to be transformed *from*
    Function    list all pairs of mappings in a region
    Returntype  list of Bio::EnsEMBL::Mapper::Pair
    Exceptions  none
    Caller      Bio::EnsEMBL::Mapper

=cut
*/

MapperPairSet *Mapper_listPairs(Mapper *m, long id, int start, int end, CoordSystem type) {
  MapperPairSet *pairs;
  IDHash *hash;
  int from, to;
  MapperPairSet *retSet;
  int i;

  if (start > end) {
    fprintf(stderr,"ERROR: Start is greater than end for id %d, start %d, end %d\n",id,start,end);
  }

  if( Mapper_getIsSorted(m) == 0 ) {
    Mapper_sort(m);
  }

  if(type == Mapper_getTo(m)) {
    from = MAPPER_TO_IND;
    to   = MAPPER_FROM_IND;
  } else {
    from = MAPPER_FROM_IND;
    to   = MAPPER_TO_IND;
  }

  hash = Mapper_getPairHash(m, from);

  if (!hash) {
    fprintf(stderr,"ERROR: Type %d is neither to or from coordinate system\n",type);
    exit(1);
  }


  if (!IDHash_contains(hash, id)) {
    return NULL;
  }

  pairs = IDHash_getValue(hash,id);

  retSet = MapperPairSet_new();
  Mapper_dump(m,NULL);
  fprintf(stderr, "listPairs with %d %d %d\n",id,start,end);
  
  if (start == -1 && end == -1) {
    for (i=0;i<MapperPairSet_getNumPair(pairs);i++) {
      MapperPairSet_addPair(retSet,MapperPairSet_getPairAt(pairs,i));
    }
  } else {
     
    for (i=0;i<MapperPairSet_getNumPair(pairs);i++) {
      MapperPair *pair = MapperPairSet_getPairAt(pairs,i);
      MapperUnit *fromCoord   = MapperPair_getUnit(pair, from);
      fprintf(stderr," unit %d %d\n",fromCoord->start,fromCoord->end);
       
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
=head2 _dump

    Arg  1      *FileHandle $fh
    Function    convenience dump function
                possibly useful for debugging
    Returntype  none
    Exceptions  none
    Caller      internal

=cut
*/

void Mapper_dump(Mapper *m, FILE *fp) {
  IDHash *fromHash;
  long *keys;
  int   nKey;
  int   i;
  MapperPairSet *pairs;

  if (fp==NULL) {
    fp = stderr;
  }

  fromHash = Mapper_getPairHash(m, MAPPER_TO_IND);

  keys = IDHash_getKeys(fromHash);
  nKey = IDHash_getNumValues(fromHash);

  for (i=0;i<nKey;i++) {
    long id = keys[i];
    MapperPairSet *set = IDHash_getValue(fromHash,id);
    int j;

    fprintf(fp, "From Hash %d with %d pairs\n",id, MapperPairSet_getNumPair(set));

    for (j=0; j<MapperPairSet_getNumPair(set); j++) {
      MapperPair *pair = MapperPairSet_getPairAt(set,j);
      MapperUnit *fromCoord = MapperPair_getUnit(pair, MAPPER_FROM_IND);
      MapperUnit *toCoord   = MapperPair_getUnit(pair, MAPPER_TO_IND);

      fprintf(fp, "    %d %d:%d %d %d\n",fromCoord->start,fromCoord->end,
              toCoord->start,toCoord->end,toCoord->id);
    }
  }

  free(keys);
}

/*
=head2 _sort

    Function    sort function so that all
                mappings are sorted by
                chromosome start
    Returntype  none
    Exceptions  none
    Caller      internal

=cut
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

  Mapper_setIsSorted(m, 1);
}
