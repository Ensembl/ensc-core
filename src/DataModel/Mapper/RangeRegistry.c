/*
 * Copyright [1999-2017] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

#include "RangeRegistry.h"
#include "CoordPair.h"
/*
=head2 new

  Arg [1]    : none
  Example    : my $rr = Bio::EnsEMBL::Mapper::RangeRegistry->new();
  Description: Creates a new RangeRegistry object
  Returntype : Bio::EnsEMBL::Mapper::RangeRegistry
  Exceptions : none
  Caller     : AssemblyMapperAdaptor
  Status     : Stable

=cut
*/

RangeRegistry *RangeRegistry_new() {
  RangeRegistry *registry;

  if ((registry = (RangeRegistry *)calloc(1,sizeof(RangeRegistry))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for registry\n");
    exit(1);
  }

  registry->registry = IDHash_new(IDHASH_MEDIUM);

  return registry;
}

void RangeRegistry_flush(RangeRegistry *registry) {
  IDHash_free(registry->registry, Vector_free);

  registry->registry = IDHash_new(IDHASH_MEDIUM);
}

/*
=head2 check_and_register

  Arg [1]    : string $id
               The id of the range to be checked/registered (e.g. a sequenceid)
  Arg [2]    : int $start
               The start of the range to be checked
  Arg [3]    : int $end
               The end of the range to be checked
  Arg [4]    : (optional) int $rstart
               The start of the range to be registered
               if the checked range was not fully registered
  Arg [5]    : (optional) int $rend
               The end of the range to be registerd
               if the checked range was not fully registered
  Example    : $ranges=$rr->check_and_register('X',500,600, 1,1_000_000));
  Description: Checks the range registry to see if the entire range
               denoted by ($id : $start-$end) is already registered.  If
               it already is, then undef is returned.  If it is not then
               the range specified by $rstart and $rend is registered, and
               a list of regions that were required to completely register
               $rstart-$rend is returned.  If $rstart and $rend are not
               defined they default to $start and $end respectively.

               The reason there is a single call to do both the checking and
               registering is to reduce the overhead. Much of the work to
               check if a range is registered is the same as registering a
               region around that range.
  Returntype : undef or listref of [start,end] range pairs
  Exceptions : throw if rstart is greater than start
               throw if rend is less than end
               throw if end is less than start
               throw if id, start, or end are not defined
  Caller     : AssemblyMapperAdaptor
  Status     : Stable

=cut
*/

// Perl version allows rStart and rEnd to be optional. I require them to be set - if you don't need them to be different to start and end,
// just use the same values as start and end for rStart and rEnd

// Also added a flag to indicate we actually want the gaps vector returned - quite often its not used in the caller and so would leak
// memory
Vector *RangeRegistry_checkAndRegister(RangeRegistry *registry, IDType id, long start, long end, 
                                       long rStart, long rEnd, int wantGaps) {

  // The following was commented out due to Ensembl Genomes requirements
  // for bacterial genomes.
  // The following was uncommented because I'm not caring about those requirements
  if ( start > end ) {
    fprintf(stderr, "start argument [%ld] must be less than (or equal to) end argument [%ld]\n", start, end);
    exit(1);
  }
  
  if ( rStart > rEnd ) {
    fprintf(stderr, "rStart argument [%ld] must be less than (or equal to) rEnd argument [%ld]\n", rStart, rEnd);
    exit(1);
  }

  if ( rStart > start ) {
    fprintf(stderr, "rStart argument [%ld] must be less than (or equal to) start [%ld]\n", rStart, start);
    exit(1);
  }

  if ( rEnd < end ) {
    fprintf(stderr, "rEnd argument [%ld] must be greater than (or equal to) end [%ld]\n", rEnd, end);
    exit(1);
  }

  IDHash *regReg = RangeRegistry_getRegistry(registry);
  Vector *list;
  if (IDHash_contains(regReg, id)) {
    list = IDHash_getValue(regReg, id);
  } else {
    list = Vector_new();
    IDHash_add(regReg, id, list);
  }

  
  Vector *gapPairs = NULL;
  if (wantGaps) {
    gapPairs = Vector_new();
  }

  int len = Vector_getNumElement(list);

  if (len == 0) {
    //this is the first request for this id, return a gap pair for the
    // entire range and register it as seen
    CoordPair *cp = CoordPair_new(rStart, rEnd);
    Vector_addElement(list, cp);

    return Vector_copy(list);
  }

  //####
  // loop through the list of existing ranges recording any "gaps" where
  // the existing range does not cover part of the requested range
  // 

  int startIdx = 0;
  int endIdx   = Vector_getNumElement(list)-1;
  int midIdx;
  CoordPair *range;

  // binary search the relevant pairs
  // helps if the list is big
  while ( ( endIdx - startIdx ) > 1 ) {
    midIdx = ( startIdx + endIdx ) >> 1;
    range  = Vector_getElementAt(list, midIdx);

    if ( CoordPair_getEnd(range) < rStart ) {
      startIdx = midIdx;
    } else {
      endIdx = midIdx;
    }
  }

  long gapStart;
  long gapEnd;
  int rIdx = -1;
  int rStartIdx = -1;
  int rEndIdx;

  gapStart = rStart;

  int i;
  for (i=startIdx; i < len ; i++ ) {
    CoordPair *pRange = Vector_getElementAt(list, i);
    long pStart = CoordPair_getStart(pRange);
    long pEnd   = CoordPair_getEnd(pRange);
    
    // no work needs to be done at all if we find a range pair that
    // entirely overlaps the requested region
    if ( pStart <= start && pEnd >= end ) {
      return Vector_new(); // perl returns undef, but that causes me problems
    }

    // find adjacent or overlapping regions already registered
    if ( pEnd >= ( rStart - 1 ) && pStart <= ( rEnd + 1 ) ) {
      if ( rStartIdx < 0 ) { // Not yet been set
        rStartIdx = i;
      }
      rEndIdx = i;
    }

    if ( pStart > rStart ) {
      gapEnd = ( rEnd < pStart ) ? rEnd : pStart - 1;
      if (wantGaps) {
        CoordPair *cp = CoordPair_new(gapStart, gapEnd);
        Vector_addElement(gapPairs, cp);
      }
    }

    gapStart = ( rStart > pEnd ) ? rStart : pEnd + 1;

    if ( pEnd >= rEnd && rIdx < 0 ) {
      rIdx = i;
      break;
    }
  }

  // do we have to make another gap?
  if ( gapStart <= rEnd ) {
    if (wantGaps) {
      CoordPair *cp = CoordPair_new(gapStart, rEnd);
      Vector_addElement(gapPairs, cp);
    }
  }

  // 
  // Merge the new range into the registered list
  // 
  if (rStartIdx >= 0 ) { // rStartIdx has been set to something 
    long newStart;
    long newEnd;
    CoordPair *rStartIdxRange = Vector_getElementAt(list, rStartIdx); 
    CoordPair *rEndIdxRange   = Vector_getElementAt(list, rEndIdx); 

    if ( rStart < CoordPair_getStart(rStartIdxRange)) {
      newStart = rStart;
    } else {
      newStart = CoordPair_getStart(rStartIdxRange);
    }

    if ( rEnd > CoordPair_getEnd(rEndIdxRange)) {
      newEnd = rEnd;
    } else {
      newEnd = CoordPair_getEnd(rEndIdxRange);
    }

    CoordPair *cp = CoordPair_new(newStart, newEnd);

    // Think its <=
    for (i=rStartIdx; i<=rEndIdx; i++) {
      Vector_removeElementAt(list, rStartIdx); // Always remove from rStartIdx as array is shrinking by one each time called
    }
    Vector_insertElementAt(list, rStartIdx, cp);
    //splice( @$list, $rstart_idx,
    //        $rend_idx - $rstart_idx + 1,
    //        [ $new_start, $new_end ] );

  } else if (rIdx >= 0) {
    CoordPair *cp = CoordPair_new(rStart, rEnd);
    Vector_insertElementAt(list, rIdx, cp);
    //splice( @$list, $r_idx, 0, [ $rstart, $rend ] );
  } else {
    CoordPair *cp = CoordPair_new(rStart, rEnd);
    Vector_addElement(list, cp);
  }

  // Note if wantGaps is not set then gapPairs will be NULL - but you said you didn't want it so that should be OK
  return gapPairs;
}

// overlap size is just added to make RangeRegistry generally more useful

/*
=head2 overlap_size

  Arg [1]    : string $id
  Arg [2]    : int $start
  Arg [3]    : int $end
  Example    : my $overlap_size = $rr->( "chr1", 1, 100_000_000 )
  Description: finds out how many bases of the given range are registered
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/

long RangeRegistry_overlapSize(RangeRegistry *registry, IDType id, long start, long end)  {
  long overlap = 0;

  if ( start > end ) return 0;

  IDHash *regReg = RangeRegistry_getRegistry(registry);
  Vector *list;
  if (IDHash_contains(regReg, id)) {
    list = IDHash_getValue(regReg, id);
  } else {
    return 0; // No list for this id, so can't be any overlap 
  }

  int len = Vector_getNumElement(list);

  if ( len == 0 ) {
    fprintf(stderr, "Odd have zero length list in RangeRegistry_overlapSize\n");
    return 0;
  }

  int startIdx = 0;
  int endIdx   = Vector_getNumElement(list)-1;
  int midIdx;
  CoordPair *range;

  // binary search the relevant pairs
  // helps if the list is big
  while ( ( endIdx - startIdx ) > 1 ) {
    midIdx = ( startIdx + endIdx ) >> 1;
    range   = Vector_getElementAt(list, midIdx);
    if ( CoordPair_getEnd(range) < start ) {
      startIdx = midIdx;
    } else {
      endIdx = midIdx;
    }
  }

  int i;
  for (i=startIdx; i < len ; i++ ) {
    CoordPair *pRange = Vector_getElementAt(list, i);
    long pStart = CoordPair_getStart(pRange);
    long pEnd   = CoordPair_getEnd(pRange);

    if ( pStart > end ) {
      break;
    }

    if ( pStart <= start && pEnd >= end ) {
      overlap = end - start + 1;
      break;
    }

    long mStart = ( start < pStart ? pStart : start );
    long mEnd   = ( end   < pEnd   ? end    : pEnd );

    if (mEnd - mStart >= 0) {
      overlap += ( mEnd - mStart + 1 );
    }
  }

  return overlap;
}


// low level function to access the ranges
// only use for read access
Vector *RangeRegistry_getRanges(RangeRegistry *registry, IDType id) {
  IDHash *regReg = RangeRegistry_getRegistry(registry);
  Vector *list = NULL;

  if (IDHash_contains(regReg, id)) {
    list = IDHash_getValue(regReg, id); 
  }

  return list;
}
