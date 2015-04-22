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

#include "CachingSequenceAdaptor.h"

#include "DBAdaptor.h"
#include "BaseAdaptor.h"
#include "SeqUtil.h"
#include "SliceAdaptor.h"
#include "Slice.h"
#include "StrUtil.h"
#include "LRUCache.h"


/*
=head1 DESCRIPTION

An adaptor for the retrieval of DNA sequence from a string 
*/

static long const CACHING_SEQ_CACHE_SZ  = 20; // Number of seq region sequences to cache

/*
=head2 new

  Arg [1]    : none
  Example    : my $sa = $db_adaptor->get_CachingSequenceAdaptor();
  Description: Constructor.  Calls superclass constructor and initialises
               internal cache structure.
  Returntype : Bio::EnsEMBL::DBSQL::CachingSequenceAdaptor
  Exceptions : none
  Caller     : DBAdaptor::get_CachingSequenceAdaptor
  Status     : Stable

=cut
*/
CachingSequenceAdaptor *CachingSequenceAdaptor_new(DBAdaptor *dba) {
  CachingSequenceAdaptor *csa;

  if ((csa = (CachingSequenceAdaptor *)calloc(1,sizeof(CachingSequenceAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for CachingSequenceAdaptor\n");
    return NULL;
  }
  BaseAdaptor_init((BaseAdaptor *)csa, dba, CACHINGSEQUENCE_ADAPTOR);

  // use an LRU cache to limit the size
  csa->seqCache = LRUCache_new(CACHING_SEQ_CACHE_SZ);

  return csa;
}


/*
=head2 clear_cache

  Example			: $sa->clear_cache();
  Description	: Removes all entries from the associcated sequence cache
  Returntype 	: None
  Exceptions 	: None

=cut
*/

void CachingSequenceAdaptor_clearCache(CachingSequenceAdaptor *csa) {
  //fprintf(stderr,"clearCache called\n");
  LRUCache_empty(csa->seqCache);

  return;
}


/*
=head2 fetch_by_Slice_start_end_strand

  Arg  [1]   : Bio::EnsEMBL::Slice slice
               The slice from which you want the sequence
  Arg  [2]   : (optional) int startBasePair 
               The start base pair relative to the start of the slice. Negative
               values or values greater than the length of the slice are fine.
               default = 1
  Arg  [3]   : (optional) int endBasePair
               The end base pair relative to the start of the slice. Negative
               values or values greater than the length of the slice are fine,
               but the end must be greater than or equal to the start
               count from 1
               default = the length of the slice
  Arg  [4]   : (optional) int strand 
               1, -1
               default = 1
  Example    : $dna = $seq_adptr->fetch_by_Slice_start_end_strand($slice, 1, 
                                                                  1000, -1);
  Description: retrieves from db the sequence for this slice
               uses AssemblyMapper to find the assembly
  Returntype : string 
  Exceptions : endBasePair should be less or equal to length of slice 
  Caller     : Bio::EnsEMBL::Slice::seq(), Slice::subseq() 
  Status     : Stable

=cut
*/
char *CachingSequenceAdaptor_fetchBySliceStartEndStrand(CachingSequenceAdaptor *csa,
                                                        Slice *slice, long start, long end,
                                                        int strand) {
  if (!slice) {
    fprintf(stderr,"ERROR: need a slice to work\n");
    exit(1);
  }

  if (start == POS_UNDEF) {
    start = 1;
  }

  if (BaseContig_getObjectType(slice) != CLASS_SLICE) {
    fprintf(stderr,"ERROR: slice fetch called with something that isn't a slice\n");
    exit(1);
  }

  if (end == POS_UNDEF) { // && (not $slice->is_circular) ) {
    end = Slice_getEnd(slice) - Slice_getStart(slice) + 1;
  }

  if ( start > end ) {
    fprintf(stderr,"Start must be less than or equal to end.\n");
    exit(1);
  }

  if (strand == STRAND_UNDEF) {
    strand = 1;
  }

  // get a new slice that spans the exact region to retrieve dna from
  long rightExpand  = end - Slice_getLength(slice); // negative is fine
  long leftExpand   = 1 - start; // negative is fine

  if ( Slice_getStrand(slice) != 1 ) {
    long tmp = rightExpand;
    rightExpand = leftExpand;
    leftExpand = tmp;
  }

  long seqRegStart = Slice_getStart(slice) - leftExpand;
  long seqRegEnd   = Slice_getEnd(slice) + rightExpand;

  char *seq;
  if ((seq = calloc((seqRegEnd-seqRegStart)+2, sizeof(char))) == NULL) {
    fprintf(stderr,"Failed allocating seq\n");
    exit(1);
  }
  memset(seq, 'N', seqRegEnd-seqRegStart+1);

  char idStr[1024];
  sprintf(idStr, IDFMTSTR, Slice_getSeqRegionId(slice));

  char *seqRegSeq;
  long  seqRegLen;

  //fprintf(stderr,"idStr = %s\n", idStr);
// Is this region in the cache
  if (!LRUCache_contains(csa->seqCache, idStr)) {
//   Get from underlying sequence adaptor
//   Store in cache
    Slice *seqRegSlice = Slice_getSeqRegionSlice(slice);


    seqRegSeq = Slice_getSeq(seqRegSlice);
    seqRegLen  = Slice_getLength(seqRegSlice); 

    // fprintf(stderr,"seqRegSlice name %s length %d\n", Slice_getName(seqRegSlice), seqRegLen);
    LRUCache_put(csa->seqCache, idStr, seqRegSeq, free, seqRegLen);
  //fprintf(stderr,"put with idStr = %s\n", idStr);

    Slice_free(seqRegSlice);
// else
  } else {
//   fetch sequence string (with its length) from cache
    seqRegSeq = LRUCache_get(csa->seqCache, idStr);
    seqRegLen = LRUCache_getSize(csa->seqCache, idStr);
    // fprintf(stderr,"seqRegSeq length %d\n", seqRegLen);
  //fprintf(stderr,"get with idStr = %s\n", idStr);
  }

// fetch the required region from the sequence string and put in place in output string in correct location
//
  long startOffset = 0;
  long startPos = seqRegStart - 1;
  long endPos   = seqRegEnd - 1;

  if (startPos < 0) {
    startOffset = -startPos;
  }

  if (endPos >= seqRegLen) {
    endPos = seqRegLen;
  }
  // fprintf(stderr,"startOffet = %ld startPos = %ld endPos = %ld\n", startOffset, startPos, endPos);
  memcpy(&seq[startOffset], &seqRegSeq[startPos], endPos-startPos+1);

  // if they asked for the negative slice strand revcomp the whole thing
  if (strand == -1) {
    SeqUtil_reverseComplement(seq, endPos-startPos+1);
  }

  return seq;
}

