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

#include "SequenceAdaptor.h"

#include "DBAdaptor.h"
#include "BaseAdaptor.h"
#include "MysqlUtil.h"
#include "SeqUtil.h"
#include "AssemblyMapperAdaptor.h"
#include "CoordSystemAdaptor.h"
#include "SliceAdaptor.h"
#include "BaseContig.h"
#include "Slice.h"
#include "StrUtil.h"
#include "LRUCache.h"
#include "StringHash.h"

#include "ProjectionSegment.h"
#include "StatementHandle.h"
#include "ResultRow.h"
#include <math.h>

/*
=head1 DESCRIPTION

An adaptor for the retrieval of DNA sequence from the EnsEMBL database
*/

static long const SEQ_CHUNK_PWR = 18; // 2^18 = approx. 250KB
//static long const SEQ_CHUNK_PWR = 1; // Basically means don't cache
static long const SEQ_CACHE_SZ  = 20;
static long SEQ_CACHE_MAX;

static int init = 0;

void SequenceAdaptor_initFunc() {
  SEQ_CACHE_MAX = (1 << SEQ_CHUNK_PWR) * SEQ_CACHE_SZ;
//  fprintf(stderr, "seq_cache_max = %ld\n", SEQ_CACHE_MAX);
  init  = 1;
}

char * SequenceAdaptor_fetchSeq(SequenceAdaptor *sa, IDType seqRegionId, long start, long length);
void SequenceAdaptor_rnaEdit(SequenceAdaptor *sa, Slice *slice, char **seqPP, int recLev);
void SequenceAdaptor_initialiseRnaEdit(SequenceAdaptor *sa);
/*
=head2 new

  Arg [1]    : none
  Example    : my $sa = $db_adaptor->get_SequenceAdaptor();
  Description: Constructor.  Calls superclass constructor and initialises
               internal cache structure.
  Returntype : Bio::EnsEMBL::DBSQL::SequenceAdaptor
  Exceptions : none
  Caller     : DBAdaptor::get_SequenceAdaptor
  Status     : Stable

=cut
*/
SequenceAdaptor *SequenceAdaptor_new(DBAdaptor *dba) {
  SequenceAdaptor *sa;

  if (!init) {
    SequenceAdaptor_initFunc();
  }
  if ((sa = (SequenceAdaptor *)calloc(1,sizeof(SequenceAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for SequenceAdaptor\n");
    return NULL;
  }
  BaseAdaptor_init((BaseAdaptor *)sa, dba, SEQUENCE_ADAPTOR);

  ((BaseAdaptor*)sa)->prepare = SequenceAdaptor_prepare;

  // use an LRU cache to limit the size
  sa->seqCache = LRUCache_new(SEQ_CACHE_MAX);
//  sa->seqCache = StringHash_new(STRINGHASH_MEDIUM);

//
// See if this has any seq_region_attrib of type "_rna_edit_cache" if so store these
// in a  hash.
//

  return sa;
}



StatementHandle *SequenceAdaptor_prepare(BaseAdaptor *ba, char *qStr, size_t len) {
  //printf("Query = %s len = %d\n",qStr,len);
  return DBAdaptor_prepare(ba->dba->dnadb,qStr,len);
}


/*
=head2 clear_cache

  Example			: $sa->clear_cache();
  Description	: Removes all entries from the associcated sequence cache
  Returntype 	: None
  Exceptions 	: None

=cut
*/

void SequenceAdaptor_clearCache(SequenceAdaptor *sa) {
  //fprintf(stderr,"clearCache called\n");
  LRUCache_empty(sa->seqCache);

  //StringHash_free(sa->seqCache, free);
  //sa->seqCache = StringHash_new(STRINGHASH_MEDIUM);
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

char *SequenceAdaptor_fetchBySliceStartEndStrand(SequenceAdaptor *sa,
                                                 Slice *slice, long start, long end,
                                                 int strand) {
 int recLev = 0;
 return SequenceAdaptor_fetchBySliceStartEndStrandRecursive(sa, slice, start, end, strand, &recLev);
}

char *SequenceAdaptor_fetchBySliceStartEndStrandRecursive(SequenceAdaptor *sa,
                                                          Slice *slice, long start, long end,
                                                          int strand, int *recLev) {
  (*recLev)++;

  if (!slice ) {
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
 
  /* No circular stuff
  if ( ( !defined($end) || $start > $end || $start < 0 || $end < 0 || $slice->start> $slice->end ) && $slice->is_circular ) {
         
       if ( !defined($end) || ($start > $end ) ) {
	   return $self->_fetch_by_Slice_start_end_strand_circular( $slice, $start, $end, $strand );
       }

       if ( defined($end) && ($end < 0) ) {
	   $end += $slice->seq_region_length;
       }
       
       if ($start < 0) {
           $start += $slice->seq_region_length;
       }

       if($slice->start> $slice->end) {
           return $self->_fetch_by_Slice_start_end_strand_circular( $slice, $slice->start, $slice->end, $strand );
       }
  }
  */
        
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

  int haveAllocedSlice = 0;
  if (rightExpand || leftExpand) {
    long junk;
    slice = Slice_expand(slice, leftExpand, rightExpand, 0, &junk, &junk);
    haveAllocedSlice = 1;
  }

  // retrieve normalized 'non-symlinked' slices
  // this allows us to support haplotypes and PARs
  SliceAdaptor *sliceAdaptor = (SliceAdaptor *)Slice_getAdaptor(slice);

  Vector *symProj = SliceAdaptor_fetchNormalizedSliceProjection(sliceAdaptor, slice, 0);
  Vector_setFreeFunc(symProj, ProjectionSegment_free);

  if (Vector_getNumElement(symProj) == 0) {
    fprintf(stderr,"Could not retrieve normalized Slices. Database contains incorrect assembly_exception information.\n");
    exit(1);
  }

  // call this method again with any slices that were 'symlinked' to by this
  // slice
  if (Vector_getNumElement(symProj) != 1 || ProjectionSegment_getToSlice((ProjectionSegment *)Vector_getElementAt(symProj,0)) != slice) {
    char *seq;
    if ((seq = calloc(Slice_getLength(slice)+2, sizeof(char))) == NULL) {
      fprintf(stderr,"Failed allocating seq\n");
      exit(1);
    }
    memset(seq, 'N', Slice_getLength(slice));

    int i;
    for (i=0; i<Vector_getNumElement(symProj); i++) {
      ProjectionSegment *segment = Vector_getElementAt(symProj, i);
      long segStart = ProjectionSegment_getFromStart(segment);
      long segEnd   = ProjectionSegment_getFromEnd(segment);

      Slice *symlinkSlice = ProjectionSegment_getToSlice(segment);
      // get sequence from each symlinked area
      char *tmpSeq = SequenceAdaptor_fetchBySliceStartEndStrandRecursive(sa, symlinkSlice, 1, POS_UNDEF, 1, recLev);
      // memcpy here
      memcpy(&(seq[segStart-1]), tmpSeq, segEnd-segStart+1);
      free(tmpSeq);
    }
// Will slice be the right length if there are rnaEdits???????? - disallow rnaEdits which are different length to original slice
    if (strand == -1) {
      SeqUtil_reverseComplement(seq, Slice_getLength(slice));
    }
    Vector_free(symProj);

    return seq;
  }
  Vector_free(symProj);

  // we need to project this slice onto the sequence coordinate system
  // even if the slice is in the same coord system, we want to trim out
  // flanking gaps (if the slice is past the edges of the seqregion)
  CoordSystemAdaptor *csa = DBAdaptor_getCoordSystemAdaptor(sa->dba);
  CoordSystem *seqLevelCs = CoordSystemAdaptor_fetchSeqLevel(csa);
  if (sa->rnaEditsCache == NULL && strcmp(Slice_getCoordSystemName(slice), "lrg") == 0) {
    SequenceAdaptor_initialiseRnaEdit(sa);
  }

  Vector *projection = Slice_project(slice, CoordSystem_getName(seqLevelCs), CoordSystem_getVersion(seqLevelCs));

  char *seq;
  if ((seq = calloc(Slice_getLength(slice)+2, sizeof(char))) == NULL) {
    fprintf(stderr,"Failed allocating seq\n");
    exit(1);
  }
  memset(seq, 'N', Slice_getLength(slice));

//  my $total = 0;
//  my $tmp_seq;

  // fetch sequence from each of the sequence regions projected onto
  int i;
  for (i=0; i<Vector_getNumElement(projection); i++) {
    ProjectionSegment *segment = Vector_getElementAt(projection, i);

    long start    = ProjectionSegment_getFromStart(segment);
    long end      = ProjectionSegment_getFromEnd(segment);
    Slice * seqSlice = ProjectionSegment_getToSlice(segment);
    long lenSegment = end - start + 1;

/*
    // check for gaps between segments and pad them with Ns
    my $gap = $start - $total - 1;
    if($gap) {
      $seq .= 'N' x $gap;
    }
*/

    IDType seqRegionId = SliceAdaptor_getSeqRegionId(sliceAdaptor, seqSlice);

    if (seqRegionId) {
      char *tmpSeq = SequenceAdaptor_fetchSeq(sa, seqRegionId,
                                              Slice_getStart(seqSlice), Slice_getLength(seqSlice));

      // reverse complement on negatively oriented slices
      if (Slice_getStrand(seqSlice) == -1) {
        SeqUtil_reverseComplement(tmpSeq, lenSegment);
      }

      if (start+lenSegment-1 > Slice_getLength(slice)) { 
        fprintf(stderr," segment off end of slice start = %ld lenSegment = %ld slice length = %ld len tmpSeq (%ld)\n", 
                start, lenSegment, Slice_getLength(slice),strlen(tmpSeq));
      }
      // Do with memcpy rather than strcpy
      // Think projection segment from coordinates are slice coordinates so relative to slice and starting at 1
      memcpy(&(seq[start-1]), tmpSeq, lenSegment);
      free(tmpSeq);
    } else {
      fprintf(stderr, "Error getting sequence region ID for slice");
    }

// Hopefully OK to free this here
    if (seqSlice != slice) {
      Slice_free(seqSlice);
    }

//    $seq .= $tmp_seq;
//    $total = $end;
  }
  Vector_setFreeFunc(projection, ProjectionSegment_free);
  Vector_free(projection);

/*
  // check for any remaining gaps at the end
  my $gap = $slice->length - $total;
  if ($gap) {
    $seq .= 'N' x $gap;
  }

  // if the sequence is too short it is because we came in with a seqlevel
  // slice that was partially off of the seq_region.  Pad the end with Ns
  // to make long enough
  if(length($seq) != $slice->length()) {
    $seq .= 'N' x ($slice->length() - length($seq));
  }
*/

  if (sa->rnaEditsCache != NULL && IDHash_contains(sa->rnaEditsCache, Slice_getSeqRegionId(slice))) {
    SequenceAdaptor_rnaEdit(sa, slice, &seq, *recLev);
  }

  // if they asked for the negative slice strand revcomp the whole thing
  if (strand == -1) {
    SeqUtil_reverseComplement(seq, Slice_getLength(slice));
  }

  if (haveAllocedSlice == 1) {
    Slice_free(slice);
  }

  return seq;
}


/* Not implementing circular!
sub _fetch_by_Slice_start_end_strand_circular {
  my ( $self, $slice, $start, $end, $strand ) = @_;

  assert_ref( $slice, 'Bio::EnsEMBL::Slice' );
  
  $strand ||= 1;
  if ( !defined($start) ) {
    $start ||= 1;
  }

  if ( !defined($end) ) {
      $end = $slice->end() - $slice->start() + 1;
  }

  if ( $start > $end && $slice->is_circular() ) {
    my ($seq, $seq1, $seq2);

    my $midpoint = $slice->seq_region_length - $slice->start + 1;
    $seq1 = ${ $self->_fetch_by_Slice_start_end_strand_circular( $slice, 1,  $midpoint, 1 )};
    $seq2 = ${ $self->_fetch_by_Slice_start_end_strand_circular( $slice, $midpoint + 1, $slice->length(), 1 )};

    $seq = $slice->strand > 0 ? "$seq1$seq2" : "$seq2$seq1";

    reverse_comp( \$seq ) if ( $strand == -1 );

    return \$seq;
  }



  # Get a new slice that spans the exact region to retrieve dna from
  my $right_expand = $end - $slice->length();    #negative is fine
  my $left_expand  = 1 - $start;                 #negative is fine

  if ( $right_expand || $left_expand ) {
    $slice =
        $slice->strand > 0
      ? $slice->expand( $left_expand,  $right_expand )
      : $slice->expand( $right_expand, $left_expand );
  }

  # Retrieve normalized 'non-symlinked' slices.  This allows us to
  # support haplotypes and PARs.
  my $slice_adaptor = $slice->adaptor();
  my @symproj =
    @{ $slice_adaptor->fetch_normalized_slice_projection($slice) };

  if ( @symproj == 0 ) {
    throw(   'Could not retrieve normalized Slices. Database contains '
           . 'incorrect assembly_exception information.' );
  }

  # Call this method again with any slices that were 'symlinked' to by
  # this slice.
  if ( @symproj != 1 || $symproj[0]->[2] != $slice ) {
    my $seq;
    foreach my $segment (@symproj) {
      my $symlink_slice = $segment->[2];

      # Get sequence from each symlinked area.
      $seq .= ${
        $self->fetch_by_Slice_start_end_strand( $symlink_slice, 1,
                                                undef, 1 ) };
    }
    if ( $strand == -1 ) {
      reverse_comp( \$seq );
    }

    return \$seq;
  }

  # We need to project this slice onto the sequence coordinate system
  # even if the slice is in the same coord system, we want to trim out
  # flanking gaps (if the slice is past the edges of the seqregion).
  my $csa      = $self->db->get_CoordSystemAdaptor();
  my $seqlevel = $csa->fetch_sequence_level();

  my @projection =
    @{ $slice->project( $seqlevel->name(), $seqlevel->version() ) };

  my $seq   = '';
  my $total = 0;
  my $tmp_seq;

  # Fetch sequence from each of the sequence regions projected onto.
  foreach my $segment (@projection) {
    my ( $start, $end, $seq_slice ) = @{$segment};

    # Check for gaps between segments and pad them with Ns
    my $gap = $start - $total - 1;
    if ($gap) {
      $seq .= 'N' x $gap;
    }

    my $seq_region_id = $slice_adaptor->get_seq_region_id($seq_slice);

    $tmp_seq = ${
      $self->_fetch_seq( $seq_region_id, $seq_slice->start(),
                         $seq_slice->length() ) };

    # Reverse compliment on negatively oriented slices.
    if ( $seq_slice->strand == -1 ) {
      reverse_comp( \$tmp_seq );
    }

    $seq .= $tmp_seq;

    $total = $end;
  }

  # Check for any remaining gaps at the end.
  my $gap = $slice->length() - $total;

  if ($gap) {
    $seq .= 'N' x $gap;
  }

  # If the sequence is too short it is because we came in with a
  # seqlevel slice that was partially off of the seq_region.  Pad the
  # end with Ns to make long enough
  if ( length($seq) != $slice->length() ) {
    $seq .= 'N' x ( $slice->length() - length($seq) );
  }

  if ( defined( $self->{_rna_edits_cache} )
       && defined(
            $self->{_rna_edits_cache}->{ $slice->get_seq_region_id } ) )
  {
    $self->_rna_edit( $slice, \$seq );
  }

  return \$seq;
} ## end sub _fetch_by_Slice_start_end_strand_circular

*/




void SequenceAdaptor_rnaEdit(SequenceAdaptor *sa, Slice *slice, char **seqPP, int recLev) {
  long sStart = Slice_getStart(slice);   // substr start at 0 , but seq starts at 1 (so no -1 here)

// Note: In C use Slice length rather than length($$seq) because length should be slice length and its quicker than doing a strlen on a huge string
//       ???? Maybe this isn't true if edits have been applied at a higher recursion level ????
//       Edits may not be different length so is true
  long sEnd   = sStart + Slice_getLength(slice); // length($$seq);

  //long lenString = Slice_getLength(slice);

  Vector *editsVec = IDHash_getValue(sa->rnaEditsCache, Slice_getSeqRegionId(slice));

  int i;
  for (i=0; i<Vector_getNumElement(editsVec); i++) {
    char *edit = Vector_getElementAt(editsVec, i);

    long start;
    long end;
    char *txt;
    
    char **tokens;
    int nTok;
    StrUtil_tokenize(&tokens, &nTok, edit);
    if (nTok != 3) {
      fprintf(stderr,"Wrong number of tokens parsing rnaEdit %s\n", edit);
      exit(1);
    }

    if (!StrUtil_isLongInteger(&start, tokens[0])) {
      fprintf(stderr,"Token 0 isn't an integer parsing rnaEdit %s\n", edit);
      exit(1);
    }

    if (!StrUtil_isLongInteger(&end, tokens[1])) {
      fprintf(stderr,"Token 1 isn't an integer parsing rnaEdit %s\n", edit);
      exit(1);
    }

    txt = tokens[2];

// check that RNA edit is not outside the requested region : happens quite often with LRG regions
    if (end < sStart) continue;
    if (sEnd < start) continue;

    int lenTxt = strlen(txt);
    int lenRange = (end-start+1);

    if (lenRange != lenTxt ) {
      fprintf(stderr, "ensc-core does not support rna edits on sequences where the replacement string is a different length to the original\n");
      exit(1);
    }

    // substr($$seq,$start-$s_start, ($end-$start)+1, $txt);
    // For equal length replacement string do in place edit
    if (lenTxt == lenRange) {
      memcpy(&((*seqPP)[start - sStart]), txt, lenTxt);
    } 
/* Unfinished support for unequal length ranges, but probably don't want to support this anyway!
    else { // Not equal length replacement string - if bigger need to realloc,  and then we need to move stuff around
      if (lenTxt > lenRange) {
        if ((*seqPP = realloc(*seqPP, )) == NULL) {
          fprintf(stderr,"Failed reallocating sequence string\n");
          exit(1);
        }
      }
      int lenDiff = lenTxt-lenRange;
      int indAfterEdit = end - sStart + 1;
      // Move the region after the edit
      memmove(&((*seqPP)[indAfterEdit + lenDiff]),  &((*seqPP)[indAfterEdit]), lenString-indAfterEdit+1);
      memcpy(&((*seqPP)[start - sStart]), txt, lenTxt);

      lenString += lenDiff;
    }
*/


    int j;
    for (j=0; j<nTok; j++) {
      free(tokens[j]);
    }
    free(tokens);
  }

  return;
}


char * SequenceAdaptor_fetchSeq(SequenceAdaptor *sa, IDType seqRegionId, long start, long length) {
  int status = 0;

  if (length < SEQ_CACHE_MAX) {
    long chunkMin = (start-1) >> SEQ_CHUNK_PWR;
    long chunkMax = (start + length - 1) >> SEQ_CHUNK_PWR;

    long minChunkMin = (chunkMin << SEQ_CHUNK_PWR) + 1;

    // piece together sequence from cached component parts

    char *entireSeq = NULL;
//    if ((entireSeq = calloc(((chunkMax-chunkMin+1) * (1<<SEQ_CHUNK_PWR))+1, sizeof(char))) == NULL) {
    if ((entireSeq = malloc((((chunkMax-chunkMin+1) * (1<<SEQ_CHUNK_PWR))+1) * sizeof(char))) == NULL) {
      fprintf(stderr,"Failed allocating entireSeq\n");
      exit(1);
    }

    int i;
    for (i = chunkMin; i <= chunkMax; i++) {
      char chunkKey[1024];
      sprintf(chunkKey, IDFMTSTR":%d",seqRegionId,i);

      long min = (i << SEQ_CHUNK_PWR) + 1;

      
      if (LRUCache_contains(sa->seqCache, chunkKey)) {
      //if (StringHash_contains(sa->seqCache, chunkKey)) {
 // What happens for length of last chunk???
        memcpy(&(entireSeq[min-minChunkMin]), LRUCache_get(sa->seqCache, chunkKey), LRUCache_getSize(sa->seqCache, chunkKey));
        //memcpy(&(entireSeq[min-minChunkMin]), StringHash_getValue(sa->seqCache, chunkKey), 1<<SEQ_CHUNK_PWR); 
        
      } else {
        // retrieve uncached portions of the sequence
        char qStr[1024];
        // Modified from perl to also return the length of the substring
        sprintf(qStr,"SELECT SUBSTRING(d.sequence, %ld, %ld) "
                     "FROM dna d "
                     "WHERE d.seq_region_id = "IDFMTSTR, min, 1L<<SEQ_CHUNK_PWR, seqRegionId);
//        sprintf(qStr,"SELECT SUBSTRING(d.sequence, %ld, %ld), LENGTH(SUBSTRING(d.sequence, %ld, %ld)) "
//                     "FROM dna d "
//                     "WHERE d.seq_region_id = "IDFMTSTR, min, 1L<<SEQ_CHUNK_PWR, min, 1L<<SEQ_CHUNK_PWR, seqRegionId);
//        sprintf(qStr,"SELECT SUBSTRING(d.sequence, %ld, %ld), if (sr.length-%ld+1 >  %ld, %ld, sr.length-%ld+1) "
//                     "FROM dna d, seq_region sr "
//                     "WHERE sr.seq_region_id = d.seq_region_id and d.seq_region_id = "IDFMTSTR, min, 1L<<SEQ_CHUNK_PWR, min, 1L<<SEQ_CHUNK_PWR, 1L<<SEQ_CHUNK_PWR, min, seqRegionId);

        StatementHandle *sth = sa->prepare((BaseAdaptor *)sa,qStr,strlen(qStr));

        sth->execute(sth);
        ResultRow *row = sth->fetchRow(sth);

        if (row) {
            char *tmpSeq   = row->getStringCopyAt(row, 0);
            //long lenTmpSeq = row->getLongAt(row, 1);
            long lenTmpSeq = strlen(tmpSeq);

            // always give back uppercased sequence so it can be properly softmasked
            //StrUtil_strupr(tmpSeq);

            memcpy(&(entireSeq[min-minChunkMin]), tmpSeq, lenTmpSeq);
            //StrUtil_appendString(entireSeq,tmpSeq);
            LRUCache_put(sa->seqCache, chunkKey, tmpSeq, free, lenTmpSeq);
            //StringHash_add(sa->seqCache, chunkKey, tmpSeq);
        }
        else {
          fprintf(stderr,"Failed to fetch sequence: SQL results empty");
          status = 1;
        }

        sth->finish(sth);
      }
    }

    if (status == 0) {
        // return only the requested portion of the entire sequence
        long min = ( chunkMin << SEQ_CHUNK_PWR ) + 1;
        //# my $max = ( $chunk_max + 1 ) << $SEQ_CHUNK_PWR;

        //seq = substr( $entire_seq, $start - $min, $length );
        // memmove it down and set '\0' at position length
        if (start-min != 0) {
          //memmove(entireSeq, &entireSeq[start-min], length);
          bcopy(&entireSeq[start-min], entireSeq, length);
        }
        entireSeq[length] = '\0';
    }

    return entireSeq;
  } else {
    // do not do any caching for requests of very large sequences

    char qStr[1024];
    sprintf(qStr,"SELECT SUBSTRING(d.sequence, %ld, %ld) "
                 "FROM dna d "
                 "WHERE d.seq_region_id = "IDFMTSTR, start, length, seqRegionId);

    StatementHandle *sth = sa->prepare((BaseAdaptor *)sa,qStr,strlen(qStr));

    sth->execute(sth);
    ResultRow *row = sth->fetchRow(sth);
    char *tmpSeq = NULL;
    
    if (row)
      tmpSeq = row->getStringCopyAt(row, 0);

    sth->finish(sth);
    // always give back uppercased sequence so it can be properly softmasked
    //StrUtil_strupr(tmpSeq);

    return tmpSeq;
  }
}


/*
=head2 store

  Arg [1]    : int $seq_region_id the id of the sequence region this dna
               will be associated with.
  Arg [2]    : string $sequence the dna sequence to be stored 
               in the database.  Note that the sequence passed in will be
               converted to uppercase.
  Example    : $seq_adaptor->store(11, 'ACTGGGTACCAAACAAACACAACA');
  Description: stores a dna sequence in the databases dna table and returns the
               database identifier for the new record.
  Returntype : none
  Exceptions : throw if the database insert fails
  Caller     : sequence loading scripts
  Status     : Stable

=cut

sub store {
  my ($self, $seq_region_id, $sequence) = @_;

  if(!$seq_region_id) {
    throw('seq_region_id is required');
  }

  $sequence = uc($sequence);

  my $statement = 
    $self->prepare("INSERT INTO dna(seq_region_id, sequence) VALUES(?,?)");

  $statement->bind_param(1,$seq_region_id,SQL_INTEGER);
  $statement->bind_param(2,$sequence,SQL_LONGVARCHAR);
  $statement->execute();

  $statement->finish();

  return;
}

*/


/*
=head2 fetch_by_assembly_location

  Description: DEPRECATED use fetch_by_Slice_start_end_strand() instead.

=cut

sub fetch_by_assembly_location {
   my ( $self, $chrStart, $chrEnd, 
        $strand, $chrName, $assemblyType ) = @_;

   deprecate('Use fetch_by_Slice_start_end_strand() instead');

   my $csa = $self->db->get_CoordSystem();
   my $top_cs = @{$csa->fetch_all};

   my $slice_adaptor = $self->db->get_SliceAdaptor();
   my $slice = $slice_adaptor->fetch_by_region($top_cs->name(), $chrName,
                                               $chrStart, $chrEnd,
                                               $strand, $top_cs->version);

   return $self->fetch_by_Slice_start_end_strand($slice,1, $slice->length,1);
}


=head2 fetch_by_RawContig_start_end_strand

  Description: DEPRECATED use fetch_by_Slice_start_end_strand instead

=cut

sub fetch_by_RawContig_start_end_strand {
  deprecate('Use fetch_by_Slice_start_end_strand instead.');
  fetch_by_Slice_start_end_strand(@_);
}


*/

void SequenceAdaptor_initialiseRnaEdit(SequenceAdaptor *sa) {
  char *qStr = "SELECT sra.seq_region_id, sra.value "
               "FROM seq_region_attrib sra, attrib_type at "
               "WHERE sra.attrib_type_id = at.attrib_type_id "
               "AND code like '_rna_edit'";

  StatementHandle *sth = sa->prepare((BaseAdaptor *)sa,qStr,strlen(qStr));

  sth->execute(sth);

  IDHash *edits = IDHash_new(IDHASH_SMALL);

  int count = 0;
  ResultRow *row;
  while ((row = sth->fetchRow(sth))){
    IDType seqRegionId = row->getLongLongAt(row, 0);
    char *value        = row->getStringAt(row, 1);

    count++;
    if (! IDHash_contains(edits, seqRegionId)) {
      IDHash_add(edits, seqRegionId, Vector_new());
    }
    Vector *vec = IDHash_getValue(edits, seqRegionId);

    char *tmp;
    Vector_addElement(vec, StrUtil_copyString(&tmp, value, 0));
  }
  sth->finish(sth);

  if (count) {
    sa->rnaEditsCache = edits;
  } else {
    IDHash_free(edits, NULL);
  }
}
