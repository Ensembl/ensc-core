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

#define __GENOMICALIGN_MAIN__
#include "GenomicAlign.h"
#undef __GENOMICALIGN_MAIN__
#include <stdlib.h>

#include "StrUtil.h"
#include "SeqUtil.h"
#include "CigarStrUtil.h"
#include "Slice.h"

GenomicAlign *GenomicAlign_new() {
  GenomicAlign *ga;

  if ((ga = (GenomicAlign *)calloc(1,sizeof(GenomicAlign))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for ga\n");
    return NULL;
  }

  ga->objectType = CLASS_GENOMICALIGN;

  ga->funcs = &genomicAlignFuncs;

  Object_incRefCount(ga);
  return ga;
}

char *GenomicAlign_setCigarString(GenomicAlign *ga, char *cigStr) {
  StrUtil_copyString(&(ga->cigarString), cigStr, 0);

  return ga->cigarString;
}

ECOSTRING GenomicAlign_setAlignmentType(GenomicAlign *ga, char *alType) {
  EcoString_copyStr(ecoSTable, &(ga->alignmentType), alType, 0);

  return ga->alignmentType;
}

/*
=head2 sequence_align_string

  Arg  1     : Bio::EnsEMBL::Slice $consensus_slice
               A slice covering the conensus area of this alignment
  Arg  2     : Bio::EnsEMBL::Slice $query_slice
               A slice covering the query_area of this alignment
  Arg  3..   : list of String $flags
               FRAG_SLICES = slices cover the dna frags
               ALIGN_SLICES = slices cover just the aligned area
               FIX_CONSENSUS = dont put dashes in consensus sequence on 
               alignment printout
               FIX_QUERY = dont put dashes in query sequence on alignment
               CONSENSUS = return the consensus aligned sequence
               QUERY = return the query aligned sequence
  Example    : none
  Description: returns representations of the aligned sequences according to
               the flags.
  Returntype : String
  Exceptions : none
  Caller     : general

=cut
*/

char * GenomicAlign_getSequenceAlignString(GenomicAlign *ga, Slice *consensusSlice, Slice *querySlice, unsigned int flags) {
  Vector *cig = CigarStrUtil_getPieces(GenomicAlign_getCigarString(ga));
  char *cSeq, *qSeq;
  char *rSeq = NULL;
  int cPos = 0;
  int qPos = 0;
  int cigCount;
  char cigType;
  int i;
  int ok = 1;

  // here fill cseq and qseq with the aligned area sequence

  if (flags & GENOMICALIGN_ALIGNSLICES) {
    cSeq = Slice_getSeq(consensusSlice);
    qSeq = Slice_getSeq(querySlice);
  } else {
    cSeq = Slice_getSubSeq(consensusSlice, 
                           GenomicAlign_getConsensusStart(ga),
                           GenomicAlign_getConsensusEnd(ga), 
                           1);
    qSeq = Slice_getSubSeq(querySlice,
                           GenomicAlign_getQueryStart(ga),
                           GenomicAlign_getQueryEnd(ga), 
                           GenomicAlign_getQueryStrand(ga));
  }

  // rSeq - result sequence
  rSeq = StrUtil_copyString(&rSeq, "", 0);

  for (i=0; i<Vector_getNumElement(cig); i++) {
    char *cigElem = Vector_getElementAt(cig,i);
    int cigElemLen = strlen(cigElem);

    if (!cigElemLen) {
      fprintf(stderr,"Error: Empty cig element in string %s\n",GenomicAlign_getCigarString(ga));
      ok = 0;
      break;
    }

    cigType = cigElem[cigElemLen-1];
    cigElem[cigElemLen-1] = '\0';
    if (cigElemLen > 1) {
      cigCount = atoi(cigElem);
    } else {
      cigCount = 1;
    }

    switch (cigType) {
      case 'M':
        if (flags & GENOMICALIGN_CONSENSUS) {
          rSeq = StrUtil_appendNString(rSeq,&(cSeq[cPos]),cigCount);
        } else {
          rSeq = StrUtil_appendNString(rSeq,&(qSeq[qPos]),cigCount);
        }
        cPos += cigCount;
        qPos += cigCount;
        break;

      case 'D':
        if (flags & GENOMICALIGN_CONSENSUS) {
          if (!(flags & GENOMICALIGN_FIXCONSENSUS)) {
            rSeq = SeqUtil_addGaps(rSeq, cigCount);
          }
        } else {
          if (!(flags & GENOMICALIGN_FIXCONSENSUS)) {
            rSeq = StrUtil_appendNString(rSeq,&(qSeq[qPos]),cigCount);
          }
        }
  
        qPos += cigCount;
        break;

      case 'I':
        if (flags & GENOMICALIGN_CONSENSUS) {
          if (!(flags & GENOMICALIGN_FIXQUERY)) {
            rSeq = StrUtil_appendNString(rSeq,&(cSeq[cPos]),cigCount);
          }
        } else {
          if (!(flags & GENOMICALIGN_FIXQUERY)) {
            rSeq = SeqUtil_addGaps(rSeq, cigCount);
          }
        }
       
        cPos += cigCount;
        break;

      default:
        fprintf(stderr,"Error: Unknown cigar elem type %c\n",cigType);
        ok = 0;
        break;
    }
  }     
  Vector_free(cig);

  if (!ok && rSeq) {
    free(rSeq);
    rSeq = NULL;
  }
    

  return rSeq;
}

void GenomicAlign_free(GenomicAlign *ga) {
  Object_decRefCount(ga);

  if (Object_getRefCount(ga) > 0) {
    return;
  } else if (Object_getRefCount(ga) < 0) {
    fprintf(stderr,"Error: Negative reference count for GenomicAlign\n"
                   "       Freeing it anyway\n");
  }

  free(ga->cigarString);

  if (ga->consensusDNAFrag) DNAFrag_free(ga->consensusDNAFrag);
  if (ga->queryDNAFrag)     DNAFrag_free(ga->queryDNAFrag);

  if (ga->alignmentType) EcoString_freeStr(ecoSTable, ga->alignmentType);

  free(ga);
}
