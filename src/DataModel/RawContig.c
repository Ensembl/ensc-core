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

#define __RAWCONTIG_MAIN__
#include "RawContig.h"
#undef __RAWCONTIG_MAIN__

#include "DNAAlignFeatureAdaptor.h"
#include "GeneAdaptor.h"
#include "PredictionTranscriptAdaptor.h"
#include "ProteinAlignFeatureAdaptor.h"
#include "RawContigAdaptor.h"
#include "RepeatFeatureAdaptor.h"
#include "SimpleFeatureAdaptor.h"
#include "SequenceAdaptor.h"
#include "DBAdaptor.h"

#include "StrUtil.h"
#include "SeqUtil.h"

RawContig *RawContig_new() {
  RawContig *rc;

  if ((rc = (RawContig *)calloc(1,sizeof(RawContig))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for rc\n");
    return NULL;
  }

  rc->length = rc->emblOffset = rc->cloneId = -1;

  rc->objectType = CLASS_RAWCONTIG;
  Object_incRefCount(rc);

  rc->funcs = &rawContigFuncs;

  return rc;
}

ECOSTRING RawContig_getName(RawContig *rc) {
  RawContigAdaptor *rca = (RawContigAdaptor *)RawContig_getAdaptor(rc);

  if (rc->name == NULL && rca) {
    RawContigAdaptor_fetchAttributes(rca,rc);
  }
  return rc->name;
}

ECOSTRING RawContig_setName(RawContig *rc, char *name) {
  if (!EcoString_copyStr(ecoSTable,&(rc->name),name,0)) {
    fprintf(stderr,"ERROR: Failed allocating space for contig name\n");
    return NULL;
  }

  return rc->name;
}

int RawContig_getLength(RawContig *rc) {
  RawContigAdaptor *rca = (RawContigAdaptor *)RawContig_getAdaptor(rc);

  if (rc->length == -1 && rca) {
    RawContigAdaptor_fetchAttributes(rca,rc);
  }
  return rc->length;
}

int RawContig_getEMBLOffset(RawContig *rc) {
  RawContigAdaptor *rca = (RawContigAdaptor *)RawContig_getAdaptor(rc);

  if (rc->emblOffset == -1 && rca) {
    RawContigAdaptor_fetchAttributes(rca,rc);
  }
  return rc->emblOffset;
}

long RawContig_getCloneID(RawContig *rc) {
  RawContigAdaptor *rca = (RawContigAdaptor *)RawContig_getAdaptor(rc);

  if (rc->cloneId == -1 && rca) {
    RawContigAdaptor_fetchAttributes(rca,rc);
  }
  return rc->cloneId;
}

#ifdef OLD
Vector *RawContig_getAllSimpleFeatures(RawContig *rc, char *logicName, double *scoreP) {
  RawContigAdaptor *rca = (RawContigAdaptor *)RawContig_getAdaptor(rc);
  SimpleFeatureAdaptor *sfa;

  if (!rca) {
    fprintf(stderr, "Warning: Contig has no adaptor - cannot retrieve simple features\n");
    return emptyVector;
  }
  sfa = DBAdaptor_getSimpleFeatureAdaptor(rca->dba);
  return SimpleFeatureAdaptor_fetchAllByRawContigAndScore(sfa,rc,scoreP,logicName);
}

Vector *RawContig_getAllPredictionTranscripts(RawContig *rc, char *logicName) {
  RawContigAdaptor *rca = (RawContigAdaptor *)RawContig_getAdaptor(rc);
  PredictionTranscriptAdaptor *pta;

  if (!rca) {
    fprintf(stderr, "Warning: Contig has no adaptor - cannot retrieve prediction transcripts\n");
    return emptyVector;
  }
  pta = DBAdaptor_getPredictionTranscriptAdaptor(rca->dba);
  return PredictionTranscriptAdaptor_fetchAllByRawContig(pta,rc,logicName);
}

Vector *RawContig_getAllRepeatFeatures(RawContig *rc, char *logicName) {
  RawContigAdaptor *rca = (RawContigAdaptor *)RawContig_getAdaptor(rc);
  RepeatFeatureAdaptor *rfa;

  if (!rca) {
    fprintf(stderr, "Warning: Contig has no adaptor - cannot retrieve repeat features\n");
    return emptyVector;
  }
  rfa = DBAdaptor_getRepeatFeatureAdaptor(rca->dba);
  return RepeatFeatureAdaptor_fetchAllByRawContig(rfa,rc,logicName);
}


Vector *RawContig_getAllDNAAlignFeatures(RawContig *rc, char *logicName, double *scoreP) {
  RawContigAdaptor *rca = (RawContigAdaptor *)RawContig_getAdaptor(rc);
  DNAAlignFeatureAdaptor *dafa;

  if (!rca) {
    fprintf(stderr, "Warning: Contig has no adaptor - cannot retrieve dna align features\n");
    return emptyVector;
  }
  dafa = DBAdaptor_getDNAAlignFeatureAdaptor(rca->dba);
  return DNAAlignFeatureAdaptor_fetchAllByRawContigAndScore(dafa,rc,scoreP,logicName);
}

Vector *RawContig_getAllProteinAlignFeatures(RawContig *rc, char *logicName, double *scoreP) {
  RawContigAdaptor *rca = (RawContigAdaptor *)RawContig_getAdaptor(rc);
  ProteinAlignFeatureAdaptor *pafa;

  if (!rca) {
    fprintf(stderr, "Warning: Contig has no adaptor - cannot retrieve pep align features\n");
    return emptyVector;
  }
  pafa = DBAdaptor_getProteinAlignFeatureAdaptor(rca->dba);
  return ProteinAlignFeatureAdaptor_fetchAllByRawContigAndScore(pafa,rc,scoreP,logicName);
}
#endif

char *RawContig_setSeq(RawContig *contig, char *seq) {
  // Sequence can be set manually
  StrUtil_copyString(&(contig->seq),seq,0);

  return contig->seq;
}

char *RawContig_getSeq(RawContig *contig) {

  if (contig->seq) {
    return contig->seq;
  }

  //or retrieved from the database
  if (RawContig_getAdaptor(contig)) {
    DBAdaptor *dnadba;
    SequenceAdaptor *sa;

    
    dnadba = RawContig_getAdaptor(contig)->dba->dnadb;
    sa = DBAdaptor_getSequenceAdaptor(dnadba);
    //return SequenceAdaptor_fetchByRawContigStartEndStrand(sa, RawContig_getDbID(contig), 1, -1, 1);
    return NULL;
  }

  fprintf(stderr,"Warning: RawContig seq not set, and no db is available\n");

  return emptyString;
}

char *RawContig_getSubSeq(RawContig *contig, int start, int end, int strand) { 
  SequenceAdaptor *sa;

  if (end < start) {
    fprintf(stderr, "Error: End coord is less then start coord to call on RawContig subseq.");
    exit(1);
  }

  if (strand != -1 && strand != 1) {
    strand = 1;
  }

  // if the sequence of this contig has been manually set retrieve its substring
  if (contig->seq) {
    char *str = StrUtil_substr(contig->seq, start-1, end - start + 1);

    if (strand == -1) {
      SeqUtil_reverseComplement(str,strlen(str));
    }
    return str;
  }

  if (!RawContig_getAdaptor(contig)) {
    fprintf(stderr, "Error: RawContig subseq no sequence set and no db available\n");
    return emptyString;
  }

  sa = DBAdaptor_getSequenceAdaptor(RawContig_getAdaptor(contig)->dba->dnadb);

  //return SequenceAdaptor_fetchByRawContigStartEndStrand(sa, RawContig_getDbID(contig), start, end, strand);
  return NULL;

}

void RawContig_free(RawContig *rc) {
  Object_decRefCount(rc);

  if (Object_getRefCount(rc) > 0) {
    return;
  } else if (Object_getRefCount(rc) < 0) {
    fprintf(stderr,"Error: Negative reference count for RawContig\n"
                   "       Freeing it anyway\n");
  }

  printf("Freeing rawcontig\n");
  BaseContig_freePtrs((BaseContig*)rc);

  free(rc);
}

