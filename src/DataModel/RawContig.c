#include "RawContig.h"

#include "DNAAlignFeatureAdaptor.h"
#include "GeneAdaptor.h"
#include "PredictionTranscriptAdaptor.h"
#include "ProteinAlignFeatureAdaptor.h"
#include "RawContigAdaptor.h"
#include "RepeatFeatureAdaptor.h"
#include "SimpleFeatureAdaptor.h"

RawContig *RawContig_new() {
  RawContig *rc;

  if ((rc = (RawContig *)calloc(1,sizeof(RawContig))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for rc\n");
    return NULL;
  }

  rc->length = rc->emblOffset = rc->cloneId = -1;

  BaseContig_setContigType(rc,RAWCONTIG);

  return rc;
}

char *RawContig_setName(RawContig *rc, char *name) {
  if ((rc->name = (char *)malloc(strlen(name)+1)) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for seqfeature name\n");
    return NULL;
  }

  strcpy(rc->name,name);

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

char *RawContig_getName(RawContig *rc) {
  RawContigAdaptor *rca = (RawContigAdaptor *)RawContig_getAdaptor(rc);

  if (rc->name == NULL && rca) {
    RawContigAdaptor_fetchAttributes(rca,rc);
  }
  return rc->name;
}

Set *RawContig_getAllSimpleFeatures(RawContig *rc, char *logicName, double *scoreP) {
  RawContigAdaptor *rca = (RawContigAdaptor *)RawContig_getAdaptor(rc);
  SimpleFeatureAdaptor *sfa;

  if (!rca) {
    fprintf(stderr, "Warning: Contig has no adaptor - cannot retrieve simple features\n");
    return emptySet;
  }
  sfa = DBAdaptor_getSimpleFeatureAdaptor(rca->dba);
  return SimpleFeatureAdaptor_fetchAllByRawContigAndScore(sfa,rc,scoreP,logicName);
}

Set *RawContig_getAllPredictionTranscripts(RawContig *rc, char *logicName) {
  RawContigAdaptor *rca = (RawContigAdaptor *)RawContig_getAdaptor(rc);
  PredictionTranscriptAdaptor *pta;

  if (!rca) {
    fprintf(stderr, "Warning: Contig has no adaptor - cannot retrieve prediction transcripts\n");
    return emptySet;
  }
  pta = DBAdaptor_getPredictionTranscriptAdaptor(rca->dba);
  return PredictionTranscriptAdaptor_fetchAllByRawContig(pta,rc,logicName);
}

Set *RawContig_getAllRepeatFeatures(RawContig *rc, char *logicName) {
  RawContigAdaptor *rca = (RawContigAdaptor *)RawContig_getAdaptor(rc);
  RepeatFeatureAdaptor *rfa;

  if (!rca) {
    fprintf(stderr, "Warning: Contig has no adaptor - cannot retrieve repeat features\n");
    return emptySet;
  }
  rfa = DBAdaptor_getRepeatFeatureAdaptor(rca->dba);
  return RepeatFeatureAdaptor_fetchAllByRawContig(rfa,rc,logicName);
}


Set *RawContig_getAllDNAAlignFeatures(RawContig *rc, char *logicName, double *scoreP) {
  RawContigAdaptor *rca = (RawContigAdaptor *)RawContig_getAdaptor(rc);
  DNAAlignFeatureAdaptor *dafa;

  if (!rca) {
    fprintf(stderr, "Warning: Contig has no adaptor - cannot retrieve dna align features\n");
    return emptySet;
  }
  dafa = DBAdaptor_getDNAAlignFeatureAdaptor(rca->dba);
  return DNAAlignFeatureAdaptor_fetchAllByRawContigAndScore(dafa,rc,scoreP,logicName);
}

Set *RawContig_getAllProteinAlignFeatures(RawContig *rc, char *logicName, double *scoreP) {
  RawContigAdaptor *rca = (RawContigAdaptor *)RawContig_getAdaptor(rc);
  ProteinAlignFeatureAdaptor *pafa;

  if (!rca) {
    fprintf(stderr, "Warning: Contig has no adaptor - cannot retrieve pep align features\n");
    return emptySet;
  }
  pafa = DBAdaptor_getProteinAlignFeatureAdaptor(rca->dba);
  return ProteinAlignFeatureAdaptor_fetchAllByRawContigAndScore(pafa,rc,scoreP,logicName);
}
