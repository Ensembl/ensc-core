#ifndef __REFINESOLEXAGENES_H__
#define __REFINESOLEXAGENES_H__

#include "DBAdaptor.h"
#include "Vector.h"
#include "Slice.h"
#include "StringHash.h"
#include "SliceAdaptor.h"
#include "Analysis.h"
#include "Transcript.h"

#include "sam.h"
#include "bam.h"

#define RSGEXON_RETAINED 1<<1
#define RSGEXON_EXTRA    1<<2

typedef struct RefineSolexaGenesStruct {
  char *badModelsType;
  char *bestScoreType;
  char *inputId;
  char *intronDb;
  char *modelDb;
  char *modelLogicName;
  char *outputDb;
  char *otherIsoformsType;
  char *singleExonModelType;

  int max3PrimeExons;
  int max5PrimeExons;
  int max3PrimeLength;
  int max5PrimeLength;
  int maxIntronSize;
  int maxNum;
  int maxRecursions;
  int minIntronSize;
  int minSingleExonCDSLength;
  int minSingleExonLength;
  int otherNum;
  int recursiveLimit;
  int strictInternalSpliceSites;
  int strictInternalEndSpliceSites;
  int trimUtr;
  int writeIntrons;

  SliceAdaptor *geneSliceAdaptor;
  SliceAdaptor *intronSliceAdaptor;

  Slice *chrSlice;

  StringHash *extraExons;

  Vector *intronBamFiles;
  Vector *intronFeatures;
  Vector *logicNames;
  Vector *output;
  Vector *prelimGenes;
  
  double filterOnOverlapThreshold;
  double rejectIntronCutoff;
  double retainedIntronPenalty;

  StringHash *adaptorAliasHash;

  Analysis *analysis;

  DBAdaptor *db;
} RefineSolexaGenes;

typedef struct ModelClusterStruct {
  long start;
  long end;
  int strand;
  Vector *models; // A vector of Model objects
  Vector *finalModels; // A vector of Gene objects
} ModelCluster;


typedef struct IntronBamConfigStruct {
  int depth;
  int mixedBam;
  Vector *groupNames;
  char *fileName;
} IntronBamConfig;

typedef struct ORFRangeStruct {
  long length;
  long start;
  long end;
} ORFRange;



RefineSolexaGenes *RefineSolexaGenes_new(char *configFile);
DBAdaptor *RefineSolexaGenes_getDbAdaptor(RefineSolexaGenes *rsg, char *alias);
void RefineSolexaGenes_fetchInput(RefineSolexaGenes *rsg);
void  RefineSolexaGenes_run(RefineSolexaGenes *rsg);
void RefineSolexaGenes_refineGenes(RefineSolexaGenes *rsg);
Analysis *RefineSolexaGenes_createAnalysisObject(RefineSolexaGenes *rsg, char *logicName);
Vector *RefineSolexaGenes_reclusterModels(RefineSolexaGenes *rsg, Vector *clusters, Vector **retNewClusters);
ModelCluster *RefineSolexaGenes_recalculateCluster(RefineSolexaGenes *rsg, Vector *genes);
void RefineSolexaGenes_filterModels(RefineSolexaGenes *rsg, Vector *clusters);
Vector *RefineSolexaGenes_makeModels(RefineSolexaGenes *rsg, StringHash *paths, int strand, Vector *exons, Gene *gene, StringHash *intronHash);
Transcript *RefineSolexaGenes_modifyTranscript(RefineSolexaGenes *rsg, Transcript *tran, Vector *exons);
void RefineSolexaGenes_writeOutput(RefineSolexaGenes *rsg);
int RefineSolexaGenes_processTree(RefineSolexaGenes *rsg, StringHash *hashref, char *index, char *soFar, StringHash *paths);
StringHash *RefineSolexaGenes_processPaths(RefineSolexaGenes *rsg, Vector *exons, Vector *exonIntron, StringHash *intronExon, int strict, int *giveUpFlag);
Vector *RefineSolexaGenes_makeModelClusters(RefineSolexaGenes *rsg, Vector *models, int strand);
Vector *RefineSolexaGenes_mergeExons(RefineSolexaGenes *rsg, Gene *gene, int strand);
Exon *RefineSolexaGenes_binSearchForOverlap(RefineSolexaGenes *rsg, Vector *exons, int pos);
void RefineSolexaGenes_bamToIntronFeatures(RefineSolexaGenes *rsg, IntronBamConfig *intronBamConf, samfile_t *sam, bam_index_t *idx, int ref, int begRange, int endRange);
Vector *RefineSolexaGenes_getUngappedFeatures(RefineSolexaGenes *rsg, bam1_t *b);
void RefineSolexaGenes_dnaToIntronFeatures(RefineSolexaGenes *rsg, long start, long end);
Vector *RefineSolexaGenes_fetchIntronFeatures(RefineSolexaGenes *rsg, long start, long end, long *offsetP);
Exon *RefineSolexaGenes_makeExon(RefineSolexaGenes *rsg, long start, long end, double score, char *diplayId);
void RefineSolexaGenes_setRecursiveLimit(RefineSolexaGenes *rsg, int limit);
int RefineSolexaGenes_getRecursiveLimit(RefineSolexaGenes *rsg);
SliceAdaptor *RefineSolexaGenes_getGeneSliceAdaptor(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setGeneSliceAdaptor(RefineSolexaGenes *rsg, SliceAdaptor *sa);
SliceAdaptor *RefineSolexaGenes_getIntronSliceAdaptor(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setIntronSliceAdaptor(RefineSolexaGenes *rsg, SliceAdaptor *sa);
Vector *RefineSolexaGenes_getPrelimGenes(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setPrelimGenes(RefineSolexaGenes *rsg, Vector *genes);
Slice *RefineSolexaGenes_getChrSlice(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setChrSlice(RefineSolexaGenes *rsg, Slice *slice);
Vector *RefineSolexaGenes_getIntronFeatures(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setIntronFeatures(RefineSolexaGenes *rsg, Vector *features);
StringHash *RefineSolexaGenes_getExtraExons(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setExtraExons(RefineSolexaGenes *rsg, StringHash *extraExons);
void RefineSolexaGenes_setIntronDb(RefineSolexaGenes *rsg, char *intronDb);
char *RefineSolexaGenes_getIntronDb(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setOutputDb(RefineSolexaGenes *rsg, char *outputDb);
char *RefineSolexaGenes_getOutputDb(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setModelDb(RefineSolexaGenes *rsg, char *modelDb);
char *RefineSolexaGenes_getModelDb(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setLogicNames(RefineSolexaGenes *rsg, Vector *logicNames);
Vector *RefineSolexaGenes_getLogicNames(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setRetainedIntronPenalty(RefineSolexaGenes *rsg, double retainedIntronPenalty);
double RefineSolexaGenes_getRetainedIntronPenalty(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setMinIntronSize(RefineSolexaGenes *rsg, int minIntronSize);
int RefineSolexaGenes_getMinIntronSize(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setMaxIntronSize(RefineSolexaGenes *rsg, int maxIntronSize);
int RefineSolexaGenes_getMaxIntronSize(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setBestScoreType(RefineSolexaGenes *rsg, char *bestScoreType);
char *RefineSolexaGenes_getBestScoreType(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setOtherNum(RefineSolexaGenes *rsg, int otherNum);
int RefineSolexaGenes_getOtherNum(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setOtherIsoformsType(RefineSolexaGenes *rsg, char *otherIsoformsType);
char *RefineSolexaGenes_getOtherIsoformsType(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setModelLogicName(RefineSolexaGenes *rsg, char *modelLN);
char *RefineSolexaGenes_getModelLogicName(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setBadModelsType(RefineSolexaGenes *rsg, char *badModelsType);
char *RefineSolexaGenes_getBadModelsType(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setMaxNum(RefineSolexaGenes *rsg, int maxNum);
int RefineSolexaGenes_getMaxNum(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setMaxRecursions(RefineSolexaGenes *rsg, int maxRecursions);
int RefineSolexaGenes_getMaxRecursions(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setMinSingleExonLength(RefineSolexaGenes *rsg, int minSingleExonLength);
int RefineSolexaGenes_getMinSingleExonLength(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setMinSingleExonCDSLength(RefineSolexaGenes *rsg, int minSingleExonCDSLength);
int RefineSolexaGenes_getMinSingleExonCDSLength(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setSingleExonModelType(RefineSolexaGenes *rsg, char *singleExonModelType);
char *RefineSolexaGenes_getSingleExonModelType(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setStrictInternalSpliceSites(RefineSolexaGenes *rsg, int strictInternalSpliceSites);
int RefineSolexaGenes_strictInternalSpliceSites(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setStrictInternalEndSpliceSites(RefineSolexaGenes *rsg, int strictInternalEndSpliceSites);
int RefineSolexaGenes_strictInternalEndSpliceSites(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setIntronBamFiles(RefineSolexaGenes *rsg, Vector *intronBamFiles);
Vector *RefineSolexaGenes_getIntronBamFiles(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setWriteIntrons(RefineSolexaGenes *rsg, int writeIntrons);
int RefineSolexaGenes_writeIntrons(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setTrimUtr(RefineSolexaGenes *rsg, int trimUtr);
int RefineSolexaGenes_trimUtr(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setMax3PrimeExons(RefineSolexaGenes *rsg, int max3PrimeExons);
int RefineSolexaGenes_getMax3PrimeExons(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setMax3PrimeLength(RefineSolexaGenes *rsg, int max3PrimeLength);
int RefineSolexaGenes_getMax3PrimeLength(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setMax5PrimeExons(RefineSolexaGenes *rsg, int max5PrimeExons);
int RefineSolexaGenes_getMax5PrimeExons(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setMax5PrimeLength(RefineSolexaGenes *rsg, int max5PrimeLength);
int RefineSolexaGenes_getMax5PrimeLength(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setFilterOnOverlapThreshold(RefineSolexaGenes *rsg, int filterOnOverlapThreshold);
int RefineSolexaGenes_getFilterOnOverlapThreshold(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setRejectIntronCutoff(RefineSolexaGenes *rsg, double rejectIntronCutoff);
double RefineSolexaGenes_getRejectIntronCutoff(RefineSolexaGenes *rsg);


// To move
  Transcript *TranslationUtils_addORFToTranscript(ORFRange *orf, Transcript *transcript);
  Transcript *TranslationUtils_computeTranslation(Transcript *tran);
  Gene *TranscriptUtils_convertToGene(Transcript *t, Analysis *analysis, char *biotype);
  Exon *ExonUtils_createExon(long start, long end, int phase, int endPhase, int strand, Analysis *analysis, Vector *supportingFeatures, IDType dbId, Slice *slice, char *stableId, int version);
  Exon *ExonUtils_cloneExon(Exon *exon);
  BaseAlignFeature *EvidenceUtils_cloneEvidence(BaseAlignFeature *feature);
  Slice *RefineSolexaGenes_fetchSequence(RefineSolexaGenes *rsg, char *name, DBAdaptor *db, Vector *repeatMaskTypes, int softMask);
  void RefineSolexaGenes_setInputId(RefineSolexaGenes *rsg, char *inputId);
  char *RefineSolexaGenes_getInputId(RefineSolexaGenes *rsg);
  Vector *RefineSolexaGenes_getOutput(RefineSolexaGenes *rsg);
  void RefineSolexaGenes_addToOutput(RefineSolexaGenes *rsg, Gene *gene);
  void RefineSolexaGenes_setDb(RefineSolexaGenes *rsg, DBAdaptor *db);
  DBAdaptor *RefineSolexaGenes_getDb(RefineSolexaGenes *rsg);
  Vector *TranslationUtils_generateORFRanges(Transcript *transcript, int requireMet, int minLength);

  Analysis *RefineSolexaGenes_setAnalysis(RefineSolexaGenes *rsg, Analysis *analysis);
  Analysis *RefineSolexaGenes_getAnalysis(RefineSolexaGenes *rsg);

#endif
