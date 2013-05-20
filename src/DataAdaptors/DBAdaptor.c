#include "DBAdaptor.h"
#include "MetaContainer.h"

#include "AnalysisAdaptor.h"
#include "AssemblyMapperAdaptor.h"
#include "ChromosomeAdaptor.h"
#include "CoordSystemAdaptor.h"
#include "CloneAdaptor.h"
#include "DBEntryAdaptor.h"
#include "DNAAlignFeatureAdaptor.h"
#include "ExonAdaptor.h"
#include "GeneAdaptor.h"
#include "PredictionTranscriptAdaptor.h"
#include "ProteinAlignFeatureAdaptor.h"
#include "RawContigAdaptor.h"
#include "RepeatConsensusAdaptor.h"
#include "RepeatFeatureAdaptor.h"
#include "SequenceAdaptor.h"
#include "SimpleFeatureAdaptor.h"
#include "SliceAdaptor.h"
#include "SupportingFeatureAdaptor.h"
#include "TranscriptAdaptor.h"
#include "TranslationAdaptor.h"
#include "MetaCoordContainer.h"
#include "TranscriptSupportingFeatureAdaptor.h"

#include "StrUtil.h"
#include "SeqRegionCacheEntry.h"

#include "ProcUtil.h"

DBAdaptor *DBAdaptor_new(char *host, char *user, char *pass, char *dbname,
                         unsigned int port, DBAdaptor *dnadb) {
  DBAdaptor *dba;

  if ((dba = (DBAdaptor *)calloc(1,sizeof(DBAdaptor))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating DBAdaptor\n");
    return NULL;
  }

  dba->dbc = DBConnection_new(host,user,pass,dbname,port);

  if (dnadb) {
    dba->dnadb = dnadb;
  } else {
    dba->dnadb = dba;
  }

  dba->srIdCache = IDHash_new(IDHASH_MEDIUM);
  dba->srNameCache = StringHash_new(STRINGHASH_MEDIUM);
  return dba;
}

void DBAdaptor_addToSrCaches(DBAdaptor *dba, IDType regionId, char *regionName, IDType csId, long regionLength) {
  char key[1024];
  SeqRegionCacheEntry *cacheData;

  // Do a quick sanity check
  if (IDHash_contains(dba->srIdCache, regionId)) {
    //fprintf(stderr,"Hmm - seq region already in id cache - odd\n");
    //ProcUtil_showBacktrace(EnsC_progName);
    return;
  }

  cacheData = SeqRegionCacheEntry_new(regionId, regionName, csId, regionLength);

  sprintf(key,"%s:"IDFMTSTR, regionName, csId);

  //if (StringHash_contains(dba->srNameCache, key)) {
    //fprintf(stderr,"Hmm - seq region already in name cache - odd\n");
  //} else {
    StringHash_add(dba->srNameCache, key, cacheData);
  //}

  //if (IDHash_contains(dba->srIdCache, regionId)) {
    //fprintf(stderr,"Hmm - seq region already in id cache - odd\n");
  //} else {
    IDHash_add(dba->srIdCache, regionId, cacheData);
  //}

  return;
}
char *DBAdaptor_setAssemblyType(DBAdaptor *dba, char *type) {
  StrUtil_copyString(&(dba->assemblyType),type,0);

  return dba->assemblyType;
}

char *DBAdaptor_getAssemblyType(DBAdaptor *dba) {
  if (!dba->assemblyType) {
    dba->assemblyType = MetaContainer_getDefaultAssembly(DBAdaptor_getMetaContainer(dba->dnadb));
    if (!dba->assemblyType) {
      fprintf(stderr,"ERROR: Didn't get default assembly type\n");
      return NULL;
    }
  }
  return dba->assemblyType;
}

MetaContainer *DBAdaptor_getMetaContainer(DBAdaptor *dba) {
  if (!dba->metaContainer) {
    dba->metaContainer = MetaContainer_new(dba);
  }
  return dba->metaContainer;
}

MetaCoordContainer *DBAdaptor_getMetaCoordContainer(DBAdaptor *dba) {
  if (!dba->metaCoordContainer) {
    dba->metaCoordContainer = MetaCoordContainer_new(dba);
  }
  return dba->metaCoordContainer;
}

GeneAdaptor *DBAdaptor_getGeneAdaptor(DBAdaptor *dba) {
  if (!DBConnection_getAdaptor(dba->dbc,GENE_ADAPTOR)) {
    DBConnection_addAdaptor(dba->dbc,
                            (BaseAdaptor *)GeneAdaptor_new(dba));
  }
  return (GeneAdaptor *)DBConnection_getAdaptor(dba->dbc,GENE_ADAPTOR);
}

AnalysisAdaptor *DBAdaptor_getAnalysisAdaptor(DBAdaptor *dba) {
  if (!DBConnection_getAdaptor(dba->dbc,ANALYSIS_ADAPTOR)) {
    DBConnection_addAdaptor(dba->dbc,
                            (BaseAdaptor *)AnalysisAdaptor_new(dba));
  }
  return (AnalysisAdaptor *)DBConnection_getAdaptor(dba->dbc,ANALYSIS_ADAPTOR);
}

SimpleFeatureAdaptor *DBAdaptor_getSimpleFeatureAdaptor(DBAdaptor *dba) {
  if (!DBConnection_getAdaptor(dba->dbc,SIMPLEFEATURE_ADAPTOR)) {
    DBConnection_addAdaptor(dba->dbc,
                            (BaseAdaptor *)SimpleFeatureAdaptor_new(dba));
  }
  return (SimpleFeatureAdaptor *)DBConnection_getAdaptor(dba->dbc,SIMPLEFEATURE_ADAPTOR);
}

SupportingFeatureAdaptor *DBAdaptor_getSupportingFeatureAdaptor(DBAdaptor *dba) {
  if (!DBConnection_getAdaptor(dba->dbc,SUPPORTINGFEATURE_ADAPTOR)) {
    DBConnection_addAdaptor(dba->dbc,
                            (BaseAdaptor *)SupportingFeatureAdaptor_new(dba));
  }
  return (SupportingFeatureAdaptor *)DBConnection_getAdaptor(dba->dbc,SUPPORTINGFEATURE_ADAPTOR);
}

TranscriptSupportingFeatureAdaptor *DBAdaptor_getTranscriptSupportingFeatureAdaptor(DBAdaptor *dba) {
  if (!DBConnection_getAdaptor(dba->dbc,TRANSCRIPTSUPPORTINGFEATURE_ADAPTOR)) {
    DBConnection_addAdaptor(dba->dbc,
                            (BaseAdaptor *)TranscriptSupportingFeatureAdaptor_new(dba));
  }
  return (TranscriptSupportingFeatureAdaptor *)DBConnection_getAdaptor(dba->dbc,TRANSCRIPTSUPPORTINGFEATURE_ADAPTOR);
}

DNAAlignFeatureAdaptor *DBAdaptor_getDNAAlignFeatureAdaptor(DBAdaptor *dba) {
  if (!DBConnection_getAdaptor(dba->dbc,DNAALIGNFEATURE_ADAPTOR)) {
    DBConnection_addAdaptor(dba->dbc,
                            (BaseAdaptor *)DNAAlignFeatureAdaptor_new(dba));
  }
  return (DNAAlignFeatureAdaptor *)DBConnection_getAdaptor(dba->dbc,DNAALIGNFEATURE_ADAPTOR);
}

ProteinAlignFeatureAdaptor *DBAdaptor_getProteinAlignFeatureAdaptor(DBAdaptor *dba) {
  if (!DBConnection_getAdaptor(dba->dbc,PROTEINALIGNFEATURE_ADAPTOR)) {
    DBConnection_addAdaptor(dba->dbc,
                            (BaseAdaptor *)ProteinAlignFeatureAdaptor_new(dba));
  }
  return (ProteinAlignFeatureAdaptor *)DBConnection_getAdaptor(dba->dbc,PROTEINALIGNFEATURE_ADAPTOR);
}

PredictionTranscriptAdaptor *DBAdaptor_getPredictionTranscriptAdaptor(DBAdaptor *dba) {
  if (!DBConnection_getAdaptor(dba->dbc,PREDICTIONTRANSCRIPT_ADAPTOR)) {
    DBConnection_addAdaptor(dba->dbc,
                            (BaseAdaptor *)PredictionTranscriptAdaptor_new(dba));
  }
  return (PredictionTranscriptAdaptor *)DBConnection_getAdaptor(dba->dbc,PREDICTIONTRANSCRIPT_ADAPTOR);
}

DBEntryAdaptor *DBAdaptor_getDBEntryAdaptor(DBAdaptor *dba) {
  if (!DBConnection_getAdaptor(dba->dbc,DBENTRY_ADAPTOR)) {
    DBConnection_addAdaptor(dba->dbc,
                            (BaseAdaptor *)DBEntryAdaptor_new(dba));
  }
  return (DBEntryAdaptor *)DBConnection_getAdaptor(dba->dbc,DBENTRY_ADAPTOR);
}

RepeatFeatureAdaptor *DBAdaptor_getRepeatFeatureAdaptor(DBAdaptor *dba) {
  if (!DBConnection_getAdaptor(dba->dbc,REPEATFEATURE_ADAPTOR)) {
    DBConnection_addAdaptor(dba->dbc,
                            (BaseAdaptor *)RepeatFeatureAdaptor_new(dba));
  }
  return (RepeatFeatureAdaptor *)DBConnection_getAdaptor(dba->dbc,REPEATFEATURE_ADAPTOR);
}

RepeatConsensusAdaptor *DBAdaptor_getRepeatConsensusAdaptor(DBAdaptor *dba) {
  if (!DBConnection_getAdaptor(dba->dbc,REPEATCONSENSUS_ADAPTOR)) {
    DBConnection_addAdaptor(dba->dbc,
                            (BaseAdaptor *)RepeatConsensusAdaptor_new(dba));
  }
  return (RepeatConsensusAdaptor *)DBConnection_getAdaptor(dba->dbc,REPEATCONSENSUS_ADAPTOR);
}

TranscriptAdaptor *DBAdaptor_getTranscriptAdaptor(DBAdaptor *dba) {
  if (!DBConnection_getAdaptor(dba->dbc,TRANSCRIPT_ADAPTOR)) {
    DBConnection_addAdaptor(dba->dbc,
                            (BaseAdaptor *)TranscriptAdaptor_new(dba));
  }
  return (TranscriptAdaptor *)DBConnection_getAdaptor(dba->dbc,TRANSCRIPT_ADAPTOR);
}

ExonAdaptor *DBAdaptor_getExonAdaptor(DBAdaptor *dba) {
  if (!DBConnection_getAdaptor(dba->dbc,EXON_ADAPTOR)) {
    DBConnection_addAdaptor(dba->dbc,
                            (BaseAdaptor *)ExonAdaptor_new(dba));
  }
  return (ExonAdaptor *)DBConnection_getAdaptor(dba->dbc,EXON_ADAPTOR);
}

RawContigAdaptor *DBAdaptor_getRawContigAdaptor(DBAdaptor *dba) {
  if (!DBConnection_getAdaptor(dba->dbc,RAWCONTIG_ADAPTOR)) {
    DBConnection_addAdaptor(dba->dbc,
                            (BaseAdaptor *)RawContigAdaptor_new(dba));
  }
  return (RawContigAdaptor *)DBConnection_getAdaptor(dba->dbc,RAWCONTIG_ADAPTOR);
}

SliceAdaptor *DBAdaptor_getSliceAdaptor(DBAdaptor *dba) {
  if (!DBConnection_getAdaptor(dba->dbc,SLICE_ADAPTOR)) {
    DBConnection_addAdaptor(dba->dbc,
                            (BaseAdaptor *)SliceAdaptor_new(dba));
  }
  return (SliceAdaptor *)DBConnection_getAdaptor(dba->dbc,SLICE_ADAPTOR);
}

SequenceAdaptor *DBAdaptor_getSequenceAdaptor(DBAdaptor *dba) {
  if (!DBConnection_getAdaptor(dba->dbc,SEQUENCE_ADAPTOR)) {
    DBConnection_addAdaptor(dba->dbc,
                            (BaseAdaptor *)SequenceAdaptor_new(dba));
  }
  return (SequenceAdaptor *)DBConnection_getAdaptor(dba->dbc,SEQUENCE_ADAPTOR);
}

TranslationAdaptor *DBAdaptor_getTranslationAdaptor(DBAdaptor *dba) {
  if (!DBConnection_getAdaptor(dba->dbc,TRANSLATION_ADAPTOR)) {
    DBConnection_addAdaptor(dba->dbc,
                            (BaseAdaptor *)TranslationAdaptor_new(dba));
  }
  return (TranslationAdaptor *)DBConnection_getAdaptor(dba->dbc,TRANSLATION_ADAPTOR);
}

ChromosomeAdaptor *DBAdaptor_getChromosomeAdaptor(DBAdaptor *dba) {
  if (!DBConnection_getAdaptor(dba->dbc,CHROMOSOME_ADAPTOR)) {
    DBConnection_addAdaptor(dba->dbc,
                            (BaseAdaptor *)ChromosomeAdaptor_new(dba));
  }
  return (ChromosomeAdaptor *)DBConnection_getAdaptor(dba->dbc,CHROMOSOME_ADAPTOR);
}

CloneAdaptor *DBAdaptor_getCloneAdaptor(DBAdaptor *dba) {
  if (!DBConnection_getAdaptor(dba->dbc,CLONE_ADAPTOR)) {
    DBConnection_addAdaptor(dba->dbc,
                            (BaseAdaptor *)CloneAdaptor_new(dba));
  }
  return (CloneAdaptor *)DBConnection_getAdaptor(dba->dbc,CLONE_ADAPTOR);
}

AssemblyMapperAdaptor *DBAdaptor_getAssemblyMapperAdaptor(DBAdaptor *dba) {
  if (!DBConnection_getAdaptor(dba->dbc,ASSEMBLYMAPPER_ADAPTOR)) {
    DBConnection_addAdaptor(dba->dbc,
                            (BaseAdaptor *)AssemblyMapperAdaptor_new(dba));
  }
  return (AssemblyMapperAdaptor *)DBConnection_getAdaptor(dba->dbc,ASSEMBLYMAPPER_ADAPTOR);
}

CoordSystemAdaptor *DBAdaptor_getCoordSystemAdaptor(DBAdaptor *dba) {
  if (!DBConnection_getAdaptor(dba->dbc,COORDSYSTEM_ADAPTOR)) {
    DBConnection_addAdaptor(dba->dbc,
                            (BaseAdaptor *)CoordSystemAdaptor_new(dba));
  }
  return (CoordSystemAdaptor *)DBConnection_getAdaptor(dba->dbc,COORDSYSTEM_ADAPTOR);
}
