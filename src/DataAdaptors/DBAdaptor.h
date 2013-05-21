#ifndef __DBADAPTOR_H__
#define __DBADAPTOR_H__

#include "BaseDBAdaptor.h"
#include "DBConnection.h"
#include "AdaptorTypes.h"
#include "EnsC.h"
#include "IDHash.h"
#include "StringHash.h"

struct DBAdaptorStruct {
  BASEDBADAPTOR_DATA
  DBAdaptor     *dnadb;
  char          *assemblyType;
  IDHash        *srIdCache;
  StringHash    *srNameCache;
  int            noCache;
  int            speciesId;
};

DBAdaptor *DBAdaptor_new(char *host, char *user, char *pass, char *dbname,
                         unsigned int port, DBAdaptor *dnadb);

char *DBAdaptor_setAssemblyType(DBAdaptor *dba, char *type);
char *DBAdaptor_getAssemblyType(DBAdaptor *dba);

AnalysisAdaptor             *DBAdaptor_getAnalysisAdaptor(DBAdaptor *dba);
AssemblyMapperAdaptor       *DBAdaptor_getAssemblyMapperAdaptor(DBAdaptor *dba);
AttributeAdaptor            *DBAdaptor_getAttributeAdaptor(DBAdaptor *dba);
ChromosomeAdaptor           *DBAdaptor_getChromosomeAdaptor(DBAdaptor *dba);
CoordSystemAdaptor          *DBAdaptor_getCoordSystemAdaptor(DBAdaptor *dba);
DBEntryAdaptor              *DBAdaptor_getDBEntryAdaptor(DBAdaptor *dba);
DNAAlignFeatureAdaptor      *DBAdaptor_getDNAAlignFeatureAdaptor(DBAdaptor *dba);
ExonAdaptor                 *DBAdaptor_getExonAdaptor(DBAdaptor *dba);
GeneAdaptor                 *DBAdaptor_getGeneAdaptor(DBAdaptor *dba);
IntronSupportingEvidenceAdaptor *DBAdaptor_getIntronSupportingEvidenceAdaptor(DBAdaptor *dba);
MetaContainer               *DBAdaptor_getMetaContainer(DBAdaptor *dba);
MetaCoordContainer          *DBAdaptor_getMetaCoordContainer(DBAdaptor *dba);
PredictionTranscriptAdaptor *DBAdaptor_getPredictionTranscriptAdaptor(DBAdaptor *dba);
ProteinAlignFeatureAdaptor  *DBAdaptor_getProteinAlignFeatureAdaptor(DBAdaptor *dba);
RawContigAdaptor            *DBAdaptor_getRawContigAdaptor(DBAdaptor *dba);
RepeatConsensusAdaptor      *DBAdaptor_getRepeatConsensusAdaptor(DBAdaptor *dba);
RepeatFeatureAdaptor        *DBAdaptor_getRepeatFeatureAdaptor(DBAdaptor *dba);
SliceAdaptor                *DBAdaptor_getSliceAdaptor(DBAdaptor *dba); 
SimpleFeatureAdaptor        *DBAdaptor_getSimpleFeatureAdaptor(DBAdaptor *dba); 
SupportingFeatureAdaptor    *DBAdaptor_getSupportingFeatureAdaptor(DBAdaptor *dba);
SequenceAdaptor             *DBAdaptor_getSequenceAdaptor(DBAdaptor *dba); 
TranscriptAdaptor           *DBAdaptor_getTranscriptAdaptor(DBAdaptor *dba);
TranscriptSupportingFeatureAdaptor    *DBAdaptor_getTranscriptSupportingFeatureAdaptor(DBAdaptor *dba);
TranslationAdaptor          *DBAdaptor_getTranslationAdaptor(DBAdaptor *dba);


#define DBAdaptor_getDNADBAdaptor(dba) (dba)->dnadb
#define DBAdaptor_setDNADBAdaptor(dba, ddb) (dba)->dnadb = ddb

#define DBAdaptor_getSeqRegionIdCache(dba) (dba)->srIdCache
#define DBAdaptor_getSeqRegionNameCache(dba) (dba)->srNameCache

#define DBAdaptor_setNoCache(dba, val) (dba)->noCache = (val)
#define DBAdaptor_noCache(dba) (dba)->noCache

#define DBAdaptor_setSpeciesId(dba, val) (dba)->speciesId = (val)
#define DBAdaptor_getSpeciesId(dba) (dba)->speciesId

#define DBAdaptor_prepare(dba,qStr,qLen) BaseDBAdaptor_prepare((dba),(qStr),(qLen))


#endif
