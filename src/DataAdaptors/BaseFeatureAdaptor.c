#include "BaseFeatureAdaptor.h"

#include "AnalysisAdaptor.h"
#include "AssemblyMapperAdaptor.h"
#include "StrUtil.h"
#include "Error.h"

int SLICE_FEATURE_CACHE_SIZE = 4;


#define NAME 0
#define SYN  1

void BaseFeatureAdaptor_init(BaseFeatureAdaptor *bfa, DBAdaptor *dba, int adaptorType) {
  BaseAdaptor_init((BaseAdaptor *)bfa,dba,adaptorType);

  bfa->sliceFeatureCache = Cache_new(SLICE_FEATURE_CACHE_SIZE);

  bfa->objectsFromStatementHandle = BaseFeatureAdaptor_objectsFromStatementHandle;
  bfa->getTables                  = BaseFeatureAdaptor_getTables;
  bfa->getColumns                 = BaseFeatureAdaptor_getColumns;
  bfa->finalClause                = BaseFeatureAdaptor_finalClause;
  bfa->leftJoin                   = BaseFeatureAdaptor_leftJoin;
  bfa->defaultWhereClause         = BaseFeatureAdaptor_defaultWhereClause;
  bfa->store                      = BaseFeatureAdaptor_store;

  return;
}

Vector *BaseFeatureAdaptor_genericFetch(BaseFeatureAdaptor *bfa, char *constraint,
                                     char *logicName, AssemblyMapper *mapper, Slice *slice) {
/* HACK HACK HACK */
  char qStr[65500]; 
  NameTableType *tables = bfa->getTables();
  char *columns = bfa->getColumns();
  StatementHandle *sth;
  Vector *features;
/* HACK HACK HACK */
  char allConstraints[65500];
  char tableNamesStr[512]; 
  char leftJoinStr[512]; 
  char **lj;
/* HACK HACK HACK */
  char tmpStr[65500];
  int i;
  
  allConstraints[0] = '\0';
  if (constraint[0]) strcpy(allConstraints,constraint);
  
  if (logicName[0]) {
    AnalysisAdaptor *aa = DBAdaptor_getAnalysisAdaptor(bfa->dba);
    Analysis *analysis;
    char *syn;
    IDType analysisId;

    //determine the analysis id via the logic_name
    analysis = AnalysisAdaptor_fetchByLogicName(aa, logicName);

    if (!analysis || !Analysis_getDbID(analysis) ) {
      fprintf(stderr,"No analysis for logic name %s exists\n",logicName);
      return emptyVector;
    }
    
    analysisId = Analysis_getDbID(analysis);

    // get the synonym for the primary table
    syn = (*tables)[0][SYN];

    if (constraint[0]) {
      sprintf(allConstraints,"%s  AND %s.analysis_id = " IDFMTSTR, constraint, syn, analysisId);
    } else {
      sprintf(allConstraints," %s.analysis_id = " IDFMTSTR, syn, analysisId);
    }
  } 

  //
  // Construct a left join statement if one was defined, and remove the
  // left-joined table from the table list
  // 

  leftJoinStr[0]   = '\0';
  tableNamesStr[0] = '\0';

  lj = bfa->leftJoin();

  i = 0;
  while ((*tables)[i][0]) {
    char **t = (*tables)[i];
    if (lj!=NULL && !strcmp(lj[0],t[0])) {
      sprintf(leftJoinStr,"LEFT JOIN %s %s %s",lj[0], t[SYN], lj[1]);
    } else {
      if (tableNamesStr[0]) {
        sprintf(tmpStr,", %s %s",t[NAME],t[SYN]);
      } else {
        sprintf(tmpStr,"%s %s",t[NAME],t[SYN]);
      }
      strcat(tableNamesStr,tmpStr);
    }
    i++;
  }
      
  sprintf(qStr,"SELECT %s FROM %s %s", columns, tableNamesStr, leftJoinStr);

  //append a where clause if it was defined
  if (allConstraints[0]) { 
    // printf("qStr = %s allConstraints = %s\n",qStr,allConstraints);
    sprintf(tmpStr," where %s", allConstraints);
    strcat(qStr,tmpStr);
    if ((bfa->defaultWhereClause())[0]) {
      sprintf(tmpStr," and %s", bfa->defaultWhereClause());
      strcat(qStr,tmpStr);
    }
  } else if ((bfa->defaultWhereClause())[0]) {
    sprintf(tmpStr," where %s", bfa->defaultWhereClause());
    strcat(qStr,tmpStr);
  }

  //append additional clauses which may have been defined
  if ((bfa->finalClause())[0]) strcat(qStr, bfa->finalClause());

  sth = bfa->prepare((BaseAdaptor *)bfa,qStr,strlen(qStr));
  sth->execute(sth);  

  features = bfa->objectsFromStatementHandle(bfa, sth, mapper, slice);
  sth->finish(sth);

  return features;
}

SeqFeature *BaseFeatureAdaptor_fetchByDbID(BaseFeatureAdaptor *bfa, IDType dbID) {
  Vector *features;
  SeqFeature *sf;
  char constraintStr[256];
  NameTableType *tables = bfa->getTables();

  //construct a constraint like 't1.table1_id = 1'
  sprintf(constraintStr,"%s.%s_id = " IDFMTSTR, (*tables)[0][SYN], (*tables)[0][NAME], dbID);

  //return first element of _generic_fetch list
  features = BaseFeatureAdaptor_genericFetch(bfa, constraintStr, "", NULL, NULL);
  sf = Vector_getElementAt(features, 0);

  if (sf) Object_incRefCount(sf);

// NIY free func
  if (Vector_getNumElement(features) > 1) {
    fprintf(stderr,
            "Error: Memory leak because got more than one feature with id " 
            IDFMTSTR, dbID);
  }
  Vector_free(features);

  return (SeqFeature *)sf;
}

Vector *BaseFeatureAdaptor_fetchAllByRawContigConstraint(BaseFeatureAdaptor *bfa, RawContig *contig,
                                                      char *constraint, char *logicName)  {
  IDType cid;
  char allConstraints[256];
  NameTableType *tables = bfa->getTables();

  if (contig == NULL) {
    fprintf(stderr,"ERROR: fetch_by_Contig_constraint must have an contig\n");
    exit(1);
  }

  cid = RawContig_getDbID(contig);

  if (constraint[0]) {
    sprintf(allConstraints,"%s AND %s.contig_id = " IDFMTSTR, constraint, (*tables)[0][SYN], cid);
  } else {
    sprintf(allConstraints,"%s.contig_id = " IDFMTSTR, (*tables)[0][SYN], cid);
  }

  return BaseFeatureAdaptor_genericFetch(bfa, allConstraints, logicName, NULL, NULL);
}

Vector *BaseFeatureAdaptor_fetchAllByRawContig(BaseFeatureAdaptor *bfa, RawContig *contig,
                                            char *logicName) {
  return BaseFeatureAdaptor_fetchAllByRawContigConstraint(bfa,contig,"",logicName);
}

Vector *BaseFeatureAdaptor_fetchAllByRawContigAndScore(BaseFeatureAdaptor *bfa, RawContig *contig,
                                                    double *scoreP, char *logicName) {
  char constraintStr[256];
  NameTableType *tables = bfa->getTables();

  constraintStr[0] = '\0';
// Perl does a defined check on score
  if (scoreP) {
    sprintf(constraintStr,"%s.score > %f",(*tables)[0][SYN], *scoreP);
  }
    
  return BaseFeatureAdaptor_fetchAllByRawContigConstraint(bfa, contig, constraintStr, logicName);
}

Vector *BaseFeatureAdaptor_fetchAllBySlice(BaseFeatureAdaptor *bfa, Slice *slice,
                                        char *logicName) {
  return BaseFeatureAdaptor_fetchAllBySliceConstraint(bfa, slice, "", logicName);
}

Vector *BaseFeatureAdaptor_fetchAllBySliceAndScore(BaseFeatureAdaptor *bfa, Slice *slice,
                                                double *scoreP, char *logicName) {
  char constraintStr[256];
  NameTableType *tables = bfa->getTables();

  constraintStr[0] = '\0';
// Perl does a defined check on score
  if (scoreP) {
    sprintf(constraintStr,"%s.score > %f",(*tables)[0][SYN], *scoreP);
  }

  return BaseFeatureAdaptor_fetchAllBySliceConstraint(bfa, slice, constraintStr, logicName);
}  

Vector *BaseFeatureAdaptor_fetchAllBySliceConstraint(BaseFeatureAdaptor *bfa, Slice *slice,
                                                  char *constraint, char *logicName) {

  char cacheKey[EXTREMELEN];
  void *val;
  Vector *features;
  Vector *out;
  char *allConstraints;
  int sliceChrId;
  int sliceEnd;
  int sliceStart;
  int sliceStrand;
  NameTableType *tables = bfa->getTables();
  int nContigId;
  IDType *contigIds;
  char tmpStr[512];
  AssemblyMapper *assMapper;
  AssemblyMapperAdaptor *ama;
  int i;

  // check the cache and return if we have already done this query

  sprintf(cacheKey,"%s%s%s",Slice_getName(slice), constraint, logicName);
  //printf("cacheKey = %s\n",cacheKey);
  StrUtil_strupr(cacheKey);

  if ((val = Cache_findElem(bfa->sliceFeatureCache, cacheKey)) != NULL) {
    return (Vector *)val;
  }
    
  sliceChrId = Slice_getChrId(slice);
  sliceStart = Slice_getChrStart(slice);
  sliceEnd   = Slice_getChrEnd(slice);
  sliceStrand= Slice_getStrand(slice);

/*
  ama = DBAdaptor_getAssemblyMapperAdaptor(bfa->dba);
  assMapper = AssemblyMapperAdaptor_fetchByType(ama,Slice_getAssemblyType(slice));

  nContigId = AssemblyMapper_listContigIds(assMapper,
                                           sliceChrId,
                                           sliceStart,
                                           sliceEnd,
                                           &contigIds);


  if (!nContigId) {
    return emptyVector;
  }

  //construct the SQL constraint for the contig ids 
  if (constraint[0]) {
    sprintf(tmpStr,"%s AND %s.contig_id IN (", constraint, (*tables)[0][SYN]);
  } else {
    sprintf(tmpStr,"%s.contig_id IN (", (*tables)[0][SYN]);
  }

  allConstraints = StrUtil_copyString(&allConstraints,tmpStr,0);

  if (!allConstraints) {
    Error_trace("fetch_all_by_Slice",NULL);
    return emptyVector;
  }

  for (i=0; i<nContigId; i++) {
    char numStr[256];
    if (i!=(nContigId-1)) {
      sprintf(numStr,IDFMTSTR ",",contigIds[i]);
    } else {
      sprintf(numStr,IDFMTSTR,contigIds[i]);
    }
    allConstraints = StrUtil_appendString(allConstraints, numStr);
  }


  allConstraints = StrUtil_appendString(allConstraints,")");
*/

  if (constraint[0]) {
    sprintf(tmpStr,"%s AND %s.seq_region_id = sr.seq_region_id and sr.name = '%s' and %s.seq_region_end >= %d and %s.seq_region_start <= %d", constraint, (*tables)[0][SYN],Slice_getChrName(slice),(*tables)[0][SYN],Slice_getChrStart(slice),(*tables)[0][SYN],Slice_getChrEnd(slice));
  } else {
    sprintf(tmpStr,"%s.seq_region_id = sr.seq_region_id and sr.name = '%s' and %s.seq_region_end >= %d and %s.seq_region_start <= %d", (*tables)[0][SYN],Slice_getChrName(slice),(*tables)[0][SYN],Slice_getChrStart(slice),(*tables)[0][SYN],Slice_getChrEnd(slice));
  }
  allConstraints = StrUtil_copyString(&allConstraints,tmpStr,0);

  // for speed the remapping to slice may be done at the time of object creation
  //printf("allConstraints = %s\n",allConstraints);
  features = 
    BaseFeatureAdaptor_genericFetch(bfa, allConstraints, logicName, assMapper, slice); 
  
  if (Vector_getNumElement(features)) {
// Can't easily do this in C     && (!$features->[0]->can('contig') || 
// Its for the PredictionTranscriptAdaptor so HACK HACK HACK
    SeqFeature *sf = (SeqFeature *)Vector_getElementAt(features,0);
    if (bfa->adaptorType == PREDICTIONTRANSCRIPT_ADAPTOR || SeqFeature_getContig(sf) == (BaseContig *)slice) {
      // features have been converted to slice coords already, cache and return
      Cache_addElement(bfa->sliceFeatureCache, cacheKey, features, NULL);
      return features;
    }
  } 

  //remapping has not been done, we have to do our own conversion from
  //raw contig coords to slice coords
  out = Vector_new();

    
  for (i=0;i<Vector_getNumElement(features); i++) {
    //since feats were obtained in contig coords, attached seq is a contig
    SeqFeature *f = Vector_getElementAt(features, i);
    IDType contigId = RawContig_getDbID(SeqFeature_getContig(f));
    MapperCoordinate fRange;
  
    int mapSucceeded = AssemblyMapper_fastToAssembly(assMapper, contigId, 
                                               SeqFeature_getStart(f), 
                                               SeqFeature_getEnd(f), 
                                               SeqFeature_getStrand(f), 
  				               &fRange);
  
    // undefined start means gap
    if (!mapSucceeded) continue;
  
    // maps to region outside desired area 
    if (fRange.start > sliceEnd || fRange.end < sliceStart) continue;
      
    // shift the feature start, end and strand in one call
    // In C I can't be arsed to write this call - it should be quick enough
    if(sliceStrand == -1) {
      SeqFeature_setStart (f, sliceEnd - fRange.end + 1);
      SeqFeature_setEnd   (f, sliceEnd - fRange.start + 1);
      SeqFeature_setStrand(f, fRange.strand * -1 );
    } else {
      SeqFeature_setStart (f, fRange.start - sliceStart + 1);
      SeqFeature_setEnd   (f, fRange.end - sliceStart + 1);
      SeqFeature_setStrand(f, fRange.strand);
    }
      
    SeqFeature_setContig(f,slice);
      
    Vector_addElement(out,f);
  }
    
  //update the cache
  Cache_addElement(bfa->sliceFeatureCache, cacheKey, out, NULL);
  return out;
}

int BaseFeatureAdaptor_store(BaseFeatureAdaptor *bfa, Vector *features) {
  fprintf(stderr,"ERROR: Abstract method store not defined by implementing subclass\n");
  exit(1);

// Please the compiler
  return 0;
}

int BaseFeatureAdaptor_remove(BaseFeatureAdaptor *bfa, SeqFeature *feature) {
  char qStr[256];
  NameTableType *tables = bfa->getTables();
  char *tableName = (*tables)[0][0];
  StatementHandle *sth;
  
  if (!SeqFeature_getDbID(feature)) {
    fprintf(stderr, "BaseFeatureAdaptor_remove - dbID not defined - "
                    "feature could not be removed\n");
  }


  sprintf(qStr,"DELETE FROM %s WHERE %s_id = " IDFMTSTR,tableName,tableName,SeqFeature_getDbID(feature));

  sth = bfa->prepare((BaseAdaptor *)bfa, qStr,strlen(qStr));
  sth->execute(sth);
  sth->finish(sth);

  //unset the feature dbID
  SeqFeature_setDbID(feature, 0);
  
  return 0;
}

int BaseFeatureAdaptor_removeByRawContig(BaseFeatureAdaptor *bfa, RawContig *contig) {
  char qStr[256];
  NameTableType *tables = bfa->getTables();
  char *tableName; 
  StatementHandle *sth;
  tableName = (*tables)[0][0];

  if (contig == NULL) {
    fprintf(stderr,"BaseFeatureAdaptor_removeByRawContig - no contig supplied: "
		   "Deletion of features failed.\n");
    return 0;
  }


  sprintf(qStr, "DELETE FROM %s WHERE contig_id = " IDFMTSTR, tableName, RawContig_getDbID(contig));

  sth = bfa->prepare((BaseAdaptor *)bfa,qStr,strlen(qStr));
  sth->execute(sth);
  sth->finish(sth);

  return 1;
}

NameTableType *BaseFeatureAdaptor_getTables(void) {
  fprintf(stderr,"ERROR: Abstract method getTables not defined by implementing subclass\n");
  exit(1);

// Please the compiler
  return NULL;
}

char *BaseFeatureAdaptor_getColumns(void) {
  fprintf(stderr,"ERROR: Abstract method getColumns not defined by implementing subclass\n");
  exit(1);

// Please the compiler
  return NULL;
}

/* Actually added to default where */
char *BaseFeatureAdaptor_defaultWhereClause(void) {
  return "";
}

/* Overridden in Markers and Qtls (not implemented in C) */
char **BaseFeatureAdaptor_leftJoin(void) {
  return NULL;
}

/* Overridden in PredictionTranscriptAdaptor */
char *BaseFeatureAdaptor_finalClause(void) {
  return "";
}

Vector *BaseFeatureAdaptor_objectsFromStatementHandle(BaseFeatureAdaptor *bfa, StatementHandle *sth, 
                                                   AssemblyMapper *mapper, Slice *slice) {
  fprintf(stderr,"ERROR: Abstract method objectsFromStatementHandle not defined by implementing subclass\n");
  exit(1);

// Please the compiler
  return NULL;
} 
