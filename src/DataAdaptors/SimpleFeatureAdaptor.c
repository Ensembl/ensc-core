#include "SimpleFeatureAdaptor.h"

#include "RawContigAdaptor.h"
#include "AnalysisAdaptor.h"

NameTableType SimpleFeatureAdaptor_tableNames = {{"simple_feature","sf"},{NULL,NULL}};

SimpleFeatureAdaptor *SimpleFeatureAdaptor_new(DBAdaptor *dba) {
  SimpleFeatureAdaptor *sfa;

  if ((sfa = (SimpleFeatureAdaptor *)calloc(1,sizeof(SimpleFeatureAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for SimpleFeatureAdaptor\n");
    return NULL;
  }
  BaseFeatureAdaptor_init((BaseFeatureAdaptor *)sfa, dba, SIMPLEFEATURE_ADAPTOR);

  sfa->getTables = SimpleFeatureAdaptor_getTables;
  sfa->getColumns = SimpleFeatureAdaptor_getColumns;
  sfa->store = SimpleFeatureAdaptor_store;
  sfa->objectsFromStatementHandle = SimpleFeatureAdaptor_objectsFromStatementHandle;


  return sfa;
}

int SimpleFeatureAdaptor_store(BaseFeatureAdaptor *bfa, Vector *features) {
  char qStr[512];
  StatementHandle *sth;
  int i;
  
  if (!Vector_getNumElement(features)) {
    fprintf(stderr, "Warning: SimpleFeatureAdaptor_store called with no features\n");
    return 0;
  }
  
  sprintf(qStr,"INSERT INTO simple_feature (contig_id, contig_start,"
                                           " contig_end, contig_strand,"
                                           "display_label, analysis_id,"
                                           "score) VALUES (%"
                                           IDFMTSTR ",%%d,%%d,%%d,\'%%s\',%"
                                           IDFMTSTR ",%%f)");
  sth = bfa->prepare((BaseAdaptor *)bfa, qStr,strlen(qStr));
  printf("%s\n",qStr);

  for (i=0; i<Vector_getNumElement(features); i++) {
    SimpleFeature *sf = Vector_getElementAt(features, i);
    Analysis *analysis = SimpleFeature_getAnalysis(sf);
    RawContig *contig;

/* NIY
    if( !ref $sf || !$sf->isa("Bio::EnsEMBL::SimpleFeature") ) {
      $self->throw("Simple feature must be an Ensembl SimpleFeature, " .
		   "not a [$sf]");
    }
*/
    
    if (!analysis) {
      fprintf(stderr,"Cannot store sequence features without analysis");
      exit(1);
    }
    if( !Analysis_getDbID(analysis)) {
      fprintf(stderr,"I think we should always have an analysis object "
		     "which has originated from the database. No dbID, "
		     "not putting in!");
      exit(1);
    }
    
    if (!SimpleFeature_getContig(sf)->objectType != CLASS_RAWCONTIG) {
      fprintf(stderr,"Error: contig isn't raw contig when trying to store\n");
      exit(1);
    }


    contig = (RawContig *)SimpleFeature_getContig(sf);
/* NIY
    my $contig = $sf->entire_seq();
    unless(defined $contig && $contig->isa("Bio::EnsEMBL::RawContig")) {
      $self->throw("Cannot store feature without a Contig object attached via "
		   . "attach_seq\n");
    }
*/

    sth->execute(sth, (IDType)(RawContig_getDbID(contig)), SimpleFeature_getStart(sf), SimpleFeature_getEnd(sf), 
                 SimpleFeature_getStrand(sf),SimpleFeature_getDisplayLabel(sf),
                 (IDType)(Analysis_getDbID(analysis)), SimpleFeature_getScore(sf));
  } 
  return 1;
}

NameTableType *SimpleFeatureAdaptor_getTables(void) {
  return &SimpleFeatureAdaptor_tableNames;
}

char *SimpleFeatureAdaptor_getColumns(void) {
  return "sf.simple_feature_id," 
	 "sf.contig_id,"
         "sf.contig_start,"
         "sf.contig_end,"
         "sf.contig_strand,"
	 "sf.display_label,"
         "sf.analysis_id,"
         "sf.score";
}

Vector *SimpleFeatureAdaptor_objectsFromStatementHandle(BaseFeatureAdaptor *bfa,
                                                     StatementHandle *sth,
                                                     AssemblyMapper *mapper,
                                                     Slice *slice) {
  AnalysisAdaptor *aa;
  RawContigAdaptor *rca;
  Vector *features;
  ResultRow *row;

  aa = DBAdaptor_getAnalysisAdaptor(bfa->dba);
  rca = DBAdaptor_getRawContigAdaptor(bfa->dba);

  features = Vector_new();
  
  while ((row = sth->fetchRow(sth))) {
    RawContig *contig = RawContigAdaptor_fetchByDbID(rca, row->getLongLongAt(row,1));
    Analysis  *analysis = AnalysisAdaptor_fetchByDbID(aa, row->getLongLongAt(row,6));

    SimpleFeature *sf = SimpleFeature_new();
    SimpleFeature_setStart(sf,row->getIntAt(row,2));
    SimpleFeature_setEnd(sf,row->getIntAt(row,3));
    SimpleFeature_setStrand(sf,row->getIntAt(row,4));
    SimpleFeature_setAnalysis(sf,analysis);
    SimpleFeature_setDisplayLabel(sf,row->getStringAt(row,5));
    SimpleFeature_setContig(sf,contig); 

    if (row->col(row,7)) {
      SimpleFeature_setScore(sf,row->getDoubleAt(row,7));
    }
    
    SimpleFeature_setDbID(sf,row->getLongLongAt(row,0));

    Vector_addElement(features,sf);
  }

  return features;
}
