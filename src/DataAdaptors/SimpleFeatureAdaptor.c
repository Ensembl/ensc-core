#include "SimpleFeatureAdaptor.h"

#include "RawContigAdaptor.h"
#include "AnalysisAdaptor.h"

char *SimpleFeatureAdaptor_tableNames[][2] = {{"simple_feature","sf"}};

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

int SimpleFeatureAdaptor_store(BaseFeatureAdaptor *baf, Set *features) {
/*
  my ($self,@sf) = @_;
  
  if( scalar(@sf) == 0 ) {
    $self->throw("Must call store with list of sequence features");
  }
  
  my $sth = 
    $self->prepare("INSERT INTO simple_feature (contig_id, contig_start,
                                                contig_end, contig_strand,
                                                display_label, analysis_id,
                                                score) 
                    VALUES (?,?,?,?,?,?,?)");

  foreach my $sf ( @sf ) {
    if( !ref $sf || !$sf->isa("Bio::EnsEMBL::SimpleFeature") ) {
      $self->throw("Simple feature must be an Ensembl SimpleFeature, " .
		   "not a [$sf]");
    }
    
    if( !defined $sf->analysis ) {
      $self->throw("Cannot store sequence features without analysis");
    }
    if( !defined $sf->analysis->dbID ) {
      $self->throw("I think we should always have an analysis object " .
		   "which has originated from the database. No dbID, " .
		   "not putting in!");
    }
    
    my $contig = $sf->entire_seq();
    unless(defined $contig && $contig->isa("Bio::EnsEMBL::RawContig")) {
      $self->throw("Cannot store feature without a Contig object attached via "
		   . "attach_seq\n");
    }

    $sth->execute($contig->dbID(), $sf->start, $sf->end, $sf->strand,
		  $sf->display_label, $sf->analysis->dbID, $sf->score);
  } 
*/
}

char ***SimpleFeatureAdaptor_getTables() {
  return SimpleFeatureAdaptor_tableNames;
}

char *SimpleFeatureAdaptor_getColumns() {
  return "sf.simple_feature_id," 
	 "sf.contig_id,"
         "sf.contig_start,"
         "sf.contig_end,"
         "sf.contig_strand,"
	 "sf.display_label,"
         "sf.analysis_id,"
         "sf.score";
}

Set *SimpleFeatureAdaptor_objectsFromStatementHandle(BaseFeatureAdaptor *baf,
                                                     StatementHandle *sth,
                                                     AssemblyMapper *mapper,
                                                     Slice *slice) {
  AnalysisAdaptor *aa;
  RawContigAdaptor *rca;
  Set *features;
  ResultRow *row;

  aa = DBAdaptor_getAnalysisAdaptor(baf->dba);
  rca = DBAdaptor_getRawContigAdaptor(baf->dba);

  features = Set_new();
  
  while(row = sth->fetchRow(sth)) {
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

    Set_addElement(features,sf);
  }

  return features;
}
