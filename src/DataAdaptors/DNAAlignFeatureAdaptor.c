#include "DNAAlignFeatureAdaptor.h"

#include "RawContigAdaptor.h"
#include "AnalysisAdaptor.h"
#include "DNAAlignFeature.h"

NameTableType DNAAlignFeatureAdaptor_tableNames = {{"dna_align_feature","daf"},{NULL, NULL}};

DNAAlignFeatureAdaptor *DNAAlignFeatureAdaptor_new(DBAdaptor *dba) {
  DNAAlignFeatureAdaptor *dafa;

  if ((dafa = (DNAAlignFeatureAdaptor *)calloc(1,sizeof(DNAAlignFeatureAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for DNAAlignFeatureAdaptor\n");
    return NULL;
  }
  BaseFeatureAdaptor_init((BaseFeatureAdaptor *)dafa, dba, DNAALIGNFEATURE_ADAPTOR);

  dafa->getTables = DNAAlignFeatureAdaptor_getTables;
  dafa->getColumns = DNAAlignFeatureAdaptor_getColumns;
  dafa->store = DNAAlignFeatureAdaptor_store;
  dafa->objectsFromStatementHandle = DNAAlignFeatureAdaptor_objectsFromStatementHandle;

  return dafa;
}

int DNAAlignFeatureAdaptor_store(BaseFeatureAdaptor *bfa, Vector *features) {
  char qStr[512];
  StatementHandle *sth;
  int i;
  
  if (!Vector_getNumElement(features)) {
    fprintf(stderr, "Warning: ProteinAlignFeatureAdaptor_store called with no features\n");
    return 0;
  }
  
  sprintf(qStr,"INSERT INTO dna_align_feature(contig_id, contig_start, contig_end,"
                       "contig_strand, hit_start, hit_end, hit_strand, hit_name," 
                       "cigar_line, analysis_id,score, evalue, perc_ident) "
                       "VALUES (%" IDFMTSTR ",%%d,%%d,%%d,%%d,%%d,%%d,'%%s','%%s',%" 
                                IDFMTSTR ",%%f,%%f,%%f)");

  sth = bfa->prepare((BaseAdaptor *)bfa, qStr,strlen(qStr));
  printf("%s\n",qStr);

  for (i=0; i<Vector_getNumElement(features); i++) {
    DNAAlignFeature *sf = Vector_getElementAt(features, i);
    Analysis *analysis = DNAAlignFeature_getAnalysis(sf);
    AnalysisAdaptor *aa = DBAdaptor_getAnalysisAdaptor(bfa->dba);
    RawContig *contig;

/* NIY
    if( !ref $sf || !$sf->isa("Bio::EnsEMBL::DNADnaAlignFeature") ) {
      $self->throw("Feature must be a Bio::EnsEMBL::DNAAlignFeature, " .
                   "not a [$sf]");
    }
*/
     
    if (!analysis) {
      fprintf(stderr,"Cannot store sequence features without analysis");
      exit(1);
    }

    // will only store if object is not already stored in this database
    AnalysisAdaptor_store(aa,analysis);

    if (DNAAlignFeature_getContig(sf)->objectType != CLASS_RAWCONTIG) {
      fprintf(stderr,"Error: contig isn't raw contig when trying to store\n");
      exit(1);
    }

    contig = (RawContig *)DNAAlignFeature_getContig(sf);

/* NIY
     unless(defined $contig && $contig->isa("Bio::EnsEMBL::RawContig")) { 
       $self->throw("Cannot store feature without Contig attached via " .
                    "attach_seq\n");
     }   
*/
     
    sth->execute(sth, (IDType)RawContig_getDbID(contig), 
                      DNAAlignFeature_getStart(sf), 
                      DNAAlignFeature_getEnd(sf), 
                      DNAAlignFeature_getStrand(sf), 
                      DNAAlignFeature_getHitStart(sf), 
                      DNAAlignFeature_getHitEnd(sf),
                      DNAAlignFeature_getHitStrand(sf),
                      DNAAlignFeature_getHitSeqName(sf), 
                      DNAAlignFeature_getCigarString(sf),
                      (IDType)Analysis_getDbID(analysis),
                      DNAAlignFeature_getScore(sf), 
                      DNAAlignFeature_getpValue(sf), 
                      DNAAlignFeature_getPercId(sf));
     
    DNAAlignFeature_setDbID(sf,sth->getInsertId(sth));
  }
  sth->finish(sth);
  return 1;
}

NameTableType *DNAAlignFeatureAdaptor_getTables(void) {
  return &DNAAlignFeatureAdaptor_tableNames;
}

char *DNAAlignFeatureAdaptor_getColumns(void) {
  return "daf.dna_align_feature_id,"
         "daf.contig_id,"
         "daf.analysis_id,"
         "daf.contig_start,"
         "daf.contig_end,"
         "daf.contig_strand,"
         "daf.hit_start,"
         "daf.hit_end,"
         "daf.hit_name,"
         "daf.hit_strand,"
         "daf.cigar_line,"
         "daf.evalue,"
         "daf.perc_ident,"
         "daf.score";
}

Vector *DNAAlignFeatureAdaptor_objectsFromStatementHandle(BaseFeatureAdaptor *bfa,
                                                       StatementHandle *sth,
                                                       AssemblyMapper *assMapper,
                                                       Slice *slice) {

  AnalysisAdaptor *aa;
  RawContigAdaptor *rca;
  Vector *features;
  ResultRow *row;

  aa = DBAdaptor_getAnalysisAdaptor(bfa->dba);
  rca = DBAdaptor_getRawContigAdaptor(bfa->dba);

  features = Vector_new();

  if (slice) {
    int featStart, featEnd, featStrand;
    int sliceChrId;
    int sliceEnd;
    int sliceStart;
    int sliceStrand;


    sliceChrId = Slice_getChrId(slice);
    sliceStart = Slice_getChrStart(slice);
    sliceEnd   = Slice_getChrEnd(slice);
    sliceStrand= Slice_getStrand(slice);

    // Does this really need to be set ??? my $slice_name   = $slice->name();


    while ((row = sth->fetchRow(sth))) {
      DNAAlignFeature *daf;
      int contigId    = row->getLongLongAt(row,1);
      int contigStart = row->getIntAt(row,3);
      int contigEnd   = row->getIntAt(row,4);
      int contigStrand= row->getIntAt(row,5);


// Perl has a cache for analysis types but the analysis adaptor should have one
      Analysis  *analysis = AnalysisAdaptor_fetchByDbID(aa, row->getLongLongAt(row,2));
      MapperCoordinate fRange;


      //convert contig coordinates to assembly coordinates
      int mapSucceeded = AssemblyMapper_fastToAssembly(assMapper, contigId, 
                                               contigStart, 
                                               contigEnd, 
                                               contigStrand, 
                                               &fRange);

      
      // undefined start means gap
      if (!mapSucceeded) continue;
  
      // maps to region outside desired area 
      if (fRange.start > sliceEnd || fRange.end < sliceStart) continue;


      // convert assembly coordinates to slice coordinates
      if(sliceStrand == -1) {
        featStart  = sliceEnd - fRange.end + 1;
        featEnd    = sliceEnd - fRange.start + 1;
        featStrand = fRange.strand * -1 ;
      } else {
        featStart  = fRange.start - sliceStart + 1;
        featEnd    = fRange.end - sliceStart + 1;
        featStrand = fRange.strand;
      }

      daf = DNAAlignFeature_new();

      DNAAlignFeature_setDbID(daf,row->getLongLongAt(row,0));
      DNAAlignFeature_setContig(daf,slice); 
      DNAAlignFeature_setAnalysis(daf,analysis);

      DNAAlignFeature_setStart(daf,featStart);
      DNAAlignFeature_setEnd(daf,featEnd);
      DNAAlignFeature_setStrand(daf,featStrand);

      DNAAlignFeature_setHitStart(daf,row->getIntAt(row,6));
      DNAAlignFeature_setHitEnd(daf,row->getIntAt(row,7));
      DNAAlignFeature_setHitSeqName(daf,row->getStringAt(row,8));
      DNAAlignFeature_setHitStrand(daf,row->getIntAt(row,9));

      DNAAlignFeature_setCigarString(daf,row->getStringAt(row,10));
  
      if (row->col(row,11)) DNAAlignFeature_setpValue(daf,row->getDoubleAt(row,11));
      if (row->col(row,12)) DNAAlignFeature_setPercId(daf,row->getDoubleAt(row,12));
      if (row->col(row,13)) DNAAlignFeature_setScore(daf,row->getDoubleAt(row,13));

      Vector_addElement(features,daf);
    }
  } else { // No slice

    while ((row = sth->fetchRow(sth))) {
      DNAAlignFeature *daf;
      
// Perl has a cache for analysis types but the analysis adaptor should have one
      Analysis  *analysis = AnalysisAdaptor_fetchByDbID(aa, row->getLongLongAt(row,2));
// Perl has a cache for contigs - maybe important
      RawContig *contig = RawContigAdaptor_fetchByDbID(rca, row->getLongLongAt(row,1));

      daf = DNAAlignFeature_new();

      DNAAlignFeature_setDbID(daf,row->getLongLongAt(row,0));
      DNAAlignFeature_setContig(daf,contig); 
      DNAAlignFeature_setAnalysis(daf,analysis);

      DNAAlignFeature_setStart(daf,row->getIntAt(row,3));
      DNAAlignFeature_setEnd(daf,row->getIntAt(row,4));
      DNAAlignFeature_setStrand(daf,row->getIntAt(row,5));

      DNAAlignFeature_setHitStart(daf,row->getIntAt(row,6));
      DNAAlignFeature_setHitEnd(daf,row->getIntAt(row,7));
      DNAAlignFeature_setHitSeqName(daf,row->getStringAt(row,8));
      DNAAlignFeature_setHitStrand(daf,row->getIntAt(row,9));

      DNAAlignFeature_setCigarString(daf,row->getStringAt(row,10));
  
      if (row->col(row,11)) DNAAlignFeature_setpValue(daf,row->getDoubleAt(row,11));
      if (row->col(row,12)) DNAAlignFeature_setPercId(daf,row->getDoubleAt(row,12));
      if (row->col(row,13)) DNAAlignFeature_setScore(daf,row->getDoubleAt(row,13));

      Vector_addElement(features,daf);
    }
  }
  
  return features;
}
