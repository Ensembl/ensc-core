#include "ProteinAlignFeatureAdaptor.h"

#include "RawContigAdaptor.h"
#include "AnalysisAdaptor.h"

#include "DNAPepAlignFeature.h"

NameTableType ProteinAlignFeatureAdaptor_tableNames = {{"protein_align_feature","paf"},{NULL,NULL}};

ProteinAlignFeatureAdaptor *ProteinAlignFeatureAdaptor_new(DBAdaptor *dba) {
  ProteinAlignFeatureAdaptor *pafa;

  if ((pafa = (ProteinAlignFeatureAdaptor *)calloc(1,sizeof(ProteinAlignFeatureAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for ProteinAlignFeatureAdaptor\n");
    return NULL;
  }
  BaseFeatureAdaptor_init((BaseFeatureAdaptor *)pafa, dba, PROTEINALIGNFEATURE_ADAPTOR);

  pafa->getTables = ProteinAlignFeatureAdaptor_getTables;
  pafa->getColumns = ProteinAlignFeatureAdaptor_getColumns;
  pafa->store = ProteinAlignFeatureAdaptor_store;
  pafa->objectsFromStatementHandle = ProteinAlignFeatureAdaptor_objectsFromStatementHandle;

  return pafa;
}

int ProteinAlignFeatureAdaptor_store(BaseFeatureAdaptor *bfa, Set *features) {
/*
   my ($self ,@sf) = @_;

   if( scalar(@sf) == 0 ) {
       $self->throw("Must call store with list of sequence features");
   }

   my $sth = $self->prepare(
        "INSERT INTO protein_align_feature(contig_id, contig_start, contig_end,
                                           contig_strand, hit_start, hit_end,
                                           hit_name, cigar_line, analysis_id,
                                           score, evalue, perc_ident) 
         VALUES (?,?,?,?,?,?,?,?,?,?, ?, ?)");

   foreach my $sf ( @sf ) {
     if( !ref $sf || !$sf->isa("Bio::EnsEMBL::DnaPepAlignFeature") ) {
       $self->throw("Feature must be a Bio::EnsEMBL::DnaPepAlignFeature, " .
                    "not a [$sf]");
     }
     
     if( !defined $sf->analysis ) {
       $self->throw("Cannot store sequence features without analysis");
     }

     # will only store if object is not already stored in this database
     $self->db()->get_AnalysisAdaptor()->store( $sf->analysis() );

     my $contig = $sf->entire_seq();
     #print STDERR $contig."\n";
     unless(defined $contig && $contig->isa("Bio::EnsEMBL::RawContig")) { 
       $self->throw("Cannot store feature without Contig attached via " .
                    "attach_seq\n");
     }   
     
     $sth->execute($contig->dbID(), $sf->start, $sf->end,
                   $sf->strand, $sf->hstart, $sf->hend,
                   $sf->hseqname, $sf->cigar_string, $sf->analysis->dbID,
                   $sf->score, $sf->p_value, $sf->percent_id);
     
     $sf->dbID($sth->{'mysql_insertid'});
   }
*/

}


NameTableType *ProteinAlignFeatureAdaptor_getTables(void) {
  return &ProteinAlignFeatureAdaptor_tableNames;
}

char *ProteinAlignFeatureAdaptor_getColumns(void) {

  return "paf.protein_align_feature_id,"
         "paf.contig_id,"
         "paf.analysis_id,"
         "paf.contig_start,"
         "paf.contig_end,"
         "paf.contig_strand,"
         "paf.hit_start,"
         "paf.hit_end,"
         "paf.hit_name,"
         "paf.cigar_line,"
         "paf.evalue,"
         "paf.perc_ident,"
         "paf.score";
}

Set *ProteinAlignFeatureAdaptor_objectsFromStatementHandle(BaseFeatureAdaptor *bfa,
                                                           StatementHandle *sth,
                                                           AssemblyMapper *assMapper,
                                                           Slice *slice) {
  AnalysisAdaptor *aa;
  RawContigAdaptor *rca;
  Set *features;
  ResultRow *row;
  int i;

  aa = DBAdaptor_getAnalysisAdaptor(bfa->dba);
  rca = DBAdaptor_getRawContigAdaptor(bfa->dba);

  features = Set_new();

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


    while (row = sth->fetchRow(sth)) {
      DNAPepAlignFeature *dpaf;
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

      dpaf = DNAPepAlignFeature_new();

      DNAPepAlignFeature_setDbID(dpaf,row->getLongLongAt(row,0));
      DNAPepAlignFeature_setContig(dpaf,slice); 
      DNAPepAlignFeature_setAnalysis(dpaf,analysis);

      DNAPepAlignFeature_setStart(dpaf,featStart);
      DNAPepAlignFeature_setEnd(dpaf,featEnd);
      DNAPepAlignFeature_setStrand(dpaf,featStrand);

      DNAPepAlignFeature_setHitStart(dpaf,row->getIntAt(row,6));
      DNAPepAlignFeature_setHitEnd(dpaf,row->getIntAt(row,7));
      DNAPepAlignFeature_setHitId(dpaf,row->getStringAt(row,8));
      DNAPepAlignFeature_setHitStrand(dpaf,1);

      DNAPepAlignFeature_setCigarString(dpaf,row->getStringAt(row,9));
  
      if (row->col(row,10)) DNAPepAlignFeature_setEValue(dpaf,row->getDoubleAt(row,10));
      if (row->col(row,11)) DNAPepAlignFeature_setPercId(dpaf,row->getDoubleAt(row,11));
      if (row->col(row,12)) DNAPepAlignFeature_setScore(dpaf,row->getDoubleAt(row,12));

      Set_addElement(features,dpaf);
    }
  } else { // No slice

    while(row = sth->fetchRow(sth)) {
      DNAPepAlignFeature *dpaf;
      
// Perl has a cache for analysis types but the analysis adaptor should have one
      Analysis  *analysis = AnalysisAdaptor_fetchByDbID(aa, row->getLongLongAt(row,2));
// Perl has a cache for contigs - maybe important
      RawContig *contig = RawContigAdaptor_fetchByDbID(rca, row->getLongLongAt(row,1));

      dpaf = DNAPepAlignFeature_new();

      DNAPepAlignFeature_setDbID(dpaf,row->getLongLongAt(row,0));
      DNAPepAlignFeature_setContig(dpaf,contig); 
      DNAPepAlignFeature_setAnalysis(dpaf,analysis);

      DNAPepAlignFeature_setStart(dpaf,row->getIntAt(row,3));
      DNAPepAlignFeature_setEnd(dpaf,row->getIntAt(row,4));
      DNAPepAlignFeature_setStrand(dpaf,row->getIntAt(row,5));

      DNAPepAlignFeature_setHitStart(dpaf,row->getIntAt(row,6));
      DNAPepAlignFeature_setHitEnd(dpaf,row->getIntAt(row,7));
      DNAPepAlignFeature_setHitId(dpaf,row->getStringAt(row,8));
      DNAPepAlignFeature_setHitStrand(dpaf,1);

      DNAPepAlignFeature_setCigarString(dpaf,row->getStringAt(row,9));
  
      if (row->col(row,10)) DNAPepAlignFeature_setEValue(dpaf,row->getDoubleAt(row,10));
      if (row->col(row,11)) DNAPepAlignFeature_setPercId(dpaf,row->getDoubleAt(row,11));
      if (row->col(row,12)) DNAPepAlignFeature_setScore(dpaf,row->getDoubleAt(row,12));

      Set_addElement(features,dpaf);
    }
  }
  
  return features;
}
