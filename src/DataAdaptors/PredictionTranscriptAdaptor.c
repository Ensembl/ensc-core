#include "PredictionTranscriptAdaptor.h"

#include "PredictionTranscript.h"
#include "AnalysisAdaptor.h"
#include "RawContigAdaptor.h"

NameTableType PredictionTranscriptAdaptor_tableNames = {{"prediction_transcript","p"},
                                                        {NULL,NULL}};

PredictionTranscriptAdaptor *PredictionTranscriptAdaptor_new(DBAdaptor *dba) {
  PredictionTranscriptAdaptor *pta;

  if ((pta = (PredictionTranscriptAdaptor *)calloc(1,sizeof(PredictionTranscriptAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for PredictionTranscriptAdaptor\n");
    return NULL;
  }
  BaseFeatureAdaptor_init((BaseFeatureAdaptor *)pta, dba, PREDICTIONTRANSCRIPT_ADAPTOR);

  pta->getTables = PredictionTranscriptAdaptor_getTables;
  pta->getColumns = PredictionTranscriptAdaptor_getColumns;
  pta->store = PredictionTranscriptAdaptor_store;
  pta->objectsFromStatementHandle = PredictionTranscriptAdaptor_objectsFromStatementHandle;
  pta->finalClause = PredictionTranscriptAdaptor_finalClause;

  return pta;
}

int PredictionTranscriptAdaptor_store(BaseFeatureAdaptor *bfa, Set *features) {
/*
  my ( $self, @pre_transcripts ) = @_;

  my $exon_sql = q{
      INSERT INTO prediction_transcript ( prediction_transcript_id, exon_rank,
                                          contig_id, contig_start, contig_end,
                                          contig_strand, start_phase, score,
                                          p_value, analysis_id, exon_count )
        VALUES ( ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ? )
      };

  my $exonst = $self->prepare($exon_sql);

  foreach my $pre_trans (@pre_transcripts) {
    if( ! $pre_trans->isa('Bio::EnsEMBL::PredictionTranscript') ) {
      $self->throw("$pre_trans is not a EnsEMBL PredictionTranscript "
                   . "- not dumping!");
    }

    if( $pre_trans->dbID && $pre_trans->adaptor == $self ) {
      $self->warn("Already stored");
    }

    my $exonId = undef;
    my $exons = $pre_trans->get_all_Exons();
    my $dbID = undef;
    my $rank = 1;

    for my $exon ( @$exons ) {
      if( ! defined $exon ) { $rank++; next; }

      my $contig_id = $exon->contig->dbID();
      my $contig_start = $exon->start();
      my $contig_end = $exon->end();
      my $contig_strand = $exon->strand();

      my $start_phase = $exon->phase();
      my $end_phase = $exon->end_phase();

      # this is only in PredictionExon
      my $score = $exon->score();
      my $p_value = $exon->p_value();

      my $analysis = $pre_trans->analysis->dbID;

      if( $rank == 1 ) {
        $exonst->execute( undef, 1, $contig_id, $contig_start,
                          $contig_end, $contig_strand,
                          $start_phase, $score, $p_value, $analysis,
                          scalar( @{$exons} ));
        $dbID = $exonst->{'mysql_insertid'};
      } else {
        $exonst->execute( $dbID, $rank, $contig_id, $contig_start,
                          $contig_end, $contig_strand,
                          $start_phase, $score, $p_value, $analysis,
                         scalar( @{$exons} ) );
      }
      $rank++;
    }

    $pre_trans->dbID( $dbID );
    $pre_trans->adaptor( $self );
  }

  $exonst->finish;
*/
}


NameTableType *PredictionTranscriptAdaptor_getTables(void) {
  return &PredictionTranscriptAdaptor_tableNames;
}

char *PredictionTranscriptAdaptor_getColumns(void) {
  return "p.prediction_transcript_id,"
         "p.contig_id,"
         "p.analysis_id,"
         "p.contig_start,"
         "p.contig_end,"
         "p.contig_strand,"
         "p.start_phase,"
         "p.exon_rank,"
         "p.score,"
         "p.p_value,"
         "p.exon_count";
}

Set *PredictionTranscriptAdaptor_objectsFromStatementHandle(BaseFeatureAdaptor *bfa,
                                                       StatementHandle *sth,
                                                       AssemblyMapper *assMapper,
                                                       Slice *slice) {
  AnalysisAdaptor *aa;
  RawContigAdaptor *rca;
  Set *out;
  ResultRow *row;
  int i;
  int sliceChrId;
  int sliceEnd;
  int sliceStart;
  int sliceStrand;
  int sliceLen;
  int64 ptId = -1;
  PredictionTranscript *predTrans = NULL;
  int transcriptSliceStart;
  int transcriptSliceEnd;
  int nExon = 0;
  int stableStart;
  int stableEnd;
  char *stableCtg;
  Analysis *analysis;
  RawContig *contig;

  out = Set_new();

  aa = DBAdaptor_getAnalysisAdaptor(bfa->dba);
  rca = DBAdaptor_getRawContigAdaptor(bfa->dba);

  if (slice) {
    sliceChrId = Slice_getChrId(slice);
    sliceStart = Slice_getChrStart(slice);
    sliceEnd   = Slice_getChrEnd(slice);
    sliceStrand= Slice_getStrand(slice);
    sliceLen   = Slice_getLength(slice);
  }

  while(row = sth->fetchRow(sth)) {
    int64 predictionTranscriptId = row->getLongLongAt(row,0);
    int contigId    = row->getLongLongAt(row,1);
    int64 analysisId= row->getLongLongAt(row,2);
    int contigStart = row->getIntAt(row,3);
    int contigEnd   = row->getIntAt(row,4);
    int contigStrand= row->getIntAt(row,5);
    int startPhase  = row->getIntAt(row,6);
    int exonRank    = row->getIntAt(row,7);
    double score    = row->getDoubleAt(row,8);
    double pValue   = row->getDoubleAt(row,9);
    int exonCount   = row->getIntAt(row,10);
    Exon *exon = NULL;
    int exonStart;
    int exonEnd;
    int exonStrand;
    
    //create a new transcript for each new prediction transcript id
    if (!predTrans || !(ptId == predictionTranscriptId)) {
 
      // throw away last pt if no exons or introns were on the slice
      if (predTrans) {
        if (slice && (transcriptSliceEnd < 1 ||
                      transcriptSliceStart > sliceLen)) {
          PredictionTranscript_free(predTrans);
        } else {
          // set the stable_id of the previous prediction
          char tmpStr[256]; 
          sprintf(tmpStr,"%s.%d.%d\n",stableCtg,stableStart,stableEnd);
          PredictionTranscript_setStableId(predTrans,tmpStr);
          Set_addElement(out,predTrans);
        }
      }

      // Now we've stored the old (if it was any good) start a new prediction transcript
      predTrans = PredictionTranscript_new();

      ptId = predictionTranscriptId;
      PredictionTranscript_setDbID(predTrans,ptId);

      analysis = AnalysisAdaptor_fetchByDbID(aa, analysisId);

      PredictionTranscript_setAnalysis(predTrans,analysis);
      PredictionTranscript_setExonCount(predTrans,exonCount);

      //reset values used for last predtrans
      stableStart = -1;
      stableEnd   = -1;
      stableCtg   = NULL;

      nExon = 0;
    }

    //recalculate stable id values
    if(stableStart == -1 || contigStart < stableStart) {
      stableStart = contigStart;
    }
    if(contigEnd > stableEnd) {
      stableEnd = contigEnd;
    }

    contig = RawContigAdaptor_fetchByDbID(rca, contigId);

    stableCtg = RawContig_getName(contig);

    if (slice) {
      //a slice was passed in so we want slice coords
      MapperCoordinate fRange;

      //convert contig coordinates to assembly coordinates
      int mapSucceeded = AssemblyMapper_fastToAssembly(assMapper, contigId, 
                                               contigStart, 
                                               contigEnd, 
                                               contigStrand, 
                                               &fRange);

      // undefined start means gap
      if (!mapSucceeded) continue;
  
      // convert assembly coordinates to slice coordinates
      if(sliceStrand == -1) {
        exonStart  = sliceEnd - fRange.end + 1;
        exonEnd    = sliceEnd - fRange.start + 1;
        exonStrand = fRange.strand * -1 ;
      } else {
        exonStart  = fRange.start - sliceStart + 1;
        exonEnd    = fRange.end - sliceStart + 1;
        exonStrand = fRange.strand;
      }

      if( !nExon || transcriptSliceStart > exonStart ) {
        transcriptSliceStart = exonStart;
      }

      if( !nExon || transcriptSliceEnd < exonEnd ) {
        transcriptSliceEnd = exonEnd;
      }
    } else {
      // we just want plain old contig coords
      exonStart =  contigStart;
      exonEnd   =  contigEnd;
      exonStrand = contigStrand;
    }

    // create an exon and add it to the prediction transcript
    exon = Exon_new();

    Exon_setStart(exon, exonStart);
    Exon_setEnd(exon, exonEnd);
    Exon_setStrand(exon, exonStrand);
    if (slice) {
      Exon_setContig(exon, slice);
    } else {
      Exon_setContig(exon, contig);
    }
    Exon_setPhase(exon, startPhase);
    Exon_setEndPhase(exon, (exonEnd - exonStart + 1 + startPhase) % 3 );
    Exon_setScore(exon, score);
    Exon_setpValue(exon, pValue);

    //PredictionTranscript_addExon(predTrans, exon, exonRank);
    PredictionTranscript_addExon(predTrans, exon);
    nExon++;
  }

  // throw away last pt if no exons or introns were on the slice
  if(predTrans) {
    if (slice && (transcriptSliceEnd < 1 ||
                  transcriptSliceStart > sliceLen)) {
      PredictionTranscript_free(predTrans);
    } else {
      // set the stable_id of the previous prediction
      char tmpStr[256]; 
      sprintf(tmpStr,"%s.%d.%d\n",stableCtg,stableStart,stableEnd);
      PredictionTranscript_setStableId(predTrans,tmpStr);
      Set_addElement(out,predTrans);
    }
  }

  return out;
}

char *PredictionTranscriptAdaptor_finalClause(void) {
  return  " order by p.prediction_transcript_id, p.exon_rank";
}

