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

int PredictionTranscriptAdaptor_store(BaseFeatureAdaptor *bfa, Vector *features) {
  StatementHandle *sth;
  char qStr[1024];
  int i;

  sprintf(qStr,
      "INSERT INTO prediction_transcript ( prediction_transcript_id, exon_rank, "
                                         "contig_id, contig_start, contig_end, "
                                         "contig_strand, start_phase, score, "
                                         "p_value, analysis_id, exon_count ) "
        " VALUES ( %" IDFMTSTR ", %%d, %" IDFMTSTR 
                  ", %%d, %%d, %%d, %%d, %%f, %%f, %" IDFMTSTR ", %%d)"); 

  sth = bfa->prepare((BaseAdaptor *)bfa,qStr,strlen(qStr));

  for (i=0; i<Vector_getNumElement(features); i++) {
    PredictionTranscript *predTrans = Vector_getElementAt(features, i);
    int j;
    IDType dbID = 0;
    int rank = 1;
    IDType analysisId = Analysis_getDbID(PredictionTranscript_getAnalysis(predTrans));
    int nExon = PredictionTranscript_getExonCount(predTrans);

    if (predTrans->objectType != CLASS_PREDICTIONTRANSCRIPT) {
      fprintf(stderr,"Error: Object is not a EnsEMBL PredictionTranscript "
                     "- not dumping!");
      exit(1);
    }

    if (PredictionTranscript_getDbID(predTrans) && 
        (BaseFeatureAdaptor *)PredictionTranscript_getAdaptor(predTrans) == bfa ) {
      fprintf(stderr,"PT already stored\n");
    }


    for (j=0; j<nExon; j++) {
      Exon *exon = PredictionTranscript_getExonAt(predTrans, j);
      IDType contigId;
      int contigStart;
      int contigEnd;
      int contigStrand;
      int startPhase;
      int endPhase;
      double score;
      double pValue;

      if (!exon) { rank++; continue; }

      contigId = RawContig_getDbID(Exon_getContig(exon));
      contigStart = Exon_getStart(exon);
      contigEnd = Exon_getEnd(exon);
      contigStrand = Exon_getStrand(exon);

      startPhase = Exon_getPhase(exon);
      endPhase = Exon_getEndPhase(exon);

      score = Exon_getScore(exon);
      pValue = Exon_getpValue(exon);

      if (rank == 1) {
        sth->execute(sth, NULL, 1, contigId, contigStart,
                          contigEnd, contigStrand,
                          startPhase, score, pValue, analysisId,
                          nExon);
        dbID = sth->getInsertId(sth);
      } else {
        sth->execute(sth, dbID, rank, contigId, contigStart,
                      contigEnd, contigStrand,
                      startPhase, score, pValue, analysisId,
                      nExon);
      }
      rank++;
    }

    PredictionTranscript_setDbID(predTrans, dbID);
    PredictionTranscript_setAdaptor(predTrans,(BaseAdaptor *)bfa);
  }

  sth->finish(sth);
  return 1;
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

Vector *PredictionTranscriptAdaptor_objectsFromStatementHandle(BaseFeatureAdaptor *bfa,
                                                       StatementHandle *sth,
                                                       AssemblyMapper *assMapper,
                                                       Slice *slice) {
  AnalysisAdaptor *aa;
  RawContigAdaptor *rca;
  Vector *out;
  ResultRow *row;
  int sliceChrId;
  int sliceEnd;
  int sliceStart;
  int sliceStrand;
  int sliceLen;
  IDType ptId = -1;
  PredictionTranscript *predTrans = NULL;
  int transcriptSliceStart;
  int transcriptSliceEnd;
  int nExon = 0;
  int stableStart;
  int stableEnd;
  char *stableCtg;
  Analysis *analysis;
  RawContig *contig;

  out = Vector_new();

  aa = DBAdaptor_getAnalysisAdaptor(bfa->dba);
  rca = DBAdaptor_getRawContigAdaptor(bfa->dba);

  if (slice) {
    sliceChrId = Slice_getChrId(slice);
    sliceStart = Slice_getChrStart(slice);
    sliceEnd   = Slice_getChrEnd(slice);
    sliceStrand= Slice_getStrand(slice);
    sliceLen   = Slice_getLength(slice);
  }

  while ((row = sth->fetchRow(sth))) {
    IDType predictionTranscriptId = row->getLongLongAt(row,0);
    int contigId    = row->getLongLongAt(row,1);
    IDType analysisId= row->getLongLongAt(row,2);
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
          Vector_addElement(out,predTrans);
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

    PredictionTranscript_addExon(predTrans, exon, &exonRank);
    //PredictionTranscript_addExon(predTrans, exon);
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
      Vector_addElement(out,predTrans);
    }
  }

  return out;
}

char *PredictionTranscriptAdaptor_finalClause(void) {
  return  " order by p.prediction_transcript_id, p.exon_rank";
}

