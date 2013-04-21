#include "GeneAdaptor.h"
#include "AnalysisAdaptor.h"
#include "AssemblyMapperAdaptor.h"
#include "DBAdaptor.h"
#include "DBEntryAdaptor.h"
#include "StrUtil.h"
#include "BaseAdaptor.h"
#include "MysqlUtil.h"
#include "Exon.h"
#include "IDHash.h"
#include "Transcript.h"
#include "Slice.h"
#include "AssemblyMapper.h"

#include "ExonAdaptor.h"
#include "TranscriptAdaptor.h"

#include "StatementHandle.h"
#include "ResultRow.h"

#include "Error.h"

GeneAdaptor *GeneAdaptor_new(DBAdaptor *dba) {
  GeneAdaptor *ga;

  if ((ga = (GeneAdaptor *)calloc(1,sizeof(GeneAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for GeneAdaptor\n");
    return NULL;
  }
  BaseAdaptor_init((BaseAdaptor *)ga, dba, GENE_ADAPTOR);
  ga->sliceGeneCache = StringHash_new(STRINGHASH_MEDIUM);

  return ga;
}

int GeneAdaptor_listGeneIds(GeneAdaptor *ga, IDType **geneIds) {
  char *qStr = "SELECT gene_id from gene";
  StatementHandle *sth = ga->prepare((BaseAdaptor *)ga,qStr,strlen(qStr));
  ResultRow *row;
  int idCount = 0;

  sth->execute(sth);

  *geneIds = NULL;
  
  while ((row = sth->fetchRow(sth))) {
    if (!idCount || !(idCount%10)) {
      if ((*geneIds = (IDType *)realloc(*geneIds,(idCount+10)*sizeof(IDType))) == NULL) {
        fprintf(stderr,"ERROR: Failed allocating geneIds\n");
        return -1;
      }
    }
    (*geneIds)[idCount++] = row->getLongLongAt(row,0);
  }
  sth->finish(sth);
  printf("Finished\n");
  return idCount;
}

#define MAXEXON 1024

typedef struct TranscriptExonsStruct {
  int nExon;
  IDType exonIds[MAXEXON];
} TranscriptExons;

Gene *GeneAdaptor_fetchByDbID(GeneAdaptor *ga, IDType geneId, int chrCoords) {
  Gene *gene = NULL;
  ExonAdaptor *ea;
  TranscriptAdaptor *ta;
  AnalysisAdaptor *aa;
  DBEntryAdaptor *dbea;
  Exon **exons;
  int nExon;
  int first = 1;
  Analysis *ana = NULL;
  IDHash *exonHash;
  IDHash *translationHash;
  IDHash *transcriptExonsHash;
  IDHash *transcriptHash;
  IDType *transcriptIds;
  int nTranscriptId;
  int i;
  char qStr[1024];
  StatementHandle *sth;
  ResultRow *row;

  ea = DBAdaptor_getExonAdaptor(ga->dba);
  ta = DBAdaptor_getTranscriptAdaptor(ga->dba);
  aa = DBAdaptor_getAnalysisAdaptor(ga->dba);
  dbea = DBAdaptor_getDBEntryAdaptor(ga->dba);

  nExon = ExonAdaptor_fetchAllByGeneId( ea, geneId, &exons);

  if( !nExon ) {
    fprintf(stderr,"ERROR: No exons for gene " IDFMTSTR ", assumming no gene\n",geneId);
    return NULL;
  }

  
/*
  fetching all exons by gene
  fetching all transcripts
  adding the exons
  adding the transcripts
*/

  sprintf(qStr,
    "SELECT tscript.gene_id"
    "  , tscript.transcript_id"
    "  , e_t.exon_id"
    "  , e_t.rank"
    "  , gene.analysis_id"
    "  , gene.biotype"
    "  , tscript.canonical_translation_id"
    "  , tscript.seq_region_start"
    "  , tscript.seq_region_end"
    "  , tscript.seq_region_strand"
    "  , tscript.biotype"
    "  , gene.display_xref_id"
    " FROM gene"
    "  , transcript tscript"
    "  , exon_transcript e_t"
    " WHERE gene.gene_id = tscript.gene_id"
    "  AND tscript.transcript_id = e_t.transcript_id"
    "  AND gene.gene_id = "
    IDFMTSTR
    " ORDER BY tscript.gene_id"
    "  , tscript.transcript_id"
    "  , e_t.rank",geneId);

  //printf("Query = %s\n",qStr);
  sth = ga->prepare((BaseAdaptor *)ga, qStr, strlen(qStr) );
  sth->execute(sth);

  transcriptExonsHash = IDHash_new(IDHASH_SMALL);
  transcriptHash = IDHash_new(IDHASH_SMALL);
  translationHash = IDHash_new(IDHASH_SMALL);

  while ((row = sth->fetchRow(sth))) {
    // building a gene
    TranscriptExons *tes = NULL;
    IDType transcriptId;

    if( first ) {
      gene = Gene_new();
      Gene_setAdaptor(gene,(BaseAdaptor *)ga);
      Gene_setDbID(gene, geneId );
      ana = AnalysisAdaptor_fetchByDbID(aa,row->getLongLongAt(row,4));
      Gene_setAnalysis(gene,ana);
// Allocating twice???
      Gene_setType(gene,row->getStringAt(row,5));

      IDType displayXrefId = row->getLongLongAt(row,11);
      if (displayXrefId) {
        DBEntry *displayXref = DBEntryAdaptor_fetchByDbID(dbea,displayXrefId);
        Gene_setDisplayXref(gene, displayXref);
      }
      first = 0;
    }

    // store an array of exon ids for each transcript
    transcriptId = row->getLongLongAt(row,1);
    if( !IDHash_contains(transcriptExonsHash,transcriptId)) {

      //printf("Adding new trans %d\n",transcriptId);
      if ((tes = (TranscriptExons *)calloc(1,sizeof(TranscriptExons))) == NULL) {
        fprintf(stderr,"ERROR: Failed allocating TranscriptExons\n");
        exit(1);
      }
      IDHash_add(transcriptExonsHash,transcriptId,tes);

      Transcript *transcript = Transcript_new();
      Transcript_setStart(transcript, row->getIntAt(row,7));
      Transcript_setEnd(transcript, row->getIntAt(row,8));
      Transcript_setStrand(transcript, row->getIntAt(row,9));
      Transcript_setType(transcript, row->getStringAt(row,10));
      IDHash_add(transcriptHash,transcriptId,transcript);
    }

    tes = IDHash_getValue(transcriptExonsHash,transcriptId);

    if (tes->nExon >= MAXEXON) {
      fprintf(stderr,"ERROR: MAXEXON exceeded\n");
      exit(1);
    }

    tes->exonIds[tes->nExon++] = row->getLongLongAt(row,2);

// Note using string because its allocated
    if (!IDHash_contains(translationHash,row->getLongLongAt(row,1))) {
      IDHash_add(translationHash, row->getLongLongAt(row,1), row->getStringAt(row,6));
    }
  }

  sth->finish(sth);

  if ( first ) {
    fprintf(stderr,"ERROR: Weird no first error\n");
    IDHash_free(translationHash, NULL);
// Should free exons
    return NULL;
  }
  
  //discard duplicate exons, add analysis object to exons
  exonHash = IDHash_new(IDHASH_SMALL);
  for (i=0; i<nExon; i++) {
    Exon_setAnalysis(exons[i],ana);
    IDHash_add(exonHash,Exon_getDbID(exons[i]),exons[i]);
  }

  // Exons array not needed anymore
  free(exons);

  transcriptIds = IDHash_getKeys(transcriptExonsHash);
  nTranscriptId = IDHash_getNumValues(transcriptExonsHash);
  
  for (i=0;i<nTranscriptId;i++) {
    IDType transcriptId = transcriptIds[i];
    //printf("transcript id = %ld\n",transcriptId);
    Transcript *transcript = IDHash_getValue(transcriptHash,transcriptId);
    IDType translationId;
    TranscriptExons *tes = IDHash_getValue(transcriptExonsHash,transcriptId);
    int i;
    sscanf((char *)IDHash_getValue(translationHash,transcriptId),IDFMTSTR,&translationId);

    Transcript_setDbID( transcript, transcriptId );
    Transcript_setAdaptor(transcript, (BaseAdaptor *)ta);
    
    // Test for truth because translation_id will be 0 if not set
    // because column is NOT NULL in schema.
    
    if (translationId) {
      Transcript_setTranslationId(transcript,translationId);
    }

    for (i=0;i<tes->nExon;i++) {
      Exon *exon = IDHash_getValue(exonHash, tes->exonIds[i]);
      //printf(" Exon from tes = " IDFMTSTR "\n",Exon_getDbID(exon));
      Transcript_addExon(transcript, exon);
    }

    //Transcript_setType(transcript, Gene_getType(gene));
    Gene_addTranscript(gene, transcript);
  }

  free(transcriptIds);
  
  // if chromosomal coordinates are needed, transform with empty slice

  if( chrCoords ) {
    SliceAdaptor *sa = DBAdaptor_getSliceAdaptor(ga->dba);
    //Slice *emptySlice = Slice_new(NULL,0,0,0,NULL,sa,0,TRUE);
    Slice_setAdaptor(emptySlice, (BaseAdaptor *)sa);
    //Slice_setEmptyFlag(emptySlice,1);
    // NIY $gene->transform( $empty_slice );
  }

  IDHash_free(exonHash,NULL);
  IDHash_free(transcriptHash,NULL);
  IDHash_free(translationHash,free);
  IDHash_free(transcriptExonsHash,free);

  return gene;
}

int GeneAdaptor_getStableEntryInfo(GeneAdaptor *ga, Gene *gene) {
  char qStr[256];
  StatementHandle *sth;
  ResultRow *row;

  if( !gene ) {
    fprintf(stderr, "ERROR: GeneAdaptor_getStableEntryInfo needs a gene object\n");
    exit(1); 
  }

  sprintf(qStr,
          "SELECT stable_id, UNIX_TIMESTAMP(created_date),"
          "                  UNIX_TIMESTAMP(modified_date), version"
          " FROM gene"
          " WHERE gene_id = "
          IDFMTSTR, Gene_getDbID(gene));

  sth = ga->prepare((BaseAdaptor *)ga,qStr,strlen(qStr));
  sth->execute(sth);

  row = sth->fetchRow(sth);
  if( row == NULL ) {
    fprintf(stderr,"WARNING: Failed fetching stable id info\n");
    sth->finish(sth);
    return 0;
  }

  Gene_setStableId(gene,row->getStringAt(row,0));
  Gene_setCreated(gene,row->getIntAt(row,1));
  Gene_setModified(gene,row->getIntAt(row,2));
  Gene_setVersion(gene,row->getIntAt(row,3));

  sth->finish(sth);

  return 1;
}


Vector *GeneAdaptor_fetchAllBySlice(GeneAdaptor *ga, Slice *slice, char *logicName) {
  char sliceCacheKey[512];
  AssemblyMapperAdaptor *ama;
  AssemblyMapper *assMapper;
  IDType *contigIds;
  int nContigId;
  int i;
  char *qStr;
  char tmpStr[1024];
  StatementHandle *sth;
  ResultRow *row;
  Vector *geneVector;

  if (logicName) {
    sprintf(sliceCacheKey,"%s:%s",Slice_getName(slice),logicName);
  } else {
    sprintf(sliceCacheKey,"%s:",Slice_getName(slice));
  }
  StrUtil_strupr(sliceCacheKey);

  // check the cache which uses the slice name as it key
  if(StringHash_contains(ga->sliceGeneCache,sliceCacheKey)) {
    return (Vector *)StringHash_getValue(ga->sliceGeneCache,sliceCacheKey);
  }

/*
  ama = DBAdaptor_getAssemblyMapperAdaptor(ga->dba);
  assMapper = AssemblyMapperAdaptor_fetchByType(ama,Slice_getAssemblyType(slice));

  nContigId = AssemblyMapper_listContigIds(assMapper,
                                           Slice_getChrId(slice),
                                           Slice_getChrStart(slice),
                                           Slice_getChrEnd(slice),
                                           &contigIds);


  if (!nContigId) {
    return emptyVector;
  }

  qStr = StrUtil_copyString(&qStr, 
    "SELECT distinct(t.gene_id)"
    " FROM   transcript t,exon_transcript et,exon e, gene g"
    " WHERE e.contig_id in (", 0);

  if (!qStr) {
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
    qStr = StrUtil_appendString(qStr, numStr);
  }
*/
  sprintf(tmpStr,
    "SELECT distinct(t.gene_id), g.seq_region_start, g.seq_region_end, g.seq_region_strand "
    " FROM   transcript t,exon_transcript et,exon e, gene g, seq_region sr"
    " WHERE e.seq_region_id = sr.seq_region_id and sr.name='%s' and g.seq_region_start <= %d and g.seq_region_end >= %d ", 
    Slice_getChrName(slice), Slice_getChrEnd(slice), Slice_getChrStart(slice));
  qStr = StrUtil_copyString(&qStr, tmpStr, 0);

  if (!qStr) {
    Error_trace("fetch_all_by_Slice",NULL);
    return emptyVector;
  }

  qStr = StrUtil_appendString(qStr, 
              " AND   et.exon_id = e.exon_id"
              " AND   et.transcript_id = t.transcript_id"
              " AND   g.gene_id = t.gene_id");

  if (logicName) {
    // determine analysis id via logic_name
    AnalysisAdaptor *aa = DBAdaptor_getAnalysisAdaptor(ga->dba);
    Analysis *analysis = AnalysisAdaptor_fetchByLogicName(aa,logicName);
    char analStr[256];

    if (!analysis || !Analysis_getDbID(analysis)) {
      fprintf(stderr,"WARNING: No analysis for logic name %s exists\n",logicName);
      return emptyVector;
    }

    sprintf(analStr," AND g.analysis_id = " IDFMTSTR , Analysis_getDbID(analysis));
   
    qStr = StrUtil_appendString(qStr,analStr); 
  }

  geneVector = Vector_new();

  sth = ga->prepare((BaseAdaptor *)ga, qStr, strlen(qStr) );

  sth->execute(sth);

  while ((row = sth->fetchRow(sth))) {
    IDType geneId = row->getLongLongAt(row,0);
    int start = row->getIntAt(row,1);
    int end = row->getIntAt(row,2);
    int strand = row->getIntAt(row,3);

    Gene *gene  = GeneAdaptor_fetchByDbID(ga, geneId, 0 );

    Gene_setStart(gene,start);
    Gene_setEnd(gene,end);
    Gene_setStrand(gene,strand);
   // Gene *newGene = Gene_transformToSlice(gene, slice);

/* NIY
    if (Gene_getStart(newgene) <= Slice_getLength(slice) &&
        Gene_getEnd(gene) >= 1 ) {
      // only take the gene if its really overlapping the Slice
      push( @out, $newgene );
    }
*/
    //Vector_addElement(geneVector, newGene);
    Vector_addElement(geneVector, gene);
  }
  sth->finish(sth);

  free(qStr);

  StringHash_add(ga->sliceGeneCache, sliceCacheKey, geneVector);

  return geneVector;
}

IDType GeneAdaptor_store(GeneAdaptor *ga, Gene *gene) {
  TranscriptAdaptor *ta;
  AnalysisAdaptor *aa;
  IDType analysisId;
  IDType xrefId;
  IDType geneId;
  int transCount;
  char *type = NULL;
  Analysis *analysis;
  StatementHandle *sth;
  char qStr[1024];
  int i;


  ta = DBAdaptor_getTranscriptAdaptor(ga->dba);
/* NIY
   if( !defined $gene || !ref $gene || !$gene->isa('Bio::EnsEMBL::Gene') ) {
       $self->throw("Must store a gene object, not a $gene");
   }
*/

  if ( Gene_getAnalysis(gene)) {
    fprintf(stderr,"Genes must have an analysis object!");
    exit(1);
  }

  aa = DBAdaptor_getAnalysisAdaptor(ga->dba);
  analysis = Gene_getAnalysis(gene);
  analysisId = AnalysisAdaptor_analysisExists(aa, analysis);

  if (analysisId) {
    Analysis_setDbID(analysis, analysisId);
  } else {
    AnalysisAdaptor_store(aa, analysis);
    analysisId = Analysis_getDbID(analysis);
  }

  transCount = Gene_getTranscriptCount(gene);

  if (Gene_getType(gene)) {
    type = Gene_getType(gene);
  } else {
    type = emptyString; /* static "" */
  }

  // assuming that the store is used during the Genebuil process, set
  // the display_xref_id to 0.  This ought to get re-set during the protein
  // pipeline run.  This probably update to the gene table has yet to be
  // implemented.

  xrefId = 0;
  sprintf(qStr,
         "INSERT INTO gene(biotype, analysis_id, transcript_count, display_xref_id) "
         "VALUES('%s'," IDFMTSTR ", %d, " IDFMTSTR ")",
         type, (IDType)analysisId, transCount, (IDType)xrefId);

  sth = ga->prepare((BaseAdaptor *)ga,qStr,strlen(qStr));
  sth->execute(sth);

  geneId = sth->getInsertId(sth);
  sth->finish(sth);

  if (Gene_getStableId(gene)) {
    if (!Gene_getCreated(gene) ||
        !Gene_getModified(gene) ||
        Gene_getVersion(gene) == -1) {
      fprintf(stderr,"Error: Trying to store incomplete stable id information for gene\n");
      exit(1);
    }
    sprintf(qStr,"INSERT INTO gene_stable_id(gene_id," 
                              "version, stable_id, created_date, modified_date)"
                      " VALUES(" IDFMTSTR ",%d, '%s'," 
                               "FROM_UNIXTIME(%ld),"
                               "FROM_UNIXTIME(%ld))",
                  geneId, 
                  Gene_getVersion(gene),
                  Gene_getStableId(gene),
                  Gene_getCreated(gene),
                  Gene_getModified(gene));
    sth = ga->prepare((BaseAdaptor *)ga,qStr,strlen(qStr));
    sth->execute(sth);
    sth->finish(sth);
  }


  // 
  // store the dbentries associated with this gene
  //
/* NIY
  my $dbEntryAdaptor = $self->db->get_DBEntryAdaptor();

  foreach my $dbl ( @{$gene->get_all_DBLinks} ) {
    $dbEntryAdaptor->store( $dbl, $gene_dbID, "Gene" );
  }
*/

  // 
  // Update the genes display xref if it is set
  //
/* NIY
  my $display_xref = $gene->display_xref;
  if($display_xref) {
    if(my $dbe_id = $dbEntryAdaptor->exists($display_xref)) {
      my $dispxref_sth = $self->prepare('UPDATE gene SET display_xref_id = ?
                                        WHERE gene_id = ?');
      $dispxref_sth->execute($dbe_id, $gene_dbID);
    }
  }
*/


  // write exons transcripts and exon_transcript table

  for (i=0; i<Gene_getTranscriptCount(gene); i++) {
    Transcript *t = Gene_getTranscriptAt(gene,i);
    // force lazy loading of translations before new exon dbIDs are set
    Transcript_getTranslation(t);
    TranscriptAdaptor_store(ta, t, geneId);
  }

  Gene_setAdaptor(gene, (BaseAdaptor *)ga);
  Gene_setDbID(gene,geneId);

  return Gene_getDbID(gene);
}


