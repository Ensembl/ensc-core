#include "GeneAdaptor.h"
#include "AnalysisAdaptor.h"
#include "AssemblyMapperAdaptor.h"
#include "StrUtil.h"
#include "BaseAdaptor.h"
#include "MysqlUtil.h"
#include "Exon.h"
#include "IDHash.h"
#include "Transcript.h"
#include "Slice.h"
#include "AssemblyMapper.h"


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

int GeneAdaptor_listGeneIds(GeneAdaptor *ga, long **geneIds) {
  char *qStr = "SELECT gene_id from gene";
  MYSQL_RES *results = ga->prepare((BaseAdaptor *)ga,qStr,strlen(qStr));
  MYSQL_ROW row;
  int idCount = 0;

  *geneIds = NULL;
  
  while ((row = mysql_fetch_row(results))) {
    if (!idCount || !(idCount%10)) {
      if ((*geneIds = (long *)realloc(*geneIds,(idCount+10)*sizeof(long))) == NULL) {
        fprintf(stderr,"ERROR: Failed allocating geneIds\n");
        return -1;
      }
    }
    (*geneIds)[idCount++] = MysqlUtil_getLong(row,0);
  }
  return idCount;
}

#define MAXEXON 1024

typedef struct TranscriptExonsStruct {
  int nExon;
  long exonIds[MAXEXON];
} TranscriptExons;

Gene *GeneAdaptor_fetchByDbID(GeneAdaptor *ga, long geneId, int chrCoords) {
  Gene *gene;
  ExonAdaptor *ea;
  TranscriptAdaptor *ta;
  AnalysisAdaptor *aa;
  Exon **exons;
  int nExon;
  int first = 1;
  Analysis *ana;
  IDHash *exonHash;
  IDHash *translationHash;
  IDHash *transcriptExonsHash;
  long *transcriptIds;
  int nTranscriptId;
  int i;
  char qStr[512];
  MYSQL_RES *results;
  MYSQL_ROW row;

  ea = DBAdaptor_getExonAdaptor(ga->dba);
  ta = DBAdaptor_getTranscriptAdaptor(ga->dba);
  aa = DBAdaptor_getAnalysisAdaptor(ga->dba);

  nExon = ExonAdaptor_fetchAllByGeneId( ea, geneId, &exons);

  if( !nExon ) {
    fprintf(stderr,"ERROR: No exons for gene %d, assumming no gene\n",geneId);
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
    "  , e_t.exon_id, e_t.rank"
    "  , gene.analysis_id"
    "  , gene.type"
    "  , tscript.translation_id"
    " FROM gene"
    "  , transcript tscript"
    "  , exon_transcript e_t"
    " WHERE gene.gene_id = tscript.gene_id"
    "  AND tscript.transcript_id = e_t.transcript_id"
    "  AND gene.gene_id = %d"
    " ORDER BY tscript.gene_id"
    "  , tscript.transcript_id"
    "  , e_t.rank",geneId);

  results = ga->prepare((BaseAdaptor *)ga, qStr, strlen(qStr) );

  transcriptExonsHash = IDHash_new(IDHASH_SMALL);
  translationHash = IDHash_new(IDHASH_SMALL);

  while (row = mysql_fetch_row(results)) {
    // building a gene
    TranscriptExons *tes = NULL;
    long transcriptId;

    if( first ) {
      gene = Gene_new();
      Gene_setAdaptor(gene,(BaseAdaptor *)ga);
      Gene_setDbID(gene, geneId );
      ana = AnalysisAdaptor_fetchByDbID(aa,MysqlUtil_getLong(row,4));
      Gene_setAnalysis(gene,ana);
// Allocating twice???
      Gene_setType(gene,MysqlUtil_getString(row,5));
      first = 0;
    }

    // store an array of exon ids for each transcript
    transcriptId = MysqlUtil_getLong(row,1);
    printf("Transcript ID = %d\n",transcriptId);
    if( !IDHash_contains(transcriptExonsHash,transcriptId)) {

      if ((tes = (TranscriptExons *)calloc(1,sizeof(TranscriptExons))) == NULL) {
        fprintf(stderr,"ERROR: Failed allocating TranscriptExons\n");
        exit(1);
      }
      IDHash_add(transcriptExonsHash,transcriptId,tes);
    }
    tes = IDHash_getValue(transcriptExonsHash,transcriptId);

    if (tes->nExon >= MAXEXON) {
      fprintf(stderr,"ERROR: MAXEXON exceeded\n");
      exit(1);
    }

    tes->exonIds[tes->nExon++] = MysqlUtil_getLong(row,2);

// Note using string because its allocated
    if (!IDHash_contains(translationHash,MysqlUtil_getLong(row,1))) {
      IDHash_add(translationHash, MysqlUtil_getLong(row,1), MysqlUtil_getString(row,6));
    }
  }

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
    Transcript *transcript = Transcript_new();
    long transcriptId = transcriptIds[i];
    long translationId = atol((char *)IDHash_getValue(translationHash,transcriptId));
    TranscriptExons *tes = IDHash_getValue(transcriptExonsHash,transcriptId);
    int i;

    Transcript_setDbID( transcript, transcriptId );
    Transcript_setAdaptor(transcript, (BaseAdaptor *)ta);
    
    // Test for truth because translation_id will be 0 if not set
    // because column is NOT NULL in schema.
    
    if (translationId) {
      Transcript_setTranslationId(transcript,translationId);
    }

    for (i=0;i<tes->nExon;i++) {
      Transcript_addExon(transcript, (Exon *)IDHash_getValue(exonHash, tes->exonIds[i]));
    }

    Transcript_setType(transcript, Gene_getType(gene));
    Gene_addTranscript(gene, transcript);
  }

  free(transcriptIds);
  
  // if chromosomal coordinates are needed, transform with empty slice

  if( chrCoords ) {
    SliceAdaptor *sa = DBAdaptor_getSliceAdaptor(ga->dba);
    Slice *emptySlice = Slice_new(NULL,0,0,0,NULL,sa,0,TRUE);
    Slice_setAdaptor(emptySlice, (BaseAdaptor *)sa);
    Slice_setEmptyFlag(emptySlice,1);
    // NIY $gene->transform( $empty_slice );
  }

  IDHash_free(exonHash,NULL);
  IDHash_free(translationHash,free);
  IDHash_free(transcriptExonsHash,free);

  return gene;
}

int GeneAdaptor_getStableEntryInfo(GeneAdaptor *ga, Gene *gene) {
  char qStr[256];
  MYSQL_RES *results;
  MYSQL_ROW row;

  if( !gene ) {
    fprintf(stderr, "ERROR: GeneAdaptor_getStableEntryInfo needs a gene object\n");
    exit(1); 
  }

  sprintf(qStr,
          "SELECT stable_id, UNIX_TIMESTAMP(created),"
          "                  UNIX_TIMESTAMP(modified), version"
          " FROM gene_stable_id"
          " WHERE gene_id = %d",Gene_getDbID(gene));

  results = ga->prepare((BaseAdaptor *)ga,qStr,strlen(qStr));

  row = mysql_fetch_row(results);
  if( row == NULL ) {
    fprintf(stderr,"WARNING: Failed fetching stable id info\n");
    return 0;
  }

  Gene_setStableId(gene,MysqlUtil_getString(row,0));
  Gene_setCreated(gene,MysqlUtil_getInt(row,1));
  Gene_setModified(gene,MysqlUtil_getInt(row,2));
  Gene_setVersion(gene,MysqlUtil_getInt(row,3));

  return 1;
}


int GeneAdaptor_fetchAllBySlice(GeneAdaptor *ga, Slice *slice, char *logicName, Gene ***outP) {
  Gene **out = *outP;
  char sliceName[256];
  char sliceCacheKey[512];
  AssemblyMapperAdaptor *ama;
  AssemblyMapper *assMapper;
  long *contigIds;
  int nContigId;
  int i;
  char *qStr;
  MYSQL_RES *results;
  MYSQL_ROW row;

  if (logicName) {
    sprintf(sliceCacheKey,"%s:%s",Slice_getName(slice,sliceName),logicName);
  } else {
    sprintf(sliceCacheKey,"%s:",Slice_getName(slice,sliceName));
  }
  StrUtil_strupr(sliceCacheKey);

  // check the cache which uses the slice name as it key
  if(StringHash_contains(ga->sliceGeneCache,sliceCacheKey)) {
    return (Gene **)StringHash_getValue(ga->sliceGeneCache,sliceCacheKey);
  }

  ama = DBAdaptor_getAssemblyMapperAdaptor(ga->dba);
  assMapper = AssemblyMapperAdaptor_fetchByType(ama,Slice_getAssemblyType(slice));

  nContigId = AssemblyMapper_listContigIds(assMapper,
                                           Slice_getChrId(slice),
                                           Slice_getChrStart(slice),
                                           Slice_getChrEnd(slice),
                                           &contigIds);


  if (!nContigId) {
    return 0;
  }

  qStr = StrUtil_CopyString(
    "SELECT distinct(t.gene_id)"
    " FROM   transcript t,exon_transcript et,exon e, gene g"
    " WHERE e.contig_id in (");

  for (i=0; i<nContigId; i++) {
    char numStr[256];
    if (i!=(nContigId-1)) {
      sprintf(numStr,"%d,",contigIds[i]);
    } else {
      sprintf(numStr,"%d",contigIds[i]);
    }
    qStr = StrUtil_AppendString(qStr, numStr);
  }

  qStr = StrUtil_AppendString(qStr, 
              ") AND   et.exon_id = e.exon_id"
              " AND   et.transcript_id = t.transcript_id"
              " AND   g.gene_id = t.gene_id");

  if (logicName) {
    // determine analysis id via logic_name
    AnalysisAdaptor *aa = DBAdaptor_getAnalysisAdaptor(ga->dba);
    Analysis *analysis = AnalysisAdaptor_fetchByLogicName(aa,logicName);
    char analStr[256];

    if (!analysis || !Analysis_getDbID(analysis)) {
      fprintf(stderr,"WARNING: No analysis for logic name %s exists\n",logicName);
      return 0;
    }

    sprintf(analStr," AND g.analysis_id = %d", Analysis_getDbID(analysis));
   
    qStr = StrUtil_AppendString(qStr,analStr); 
  }


  results = ga->prepare((BaseAdaptor *)ga, qStr, strlen(qStr) );

  while (row = mysql_fetch_row(results)) {
    long geneId = MysqlUtil_getLong(row,0);
    Gene *gene  = GeneAdaptor_fetchByDbID(ga, geneId, NULL );
    Gene *newGene = Gene_transformToSlice(gene, slice);

/*
    if (Gene_getStart(newgene) <= Slice_getLength(slice) &&
        Gene_getEnd(gene) >= 1 ) {
      // only take the gene if its really overlapping the Slice
      push( @out, $newgene );
    }
*/
  }

  free(qStr);

  //place the results in an LRU cache
/*
  $self->{'_slice_gene_cache'}{$key} = \@out;

  return \@out;
*/
}

