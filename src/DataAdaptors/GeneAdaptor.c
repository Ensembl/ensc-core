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

#include "StatementHandle.h"
#include "ResultRow.h"

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
  IDType *transcriptIds;
  int nTranscriptId;
  int i;
  char qStr[512];
  StatementHandle *sth;
  ResultRow *row;

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
    "  AND gene.gene_id = "
    IDFMTSTR
    " ORDER BY tscript.gene_id"
    "  , tscript.transcript_id"
    "  , e_t.rank",geneId);

  sth = ga->prepare((BaseAdaptor *)ga, qStr, strlen(qStr) );
  sth->execute(sth);

  transcriptExonsHash = IDHash_new(IDHASH_SMALL);
  translationHash = IDHash_new(IDHASH_SMALL);

  while (row = sth->fetchRow(sth)) {
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
      first = 0;
    }

    // store an array of exon ids for each transcript
    transcriptId = row->getLongLongAt(row,1);
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

    tes->exonIds[tes->nExon++] = row->getLongLongAt(row,2);

// Note using string because its allocated
    if (!IDHash_contains(translationHash,row->getLongAt(row,1))) {
      IDHash_add(translationHash, row->getLongAt(row,1), row->getStringAt(row,6));
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
    Transcript *transcript = Transcript_new();
    IDType transcriptId = transcriptIds[i];
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
      Transcript_addExon(transcript, (Exon *)IDHash_getValue(exonHash, tes->exonIds[i]));
    }

    Transcript_setType(transcript, Gene_getType(gene));
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
          "SELECT stable_id, UNIX_TIMESTAMP(created),"
          "                  UNIX_TIMESTAMP(modified), version"
          " FROM gene_stable_id"
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


Set *GeneAdaptor_fetchAllBySlice(GeneAdaptor *ga, Slice *slice, char *logicName) {
  char sliceName[256];
  char sliceCacheKey[512];
  AssemblyMapperAdaptor *ama;
  AssemblyMapper *assMapper;
  IDType *contigIds;
  int nContigId;
  int i;
  char *qStr;
  StatementHandle *sth;
  ResultRow *row;
  Set *geneSet;

  if (logicName) {
    sprintf(sliceCacheKey,"%s:%s",Slice_getName(slice),logicName);
  } else {
    sprintf(sliceCacheKey,"%s:",Slice_getName(slice));
  }
  StrUtil_strupr(sliceCacheKey);

  // check the cache which uses the slice name as it key
  if(StringHash_contains(ga->sliceGeneCache,sliceCacheKey)) {
    return (Set *)StringHash_getValue(ga->sliceGeneCache,sliceCacheKey);
  }

  ama = DBAdaptor_getAssemblyMapperAdaptor(ga->dba);
  assMapper = AssemblyMapperAdaptor_fetchByType(ama,Slice_getAssemblyType(slice));

  nContigId = AssemblyMapper_listContigIds(assMapper,
                                           Slice_getChrId(slice),
                                           Slice_getChrStart(slice),
                                           Slice_getChrEnd(slice),
                                           &contigIds);


  if (!nContigId) {
    return emptySet;
  }

  qStr = StrUtil_copyString(&qStr, 
    "SELECT distinct(t.gene_id)"
    " FROM   transcript t,exon_transcript et,exon e, gene g"
    " WHERE e.contig_id in (", 0);

  if (!qStr) {
    Error_trace("fetch_all_by_Slice",NULL);
    return emptySet;
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

  qStr = StrUtil_appendString(qStr, 
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
      return emptySet;
    }

    sprintf(analStr," AND g.analysis_id = " IDFMTSTR , Analysis_getDbID(analysis));
   
    qStr = StrUtil_appendString(qStr,analStr); 
  }

  geneSet = Set_new();

  sth = ga->prepare((BaseAdaptor *)ga, qStr, strlen(qStr) );

  sth->execute(sth);

  while (row = sth->fetchRow(sth)) {
    IDType geneId = row->getLongLongAt(row,0);
    Gene *gene  = GeneAdaptor_fetchByDbID(ga, geneId, NULL );
    Gene *newGene = Gene_transformToSlice(gene, slice);

/* NIY
    if (Gene_getStart(newgene) <= Slice_getLength(slice) &&
        Gene_getEnd(gene) >= 1 ) {
      // only take the gene if its really overlapping the Slice
      push( @out, $newgene );
    }
*/
    Set_addElement(geneSet, newGene);
  }
  sth->finish(sth);

  free(qStr);

  StringHash_add(ga->sliceGeneCache, sliceCacheKey, geneSet);

  return geneSet;
}

