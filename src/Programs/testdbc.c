#include <stdio.h>

#include "EnsC.h"

#include "DBAdaptor.h"
#include "BaseAdaptor.h"
#include "GeneAdaptor.h"
#include "AnalysisAdaptor.h"
#include "ChromosomeAdaptor.h"
#include "SequenceAdaptor.h"
#include "SliceAdaptor.h"
#include "TranscriptAdaptor.h"
#include "RawContig.h"
#include "ExonAdaptor.h"

int main(int argc, char *argv[]) {
  DBAdaptor *dba;

  printf("Opening connection ...");
  //dba = DBAdaptor_new("localhost","root",NULL,"test_ensembl",3306,NULL);
  dba = DBAdaptor_new("kaka.sanger.ac.uk","anonymous",NULL,"homo_sapiens_core_12_31",3306,NULL);
  printf(" Done\n");

  printf("Assembly type %s\n",DBAdaptor_getAssemblyType(dba));

  {
    GeneAdaptor *ga = DBAdaptor_getGeneAdaptor(dba);
    long *geneIds;
    int   nGeneId;
    int   i;
    Gene *gene;

    nGeneId = GeneAdaptor_listGeneIds(ga,&geneIds);

    printf("nGeneId = %d\n",nGeneId);

    for (i=1;i<5;i++) {
      gene = GeneAdaptor_fetchByDbID(ga,i,FALSE);
      if (gene) printf("Gene %s %d\n",Gene_getStableId(gene), Gene_getDbID(gene));
    }

    
  }
  {
    AnalysisAdaptor *aa = DBAdaptor_getAnalysisAdaptor(dba);

    Analysis *anal = AnalysisAdaptor_fetchByDbID(aa,1);

    printf("Anal logic name = %s\n",anal->logicName);
    anal = AnalysisAdaptor_fetchByLogicName(aa,"ensembl");

    printf("Anal logic name = %s\n",anal->logicName);
  }
  {
    TranscriptAdaptor *ta = DBAdaptor_getTranscriptAdaptor(dba);
    Transcript *t = TranscriptAdaptor_fetchByDbID(ta,1);

    printf("Transcript %d translation id = %d\n",Transcript_getDbID(t), Transcript_getTranslationId(t));
  }
  {
    ExonAdaptor *ea = DBAdaptor_getExonAdaptor(dba);
    Exon *e = ExonAdaptor_fetchByDbID(ea,1);
    Exon **exons;
    int nExon;
    int i;

    printf("Exon %d %d %d\n",Exon_getStart(e), Exon_getEnd(e), Exon_getStrand(e));

    nExon = ExonAdaptor_fetchAllByGeneId(ea,53,&exons);

    printf("NExon = %d\n",nExon);
    for (i=0;i<nExon;i++) {
      e = exons[i];
      printf("Exon %d %d %d %d\n",Exon_getDbID(e), Exon_getStart(e), Exon_getEnd(e), Exon_getStrand(e));
    }
    
  }
  {
    RawContigAdaptor *rca = DBAdaptor_getRawContigAdaptor(dba);
  }
  {
    SliceAdaptor *sa = DBAdaptor_getSliceAdaptor(dba);
    Slice *slice = SliceAdaptor_fetchByChrStartEnd(sa,"1",1,1500000);
    Gene **genes = Slice_getAllGenes(slice,NULL);
  }
  {
    ChromosomeAdaptor *ca = DBAdaptor_getChromosomeAdaptor(dba);
    Chromosome *chrom = ChromosomeAdaptor_fetchByDbID(ca,1);
    
    printf("Chromosome name = %s\n",chrom->name);
    chrom = ChromosomeAdaptor_fetchByChrName(ca,"1");
    printf("Chromosome name = %s\n",chrom->name);
  }
  {
    SequenceAdaptor *sa = DBAdaptor_getSequenceAdaptor(dba);
    RawContig *rc = RawContig_new();
    char *rawseq;
    RawContig_setDbID(rc,24);
    rawseq = SequenceAdaptor_fetchByRawContigStartEndStrand(sa,rc,1,10,1);
    
    printf("Seq = %s\n",rawseq);
  }
}
