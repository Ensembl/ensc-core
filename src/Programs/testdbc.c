#include <stdio.h>

#include "EnsC.h"

#include "DBAdaptor.h"
#include "AnalysisAdaptor.h"
#include "BaseAdaptor.h"
#include "ChromosomeAdaptor.h"
#include "DNAAlignFeatureAdaptor.h"
#include "ExonAdaptor.h"
#include "GeneAdaptor.h"
#include "ProteinAlignFeatureAdaptor.h"
#include "RawContigAdaptor.h"
#include "RepeatFeatureAdaptor.h"
#include "SequenceAdaptor.h"
#include "SimpleFeatureAdaptor.h"
#include "SliceAdaptor.h"
#include "TranscriptAdaptor.h"

#include "DNAAlignFeature.h"
#include "DNAPepAlignFeature.h"
#include "PredictionTranscript.h"
#include "RawContig.h"
#include "RepeatFeature.h"
#include "SimpleFeature.h"

int main(int argc, char *argv[]) {
  DBAdaptor *dba;

  initEnsC();

  printf("Opening connection ...");
  //dba = DBAdaptor_new("localhost","root",NULL,"test_ensembl",3306,NULL);
  dba = DBAdaptor_new("kaka.sanger.ac.uk","anonymous",NULL,"homo_sapiens_core_12_31",3306,NULL);
  //dba = DBAdaptor_new("ecs2d.internal.sanger.ac.uk","ensro",NULL,"homo_sapiens_core_13_31",3306,NULL);
  printf(" Done\n");

  printf("Assembly type %s\n",DBAdaptor_getAssemblyType(dba));

  {
    SimpleFeatureAdaptor *sfa = DBAdaptor_getSimpleFeatureAdaptor(dba);
    SimpleFeature *sf = (SimpleFeature *)SimpleFeatureAdaptor_fetchByDbID(sfa,1);
    printf("Simple feature: %d-%d id " IDFMTSTR "\n", SimpleFeature_getStart(sf), 
           SimpleFeature_getEnd(sf), SimpleFeature_getDbID(sf));
  }
  {
    SliceAdaptor *sa = DBAdaptor_getSliceAdaptor(dba);
    Slice *slice = SliceAdaptor_fetchByChrStartEnd(sa,"1",1,2000000);
    Set *ptSet = Slice_getAllPredictionTranscripts(slice,"");
    int i;

    for (i=0;i<Set_getNumElement(ptSet);i++) {
      PredictionTranscript *pt = Set_getElementAt(ptSet,i);
      printf("PredictionTranscript feature: %s %d-%d id " IDFMTSTR "\n", 
             PredictionTranscript_getStableId(pt), 
             PredictionTranscript_getStart(pt), PredictionTranscript_getEnd(pt), 
             PredictionTranscript_getDbID(pt));
    }
  }
  {
    SliceAdaptor *sa = DBAdaptor_getSliceAdaptor(dba);
    Slice *slice = SliceAdaptor_fetchByChrStartEnd(sa,"1",1,2000000);
    Set *rfSet = Slice_getAllRepeatFeatures(slice,"");
    int i;

    for (i=0;i<Set_getNumElement(rfSet);i++) {
      RepeatFeature *rf = Set_getElementAt(rfSet,i);
      printf("Repeat feature: %d-%d id " IDFMTSTR "\n", RepeatFeature_getStart(rf), 
             RepeatFeature_getEnd(rf), RepeatFeature_getDbID(rf));
    }
  }
  {
    SliceAdaptor *sa = DBAdaptor_getSliceAdaptor(dba);
    Slice *slice = SliceAdaptor_fetchByChrStartEnd(sa,"1",1,2000000);
    Set *sSet = Slice_getAllSimpleFeatures(slice,"",NULL);
    int i;

    for (i=0;i<Set_getNumElement(sSet);i++) {
      SimpleFeature *sf = Set_getElementAt(sSet,i);
      printf("Simple feature: %d-%d id " IDFMTSTR "\n", SimpleFeature_getStart(sf), 
             SimpleFeature_getEnd(sf), SimpleFeature_getDbID(sf));
    }
  }
  {
    RawContigAdaptor *rca = DBAdaptor_getRawContigAdaptor(dba);
    RawContig *contig = RawContigAdaptor_fetchByDbID(rca,10000);
    SimpleFeatureAdaptor *sfa = DBAdaptor_getSimpleFeatureAdaptor(dba);
    Set *sSet = RawContig_getAllSimpleFeatures(contig,"",NULL);
    int i;
    for (i=0;i<Set_getNumElement(sSet);i++) {
      SimpleFeature *sf = Set_getElementAt(sSet,i);
      printf("Simple feature: %d-%d id " IDFMTSTR "\n", SimpleFeature_getStart(sf), 
             SimpleFeature_getEnd(sf), SimpleFeature_getDbID(sf));
    }
    
    SimpleFeatureAdaptor_store((BaseFeatureAdaptor *)sfa,sSet);
  }
  {
    RawContigAdaptor *rca = DBAdaptor_getRawContigAdaptor(dba);
    RawContig *contig = RawContigAdaptor_fetchByDbID(rca,10000);
    PredictionTranscriptAdaptor *pta = DBAdaptor_getPredictionTranscriptAdaptor(dba);
    Set *sSet = RawContig_getAllPredictionTranscripts(contig,"");
    int i;
    for (i=0;i<Set_getNumElement(sSet);i++) {
      PredictionTranscript *pt = Set_getElementAt(sSet,i);
      printf("Prediction transcript: %d-%d id " IDFMTSTR "\n", PredictionTranscript_getStart(pt), 
             PredictionTranscript_getEnd(pt), PredictionTranscript_getDbID(pt));
    }
  }
  {
    RawContigAdaptor *rca = DBAdaptor_getRawContigAdaptor(dba);
    RawContig *contig = RawContigAdaptor_fetchByDbID(rca,10000);
    RepeatFeatureAdaptor *rfa = DBAdaptor_getRepeatFeatureAdaptor(dba);
    Set *sSet = RawContig_getAllRepeatFeatures(contig,"");
    int i;
    for (i=0;i<Set_getNumElement(sSet);i++) {
      RepeatFeature *rf = Set_getElementAt(sSet,i);
      printf("Prediction transcrirf: %d-%d id " IDFMTSTR "\n", RepeatFeature_getStart(rf), 
             RepeatFeature_getEnd(rf), RepeatFeature_getDbID(rf));
    }
  }
  {
    SliceAdaptor *sa = DBAdaptor_getSliceAdaptor(dba);
    Slice *slice = SliceAdaptor_fetchByChrStartEnd(sa,"1",1,2000000);
    Set *dafSet = Slice_getAllDNAAlignFeatures(slice,"",NULL);
    int i;

    for (i=0;i<Set_getNumElement(dafSet);i++) {
      DNAAlignFeature *daf = Set_getElementAt(dafSet,i);
      printf("DNA align feature: %d-%d id " IDFMTSTR "\n", DNAAlignFeature_getStart(daf), 
             DNAAlignFeature_getEnd(daf), DNAAlignFeature_getDbID(daf));
    }
  }
  {
    SliceAdaptor *sa = DBAdaptor_getSliceAdaptor(dba);
    Slice *slice = SliceAdaptor_fetchByChrStartEnd(sa,"1",1,2000000);
    Set *dpafSet = Slice_getAllDNAPepAlignFeatures(slice,"",NULL);
    int i;

    for (i=0;i<Set_getNumElement(dpafSet);i++) {
      DNAPepAlignFeature *dpaf = Set_getElementAt(dpafSet,i);
      printf("Pep align feature: %d-%d id " IDFMTSTR "\n", DNAPepAlignFeature_getStart(dpaf), 
             DNAPepAlignFeature_getEnd(dpaf), DNAPepAlignFeature_getDbID(dpaf));
    }
  }
  {
    GeneAdaptor *ga = DBAdaptor_getGeneAdaptor(dba);
    IDType *geneIds;
    int   nGeneId;
    int   i;
    Gene *gene;

    nGeneId = GeneAdaptor_listGeneIds(ga,&geneIds);

    printf("nGeneId = %d\n",nGeneId);

    for (i=1;i<5;i++) {
      gene = GeneAdaptor_fetchByDbID(ga,i,FALSE);
      if (gene) printf("Gene %s " IDFMTSTR "\n",Gene_getStableId(gene), Gene_getDbID(gene));
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

    printf("Transcript " IDFMTSTR " translation id = " IDFMTSTR "\n",Transcript_getDbID(t), Transcript_getTranslationId(t));
  }
  {
    ExonAdaptor *ea = DBAdaptor_getExonAdaptor(dba);
    Exon *e = ExonAdaptor_fetchByDbID(ea,1);
    Exon **exons;
    int nExon;
    int i;

    printf("Exon %d %d %d\n",Exon_getStart(e), Exon_getEnd(e), Exon_getStrand(e));

    nExon = ExonAdaptor_fetchAllByGeneId(ea,53LL,&exons);

    printf("NExon = %d\n",nExon);
    for (i=0;i<nExon;i++) {
      e = exons[i];
      printf("Exon " IDFMTSTR " %d %d %d\n",Exon_getDbID(e), Exon_getStart(e), Exon_getEnd(e), Exon_getStrand(e));
    }
  }
  {
    SliceAdaptor *sa = DBAdaptor_getSliceAdaptor(dba);
    Slice *slice = SliceAdaptor_fetchByChrStartEnd(sa,"1",1,230000000);
    Set *geneSet = Slice_getAllGenes(slice,NULL);
    int i;

    for (i=0;i<Set_getNumElement(geneSet);i++) {
      Gene *gene = Set_getElementAt(geneSet,i);
      int j;

      fprintf(stderr,"Gene %s\n",Gene_getStableId(gene));
      for (j=0;j<Gene_getTranscriptCount(gene);j++) {
        Transcript *trans = Gene_getTranscriptAt(gene,j);
        int k;

        fprintf(stderr," Transcript %s\n",Transcript_getStableId(trans));
        for (k=0;k<Transcript_getExonCount(trans);k++) {
          Exon *exon = (Exon *)Transcript_getExonAt(trans,k);
          fprintf(stderr,"  Exon %s %d %d %d\n",
                  Exon_getStableId(exon),
                  Exon_getStart(exon),
                  Exon_getEnd(exon),
                  Exon_getStrand(exon)
                 );
        }
      }
    }
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
    SliceAdaptor *sliceAd = DBAdaptor_getSliceAdaptor(dba);
    Slice *slice = SliceAdaptor_fetchByChrStartEnd(sliceAd,"1",9900000,10000000);
    RawContig *rc = RawContig_new();
    char *rawseq;
    char *sliceSeq;

    RawContig_setDbID(rc,24);
    rawseq = SequenceAdaptor_fetchByRawContigStartEndStrand(sa,rc,1,10,1);
    
    printf("Seq = %s\n",rawseq);

    sliceSeq = SequenceAdaptor_fetchBySliceStartEndStrand(sa,slice,1,100,1);
    printf("Slice Seq = %s\n",sliceSeq);
  }
}
