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
#include "PredictionTranscriptAdaptor.h"
#include "RawContigAdaptor.h"
#include "RepeatFeatureAdaptor.h"
#include "SequenceAdaptor.h"
#include "SimpleFeatureAdaptor.h"
#include "SliceAdaptor.h"
#include "TranscriptAdaptor.h"

#include "DBEntry.h"
#include "DNAAlignFeature.h"
#include "DNAPepAlignFeature.h"
#include "PredictionTranscript.h"
#include "RawContig.h"
#include "RepeatFeature.h"
#include "SimpleFeature.h"

int main(int argc, char *argv[]) {
  DBAdaptor *dba;

  initEnsC(argc, argv);

  printf("Opening connection ...");
  //dba = DBAdaptor_new("localhost","root",NULL,"test_ensembl",3306,NULL);
  //dba = DBAdaptor_new("kaka.sanger.ac.uk","anonymous",NULL,"homo_sapiens_core_12_31",3306,NULL);
  //dba = DBAdaptor_new("localhost","root",NULL,"homo_sapiens_core_18_34",3306,NULL);
  //dba = DBAdaptor_new("ecs2d.internal.sanger.ac.uk","ensro",NULL,"homo_sapiens_core_14_31",3306,NULL);
  //dba = DBAdaptor_new("ecs2d.internal.sanger.ac.uk","ensro",NULL,"rattus_norvegicus_core_15_2",3306,NULL);
  dba = DBAdaptor_new("ens-staging.internal.sanger.ac.uk","ensro",NULL,"homo_sapiens_core_70_37",3306,NULL);
  printf(" Done\n");

  printf("Assembly type %s\n",DBAdaptor_getAssemblyType(dba));

/*
  {
    SimpleFeatureAdaptor *sfa = DBAdaptor_getSimpleFeatureAdaptor(dba);
    SimpleFeature *sf = (SimpleFeature *)SimpleFeatureAdaptor_fetchByDbID(sfa,1);
    printf("Simple feature: %d-%d id " IDFMTSTR "\n", SimpleFeature_getStart(sf), 
           SimpleFeature_getEnd(sf), SimpleFeature_getDbID(sf));
  }
  {
    SliceAdaptor *sa = DBAdaptor_getSliceAdaptor(dba);
    Slice *slice = SliceAdaptor_fetchByChrStartEnd(sa,"1",1,2000000);
    Vector *ptVector = Slice_getAllPredictionTranscripts(slice,"");
    int i;

    for (i=0;i<Vector_getNumElement(ptVector);i++) {
      PredictionTranscript *pt = Vector_getElementAt(ptVector,i);
      printf("PredictionTranscript feature: %s %d-%d id " IDFMTSTR "\n", 
             PredictionTranscript_getStableId(pt), 
             PredictionTranscript_getStart(pt), PredictionTranscript_getEnd(pt), 
             PredictionTranscript_getDbID(pt));
    }
  }
  {
    SliceAdaptor *sa = DBAdaptor_getSliceAdaptor(dba);
    Slice *slice = SliceAdaptor_fetchByChrStartEnd(sa,"1",1,2000000);
    Vector *rfVector = Slice_getAllRepeatFeatures(slice,"");
    int i;

    for (i=0;i<Vector_getNumElement(rfVector);i++) {
      RepeatFeature *rf = Vector_getElementAt(rfVector,i);
      printf("Repeat feature: %d-%d id " IDFMTSTR "\n", RepeatFeature_getStart(rf), 
             RepeatFeature_getEnd(rf), RepeatFeature_getDbID(rf));
    }
  }
  {
    SliceAdaptor *sa = DBAdaptor_getSliceAdaptor(dba);
    Slice *slice = SliceAdaptor_fetchByChrStartEnd(sa,"1",1,2000000);
    Vector *sVector = Slice_getAllSimpleFeatures(slice,"",NULL);
    int i;

    for (i=0;i<Vector_getNumElement(sVector);i++) {
      SimpleFeature *sf = Vector_getElementAt(sVector,i);
      printf("Simple feature: %d-%d id " IDFMTSTR "\n", SimpleFeature_getStart(sf), 
             SimpleFeature_getEnd(sf), SimpleFeature_getDbID(sf));
    }
  }
  {
    RawContigAdaptor *rca = DBAdaptor_getRawContigAdaptor(dba);
    RawContig *contig = RawContigAdaptor_fetchByDbID(rca,10000);
    SimpleFeatureAdaptor *sfa = DBAdaptor_getSimpleFeatureAdaptor(dba);
    Vector *sVector = RawContig_getAllSimpleFeatures(contig,"",NULL);
    int i;
    for (i=0;i<Vector_getNumElement(sVector);i++) {
      SimpleFeature *sf = Vector_getElementAt(sVector,i);
      printf("Simple feature: %d-%d id " IDFMTSTR "\n", SimpleFeature_getStart(sf), 
             SimpleFeature_getEnd(sf), SimpleFeature_getDbID(sf));
    }
    
    SimpleFeatureAdaptor_store((BaseFeatureAdaptor *)sfa,sVector);
  }
  {
    RawContigAdaptor *rca = DBAdaptor_getRawContigAdaptor(dba);
    RawContig *contig = RawContigAdaptor_fetchByDbID(rca,10000);
    PredictionTranscriptAdaptor *pta = DBAdaptor_getPredictionTranscriptAdaptor(dba);
    Vector *sVector = RawContig_getAllPredictionTranscripts(contig,"");
    int i;
    for (i=0;i<Vector_getNumElement(sVector);i++) {
      PredictionTranscript *pt = Vector_getElementAt(sVector,i);
      printf("Prediction transcript: %d-%d id " IDFMTSTR "\n", PredictionTranscript_getStart(pt), 
             PredictionTranscript_getEnd(pt), PredictionTranscript_getDbID(pt));
    }
    PredictionTranscriptAdaptor_store((BaseFeatureAdaptor *)pta,sVector);
  }
  {
    RawContigAdaptor *rca = DBAdaptor_getRawContigAdaptor(dba);
    RawContig *contig = RawContigAdaptor_fetchByDbID(rca,10000);
    RepeatFeatureAdaptor *rfa = DBAdaptor_getRepeatFeatureAdaptor(dba);
    Vector *sVector = RawContig_getAllRepeatFeatures(contig,"");
    int i;
    for (i=0;i<Vector_getNumElement(sVector);i++) {
      RepeatFeature *rf = Vector_getElementAt(sVector,i);
      printf("Repeat feature: %d-%d id " IDFMTSTR "\n", RepeatFeature_getStart(rf), 
             RepeatFeature_getEnd(rf), RepeatFeature_getDbID(rf));
    }
    RepeatFeatureAdaptor_store((BaseFeatureAdaptor *)rfa,sVector);
  }
  {
    RawContigAdaptor *rca = DBAdaptor_getRawContigAdaptor(dba);
    RawContig *contig = RawContigAdaptor_fetchByDbID(rca,10000);
    DNAAlignFeatureAdaptor *dafa = DBAdaptor_getDNAAlignFeatureAdaptor(dba);
    Vector *sVector = RawContig_getAllDNAAlignFeatures(contig,"",NULL);
    int i;
    for (i=0;i<Vector_getNumElement(sVector);i++) {
      DNAAlignFeature *daf = Vector_getElementAt(sVector,i);
      printf("DNA align feature: %d-%d id " IDFMTSTR "\n", DNAAlignFeature_getStart(daf), 
             DNAAlignFeature_getEnd(daf), DNAAlignFeature_getDbID(daf));
    }
    DNAAlignFeatureAdaptor_store((BaseFeatureAdaptor *)dafa,sVector);
  }
  {
    RawContigAdaptor *rca = DBAdaptor_getRawContigAdaptor(dba);
    RawContig *contig = RawContigAdaptor_fetchByDbID(rca,10000);
    ProteinAlignFeatureAdaptor *pafa = DBAdaptor_getProteinAlignFeatureAdaptor(dba);
    Vector *sVector = RawContig_getAllProteinAlignFeatures(contig,"",NULL);
    int i;
    for (i=0;i<Vector_getNumElement(sVector);i++) {
      DNAPepAlignFeature *paf = Vector_getElementAt(sVector,i);
      printf("Protein align feature: %d-%d id " IDFMTSTR "\n", DNAPepAlignFeature_getStart(paf), 
             DNAPepAlignFeature_getEnd(paf), DNAPepAlignFeature_getDbID(paf));
    }
    ProteinAlignFeatureAdaptor_store((BaseFeatureAdaptor *)pafa,sVector);
  }
  {
    SliceAdaptor *sa = DBAdaptor_getSliceAdaptor(dba);
    Slice *slice = SliceAdaptor_fetchByChrStartEnd(sa,"1",1,2000000);
    Vector *dafVector = Slice_getAllDNAAlignFeatures(slice,"",NULL);
    int i;

    for (i=0;i<Vector_getNumElement(dafVector);i++) {
      DNAAlignFeature *daf = Vector_getElementAt(dafVector,i);
      printf("DNA align feature: %d-%d id " IDFMTSTR "\n", DNAAlignFeature_getStart(daf), 
             DNAAlignFeature_getEnd(daf), DNAAlignFeature_getDbID(daf));
    }
  }
  {
    SliceAdaptor *sa = DBAdaptor_getSliceAdaptor(dba);
    Slice *slice = SliceAdaptor_fetchByChrStartEnd(sa,"1",1,2000000);
    Vector *dpafVector = Slice_getAllDNAPepAlignFeatures(slice,"",NULL);
    int i;

    for (i=0;i<Vector_getNumElement(dpafVector);i++) {
      DNAPepAlignFeature *dpaf = Vector_getElementAt(dpafVector,i);
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
      if (gene) {
        Vector *dblinks = Gene_getAllDBLinks(gene);
        int j;
        printf("Gene %s " IDFMTSTR "\n",Gene_getStableId(gene), Gene_getDbID(gene));
        for (j=0;j<Vector_getNumElement(dblinks);j++) {
          DBEntry *dbe = Vector_getElementAt(dblinks,j);
          printf("  dblink = %s\n",DBEntry_getDisplayId(dbe));
        }
      }
    }
  }
  {
    AnalysisAdaptor *aa = DBAdaptor_getAnalysisAdaptor(dba);

    Analysis *anal = AnalysisAdaptor_fetchByDbID(aa,2);

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
    Vector *geneVector = Slice_getAllGenes(slice,NULL);
    int i;

    for (i=0;i<Vector_getNumElement(geneVector);i++) {
      Gene *gene = Vector_getElementAt(geneVector,i);
      int j;

      fprintf(stderr,"Gene %s from %d to %d\n",
              Gene_getStableId(gene),
              Gene_getStart(gene),
              Gene_getEnd(gene));
      for (j=0;j<Gene_getTranscriptCount(gene);j++) {
        Transcript *trans = Gene_getTranscriptAt(gene,j);
        int k;

        Transcript_sort(trans);

        fprintf(stderr," Transcript %s from %d to %d\n",
                Transcript_getStableId(trans),
                Transcript_getStart(trans),
                Transcript_getEnd(trans));
        if (Transcript_getTranslation(trans)) {
          fprintf(stderr,"translation = %s\n",Transcript_translate(trans));
        } else {
          fprintf(stderr,"No translation\n");
        }
*/
/*
        for (k=0;k<Transcript_getExonCount(trans);k++) {
          Exon *exon = (Exon *)Transcript_getExonAt(trans,k);
          fprintf(stderr,"  Exon %s %d %d %d %p\n",
                  Exon_getStableId(exon),
                  Exon_getStart(exon),
                  Exon_getEnd(exon),
                  Exon_getStrand(exon),
                  Exon_getContig(exon)
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

    RawContig_setDbID(rc,10000);
    rawseq = SequenceAdaptor_fetchByRawContigStartEndStrand(sa,rc,1,10,1);
    
    printf("Seq = %s\n",rawseq);

    sliceSeq = SequenceAdaptor_fetchBySliceStartEndStrand(sa,slice,1,100,1);
    printf("Slice Seq = %s\n",sliceSeq);
  }
*/
  return 0;
}
