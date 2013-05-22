#include <stdio.h>

#include "SliceAdaptor.h"
#include "DBAdaptor.h"
#include "EnsC.h"
#include "Gene.h"
#include "Attribute.h"
#include "BaseAlignFeature.h"

#include "BaseRODBTest.h"

#include "gperftools/tcmalloc.h"
#include "ProcUtil.h"
#define UNW_LOCAL_ONLY
#include "libunwind.h"

int main(int argc, char *argv[]) {
  DBAdaptor *dba;
  GeneAdaptor *ga;
  Slice *slice = NULL;
  Vector *genes = NULL;
  int i = 0;
  int failed = 0;
  
  initEnsC(argc, argv);

//  ProcUtil_showBacktrace(EnsC_progName);

//  dba = Test_initROEnsDB();
//  slice = Test_getStandardSlice(dba);

  DBAdaptor *seqdba = DBAdaptor_new("genebuild6.internal.sanger.ac.uk","ensadmin","ensembl","steve_chicken_rnaseq_missing_reference",3306,NULL);
  dba = DBAdaptor_new("genebuild1.internal.sanger.ac.uk","ensadmin","ensembl","steve_chicken_rnaseq_missing_refined",3306,seqdba);

  ok(1, slice!=NULL);

  ga = DBAdaptor_getGeneAdaptor(dba);
  SliceAdaptor *sa = DBAdaptor_getSliceAdaptor(dba);

  ok(2, ga!=NULL);

//  slice = SliceAdaptor_fetchByRegion(sa,"chromosome","17",1000000,5000000,1,NULL,0);
// Has a seleno
  slice = SliceAdaptor_fetchByRegion(sa,"chromosome","1",1000000,27000000,1,NULL,0);
//  slice = SliceAdaptor_fetchByRegion(sa,"chromosome","MT",1,17000,1,NULL,0);
  genes =  Slice_getAllGenes(slice, NULL, NULL, 1, NULL, NULL);

  fprintf(stdout, "Have %d genes\n", Vector_getNumElement(genes));
  ok(3, genes!=NULL);
  ok(4, Vector_getNumElement(genes)!=0);

  failed = dumpGenes(genes);
  ok(5, !failed);

  exit(1);
  Vector *toplevelSlices = SliceAdaptor_fetchAll(sa, "toplevel", NULL, 0);

  for (i=0;i<Vector_getNumElement(toplevelSlices) && !failed;i++) {
    Slice *tlSlice = Vector_getElementAt(toplevelSlices, i);
    genes =  Slice_getAllGenes(tlSlice, NULL, NULL, 1, NULL, NULL);
    fprintf(stderr, "Got %d genes on %s\n", Vector_getNumElement(genes), Slice_getName(tlSlice));
    failed = dumpGenes(genes);
  }


  tc_malloc_stats();


  fprintf(stderr,"\nEcostring table stats:\n");
  EcoString_getInfo(ecoSTable);

  fprintf(stderr,"\n");
  ProcUtil_timeInfo("at end of GeneTest");

  return 0;
}

int dumpGenes(Vector *genes) {
  FILE *fp = stderr;
  int i;
  int failed = 0;
  for (i=0;i<Vector_getNumElement(genes) && !failed;i++) {
    Gene *g = Vector_getElementAt(genes,i);
    fprintf(fp,"Gene %s (%s) coords: %ld %ld %d\n",Gene_getStableId(g),(Gene_getDisplayXref(g) ? DBEntry_getDisplayId(Gene_getDisplayXref(g)) : ""),Gene_getStart(g),Gene_getEnd(g),Gene_getStrand(g));

    int j;
    for (j=0;j<Gene_getTranscriptCount(g);j++) {
      Transcript *t = Gene_getTranscriptAt(g,j);
     
      fprintf(fp," Trans %s coords: %ld %ld %d biotype: %s\n",Transcript_getStableId(t), Transcript_getStart(t),Transcript_getEnd(t),Transcript_getStrand(t),Transcript_getBiotype(t));
      Vector *support = Transcript_getAllSupportingFeatures(t);
      int k;
      for (k=0; k<Vector_getNumElement(support); k++) {
        BaseAlignFeature *baf = Vector_getElementAt(support, k);
        fprintf(fp,"   support %s coords: %ld %ld %d\n", BaseAlignFeature_getHitSeqName(baf), BaseAlignFeature_getStart(baf), BaseAlignFeature_getEnd(baf), BaseAlignFeature_getStrand(baf));
      }
      Vector *intronSupport = Transcript_getAllIntronSupportingEvidence(t);
      for (k=0; k<Vector_getNumElement(intronSupport); k++) {
        IntronSupportingEvidence *ise = Vector_getElementAt(intronSupport, k);
        fprintf(fp,"   intron support %s coords: %ld %ld %d\n", IntronSupportingEvidence_getHitName(ise), IntronSupportingEvidence_getStart(ise), IntronSupportingEvidence_getEnd(ise), IntronSupportingEvidence_getStrand(ise));
      }

      for (k=0;k<Transcript_getExonCount(t);k++) {
        Exon *e = Transcript_getExonAt(t,k);
        fprintf(fp,"  exon %s coords: %ld %ld %d\n",Exon_getStableId(e), Exon_getStart(e),Exon_getEnd(e),Exon_getStrand(e));
        Vector *support = Exon_getAllSupportingFeatures(e);
        int m;
        for (m=0; m<Vector_getNumElement(support); m++) {
          BaseAlignFeature *baf = Vector_getElementAt(support, m);
          fprintf(fp,"   support %s coords: %ld %ld %d\n", BaseAlignFeature_getHitSeqName(baf), BaseAlignFeature_getStart(baf), BaseAlignFeature_getEnd(baf), BaseAlignFeature_getStrand(baf));
        }
      }
      Translation *tln = Transcript_getTranslation(t);
      if (tln) {
 
        fprintf(fp," translation id: %s %s %d %s %d\n",Translation_getStableId(tln), 
                Exon_getStableId(Translation_getStartExon(tln)), Translation_getStart(tln),
                Exon_getStableId(Translation_getEndExon(tln)), Translation_getEnd(tln));
        char *tSeq = Transcript_translate(t);
        fprintf(fp," translation: %s\n",tSeq);
        free(tSeq);
        Vector *tlnAttribs = Translation_getAllAttributes(tln, NULL);
        if (Vector_getNumElement(tlnAttribs)) {
          fprintf(fp, " translation attributes:\n");
          int n;
          for (n=0; n<Vector_getNumElement(tlnAttribs); n++) {
            Attribute *attrib = Vector_getElementAt(tlnAttribs, n);
            fprintf(fp, "  code %s name %s desc %s value %s\n", 
                    Attribute_getCode(attrib), 
                    Attribute_getName(attrib),
                    Attribute_getDescription(attrib),
                    Attribute_getValue(attrib));
          }
        }
      }
    }
  }
  return failed;
}

