#include <stdio.h>

#include "SliceAdaptor.h"
#include "DBAdaptor.h"
#include "EnsC.h"
#include "Gene.h"

#include "BaseRODBTest.h"

#include "gperftools/tcmalloc.h"
#include "ProcUtil.h"
#define UNW_LOCAL_ONLY
#include "libunwind.h"

int main(int argc, char *argv[]) {
  DBAdaptor *dba;
  GeneAdaptor *ga;
  Slice *slice;
  Vector *genes;
  int i;
  int failed;
  
  initEnsC(argc, argv);

  ProcUtil_showBacktrace(EnsC_progName);

  dba = Test_initROEnsDB();

  slice = Test_getStandardSlice(dba);

  ok(1, slice!=NULL);

  ga = DBAdaptor_getGeneAdaptor(dba);

  ok(2, ga!=NULL);

  genes =  Slice_getAllGenes(slice, NULL, NULL, 1, NULL, NULL);

  fprintf(stdout, "Have %d genes\n", Vector_getNumElement(genes));
  ok(3, genes!=NULL);
  ok(4, Vector_getNumElement(genes)!=0);

/*
  failed = dumpGenes(genes);
  ok(5, !failed);
*/

  SliceAdaptor *sa = DBAdaptor_getSliceAdaptor(dba);
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
  int i;
  int failed = 0;
  for (i=0;i<Vector_getNumElement(genes) && !failed;i++) {
    Gene *g = Vector_getElementAt(genes,i);
    fprintf(stderr,"Gene %s (%s) coords: %ld %ld %d\n",Gene_getStableId(g),(Gene_getDisplayXref(g) ? DBEntry_getDisplayId(Gene_getDisplayXref(g)) : ""),Gene_getStart(g),Gene_getEnd(g),Gene_getStrand(g));

    int j;
    for (j=0;j<Gene_getTranscriptCount(g);j++) {
      Transcript *t = Gene_getTranscriptAt(g,j);
     
      fprintf(stderr," Trans %s coords: %ld %ld %d biotype: %s\n",Transcript_getStableId(t), Transcript_getStart(t),Transcript_getEnd(t),Transcript_getStrand(t),Transcript_getBiotype(t));
      int k;
      for (k=0;k<Transcript_getExonCount(t);k++) {
        Exon *e = Transcript_getExonAt(t,k);
        fprintf(stderr,"  exon %s coords: %ld %ld %d\n",Exon_getStableId(e), Exon_getStart(e),Exon_getEnd(e),Exon_getStrand(e));
      }
      Translation *tln = Transcript_getTranslation(t);
      if (tln) {
 
        fprintf(stderr," translation id: %s %s %d %s %d\n",Translation_getStableId(tln), 
                Exon_getStableId(Translation_getStartExon(tln)), Translation_getStart(tln),
                Exon_getStableId(Translation_getEndExon(tln)), Translation_getEnd(tln));
        char *tSeq = Transcript_translate(t);
        fprintf(stderr," translation: %s\n",tSeq);
        free(tSeq);
      }
    }
  }
  return failed;
}

