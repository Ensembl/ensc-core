#include <stdio.h>

#include "SliceAdaptor.h"
#include "DBAdaptor.h"
#include "EnsC.h"
#include "Gene.h"

#include "BaseRODBTest.h"
#include "gperftools/tcmalloc.h"

int main(int argc, char *argv[]) {
  DBAdaptor *dba;
  GeneAdaptor *ga;
  Slice *slice;
  Vector *genes;
  int i;
  int failed;
  
  initEnsC();

  dba = Test_initROEnsDB();

  slice = Test_getStandardSlice(dba);

  ok(1, slice!=NULL);

  ga = DBAdaptor_getGeneAdaptor(dba);

  ok(2, ga!=NULL);

  genes =  Slice_getAllGenes(slice, NULL, NULL, 1, NULL, NULL);

  fprintf(stdout, "Have %d genes\n", Vector_getNumElement(genes));
  ok(3, genes!=NULL);
  ok(4, Vector_getNumElement(genes)!=0);

  failed = 0;
  for (i=0;i<Vector_getNumElement(genes) && !failed;i++) {
    Gene *g = Vector_getElementAt(genes,i);
    fprintf(stderr,"Gene %s coords: %ld %ld %d\n",Gene_getStableId(g),Gene_getStart(g),Gene_getEnd(g),Gene_getStrand(g));

    int j;
    for (j=0;j<Gene_getTranscriptCount(g);j++) {
      Transcript *t = Gene_getTranscriptAt(g,j);
     
      fprintf(stderr," Trans %s coords: %ld %ld %d\n",Transcript_getStableId(t), Transcript_getStart(t),Transcript_getEnd(t),Transcript_getStrand(t));
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
        fprintf(stderr," translation: %s\n",Transcript_translate(t));
      }
    }
  }
  ok(5, !failed);
  tc_malloc_stats();

  EcoString_getInfo(ecoSTable);

  return 0;
}
