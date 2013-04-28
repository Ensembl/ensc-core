#include <stdio.h>

#include "SliceAdaptor.h"
#include "DBAdaptor.h"
#include "EnsC.h"

#include "BaseTest.h"

int main(int argc, char *argv[]) {
  DBAdaptor *dba;
  Slice *slice;
  SliceAdaptor *sa;

  initEnsC();

  //dba = DBAdaptor_new("ensembldb.ensembl.org","anonymous",NULL,"homo_sapiens_core_70_37",5306,NULL);
  dba = DBAdaptor_new("ens-livemirror.internal.sanger.ac.uk","ensro",NULL,"homo_sapiens_core_70_37",3306,NULL);

  ok(1,!strcmp("GRCh37",DBAdaptor_getAssemblyType(dba)));

  sa = DBAdaptor_getSliceAdaptor(dba);

  ok(2, sa!=NULL);

  slice = SliceAdaptor_fetchByRegion(sa,"chromosome","1",1,250000000,1,NULL,0);

  ok(3, slice!=NULL);

  printf("slice name = %s\n",Slice_getName(slice));
  ok(4,!strcmp(Slice_getName(slice),"1.1-250000000"));

  char *seq = Slice_getSeq(slice);
//  printf("slice seq = %s\n", seq);

  
  SeqUtil_writeFasta(stdout, Slice_getName(slice), seq, 60);

  fflush(stdout);


  return 0;
}
