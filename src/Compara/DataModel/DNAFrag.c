#include "DNAFrag.h"
#include "StrUtil.h"
#include "BaseContig.h"
#include "DBAdaptor.h"
#include "RawContigAdaptor.h"
#include "SliceAdaptor.h"

DNAFrag *DNAFrag_new() {
  DNAFrag *df;

  if ((df = (DNAFrag *)calloc(1,sizeof(DNAFrag))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for df\n");
    return NULL;
  }

  df->objectType = CLASS_DNAFRAG;
  return df;
}

char *DNAFrag_setName(DNAFrag *df, char *name) {
  StrUtil_copyString(&(df->name),name,0);

  return df->name;
}

char *DNAFrag_setType(DNAFrag *df, char *type) {
  StrUtil_copyString(&(df->type),type,0);

  return df->type;
}

BaseContig *DNAFrag_getContig(DNAFrag *df) {

   if (!df->contig) {
     DBAdaptor *dba = GenomeDB_getDBAdaptor(DNAFrag_getGenomeDB(df));
     if (!strcmp(DNAFrag_getType(df),"RawContig")) {
       RawContigAdaptor *rca = DBAdaptor_getRawContigAdaptor(dba);
       df->contig = (BaseContig *)RawContigAdaptor_fetchByName(rca, DNAFrag_getName(df));
     } else if (!strcmp(DNAFrag_getType(df),"VirtualContig")) {
       SliceAdaptor *sa = DBAdaptor_getSliceAdaptor(dba);
       char chrName[1024];
       int start, end;
       fprintf(stderr, "VC type not implemented \n");
       exit(1);
       //my ($chr,$start,$end) = split /\./, $self->name;
       //df->contig = $core_dbadaptor->get_SliceAdaptor->fetch_by_chr_start_end(chrName,start,end);
     } else if (!strcmp(DNAFrag_getType(df),"Chromosome")) {
       SliceAdaptor *sa = DBAdaptor_getSliceAdaptor(dba);
       df->contig = (BaseContig *)SliceAdaptor_fetchByChrName(sa, DNAFrag_getName(df));
     } else {
       fprintf(stderr, "Error: Can't fetch contig of %s with type %s\n",
               DNAFrag_getName(df), DNAFrag_getType(df));
     }
   }

   return df->contig;
}

