#include "Slice.h"
#include "Gene.h"
#include "GeneAdaptor.h"
#include "SliceAdaptor.h"
#include "ChromosomeAdaptor.h"
#include "Set.h"

Slice *Slice_new(char *chr, int start, int end, int strand, char *assemblyType,
                 SliceAdaptor *sa, long dbID, int empty) {
  Slice *slice;

  if ((slice = (Slice *)calloc(1,sizeof(Slice))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for slice\n");
    return NULL;
  }

  slice->contigType = SLICE;

  if (!empty) {
    if( !chr || !assemblyType) {
      fprintf(stderr,"ERROR: Do not have all the parameters for slice\n");
      exit(1);
    }
    Slice_setChrName(slice,chr);
    Slice_setChrStart(slice,start);
    Slice_setChrEnd(slice,end);
    Slice_setStrand(slice,strand);
  } else {
    Slice_setStrand(slice,1);
    Slice_setChrStart(slice,1);
    Slice_setEmptyFlag(slice,TRUE);

    // empty Slices are used to do mapping to chromosomal coords.
    // After the mapping, the Slice contains chr_name and is reference
    // point for the mapped object
  }

  Slice_setAssemblyType(slice,assemblyType);
  Slice_setAdaptor(slice,(BaseAdaptor *)sa);
  Slice_setDbID(slice,dbID);

  if (sa && chr) {
    ChromosomeAdaptor *ca = DBAdaptor_getChromosomeAdaptor(sa->dba);
    Chromosome *chromosome = ChromosomeAdaptor_fetchByChrName(ca, chr);
    Slice_setChrId(slice,Chromosome_getDbID(chromosome));
  }

/*
  if( defined $adaptor && !defined $type ) {
    $self->assembly_type
      ( $adaptor->db()->get_MetaContainer()->get_default_assembly());
  }
*/
  return slice;
}

char *Slice_getName(Slice *slice, char *retStr) {
  sprintf(retStr,"%s.%d-%d",Slice_getChrName(slice),
          Slice_getChrStart(slice),Slice_getChrEnd(slice));
  return retStr;
}

/* logicName is optional - set to NULL if you want them all */
Set *Slice_getAllGenes(Slice *slice, char *logicName) {
  GeneAdaptor *ga = DBAdaptor_getGeneAdaptor(Slice_getAdaptor(slice)->dba);

  return GeneAdaptor_fetchAllBySlice(ga,slice,logicName);
}

char *Slice_setChrName(Slice *sl, char *chrName) {
  if ((sl->chrName = (char *)malloc(strlen(chrName)+1)) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for chrName\n");
    exit(1);
  }

  strcpy(sl->chrName,chrName);

  return sl->chrName;
}

char *Slice_setAssemblyType(Slice *sl, char *assemblyType) {
  if ((sl->assemblyType = (char *)malloc(strlen(assemblyType)+1)) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for assemblyType\n");
    exit(1);
  }

  strcpy(sl->assemblyType,assemblyType);

  return sl->assemblyType;
}
