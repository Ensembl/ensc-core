#define __EXON_MAIN__
#include "Exon.h"
#undef __EXON_MAIN__
#include "ExonAdaptor.h"
#include "AssemblyMapperAdaptor.h"
#include "AssemblyMapper.h"
#include "DBAdaptor.h"
#include "SliceAdaptor.h"
#include "ChromosomeAdaptor.h"
#include "StickyExon.h"

Exon *Exon_new() {
  Exon *exon;

  if ((exon = (Exon *)calloc(1,sizeof(Exon))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for exon\n");
    return NULL;
  }

/* Set to empty values */
  Exon_setModified(exon,0);
  Exon_setCreated(exon,0);
  Exon_setVersion(exon,-1);

  exon->objectType = CLASS_EXON;

  exon->funcs = &exonFuncs;

  return exon;
}

char *Exon_getStableId(Exon *exon) {
  ExonAdaptor *ea = (ExonAdaptor *)Exon_getAdaptor(exon);

  if (StableIdInfo_getStableId(&(exon->si)) == NULL && ea) {
    ExonAdaptor_getStableEntryInfo(ea,exon);
  }
  return StableIdInfo_getStableId(&(exon->si));
}

time_t Exon_getCreated(Exon *exon) {
  ExonAdaptor *ea = (ExonAdaptor *)Exon_getAdaptor(exon);

  if (StableIdInfo_getCreated(&(exon->si)) == 0 && ea) {
    ExonAdaptor_getStableEntryInfo(ea,exon);
  }
  return StableIdInfo_getCreated(&(exon->si));
}

time_t Exon_getModified(Exon *exon) {
  ExonAdaptor *ea = (ExonAdaptor *)Exon_getAdaptor(exon);

  if (StableIdInfo_getModified(&(exon->si)) == 0 && ea) {
    ExonAdaptor_getStableEntryInfo(ea,exon);
  }
  return StableIdInfo_getModified(&(exon->si));
}

int Exon_getVersion(Exon *exon) {
  ExonAdaptor *ea = (ExonAdaptor *)Exon_getAdaptor(exon);

  if (StableIdInfo_getVersion(&(exon->si)) == -1 && ea) {
    ExonAdaptor_getStableEntryInfo(ea,exon);
  }
  return StableIdInfo_getVersion(&(exon->si));
}

int ExonStickyRankCompFunc(const void *a, const void *b) {
  Exon **e1 = (Exon **)a;
  Exon **e2 = (Exon **)b;

  if (Exon_getStickyRank(*e1) > Exon_getStickyRank(*e2)) {
    return 1;
  } else if (Exon_getStickyRank(*e1) < Exon_getStickyRank(*e2)) {
    return -1;
  } else {
    return 0;
  }
}

int Exon_forwardStrandCompFunc(const void *a, const void *b) {
  Exon **e1 = (Exon **)a;
  Exon **e2 = (Exon **)b;

  if (Exon_getStart(*e1) > Exon_getStart(*e2)) {
    return 1;
  } else if (Exon_getStart(*e1) < Exon_getStart(*e2)) {
    return -1;
  } else {
    return 0;
  }
}

int Exon_reverseStrandCompFunc(const void *a, const void *b) {
  Exon **e1 = (Exon **)a;
  Exon **e2 = (Exon **)b;

  if (Exon_getStart(*e1) < Exon_getStart(*e2)) {
    return 1;
  } else if (Exon_getStart(*e1) > Exon_getStart(*e2)) {
    return -1;
  } else {
    return 0;
  }
}

void Exon_sortByStickyRank(Exon *exon) {
  if (Exon_isSticky(exon)) {
    qsort(Exon_getComponents(exon), Exon_getNumComponentExon(exon), sizeof(void *), 
          ExonStickyRankCompFunc);
  }
  return; 
}

Exon *Exon_transformToSlice(Exon *exon, Slice *slice) {
  Exon *newExon;
  BaseAdaptor *adaptor;
  AssemblyMapperAdaptor *ama;
  AssemblyMapper *assMapper;
  MapperRangeSet *mapped;
  MapperCoordinate *mc;

// HACK HACK HACK
  if (Exon_isSticky(exon)) {
    return StickyExon_transformToSlice(exon,slice);
  }

  if (!Exon_getContig(exon)) {
    fprintf(stderr,"ERROR: Exon's contig must be defined to transform to Slice coords");
    exit(1);
  }
  // fprintf(stderr,"transforming exon %d from raw contig to slice coords\n",Exon_getDbID(exon));
  // fprintf(stderr,"exon %s\n",Exon_getStableId(exon));

  adaptor = Slice_getAdaptor(slice);
  if (!adaptor) {
    adaptor = BaseContig_getAdaptor((BaseContig *)Exon_getContig(exon));
  }
  
  if (!adaptor) {
    fprintf(stderr, "ERROR: Cannot transform to exon slice unless either the " 
		    "exon->contig->adaptor or slice->adaptor is defined");
    exit(1);
  }

  ama = DBAdaptor_getAssemblyMapperAdaptor(adaptor->dba);
  assMapper = AssemblyMapperAdaptor_fetchByType(ama, Slice_getAssemblyType(slice));
  
  mapped = AssemblyMapper_mapCoordinatesToAssembly(assMapper,
     BaseContig_getDbID(Exon_getContig(exon)),
     Exon_getStart(exon),
     Exon_getEnd(exon),
     Exon_getStrand(exon)
    );

  // exons should always transform so in theory no error check necessary
  // actually we could have exons inside and outside the Slice 
  // because of db design and the query that produces them
  if( !mapped || mapped->nRange == 0) {
    fprintf(stderr, "ERROR: Exon couldnt map dbID = " IDFMTSTR "\n",Exon_getDbID(exon));
    exit(1);
  }

  // should get a gap object returned if an exon lies outside of 
  // the current slice.  Simply return the exon as is - i.e. untransformed.
  // this untransformed exon will be distinguishable as it will still have
  // contig attached to it and not a slice.
  if( MapperRangeSet_getRangeAt(mapped,0)->rangeType == MAPPERRANGE_GAP) {
    fprintf(stderr,"Exon in gap\n");
    return exon;
  }

  mc = (MapperCoordinate *)MapperRangeSet_getRangeAt(mapped,0);

  // fprintf(stderr," Exon %s mapped range = %d %d\n",Exon_getStableId(exon),mc->start,mc->end);

  // the slice is an empty slice, create an enitre chromosome slice and
  // replace the empty slice with it
  if (Slice_getEmptyFlag(slice)) {
    SliceAdaptor *sa = DBAdaptor_getSliceAdaptor(adaptor->dba);
    ChromosomeAdaptor *ca = DBAdaptor_getChromosomeAdaptor(adaptor->dba);
    char *chrName = Chromosome_getName(ChromosomeAdaptor_fetchByDbID(ca,mc->id));

    slice = SliceAdaptor_fetchByChrName(sa, chrName);
// NIY free old slice (or have special empty one)???
  } 

  newExon = Exon_copy(exon, SHALLOW_DEPTH);

  if (Slice_getStrand(slice) == 1) {
    Exon_setStart(newExon,mc->start - Slice_getChrStart(slice) + 1);
    Exon_setEnd(newExon,mc->end - Slice_getChrStart(slice) + 1);
  } else {
    Exon_setStart(newExon,Slice_getChrEnd(slice) - mc->end + 1);
    Exon_setEnd(newExon,Slice_getChrEnd(slice) - mc->start + 1);
  }

  Exon_setStrand(newExon, mc->strand * Slice_getStrand(slice));
  Exon_setContig(newExon, slice);

  //copy the attached supporting features and transform them
#ifdef DONE
// NIY
  my @feats;
  if( exists $self->{_supporting_evidence} ) {
    foreach my $sf (@{$self->get_all_supporting_features()}) {
      push @feats, $sf->transform($slice);
    }
    $newexon->add_supporting_features(@feats);
  }
#endif

  // NIY free old exon (reference counting needed ???)

  // fprintf(stderr,"transformed exon %s %d %d\n",
  //        Exon_getStableId(newExon),
  //        Exon_getStart(newExon),
  //        Exon_getEnd(newExon)
  //       );
  return newExon;
}

Exon *Exon_copy(Exon *orig, CopyDepth depth) {
  Exon *newEx = Exon_new();

  if (depth != SHALLOW_DEPTH) {
    fprintf(stderr, "ERROR: Only SHALLOW copy implemented in Exon\n");
    exit(1);
  }
  Exon_setStart(newEx,Exon_getStart(orig));  
  Exon_setEnd(newEx,Exon_getEnd(orig));  
  Exon_setDbID(newEx,Exon_getEnd(orig));  
  Exon_setAdaptor(newEx,Exon_getAdaptor(orig));  
  Exon_setStrand(newEx,Exon_getStrand(orig));  
  Exon_setPhase(newEx,Exon_getPhase(orig));  
  Exon_setEndPhase(newEx,Exon_getEndPhase(orig));  
  Exon_setAnalysis(newEx,Exon_getAnalysis(orig));  

  return newEx;
}
