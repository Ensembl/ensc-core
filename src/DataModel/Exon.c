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
#include "RawContigAdaptor.h"
#include "StrUtil.h"

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

Exon *Exon_transformRawContigToSliceImpl(Exon *exon, Slice *slice) {
  Exon *newExon;
  BaseAdaptor *adaptor;
  AssemblyMapperAdaptor *ama;
  AssemblyMapper *assMapper;
  MapperRangeSet *mapped;
  MapperCoordinate *mc;

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
  if (!mapped || mapped->nRange == 0) {
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

  newExon = Exon_new();

  newExon = Exon_copy(newExon, exon, SHALLOW_DEPTH);

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

Exon *Exon_copy(Exon *copy, Exon *orig, CopyDepth depth) {
  if (depth != SHALLOW_DEPTH) {
    fprintf(stderr, "ERROR: Only SHALLOW copy implemented in Exon\n");
    exit(1);
  }
  Exon_setStart(copy,Exon_getStart(orig));  
  Exon_setEnd(copy,Exon_getEnd(orig));  
  Exon_setDbID(copy,Exon_getDbID(orig));  
  Exon_setAdaptor(copy,Exon_getAdaptor(orig));  
  Exon_setStrand(copy,Exon_getStrand(orig));  
  Exon_setPhase(copy,Exon_getPhase(orig));  
  Exon_setEndPhase(copy,Exon_getEndPhase(orig));  
  Exon_setAnalysis(copy,Exon_getAnalysis(orig));  
  Exon_setContig(copy,Exon_getContig(orig));  

  return copy;
}

void Exon_findSupportingEvidence(Exon *exon, Vector *features, int isSorted) {
  int i;
  Vector *support = Vector_new();

  for (i=0; i<Vector_getNumElement(features); i++) {
    SeqFeature *f = Vector_getElementAt(features,i);

    // return if we have a sorted feature array
    if (isSorted == 1 && SeqFeature_getStart(f) > Exon_getEnd(exon)) {
      return;
    }
/* NIY subfeatures
    if ($f->sub_SeqFeature) {
      my @subf = $f->sub_SeqFeature;

      $self->find_supporting_evidence(\@subf);
    } else {
*/
/* NIY Speedup */
      if (!strcmp(BaseContig_getName(SeqFeature_getContig(f)), 
                  BaseContig_getName(Exon_getContig(exon)))) {
        if (SeqFeature_getEnd(f) >= Exon_getStart(exon) && 
            SeqFeature_getStart(f) <= Exon_getEnd(exon) && 
            SeqFeature_getStrand(f) == Exon_getStrand(exon)) {
          Vector_addElement(support,f);
        }
      }
/* NIY subfeatures
    }
*/
  }
  Exon_addSupportingFeatures(exon, support);

  Vector_free(support,NULL);
}

Vector *Exon_getAllSupportingFeaturesImpl(Exon *exon) {
  return exon->supportingFeatures;
}

void Exon_addSupportingFeaturesImpl(Exon *exon, Vector *v) {
  Vector_append(exon->supportingFeatures,v);
}

Exon *Exon_adjustStartEndImpl(Exon *exon, int startAdjust, int endAdjust) {

  Exon *newExon = Exon_new();

// Copy - NIY won't copy support
  Exon_copy(newExon, exon, SHALLOW_DEPTH);

  // invalidate the sequence cache
  // NIY delete $new_exon->{'_seq_cache'};

  if (Exon_getStrand(exon) == 1 ) {
    Exon_setStart(newExon, Exon_getStart(exon) + startAdjust );
    Exon_setEnd(newExon, Exon_getEnd(exon) + endAdjust );
  } else {
    Exon_setStart(newExon, Exon_getStart(exon) - endAdjust );
    Exon_setEnd(newExon, Exon_getEnd(exon) - startAdjust );
  }

// NIY Delete old exon

  return newExon;
}

char *Exon_getPeptideImpl(Exon *exon, Transcript *trans) {
  char *peptide;
  MapperRangeSet *mapped;
  Vector *coords;
  int i;

  if (!trans) {
    fprintf(stderr, "Error: transcript arg non null in getPeptide\n");
    exit(1);
  }

  // convert exons coordinates to peptide coordinates
  mapped = Transcript_genomic2Pep(trans, Exon_getStart(exon), Exon_getEnd(exon),
                                  Exon_getStrand(exon),Exon_getContig(exon));

  coords = Vector_new();

  // filter out gaps
  for (i=0;i<mapped->nRange;i++) {
    MapperRange *mr = MapperRangeSet_getRangeAt(mapped,i);
    if (mr->rangeType == MAPPERRANGE_COORD) {
      Vector_addElement(coords,mr);
    }
  }

  // if this is UTR then the peptide will be empty string
  if (Vector_getNumElement(coords) > 1) {
    fprintf(stderr, "Error. Exon maps to multiple locations in peptide."
                    " Is this exon [%p] a member of this transcript [%p]?",
                    exon,trans);
    exit(1);

  } else if (Vector_getNumElement(coords) == 1) {
    MapperCoordinate *c = Vector_getElementAt(coords,0);
    int start,end;
    char *wholePeptide = Transcript_translate(trans);
    int lenPeptide = strlen(wholePeptide);

    // bioperl doesn't give back residues for incomplete codons
    // make sure we don't subseq too far...

    end = (c->end > lenPeptide) ? lenPeptide : c->end;
    start = (c->start < end) ? c->start : end;
// NIY free
// NIY check for off by one
    peptide = StrUtil_substr(wholePeptide,start,(end-start+1));
  }

  Vector_free(coords,NULL);
  MapperRangeSet_free(mapped);

  return peptide;
}

#ifdef DONE
void Exon_setSeq(Exon *exon, Sequence *seq) {
  fprintf(stderr, "Warning: Exon seq setting not supported\n");
  return;
}

sub Exon_getSeq(Exon *exon) {
  my $self = shift;
  my $arg = shift;

  if( defined $arg ) {
    $self->warn( "seq setting on Exon not supported currently" );
    $self->{'_seq_cache'} = $arg->seq();
  }

  if( defined $self->{'_seq_cache'} ) {
    return Bio::Seq->new(-seq=> $self->{'_seq_cache'});
  }

  my $seq;

  if ( ! defined $self->contig ) {
    $self->warn(" this exon doesn't have a contig you won't get a seq \n");
    return undef;
  }
  else {

    $seq = $self->contig()->subseq($self->start, $self->end);

    if($self->strand == -1){
      $seq =~ tr/ATGCatgc/TACGtacg/;
      $seq = reverse($seq);
    }

   }
  $self->{'_seq_cache'} = $seq;

  return Bio::Seq->new(-seq     => $self->{'_seq_cache'},
                       -id      => $self->stable_id,
                       -moltype => 'dna');
}
#endif
char  *Exon_getSeqStringImpl(Exon *exon) {
  char *seq;

  if (Exon_getSeqCacheString(exon)) {
    return Exon_getSeqCacheString(exon);
  }

  if (!Exon_getContig(exon)) {
    fprintf(stderr, "Warning: this exon %s doesn't have a contig you won't get a seq\n", Exon_getStableId(exon));
    return NULL;
  } else {

    seq = BaseContig_getSubSeq(Exon_getContig(exon), 
                               Exon_getStart(exon), 
                               Exon_getEnd(exon),
                               1);

    if (Exon_getStrand(exon) == -1){
      SeqUtil_reverseComplement(seq,strlen(seq));
    }
  }
  Exon_setSeqCacheString(exon, seq);

  return Exon_getSeqCacheString(exon);
}

Exon *Exon_transformSliceToRawContigImpl(Exon *exon) {
  SliceAdaptor *sa;
  Slice *slice;
  AssemblyMapperAdaptor *ama;
  AssemblyMapper *assMapper;
  int sliceChrStart;
  int sliceChrEnd;
  int exonChrStart;
  int exonChrEnd;
  RawContigAdaptor *rca;
  MapperRangeSet *mapped;

  slice = (Slice *)Exon_getContig(exon);

  if (!slice) {
    fprintf(stderr,"Error: Cannot transform exon to raw contig unless it has an attached slice\n");
    exit(1);
  }

  sa = (SliceAdaptor *)Slice_getAdaptor(slice);

  if (!sa) {
    fprintf(stderr,"Error: Cannot transform exon to raw contig unless attached slice"
                   " has adaptor defined. (i.e. exon->contig->adaptor)\n");
    exit(1);
  }

  ama = DBAdaptor_getAssemblyMapperAdaptor(sa->dba);

  assMapper     = AssemblyMapperAdaptor_fetchByType(ama, Slice_getAssemblyType(slice));
  rca           = DBAdaptor_getRawContigAdaptor(sa->dba);
  sliceChrStart = Slice_getChrStart(slice);
  sliceChrEnd   = Slice_getChrEnd(slice);


  if (Slice_getStrand(slice) == 1) {
    exonChrStart = Exon_getStart(exon) + sliceChrStart - 1;
    exonChrEnd   = Exon_getEnd(exon)   + sliceChrStart - 1;
  }
  else {
    exonChrEnd   = sliceChrEnd - Exon_getStart(exon) + 1;
    exonChrStart = sliceChrEnd - Exon_getEnd(exon)   + 1;
  }

  mapped = AssemblyMapper_mapCoordinatesToRawContig(assMapper,
     Slice_getChrId(slice),
     exonChrStart,
     exonChrEnd,
     Exon_getStrand(exon)*Slice_getStrand(slice)
    );

  if (!mapped || !mapped->nRange) {
    fprintf(stderr, "Error: Exon couldnt map" );
    return exon;
  }

  // transform the supporting features to raw contig coords (hashed on contig)
#ifdef DONE
  my %sf_hash;

  if( exists $self->{_supporting_evidence} ) {
    my $sfs = $self->get_all_supporting_features();
  SUPPORTING:foreach my $sf (@$sfs) {
      my @mapped_feats;
      eval{
        @mapped_feats = $sf->transform;
      };
      if($@){
        $self->warn("Supporting feature didn't mapped ignoring $@");
        next SUPPORTING;
      }
      foreach my $mapped_feat (@mapped_feats) {
        unless(exists $sf_hash{$mapped_feat->contig->name}) {
          $sf_hash{$mapped_feat->contig->name} = [];
        }
        push @{$sf_hash{$mapped_feat->contig->name}}, $mapped_feat;
      }
    }
  }
#endif

  if (mapped->nRange > 1 ) {
    int stickyLength = 0;
    int i;

    StickyExon *stickyExon = StickyExon_new();
    StickyExon_setPhase(stickyExon,Exon_getPhase(exon));
    StickyExon_setEndPhase(stickyExon, Exon_getEndPhase(exon));
    StickyExon_setAdaptor(stickyExon, Exon_getAdaptor(exon));
    StickyExon_setStart(stickyExon, 1);
    if (Exon_getDbID(exon)) {
      StickyExon_setDbID(stickyExon, Exon_getDbID(exon));
    }

    // and then all the component exons ...
    for (i=0; i < mapped->nRange; i++ ) {
      Exon *componentExon;
      RawContig *rawContig;
      MapperRange *mr;
      MapperCoordinate *mc;
     
      mr = MapperRangeSet_getRangeAt(mapped, i);

      if( mr->rangeType == MAPPERRANGE_GAP) {
        fprintf(stderr,"Error: exon lies on a gap cannot be mapped\n");
        exit(1);
      }

      mc = (MapperCoordinate *)mr;

      componentExon = Exon_new();

      Exon_setStart(componentExon, mc->start );
      Exon_setEnd(componentExon, mc->end );
      Exon_setStrand(componentExon, mc->strand);

      rawContig = RawContigAdaptor_fetchByDbID(rca, mc->id);
      Exon_setContig(componentExon, rawContig);

      Exon_setStickyRank(componentExon, i+1 );
      Exon_setPhase(componentExon, Exon_getPhase(exon));
      Exon_setEndPhase(componentExon, Exon_getEndPhase(exon));
      Exon_setDbID(componentExon, Exon_getDbID(exon));
      Exon_setAdaptor(componentExon, Exon_getAdaptor(exon));

      // add the supporting features on this contig to the component exon
#ifdef DONE
      if(exists $sf_hash{$rawContig->name}) {
        $componentExon->add_supporting_features(@{$sf_hash{$rawContig->name}});
      }
#endif

      StickyExon_addComponentExon(stickyExon, componentExon);
      stickyLength += ( mc->end - mc->start + 1 );
    }
    StickyExon_setEnd(stickyExon, stickyLength);
    StickyExon_setStrand(stickyExon, 1);
    if (Exon_getStableId(exon)) {
      StickyExon_setStableId(stickyExon, Exon_getStableId(exon));
    }
    if (Exon_getVersion(exon)) {
      StickyExon_setVersion(stickyExon, Exon_getVersion(exon));
    }
    if (Exon_getModified(exon)) {
      StickyExon_setModified(stickyExon, Exon_getModified(exon));
    }
    if (Exon_getCreated(exon)) {
      StickyExon_setCreated(stickyExon, Exon_getCreated(exon));
    }
    return (Exon *)stickyExon;

  } else {
    // thats a simple exon

    RawContig *rawContig;
    Exon *newExon;
    MapperRange *mr;
    MapperCoordinate *mc;
     
    mr = MapperRangeSet_getRangeAt(mapped, 0);

    if (mr->rangeType == MAPPERRANGE_GAP) {
      fprintf(stderr,"Error: exon lies on a gap cannot be mapped\n");
      exit(1);
    }

    mc = (MapperCoordinate *)mr;

    rawContig = RawContigAdaptor_fetchByDbID(rca, mc->id);
    newExon = Exon_new();

    // copy this exon
    Exon_copy(newExon, exon, SHALLOW_DEPTH);

    Exon_setStart(newExon, mc->start);
    Exon_setEnd(newExon, mc->end);
    Exon_setStrand(newExon, mc->strand);
    // attach raw contig
    Exon_setContig(newExon, rawContig);

    // replace old supporting feats with transformed supporting feats
#ifdef DONE
    $newExon->add_supporting_features(@{$sf_hash{$rawContig->name}});
#endif

    return newExon;
  }
}

void Exon_loadGenomicMapperImpl(Exon *exon, Mapper *mapper, IDType id, int start) {

// NIY Make the Exon_getContig consistent
  Mapper_addMapCoordinates( mapper, id, start, start+Exon_getLength(exon)-1,
                            Exon_getStrand(exon), Exon_getContig(exon),
                            Exon_getStart(exon),  Exon_getEnd(exon) );
}



