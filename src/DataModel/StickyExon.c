#define __STICKYEXON_MAIN__
#include "StickyExon.h"
#undef __STICKYEXON_MAIN__

#include "DBAdaptor.h"
#include "AssemblyMapperAdaptor.h"
#include "ChromosomeAdaptor.h"
#include "SliceAdaptor.h"
#include "BaseContig.h"
#include "Mapper.h"
#include "RawContigAdaptor.h"
#include "StrUtil.h"
#include "BaseAlignFeature.h"

StickyExon *StickyExon_new(void) {
  StickyExon *stickyExon;

  if ((stickyExon = (StickyExon *)calloc(1,sizeof(StickyExon))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for stickyExon\n");
    return NULL;
  }

/* Set to empty values */
  StickyExon_setModified(stickyExon,0);
  StickyExon_setCreated(stickyExon,0);
  StickyExon_setVersion(stickyExon,-1);

  stickyExon->objectType = CLASS_STICKYEXON;

  stickyExon->funcs = &stickyExonFuncs;

  return stickyExon;
}

void StickyExon_sortByStickyRank(StickyExon *exon) {
  qsort(StickyExon_getComponents(exon), StickyExon_getComponentExonCount(exon), sizeof(void *),
        StickyExon_stickyRankCompFunc);
  return;
}

int StickyExon_stickyRankCompFunc(const void *a, const void *b) {
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



Exon *StickyExon_transformRawContigToSlice(StickyExon *exon, Slice *slice)  {
  BaseContig *contig;

  contig = Exon_getContig(exon);
  
  if (contig && BaseContig_getObjectType(contig) == CLASS_RAWCONTIG) {
    AssemblyMapperAdaptor *ama;
    AssemblyMapper *assMapper;
    int mappedStart = 0;
    int mappedEnd = -1;
    int compositeExonStrand = 0;
    int compositeExonPhase = 0;
    int compositeExonEndPhase = 0;
    int i;
    Exon *firstCompEx = NULL; // Used for setting storable info on new exon
    BaseAdaptor *adaptor;
    Exon *newExon;
    Vector *supportingFeatures = Vector_new();
    Vector *mappedFeatures = Vector_new();

    adaptor = Slice_getAdaptor(slice);

    ama = DBAdaptor_getAssemblyMapperAdaptor(Slice_getAdaptor(slice)->dba);
    assMapper = AssemblyMapperAdaptor_fetchByType(ama,Slice_getAssemblyType(slice));

    // sort the component exons
    StickyExon_sortByStickyRank(exon); 
    
    for (i=0; i<StickyExon_getComponentExonCount(exon); i++) {
      Exon *compExon = StickyExon_getComponentExonAt(exon,i); 
      MapperRangeSet *mapped;
      MapperRange *mrange;
      MapperCoordinate *mcoord;

/*
      fprintf(stderr," component exon %s",Exon_getStableId(exon));
*/

      mapped = AssemblyMapper_mapCoordinatesToAssembly(assMapper,
	 BaseContig_getDbID(Exon_getContig(compExon)),
	 Exon_getStart(compExon),
	 Exon_getEnd(compExon),
	 Exon_getStrand(compExon)
	);
      
/*
      fprintf(stderr, "\nStickyExon transform method:\n"); 
*/
      
      // exons should always transform so in theory no error check
      // necessary
      if (!mapped || mapped->nRange==0) {
	fprintf(stderr, "ERROR: Component Sticky Exon couldnt map\n" );
        exit(1);
      }

      mrange = MapperRangeSet_getRangeAt(mapped,0);
    

      // should get a gap object returned if an exon lies outside of 
      // the current slice.  Simply return the exon as is - i.e. untransformed.
      // this untransformed exon will be distinguishable as it will still have
      // contig attached to it and not a slice.
      if( mrange->rangeType == MAPPERRANGE_GAP) {
        fprintf(stderr, "sticky exon mapped to gap\n");
	return (Exon *)exon;
      }

      mcoord = (MapperCoordinate *)mrange;

/*
      fprintf(stderr, "[Mapped start  %d ", mcoord->start);
      fprintf(stderr, "\tMapped end %d ", mcoord->end);
      fprintf(stderr, "\tMapped strand %d ",  mcoord->strand);
      fprintf(stderr, "]\nSticky rank : %d\n", Exon_getStickyRank(compExon));
*/

      // add the supporting features from the exons
      // each exon has the pieces of the supporting features that fall in the corresponding contig
      // they've been split before and at the moment they are not re-combined
// NIY Sticky exon vector free??
      Vector_append(supportingFeatures, Exon_getAllSupportingFeatures(compExon));
      

      // now pull out the start and end points of the newly concatenated sequence
      // if we've got the first sticky exon, store the relevant info
      
      if (Exon_getStickyRank(compExon) == 1 ) {
	// this assumes that the strand of the sticky exon is
	// set by the first one - is this correct?
        firstCompEx = compExon;

	compositeExonStrand = mcoord->strand;

        // the slice is an empty slice, create an enitre chromosome slice and
        // replace the empty slice with it
        if (Slice_getEmptyFlag(slice)) {
          SliceAdaptor *sa = DBAdaptor_getSliceAdaptor(adaptor->dba);
          ChromosomeAdaptor *ca = DBAdaptor_getChromosomeAdaptor(adaptor->dba);
          char *chrName = Chromosome_getName(ChromosomeAdaptor_fetchByDbID(ca,mcoord->id));

          slice = SliceAdaptor_fetchByChrName(sa, chrName);
        }


	if (compositeExonStrand == 1 ) {
	  // this means that in the forward strand the first component exon at the 5'end  in the slice
	  mappedStart = mcoord->start;
	  compositeExonPhase = Exon_getPhase(compExon);

	} else { // for the reverse strand case it is something different
	  
	  // this means that in the reverse strand the first component exon at the 5'end  in the slice
	  mappedEnd = mcoord->end;

	  // phase follows te 5' -> 3' sense, in contradistinction to start-end 
	  compositeExonPhase = Exon_getPhase(compExon);
	}
      }
      // now do the end point
      // keep storing as you iterate over the component exons
      // since the exons are previously sorted based on their sticky rank
      // then the last set of stored values will be the last component exon
      else {
	if (mcoord->strand == 1 ) {	  
	  mappedEnd = mcoord->end;
	  compositeExonEndPhase = Exon_getEndPhase(compExon);
	}
	else {
	  mappedStart = mcoord->start;
	  compositeExonEndPhase = Exon_getEndPhase(compExon);
	}
      }
      MapperRangeSet_free(mapped);
    }
    
    // now build the new composite exon
    newExon = Exon_new();
// NIY free old exons

    if (Slice_getStrand(slice) == 1) {
      Exon_setStart(newExon, mappedStart - Slice_getChrStart(slice) + 1);
      Exon_setEnd(newExon, mappedEnd   - Slice_getChrStart(slice) + 1);
      Exon_setStrand(newExon, compositeExonStrand);
    } else {
      Exon_setStart(newExon, Slice_getChrEnd(slice) - mappedEnd   + 1);
      Exon_setEnd(newExon, Slice_getChrEnd(slice) - mappedStart + 1);
      Exon_setStrand(newExon, compositeExonStrand * -1);
    }
    Exon_setDbID(newExon, Exon_getDbID(firstCompEx));
    Exon_setAdaptor(newExon, Exon_getAdaptor(firstCompEx));

    Exon_setContig(newExon, slice);
    
    // copy each of the supporting features and transform them
    for (i=0; i<Vector_getNumElement(supportingFeatures); i++) {
      BaseAlignFeature *baf = Vector_getElementAt(supportingFeatures,i);
      Vector_addElement(mappedFeatures, BaseAlignFeature_transformToSlice(baf,slice));
    }
    Exon_addSupportingFeatures(newExon, mappedFeatures);

    Exon_setPhase(newExon, compositeExonPhase );
    Exon_setEndPhase(newExon, compositeExonEndPhase );
/*
    fprintf(stderr, "transformed exon %s\n",Exon_getStableId(newExon));
*/
    return newExon;

  } else {
    fprintf(stderr, "ERROR: Unexpected StickyExon in Assembly coords ...\n");
    exit(1);
  }
}

void StickyExon_loadGenomicMapper(StickyExon *stickyExon, Mapper *mapper, IDType id, int start) {
  int i;
  
  for (i=0; i<StickyExon_getComponentExonCount(stickyExon); i++) {
    Exon *exon = StickyExon_getComponentExonAt(stickyExon, i);
    Mapper_addMapCoordinates(mapper, id, start, start+Exon_getLength(exon)-1,
                                  Exon_getStrand(exon), (IDType)Exon_getContig(exon),
                                  Exon_getStart(exon), Exon_getEnd(exon) );
    start += Exon_getLength(exon);
  }
}

Exon *StickyExon_adjustStartEnd(StickyExon *stickyExon, int startAdjust, int endAdjust) {
  Exon *newExon;
  int start;
  int end;
  Mapper *mapper;
  int currentStart = 1;
  int i;
  MapperRangeSet *mapped;
  RawContigAdaptor *rca = NULL;
  Exon *firstComponent;
  RawContig *firstContig;

  if (startAdjust == 0 && endAdjust == 0 ) {
    return (Exon *)stickyExon;
  }

  start = 1 + startAdjust;
  end = StickyExon_getLength(stickyExon) + endAdjust;

  mapper = Mapper_new( CDNA_COORDS, GENOMIC_COORDS );

  for (i=0; i<StickyExon_getComponentExonCount(stickyExon); i++) {
    Exon *exon = StickyExon_getComponentExonAt(stickyExon, i);
    Mapper_addMapCoordinates(mapper, (IDType)stickyExon, currentStart, currentStart+Exon_getLength(exon)-1,
                             Exon_getStrand(exon), (IDType)Exon_getContig(exon), 
                             Exon_getStart(exon), Exon_getEnd(exon) );
    currentStart += Exon_getLength(exon);
  }

  firstComponent = StickyExon_getComponentExonAt(stickyExon, 0);
  firstContig    = (RawContig *)Exon_getContig(firstComponent);

  Class_assertType(CLASS_RAWCONTIG, firstContig->objectType);

  rca = (RawContigAdaptor *)RawContig_getAdaptor(firstContig);
  
  if (!rca) {
    fprintf(stderr, "Error: No RawContigAdaptor - cannot adjust sticky exon\n");
    exit(1);
  }

  mapped = Mapper_mapCoordinates(mapper, (IDType)stickyExon, start, end, 1, CDNA_COORDS );

  if (mapped->nRange == 1 ) {
    MapperCoordinate *coord;
    RawContig *rc;

    // we can return a normal exon
    newExon = Exon_new();
    coord = (MapperCoordinate *)MapperRangeSet_getRangeAt(mapped,0);

// Copy
    Exon_copy(newExon, (Exon *)stickyExon, SHALLOW_DEPTH);

    Exon_setStart(newExon, coord->start);
    Exon_setEnd(newExon,   coord->end);
    Exon_setStrand(newExon,coord->strand);

//Is just pointer to contig
    Exon_setContig(newExon, (BaseContig *)coord->id);

  } else {
    MapperCoordinate *coord;
    int i;
    RawContig *rc;

    // make a new sticky Exon
    // cast to Exon to stop compiler complaining
    newExon = (Exon *)StickyExon_new();

// Copy
    Exon_copy(newExon, (Exon *)stickyExon, SHALLOW_DEPTH);

    StickyExon_setStart(newExon, 1);
    StickyExon_setEnd(newExon, end - start + 1);
    StickyExon_setStrand(newExon, 1);

    for (i=0; i<mapped->nRange; i++) {
      Exon *cex = Exon_new();

      coord = (MapperCoordinate *)MapperRangeSet_getRangeAt(mapped,i);

// Copy
      Exon_copy(cex, (Exon *)stickyExon, SHALLOW_DEPTH);

      Exon_setStart(cex,  coord->start);
      Exon_setEnd(cex,    coord->end);
      Exon_setStrand(cex, coord->strand);
// Is just pointer to contig
      Exon_setContig(cex, (BaseContig *)coord->id);

      StickyExon_addComponentExon((StickyExon *)newExon, cex);
    }
  }

  return newExon;
}

int StickyExon_getLength(StickyExon *stickyExon) {
  int i;
  int len = 0;

  for (i=0; i<StickyExon_getComponentExonCount(stickyExon); i++) {
    Exon *subExon = StickyExon_getComponentExonAt(stickyExon, i);
    len += Exon_getLength(subExon);
  }
  return len;
}

/* deliberately char */
#ifdef DONE
void StickyExon_setSeq(StickyExon *stickyExon, char *seq) {
  stickyExon->seq = seq;
}

sub StickyExon_getSeq(StickyExon *stickyExon) {
  int i;
  char *seqString = NULL;

  my $seqString = "";

  for (i=0; i<StickyExon_getComponentExonCount(stickyExon); i++) {
    Exon *cExon = StickyExon_getComponentExonAt(stickyExon, i);
    seqString = StrUtil_appendString(seqString, Exon_getSeqString(cExon));
  }
  $self->{'_seq'} = $seqString;

  return Bio::Seq->new( -seq => $self->{'_seq'} );
}
#endif

char *StickyExon_getSeqString(StickyExon *stickyExon) {
  int i;
  char *seqString = NULL;

  if (StickyExon_getSeqCacheString(stickyExon)) {
    return StickyExon_getSeqCacheString(stickyExon);
  } else {
    for (i=0; i<StickyExon_getComponentExonCount(stickyExon); i++) {
      Exon *cExon = StickyExon_getComponentExonAt(stickyExon, i);
      seqString = StrUtil_appendString(seqString, Exon_getSeqString(cExon));
    }
  }

  StickyExon_setSeqCacheString(stickyExon, seqString);

  return seqString;
}

char *StickyExon_getPeptide(StickyExon *stickyExon, Transcript *trans) {
  char *peptide;
  int pepStart = -1;
  int pepEnd = -1; 
  int i;

  if (!trans) {
    fprintf(stderr,"Error: transcript arg must not be null\n");
    exit(1);
  }

  for (i=0;i<StickyExon_getComponentExonCount(stickyExon); i++) {
    Exon *exon = StickyExon_getComponentExonAt(stickyExon,i);
    int j;
    MapperRangeSet *mapped;
    Vector *coords;

    mapped = Transcript_genomic2Pep(trans, Exon_getStart(exon), Exon_getEnd(exon),
                                    Exon_getStrand(exon),Exon_getContig(exon));
  
    coords = Vector_new();
  
    // filter out gaps
    for (j=0;i<mapped->nRange;i++) {
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
  
      // set the pep start to the minimum of all coords
      if(pepStart == -1 || c->start < pepStart) {
        pepStart = c->start;
      }

      // set the pep end to the maximum of all coords
      if(pepEnd == -1 || c->end > pepEnd) {
        pepEnd = c->end;
      }
    }
  
    Vector_free(coords,NULL);
    MapperRangeSet_free(mapped);
  }

  // the peptide of this sticky is the region spanned by the component exons
  if (pepStart != -1 && pepEnd != -1) {
    char *wholePeptide = Transcript_translate(trans); 
// NIYshould be $tr->translate->subseq($pep_start, $pep_end);
// NIY check for off by one
// NIY free wholePeptide
    peptide = StrUtil_substr(wholePeptide, pepStart, (pepEnd-pepStart+1));
    
  }

  return peptide;
}

Vector *StickyExon_getAllSupportingFeatures(StickyExon *stickyExon) {
  Vector *out = Vector_new();
  int i;

  for (i=0; i<StickyExon_getComponentExonCount(stickyExon); i++) {
    Exon *subExon = StickyExon_getComponentExonAt(stickyExon, i);
    int j;

    Vector_append(out, Exon_getAllSupportingFeatures(subExon));
  }

  return out;
}



void StickyExon_addSupportingFeatures(StickyExon *stickyExon, Vector *features) {
  int beenAdded = 0;
  int i;


 // check whether this feature object has been added already
  for (i=0; i<Vector_getNumElement(features); i++) {
    SeqFeature *feature = Vector_getElementAt(features,i);
    int j;

    beenAdded = 0;
    for (j=0; j<StickyExon_getComponentExonCount(stickyExon); j++) {
      Exon *cexon = StickyExon_getComponentExonAt(stickyExon,j);
      BaseContig *cexonContig = Exon_getContig(cexon);
      BaseContig *featContig = SeqFeature_getContig(feature);
      if (cexonContig && featContig &&
          EcoString_strcmp(BaseContig_getName(cexonContig), BaseContig_getName(featContig))) {
        Vector_setElementAt(singleEntryVector,0,feature);
        Exon_addSupportingFeatures(cexon, singleEntryVector);
        beenAdded = 1;
      }
    }

    if (!beenAdded) {
      fprintf(stderr, "Warning: SupportingFeature could not be added, not on same contig "
                      "as component exons.\n");
    }
  }
}

