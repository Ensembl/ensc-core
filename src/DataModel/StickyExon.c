#include "StickyExon.h"
#include "AssemblyMapperAdaptor.h"
#include "ChromosomeAdaptor.h"
#include "SliceAdaptor.h"
#include "BaseContig.h"
#include "Mapper.h"

Exon *StickyExon_transformToSlice(Exon *exon, Slice *slice)  {
  BaseContig *contig;

  contig = Exon_getContig(exon);
  
  if (contig && BaseContig_getObjectType(contig) == CLASS_RAWCONTIG) {
    AssemblyMapperAdaptor *ama;
    AssemblyMapper *assMapper;
    int mappedStart = 0;
    int mappedEnd = -1;
    int compositeExonStrand = 0;
    int compositeExonPhase;
    int compositeExonEndPhase;
    int i;
    Exon *firstCompEx; // Used for setting storable info on new exon
    BaseAdaptor *adaptor;
    Exon *newExon;

    adaptor = Slice_getAdaptor(slice);

    ama = DBAdaptor_getAssemblyMapperAdaptor(Slice_getAdaptor(slice)->dba);
    assMapper = AssemblyMapperAdaptor_fetchByType(ama,Slice_getAssemblyType(slice));


#ifdef DONE
    my @supporting_features;
#endif
    
    // sort the component exons
    Exon_sortByStickyRank(exon); 
    
    for (i=0; i<Exon_getNumComponentExon(exon); i++) {
      Exon *compExon = Exon_getComponentExonAt(exon,i); 
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
	return exon;
      }

      mcoord = (MapperCoordinate *)mrange;

/*
      fprintf(stderr, "[Mapped start  %d ", mcoord->start);
      fprintf(stderr, "\tMapped end %d ", mcoord->end);
      fprintf(stderr, "\tMapped strand %d ",  mcoord->strand);
      fprintf(stderr, "]\nSticky rank : %d\n", Exon_getStickyRank(compExon));
*/

#ifdef DONE
      // add the supporting features from the exons
      // each exon has the pieces of the supporting features that fall in the corresponding contig
      // they've been split before and at the moment they are not re-combined
      push @supporting_features, @{$compExon->get_all_supporting_features()}; 
#endif
      

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
// NIY free old slice (or have special empty one)???
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

    if(Slice_getStrand(slice) == 1) {
      Exon_setStart(newExon, mappedStart - Slice_getChrStart(slice) + 1);
      Exon_setEnd(newExon, mappedEnd   - Slice_getChrStart(slice) + 1);
      Exon_setStrand(newExon, compositeExonStrand);
    } 
    else {
      Exon_setStart(newExon, Slice_getChrEnd(slice) - mappedEnd   + 1);
      Exon_setEnd(newExon, Slice_getChrEnd(slice) - mappedStart + 1);
      Exon_setStrand(newExon, compositeExonStrand * -1);
    }
    Exon_setDbID(newExon, Exon_getDbID(firstCompEx));
    Exon_setAdaptor(newExon, Exon_getAdaptor(firstCompEx));

    Exon_setContig(newExon, slice);
    
#ifdef DONE
    // copy each of the supporting features and transform them
    my @feats;
    foreach my $sf (@supporting_features) {
      push @feats, $sf->transform($slice);
    }
    $newexon->add_supporting_features(@feats);
#endif

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

