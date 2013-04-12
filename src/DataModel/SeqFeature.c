#define __SEQFEATURE_MAIN__
#include "SeqFeature.h"
#undef __SEQFEATURE_MAIN__

#include "DBAdaptor.h"
#include "StrUtil.h"
#include "AssemblyMapperAdaptor.h"
#include "RawContigAdaptor.h"
#include "SliceAdaptor.h"
#include "SeqFeatureFactory.h"
#include "ChromosomeAdaptor.h"
#include "SimpleFeature.h"

SeqFeature *SeqFeature_new(void) {
  SeqFeature *sf;

  if ((sf = (SeqFeature *)calloc(1,sizeof(SeqFeature))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for seq feature\n");
    return NULL;
  }

  sf->funcs = &seqFeatureFuncs;

  return sf;
}


ECOSTRING SeqFeature_setSeqName(SeqFeature *sf, char *seqName) {
  EcoString_copyStr(ecoSTable, &(sf->seqName),seqName,0);

  return sf->seqName;
}

ECOSTRING SeqFeature_getSeqName(SeqFeature *sf) {
  BaseContig *contig = SeqFeature_getContig(sf);

  if (contig) {
    return BaseContig_getName(contig);
  } else if (sf->seqName) {
    return sf->seqName;
  } else {
    fprintf(stderr,"Warning: No seq name defined for feature\n");
    return emptyString;
  }
}

int SeqFeature_startCompFunc(const void *a, const void *b) {
  SeqFeature **e1 = (SeqFeature **)a;
  SeqFeature **e2 = (SeqFeature **)b;

  if (SeqFeature_getStart(*e1) > SeqFeature_getStart(*e2)) {
    return 1;
  } else if (SeqFeature_getStart(*e1) < SeqFeature_getStart(*e2)) {
    return -1;
  } else {
    return 0;
  }
}

int SeqFeature_reverseStartCompFunc(const void *a, const void *b) {
  SeqFeature **e1 = (SeqFeature **)a;
  SeqFeature **e2 = (SeqFeature **)b;

  if (SeqFeature_getStart(*e2) > SeqFeature_getStart(*e1)) {
    return 1;
  } else if (SeqFeature_getStart(*e2) < SeqFeature_getStart(*e1)) {
    return -1;
  } else {
    return 0;
  }
}

Vector *SeqFeature_transformToRawContigImpl(SeqFeature *sf) {
  BaseContig *featContig = SeqFeature_getContig(sf);

  if (featContig) {
    if (featContig->objectType == CLASS_RAWCONTIG) {
      // we are already in rawcontig coords, nothing needs to be done
      Vector_setElementAt(singleEntryVector, 0, sf);
      return singleEntryVector;
    } else if (featContig->objectType == CLASS_SLICE) {
      // transform to raw_contig coords from Slice coords
      return SeqFeature_transformSliceToRawContig(sf);
    } else {
      // Unknown contig type
      fprintf(stderr, "Error: Cannot transform unknown contig type %d\n",featContig->objectType);
      exit(1);
    }
  } else {
    // Can't convert to rawcontig coords without a contig to work with
    fprintf(stderr, "Error: Objects contig is not defined - cannot transform\n");
    exit(1);
  }
}

SeqFeature *SeqFeature_transformToSliceImpl(SeqFeature *sf, Slice *slice) {
  BaseContig *featContig = SeqFeature_getContig(sf);

  if (featContig) {
    if (featContig->objectType == CLASS_RAWCONTIG)  {
      // transform to slice coords from raw contig coords
      return SeqFeature_transformRawContigToSlice(sf, slice);
    } else if (featContig->objectType == CLASS_SLICE) {
      // transform to slice coords from other slice coords
      fprintf(stderr, "Error: Transforms between slices have not been implemented - Nag Steve\n");
      exit(1);
    } else {
      // Unknown contig type
      fprintf(stderr, "Error: Cannot transform unknown contig type %d\n",featContig->objectType);
      exit(1);
    }
  } else {
    // Can't convert to slice coords without a contig to work with
    fprintf(stderr, "Error: Objects contig is not defined - cannot transform\n");
    exit(1);
  }
}

SeqFeature *SeqFeature_transformRawContigToSliceImpl(SeqFeature *sf, Slice *slice) {
  DBAdaptor *dba;
  RawContigAdaptor *rca;
  AssemblyMapperAdaptor *ama;
  AssemblyMapper *assMapper;
  RawContig *rc;
  MapperRangeSet *mapped;
  MapperCoordinate *mc;

  if (!SeqFeature_getContig(sf)) {
    fprintf(stderr, "Error: can't transform coordinates of %p without a contig defined\n", sf);
    exit(1);
  }

  rc = (RawContig *)SeqFeature_getContig(sf);

  if (!RawContig_getAdaptor(rc)) {
    fprintf(stderr, "Error: can't transform coordinates of %p without an adaptor "
                    "attached to the feature's contig\n", sf);
    exit(1);
  }
  dba = RawContig_getAdaptor(rc)->dba;

  ama = DBAdaptor_getAssemblyMapperAdaptor(dba);
  assMapper = AssemblyMapperAdaptor_fetchByType(ama, Slice_getAssemblyType(slice));

  rca = DBAdaptor_getRawContigAdaptor(dba);

  mapped = AssemblyMapper_mapCoordinatesToAssembly(
    assMapper,
    RawContig_getDbID(rc),
    SeqFeature_getStart(sf),
    SeqFeature_getEnd(sf),
    SeqFeature_getStrand(sf)
  );

  if (!mapped || mapped->nRange == 0) {
    fprintf(stderr, "ERROR: SeqFeature couldnt map dbID = " IDFMTSTR "\n",SeqFeature_getDbID(sf));
    exit(1);
  }

  if (mapped->nRange != 1) {
    fprintf(stderr, "Error: seq feature should only map to one chromosome - "
                    "something bad has happened ...\n");
    exit(1);
  }


  if (MapperRangeSet_getRangeAt(mapped,0)->rangeType == MAPPERRANGE_GAP) {
    fprintf(stderr, "Warning: feature lies on gap\n");
    return NULL;
  }

  mc = (MapperCoordinate *)MapperRangeSet_getRangeAt(mapped,0);


  // the slice is an empty slice, create an enitre chromosome slice and
  // replace the empty slice with it
  if (Slice_getEmptyFlag(slice)) {
    SliceAdaptor *sa = DBAdaptor_getSliceAdaptor(dba);
    ChromosomeAdaptor *ca = DBAdaptor_getChromosomeAdaptor(dba);
    char *chrName = Chromosome_getName(ChromosomeAdaptor_fetchByDbID(ca,mc->id));

    slice = SliceAdaptor_fetchByChrName(sa, chrName);
  }

  // mapped coords are on chromosome - need to convert to slice
  if(Slice_getStrand(slice) == 1) {
    SeqFeature_setStart(sf, mc->start - Slice_getChrStart(slice) + 1);
    SeqFeature_setEnd(sf, mc->end   - Slice_getChrStart(slice) + 1);
    SeqFeature_setStrand(sf, mc->strand);
  } else {
    SeqFeature_setStart(sf, Slice_getChrEnd(slice) - mc->end   + 1);
    SeqFeature_setEnd(sf, Slice_getChrEnd(slice) - mc->start + 1);
    SeqFeature_setStrand(sf, mc->strand * -1);
  }

  //set the contig to the slice
  SeqFeature_setContig(sf, slice);

  return sf;
}

Vector *SeqFeature_transformSliceToRawContigImpl(SeqFeature *sf) {
  Slice *slice;
  DBAdaptor *dba;
  RawContigAdaptor *rca;
  AssemblyMapperAdaptor *ama;
  AssemblyMapper *assMapper;
  Vector *out;
  int start;
  int end;
  int strand;
  MapperRangeSet *mapped;
  MapperCoordinate *mc;

  if (!SeqFeature_getContig(sf)) {
    fprintf(stderr, "Error: can't transform coordinates of %p without a contig defined\n", sf);
    exit(1);
  }

  slice = (Slice *)SeqFeature_getContig(sf);

  if (!Slice_getAdaptor(slice)) {
    fprintf(stderr, "Error: can't transform coordinates of %p without an adaptor "
                    "attached to the feature's slice\n", sf);
    exit(1);
  }

  dba = Slice_getAdaptor(slice)->dba;

  ama = DBAdaptor_getAssemblyMapperAdaptor(dba);
  assMapper = AssemblyMapperAdaptor_fetchByType(ama, Slice_getAssemblyType(slice));

  rca = DBAdaptor_getRawContigAdaptor(dba);

  // first convert the features coordinates to assembly coordinates
  if (Slice_getStrand(slice) == 1) {
    start  = Slice_getChrStart(slice) + SeqFeature_getStart(sf) - 1;
    end    = Slice_getChrStart(slice) + SeqFeature_getEnd(sf)   - 1;
    strand = SeqFeature_getStrand(sf);
  } else {
    start  = Slice_getChrEnd(slice) - SeqFeature_getEnd(sf)   + 1;
    end    = Slice_getChrEnd(slice) - SeqFeature_getStart(sf) + 1;
    strand = SeqFeature_getStrand(sf) * -1;
  }

  // convert the assembly coordinates to RawContig coordinates
  mapped = AssemblyMapper_mapCoordinatesToRawContig(assMapper,
    Slice_getChrId(slice),
    start,
    end,
    strand
  );

  if (!mapped || mapped->nRange == 0) {
    fprintf(stderr, "ERROR: SeqFeature couldnt map dbID = " IDFMTSTR "\n",SeqFeature_getDbID(sf));
    exit(1);
  }


  if (mapped->nRange == 1) {
    RawContig *rc;

    if (MapperRangeSet_getRangeAt(mapped,0)->rangeType == MAPPERRANGE_GAP) {
      fprintf(stderr, "Warning: feature lies on gap\n");
      Vector_setElementAt(singleEntryVector,0,sf);
      return singleEntryVector;
    }

    mc = (MapperCoordinate *)MapperRangeSet_getRangeAt(mapped,0);

    rc = RawContigAdaptor_fetchByDbID(rca, mc->id);

    SeqFeature_setStart (sf, mc->start);
    SeqFeature_setEnd   (sf, mc->end);
    SeqFeature_setStrand(sf, mc->strand);
    // NIY SeqFeature_set$self->seqname   ($mapped[0]->id);
    //printf(stderr, "setting contig to be ".$mapped[0]->id."\n";
    SeqFeature_setContig(sf, rc);

    Vector_setElementAt(singleEntryVector,0,sf);
    return singleEntryVector;

  } else {
    Vector *gaps   = Vector_new();
    Vector *coords = Vector_new();
    int i;

    // more than one object returned from mapper
    // possibly more than one RawContig in region

    for (i=0; i<mapped->nRange; i++) {
      MapperRange *mr = MapperRangeSet_getRangeAt(mapped, i);
    
      if (mr->rangeType == MAPPERRANGE_GAP) {
        Vector_addElement(gaps, mr);
      } else if (mr->rangeType == MAPPERRANGE_COORD) {
        Vector_addElement(coords, mr);
      } else {
        fprintf(stderr,"Error: Unknown range type\n");
        exit(1);
      }
    }

    // case where only one RawContig maps
    if (Vector_getNumElement(coords) == 1) {
      mc = (MapperCoordinate *)Vector_getElementAt(coords,0);
 
      SeqFeature_setStart (sf, mc->start);
      SeqFeature_setEnd   (sf, mc->end);
      SeqFeature_setStrand(sf, mc->strand);
      // NIY $self->seqname($coords[0]->id);
      //print STDERR "2 setting contig to be ".$coords[0]->id."\n";
      SeqFeature_setContig(sf, RawContigAdaptor_fetchByDbID(rca, mc->id));

      fprintf(stderr, "Warning: Feature [%p] truncated as lies partially on a gap\n", sf);
      Vector_setElementAt(singleEntryVector,0,sf);
      out = singleEntryVector;

    } else {
      if (SeqFeature_getIsSplittable(sf)) {
        fprintf(stderr, "Warning: Feature spans >1 raw contig - can't split\n");
        Vector_setElementAt(singleEntryVector,0,sf);
// NIY check that this should be a return
        out = singleEntryVector;
      } else {
  
        out = Vector_new();
    
        for (i=0; i<mapped->nRange; i++) {
          SeqFeature *feat;
          MapperCoordinate *mc;
          MapperRange *mr = MapperRangeSet_getRangeAt(mapped,i);
    
          if (mr->rangeType == MAPPERRANGE_GAP) {
            fprintf(stderr, "Warning: piece of seq feature lies on gap\n");
            continue;
          }
    
          mc = (MapperCoordinate *)mr;
  
          feat = SeqFeatureFactory_newFeature(sf->objectType);
    
          SeqFeature_setStart(feat, mc->start);
          SeqFeature_setEnd(feat, mc->end);
          SeqFeature_setStrand(feat, mc->strand);
          fprintf(stderr, "3 setting contig to be " IDFMTSTR "\n",mc->id);
          SeqFeature_setContig(feat, RawContigAdaptor_fetchByDbID(rca, mc->id));
          if (SeqFeature_getAdaptor(sf)) SeqFeature_setAdaptor(sf, SeqFeature_getAdaptor(sf));
    // HACK HACK HACK
          if (Class_isDescendent(CLASS_SIMPLEFEATURE, sf->objectType)) {
            SimpleFeature_setDisplayLabel((SimpleFeature *)feat, SimpleFeature_getDisplayLabel((SimpleFeature *)sf));
          }
          SeqFeature_setAnalysis(feat,SeqFeature_getAnalysis(sf));
    
          Vector_addElement(out,feat);
        }
      }
    }
    //NIY freeing mapper and coord etc.
    Vector_free(coords);
    Vector_free(gaps);

    return out;
  }
}

void SeqFeature_freePtrs(SeqFeature *sf) {
  if (sf->seqName)  EcoString_freeStr(ecoSTable, sf->seqName);
  if (sf->analysis) Analysis_free(sf->analysis);
// NIY Is this the right thing to do
  // NIY Freeing contig if (sf->contig)   BaseContig_free(sf->contig);
}
