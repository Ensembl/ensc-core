#define __BASEALIGNFEATURE_MAIN__
#include "BaseAlignFeature.h"
#undef __BASEALIGNFEATURE_MAIN__

#include "DBAdaptor.h"
#include "StrUtil.h"
#include "AssemblyMapperAdaptor.h"
#include "RawContigAdaptor.h"
#include "SliceAdaptor.h"
#include "IDHash.h"
#include "StrUtil.h"
#include "SeqFeatureFactory.h"
#include "CigarStrUtil.h"

#define MAXCIGARPIECELEN 1024

BaseAlignFeature *BaseAlignFeature_new() {
  BaseAlignFeature *baf;

  if ((baf = (BaseAlignFeature *)calloc(1,sizeof(BaseAlignFeature))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for base align feature\n");
    return NULL;
  }

  baf->objectType = CLASS_BASEALIGNFEATURE;

  Object_incRefCount(baf);

  baf->funcs = &baseAlignFeatureFuncs;

  return baf;
}

char *BaseAlignFeature_setCigarString(BaseAlignFeature *baf, char *str) {
  baf->cigarString = StrUtil_copyString(&(baf->cigarString),str,0);

  return baf->cigarString;
}

Vector *BaseAlignFeature_getUngappedFeatures(BaseAlignFeature *baf) {
  Vector *features;

  if (BaseAlignFeature_getCigarString(baf)) {
    features = BaseAlignFeature_parseCigar(baf);
    return features;
  } else {
    fprintf(stderr,"Error: No cigar_string defined.  Can't return ungapped features\n");
    exit(1);
  }
}

void BaseAlignFeature_reverseComplementImpl(BaseAlignFeature *baf) {
  char *newCigarString;
  char *oldCigarString = BaseAlignFeature_getCigarString(baf);
  char *oldChP;
  char *newChP;
  char *newPieceStartP;
  int len = strlen(oldCigarString);


  // reverse strand in both sequences
  BaseAlignFeature_setStrand(baf, BaseAlignFeature_getStrand(baf) * -1);
  BaseAlignFeature_setHitStrand(baf, BaseAlignFeature_getHitStrand(baf) * -1);

  if (!len) {
    return;
  }

  newCigarString = CigarStrUtil_reverse(oldCigarString,len);

// NIY free old cigar string
  BaseAlignFeature_setCigarString(baf, newCigarString);
}

Vector *BaseAlignFeature_parseCigar(BaseAlignFeature *baf) {
  int queryUnit = BaseAlignFeature_getQueryUnit();
  int hitUnit   = BaseAlignFeature_getHitUnit();
  char *string  = BaseAlignFeature_getCigarString(baf);
  Vector *features;
  int strand1;
  int strand2;
  int start1;
  int start2;
  char *chP;
  char piece[MAXCIGARPIECELEN];
  char *pieceP = piece;
 
  *pieceP = '\0';


  if (!string) {
    fprintf(stderr, "Error: No cigar string defined in object.  This should be caught"
                    " by the cigar_string method and never happen\n");
    exit(1);
  }

  features = Vector_new();

  Vector_setFreeFunc(features,Object_freeImpl);

  strand1 = BaseAlignFeature_getStrand(baf);
  if (!strand1) strand1 = 1;

  strand2 = BaseAlignFeature_getHitStrand(baf);
  if (!strand2) strand2 = 1;

  if (strand1 == 1 ) {
    start1 = BaseAlignFeature_getStart(baf);
  } else {
    start1 = BaseAlignFeature_getEnd(baf);
  }

  if (strand2 == 1 ) {
    start2 = BaseAlignFeature_getHitStart(baf);
  } else {
    start2 = BaseAlignFeature_getHitEnd(baf);
  }

  chP = string;

  while (*chP != '\0') {
    if (*chP == 'M' || *chP == 'D' || *chP == 'I') {
      float mappedLength;
      int length;

      *pieceP = '\0';
      if (pieceP == piece) { // No digits
        length = 1;
      } else {
        char *retP;
        length = strtol(piece, &retP, 10);
        if (retP != pieceP) {
          fprintf(stderr, "Error: failed parsing cigar string %s\n", string);
          exit(1);
        }
      }
  
      // explicit if statements to avoid rounding problems
      // and make sure we have sane coordinate systems
  
      if (queryUnit == 1 && hitUnit == 3 ) {
        mappedLength = length*3;
      } else if (queryUnit == 3 && hitUnit == 1 ) {
        mappedLength = length/3.0;
      } else if ( queryUnit == 1 && hitUnit == 1 ) {
        mappedLength = length;
      } else {
        fprintf(stderr, "Error: Internal error %d %d, currently only allowing 1 or 3\n", queryUnit, hitUnit);
        exit(1);
      }
  
      if ((int)mappedLength > mappedLength+0.00001 ||
          (int)mappedLength < mappedLength-0.00001) {
        fprintf(stderr, "Error: Internal error with mismapped length of hit, query "
                        "%d, hit %d, length %d",queryUnit,hitUnit,length);
        exit(1);
      }
  
      if (*chP == 'M') {
        FeaturePair *fp = FeaturePair_new();
        int a, b;
  
        if (strand1 == 1 ) {
          a = start1;
          b = start1 + length - 1;
          start1 = b + 1;
        } else {
          b = start1;
          a = start1 - length + 1;
          start1 = a - 1;
        }
  
        FeaturePair_setStart(fp,a);
        FeaturePair_setEnd(fp, b);
        FeaturePair_setStrand(fp, BaseAlignFeature_getStrand(baf));
        FeaturePair_setScore(fp,BaseAlignFeature_getScore(baf));
        // NIY fp->seqname($self->seqname);
        FeaturePair_setPhase(fp,BaseAlignFeature_getPhase(baf));
        FeaturePair_setpValue(fp,BaseAlignFeature_getpValue(baf));
        FeaturePair_setPercId(fp,BaseAlignFeature_getPercId(baf));
  
        if (strand2 == 1 ) {
          a = start2;
          b = start2 + mappedLength - 1;
          start2 = b + 1;
        } else {
          b = start2;
          a = start2 - mappedLength + 1;
          start2 = a - 1;
        }
  
        FeaturePair_setHitStart(fp,a);
        FeaturePair_setHitEnd(fp,b);
        FeaturePair_setHitStrand(fp, BaseAlignFeature_getHitStrand(baf));
  // NIY memory
        FeaturePair_setHitSeqName(fp, BaseAlignFeature_getHitSeqName(baf));
  
        FeaturePair_setContig(fp,BaseAlignFeature_getContig(baf));
        FeaturePair_setAnalysis(fp,BaseAlignFeature_getAnalysis(baf));
  
        Vector_addElement(features,fp);
  
        // end M cigar bits
      } else if (*chP == 'I' ) {
        if (strand1 == 1 ) {
          start1 += length;
        } else {
          start1 -= length;
        }
      } else if (*chP == 'D' ) {
        if (strand2 == 1 ) {
          start2 += mappedLength;
        } else {
          start2 -= mappedLength;
        }
      } else {
        fprintf(stderr, "Error: Illegal cigar line %s!\n", string);
        exit(1);
      }
      chP++;
      pieceP = piece;
      *pieceP = '\0';
    } else {
      *pieceP = *chP;
      chP++;
      pieceP++;
    }
  }

  if (piece[0] != '\0') {
    fprintf(stderr, "Error: Extraneous characters at end of cigar string |%s|\n",string);
    exit(1);
  }

  // should the features be sorted ?
  // 
  return features;
}


int BaseAlignFeature_parseFeatures(BaseAlignFeature *baf, Vector *features) {
  int queryUnit = BaseAlignFeature_getQueryUnit();
  int hitUnit   = BaseAlignFeature_getHitUnit();
  int f1Start;
  int f1End;
  int f2Start;
  int f2End;
  int strand;
  FeaturePair *fp;
  int hstrand;
  BaseContig *contig;
  char *name;
  char *hname;
  double score;
  double percent;
  double pvalue;
  Analysis *analysis;
  int phase;
  int ori;
  int prev1; // where last feature q part ended
  int prev2; // where last feature s part ended
  char *string = NULL;
  FeaturePair *firstFeature;
  FeaturePair *lastFeature;
  int i;
  char piece[MAXCIGARPIECELEN];

  if (!features) {
    fprintf(stderr,"Error: features must not be null in parseFeatures\n");
    exit(1);
  }

  if (Vector_getNumElement(features) == 0) {
    fprintf(stderr, "Error: No features in the array to parse");
    exit(1);
  }

  fp = Vector_getElementAt(features,0);

// NIY check whether strand defined should be checked here
  strand = FeaturePair_getStrand(fp);

  // 
  // Sort the features on their start position
  // Ascending order on positive strand, descending on negative strand
  // 
  if (strand == 1) {
    Vector_sort(features, SeqFeature_startCompFunc);
  } else {
    Vector_sort(features, SeqFeature_reverseStartCompFunc);
  }


  hstrand = FeaturePair_getHitStrand(fp);
  contig  = FeaturePair_getContig(fp);
  name    = FeaturePair_getSeqName(fp);
  hname   = FeaturePair_getHitSeqName(fp);
  score   = FeaturePair_getScore(fp);
  percent = FeaturePair_getPercId(fp);
  pvalue  = FeaturePair_getpValue(fp);
  analysis= FeaturePair_getAnalysis(fp);
  phase   = FeaturePair_getPhase(fp);

  // implicit strand 1 for peptide sequences
// NIY not using defined but 0
/*
  ( defined $strand ) || ( $strand = 1 );
  ( defined $hstrand ) || ( $hstrand = 1 );
*/
  if (!hstrand) hstrand = 1;
  if (!strand) strand = 1;

  ori = strand * hstrand;


  // Use strandedness info of query and hit to make sure both sets of
  // start and end  coordinates are oriented the right way around.

  firstFeature = Vector_getElementAt(features,0);
  lastFeature = Vector_getLastElement(features);

  if (strand == 1) {
    f1Start = FeaturePair_getStart(firstFeature); 
    f1End = FeaturePair_getEnd(lastFeature); 
  } else {
    f1Start = FeaturePair_getStart(lastFeature); 
    f1End = FeaturePair_getEnd(firstFeature); 
  }

  if (hstrand == 1) {
    f2Start = FeaturePair_getHitStart(firstFeature); 
    f2End = FeaturePair_getHitEnd(lastFeature); 
  } else {
    f2Start = FeaturePair_getHitStart(lastFeature); 
    f2End = FeaturePair_getHitEnd(firstFeature); 
  }

  // 
  // Loop through each portion of alignment and construct cigar string
  // 
   
  string = StrUtil_copyString(&string,"",0);

  for (i=0; i<Vector_getNumElement(features); i++) {
    FeaturePair *f = Vector_getElementAt(features,i);
    int start1;
    int start2;
    int length;
    int hlength;
    float hlengthfactor;
    int insertionFlag = 0;
    int matchlength;

    //
    // Sanity checks
    //

    if (!Class_isDescendent(CLASS_FEATUREPAIR, fp->objectType)) {
      fprintf(stderr, "Error: Array element [%p] is not a FeaturePair\n",fp);
      exit(1);
    }


// NIY    if ( defined $f->hstrand() ) {
      if (FeaturePair_getHitStrand(f) != hstrand) {
        fprintf(stderr, "Error: Inconsistent hstrands in feature array (%d and %d, feature number %d)\n",hstrand, FeaturePair_getHitStrand(f),i);
        exit(1);
      }
//    }
// NIY    if ( defined $f->strand() ) {
      if (FeaturePair_getStrand(f) != strand) {
        fprintf(stderr, "Error: Inconsistent strands in feature array\n");
        exit(1);
      }
//    }
    if (strcmp(name, FeaturePair_getSeqName(f))) {
      fprintf(stderr,"Error: Inconsistent names in feature array [%s - %s]\n", 
              name, FeaturePair_getSeqName(f));
      exit(1);
    }
    if (strcmp(hname,FeaturePair_getHitSeqName(f))) {
      fprintf(stderr,"Error: Inconsistent hnames in feature array [%s - %s]\n", 
              hname, FeaturePair_getHitSeqName(f));
      exit(1);
    }
    if (score != FeaturePair_getScore(f)) {
      fprintf(stderr,"Error: Inconsistent scores in feature array [%f - %f]\n", 
              score, FeaturePair_getScore(f));
      exit(1);
    }
// Note was a defined here
    if (percent != FeaturePair_getPercId(f)) {
      fprintf(stderr,"Error: Inconsistent pids in feature array [%f - %f]\n", 
              percent, FeaturePair_getPercId(f));
      exit(1);
    }

    start1 = FeaturePair_getStart(f);       // source sequence alignment start
    start2 = FeaturePair_getHitStart(f);    // hit sequence alignment start

    //  
    // More sanity checking
    // 
    if (i!=0) {
      if (strand == 1) {
        if (FeaturePair_getStart(f) < prev1) {
          fprintf(stderr,"Error: Inconsistent coordinates feature is forward strand "
                         "hstart in current feature should be greater than "
                         "hend in previous feature %d < %d\n",
                         FeaturePair_getStart(f),prev1);
          exit(1);
        }
      } else {
        if (FeaturePair_getEnd(f) > prev1) {
          fprintf(stderr,"Error: Inconsistent coordinates in feature array feature "
                         "is reverse strand hend should be less than previous "
                         "hstart %d > %d\n", FeaturePair_getEnd(f),prev1);
          exit(1);
        }
      }
    }

    length  = (FeaturePair_getEnd(f) - FeaturePair_getStart(f) + 1); // length of source seq alignment
    hlength = (FeaturePair_getHitEnd(f) - FeaturePair_getHitStart(f) + 1); // length of hit seq alignment

    // using multiplication to avoid rounding errors, hence the
    // switch from query to hit for the ratios

    // 
    // Yet more sanity checking
    // 
    if (queryUnit > hitUnit){
      // I am going to make the assumption here that this situation will
      // only occur with DnaPepAlignFeatures, this may not be true

      int query_p_length = (int)(length/queryUnit);
      int hit_p_length   = (int)(hlength * hitUnit);

      if (query_p_length != hit_p_length) {
        fprintf(stderr, "Error: Feature lengths not comparable Lengths: %d %d Ratios: %d %d\n",
                length, hlength, queryUnit, hitUnit );
        exit(1);
      }
    } else {
      if( length * hitUnit != hlength * queryUnit ) {
        fprintf(stderr, "Error: Feature lengths not comparable Lengths: %d %d Ratios: %d %d\n",
                length, hlength, queryUnit, hitUnit );
        exit(1);
      }
    }

    hlengthfactor = (queryUnit/hitUnit);


    // 
    // Check to see if there is an I type (insertion) gap:
    //   If there is a space between the end of the last source sequence
    //   alignment and the start of this one, then this is an insertion
    // 

    insertionFlag = 0;
    if (strand == 1) {
      if (i!=0 && ( FeaturePair_getStart(f) > prev1 + 1 )) {
        // there is an insertion
        int gap = FeaturePair_getStart(f) - prev1 - 1;

        insertionFlag = 1;
        if (gap == 1 ) {
          string = StrUtil_appendString(string,"I");
        } else {
          sprintf(piece,"%dI",gap);
          string = StrUtil_appendString(string,piece);
        }
      }

      // shift our position in the source seq alignment
      prev1 = FeaturePair_getEnd(f);
    } else {

      if (i!=0 && (FeaturePair_getEnd(f)+1 < prev1)) {
        int gap = prev1 - FeaturePair_getEnd(f) - 1;

        // there is an insertion
        insertionFlag = 1;
        if (gap == 1) {
          string = StrUtil_appendString(string,"I");
        } else {
          sprintf(piece,"%dI",gap);
          string = StrUtil_appendString(string,piece);
        }
      }

      // shift our position in the source seq alignment
      prev1 = FeaturePair_getStart(f);
    }

    //
    // Check to see if there is a D type (deletion) gap
    //   There is a deletion gap if there is a space between the end of the
    //   last portion of the hit sequence alignment and this one
    // 
    if (hstrand == 1) {
      if ((i!=0) && (FeaturePair_getHitStart(f) > prev2 + 1 )) {

        // there is a deletion
        float gap = FeaturePair_getHitStart(f) - prev2 - 1;
        int gap2 = (int)(gap * hlengthfactor + 0.05 );

        if (gap2 == 1 ) {
          string = StrUtil_appendString(string,"D");
        } else {
          sprintf(piece,"%dD",gap2);
          string = StrUtil_appendString(string,piece);
        }

        // sanity check,  Should not be an insertion and deletion
        if (insertionFlag) {
          fprintf(stderr, "Warning: Should not be an deletion and insertion on the "
                          "same alignment region. cigar_line=%s\n",string);
        }
      }
      // shift our position in the hit seq alignment
      prev2 = FeaturePair_getHitEnd(f);

    } else {
      if ((i!=0) && (FeaturePair_getHitEnd(f) + 1 < prev2 )) {

        // there is a deletion
        float gap = prev2 - FeaturePair_getHitEnd(f) - 1;
        int gap2 = (int)(gap * hlengthfactor + 0.05);

        if (gap2 == 1) {
          string = StrUtil_appendString(string,"D");
        } else {
          sprintf(piece,"%dD",gap2);
          string = StrUtil_appendString(string,piece);
        }

        // sanity check,  Should not be an insertion and deletion

        if(insertionFlag) {
          fprintf(stderr, "Error: Should not be an deletion and insertion on the "
                          "same alignment region. prev2 = %d f->hend() = %d; cigar_line = %s\n",
                           prev2, FeaturePair_getHitEnd(f),string); 
          exit(1);
        }
      }
      // shift our position in the hit seq alignment
      prev2 = FeaturePair_getHitStart(f);
    }

    matchlength = FeaturePair_getEnd(f) - FeaturePair_getStart(f) + 1;
    if (matchlength == 1) {
      string = StrUtil_appendString(string,"M");
    } else {
      sprintf(piece,"%dM",matchlength);
      string = StrUtil_appendString(string,piece);
    }
  }

  if (!score) {
    score = 1;
  }

// Done as two separate seq features which are set as feature1 and feature2 
  BaseAlignFeature_setStart(baf, f1Start);
  BaseAlignFeature_setEnd(baf, f1End);
  BaseAlignFeature_setStrand(baf, strand);
  BaseAlignFeature_setScore(baf, score);
  BaseAlignFeature_setPercId(baf, percent);
  BaseAlignFeature_setpValue(baf, pvalue);

  if (contig) {
    BaseAlignFeature_setContig(baf, contig);
  } else {
    BaseAlignFeature_setSeqName(baf, name);
  }
  BaseAlignFeature_setPhase(baf, phase);
  BaseAlignFeature_setAnalysis(baf, analysis);


  BaseAlignFeature_setHitStart(baf, f2Start);
  BaseAlignFeature_setHitEnd(baf, f2End);
  BaseAlignFeature_setHitStrand(baf, hstrand);
  BaseAlignFeature_setHitSeqName(baf, hname);

  BaseAlignFeature_setCigarString(baf, string);

  return 1;
}


Vector *BaseAlignFeature_transformSliceToRawContigImpl(BaseAlignFeature *baf) {
  Slice *slice = (Slice *)BaseAlignFeature_getContig(baf);
  SliceAdaptor *sa;
  RawContigAdaptor *rca;
  Vector *out;
  Vector *features;
  Vector *mappedFeatures;
  IDType *keys;
  IDHash *rcFeatures;
  int i;

  if (!slice){
    fprintf(stderr,"Error: Can't transform coordinates of %p "
                   " without slice attached to feature\n",baf);
    exit(1);
  }

  sa = (SliceAdaptor *)Slice_getAdaptor(slice);

  if (!sa) {
    fprintf(stderr, "Error: Can't transform %p without an adaptor attached " 
                    "to the feature's slice\n", baf);
    exit(1);
  }

  rca = DBAdaptor_getRawContigAdaptor(sa->dba);
  out = Vector_new();
  rcFeatures = IDHash_new(IDHASH_SMALL);

  // parse cigarline and split this gapped feature into list of ungapped features
  features = BaseAlignFeature_parseCigar(baf);

  mappedFeatures = Vector_new();

  // transform each of the ungapped features into raw contig coordinates
  for (i=0;i<Vector_getNumElement(features); i++) {
    FeaturePair *f = Vector_getElementAt(features,i);
    Vector *featVector = BaseAlignFeature_transformFeatureSliceToRawContig(baf, f);
    Vector_append(mappedFeatures, featVector);
    Vector_free(featVector);
  }

  // sort the transformed ungapped features into contig buckets
  for (i=0;i<Vector_getNumElement(mappedFeatures); i++) {
    FeaturePair *mf = Vector_getElementAt(mappedFeatures,i);
    IDType contigId = BaseContig_getDbID(FeaturePair_getContig(mf));
    Vector *contigFeatures;

    if (!IDHash_contains(rcFeatures, contigId)) {
      IDHash_add(rcFeatures,contigId,Vector_new());
    }
    contigFeatures = IDHash_getValue(rcFeatures, contigId);
    Vector_addElement(contigFeatures, mf);
  }

  // create a single gapped feature from all the ungapped features
  // in each contig bucket
  keys = IDHash_getKeys(rcFeatures);

  for (i=0; i<IDHash_getNumValues(rcFeatures); i++) {
    RawContig *contig;
    Vector *contigFeatures = IDHash_getValue(rcFeatures, keys[i]); 

    // create a gapped feature from a list of ungapped features
    BaseAlignFeature *outputf = (BaseAlignFeature *)SeqFeatureFactory_newFeature(baf->objectType);
    BaseAlignFeature_parseFeatures(outputf, contigFeatures );
    BaseAlignFeature_setAnalysis(outputf, BaseAlignFeature_getAnalysis(baf));
    BaseAlignFeature_setScore(outputf, BaseAlignFeature_getScore(baf));
    BaseAlignFeature_setPercId(outputf, BaseAlignFeature_getPercId(baf));
    BaseAlignFeature_setpValue(outputf, BaseAlignFeature_getpValue(baf));

    contig = RawContigAdaptor_fetchByDbID(rca, keys[i]);
    BaseAlignFeature_setContig(outputf, contig);
    Vector_addElement(out, outputf);
  }
  free(keys);

// NIY proper freeing
  IDHash_free(rcFeatures, NULL);

  return out;
}

int BaseAlignFeature_getHitUnitImpl(void) {
  fprintf(stderr, "Error: Abstract method call!" );
  exit(1);
  return 0;
}


int BaseAlignFeature_getQueryUnitImpl(void) {
  fprintf(stderr, "Error: Abstract method call!" );
  exit(1);
  return 0;
}

Vector *BaseAlignFeature_transformFeatureSliceToRawContig(BaseAlignFeature *baf, FeaturePair *fp) {
  int queryUnit = BaseAlignFeature_getQueryUnit();
  int hitUnit   = BaseAlignFeature_getHitUnit();
  Slice *slice = (Slice *)FeaturePair_getContig(fp);
  Slice *bafSlice = (Slice *)BaseAlignFeature_getContig(baf);
  SliceAdaptor *sa;
  RawContigAdaptor *rca;
  AssemblyMapperAdaptor *ama;
  AssemblyMapper *assMapper;
  int sliceStart;
  int sliceEnd;
  int sliceStrand;
  int globalStart;
  int globalEnd;
  int globalStrand;
  Vector *out;
  int hitStart;
  int hitEnd;
  MapperRangeSet *mapped;

  if (!slice) {
    fprintf(stderr, "Error: can't transform coordinates of %p"
                    " without some sort of contig defined\n",fp);
    exit(1);
  }

  sa = (SliceAdaptor *)Slice_getAdaptor(slice);

  if (!sa) {
    fprintf(stderr, "Error: can't tranform coordinates of %p without " 
                    "adaptor attached to feature's slice",fp);
    exit(1);
  }

  ama = DBAdaptor_getAssemblyMapperAdaptor(sa->dba);
  assMapper = AssemblyMapperAdaptor_fetchByType(ama, Slice_getAssemblyType(bafSlice));

  rca = DBAdaptor_getRawContigAdaptor(sa->dba);

  sliceStart  = Slice_getChrStart(slice);
  sliceEnd    = Slice_getChrEnd(slice);
  sliceStrand = Slice_getStrand(slice);

  //change feature coords from slice coordinates to assembly coords
  if (sliceStrand == 1) {
    globalStart  = FeaturePair_getStart(fp) + sliceStart - 1;
    globalEnd    = FeaturePair_getEnd(fp) + sliceStart - 1;
    globalStrand = FeaturePair_getStrand(fp);
  } else {
    globalStart  = sliceEnd - FeaturePair_getEnd(fp)   + 1;
    globalEnd    = sliceEnd - FeaturePair_getStart(fp) + 1;
    globalStrand = FeaturePair_getStrand(fp) * -1;
  }


  //convert assembly coords to raw contig coords
  mapped = AssemblyMapper_mapCoordinatesToRawContig(
     assMapper,
     Slice_getChrName(slice), //Slice_getChrId(slice),
     globalStart,
     globalEnd,
     globalStrand
    );

  if (!mapped || mapped->nRange == 0) {
    fprintf(stderr, "Error: couldn't map %p\n", fp);
    exit(1);
  }

  out = Vector_new();

  if (mapped->nRange > 1) {
    // The feature needs to be mapped accross multiple contigs
    RawContig *rawContig;
    MapperCoordinate *mc;
    MapperRange *mr;
    FeaturePair *newFeature;
    int codonUnusedBases = 0;
    int i;

    if (FeaturePair_getHitStrand(fp) == 1 ) {
      hitStart = FeaturePair_getHitStart(fp);
    } else {
      hitEnd = FeaturePair_getHitEnd(fp);
    }

    // split the feature into a seperate feature for each raw contig
    for (i=0; i<mapped->nRange; i++) {
      int hitLength;
      int queryLength;
      int queryStart;
      int queryEnd;

      mr = MapperRangeSet_getRangeAt(mapped,i);

      if (mr->rangeType == MAPPERRANGE_GAP){
        fprintf(stderr,"Warning: piece of evidence lies on gap\n");
        continue;
      }

      mc = (MapperCoordinate *)mr;

      // calculate query and hit length for each portion of the split feature
      // may need to round hit length to avoid 'partial' peptide

      // decision of Michele to not cover split codons

      queryLength = (mc->end - mc->start + 1);
      queryStart = mc->start;
      queryEnd = mc->end;


      if (queryUnit == hitUnit) {

        //  DnaDna and PepPep case
        hitLength = queryLength;

      } else if (queryUnit > hitUnit ){
        int usedBases = 0;

// NIY May be wrong
        hitLength = (int)((queryLength - codonUnusedBases ) / queryUnit);

        if (codonUnusedBases) {
          if (FeaturePair_getHitStrand(fp) == 1 ) {
            hitStart++;
          } else {
            hitEnd--;
          }
        }

        usedBases = queryLength - codonUnusedBases -
          hitLength * queryUnit;

        if (mc->strand == -1 ) {
          queryEnd -= codonUnusedBases;
          queryStart += usedBases;
        } else {
          queryStart += codonUnusedBases;
          queryEnd -= usedBases;
        }


        // new rest at the end ...
        if (usedBases) {
          codonUnusedBases = 3 - usedBases;
        } else {
          codonUnusedBases = 0;
        }

      } else if (hitUnit > queryUnit){

#ifdef DONE
        //  PepDnaAlign case (rare)
        my $tmp = ($query_length*$self->_hit_unit);
        $hit_length = sprintf "%.0f", $tmp; // round value up or down
#else 
        fprintf(stderr, "Error: PepDNAAlign not implemented - Nag Steve\n");
        exit(1);
#endif
      }

      if (hitLength == 0){
        continue;
      }

      if (FeaturePair_getHitStrand(fp) == 1 ) {
        hitEnd = (hitStart + hitLength) - 1;
      } else {
        hitStart = ( hitEnd - hitLength + 1 );
      }

      rawContig = RawContigAdaptor_fetchByDbID(rca, mc->id);

      // create the new feature
      newFeature = FeaturePair_new();

      FeaturePair_setStart(newFeature, queryStart);
      FeaturePair_setEnd(newFeature, queryEnd);
      FeaturePair_setStrand(newFeature, mc->strand);
      FeaturePair_setScore(newFeature,FeaturePair_getScore(fp));
      FeaturePair_setPercId(newFeature,FeaturePair_getPercId(fp));
      FeaturePair_setpValue(newFeature,FeaturePair_getpValue(fp));
      FeaturePair_setHitStart(newFeature, hitStart);
      FeaturePair_setHitEnd(newFeature, hitEnd);
      FeaturePair_setHitStrand(newFeature,FeaturePair_getHitStrand(fp));
      FeaturePair_setHitSeqName(newFeature,FeaturePair_getHitSeqName(fp));
      FeaturePair_setAnalysis(newFeature,FeaturePair_getAnalysis(fp));
      FeaturePair_setContig(newFeature,rawContig);

      Vector_addElement(out, newFeature);

      if (FeaturePair_getHitStrand(fp) == 1 ) {
        hitStart = (hitEnd + 1);
      } else {
        hitEnd = hitStart -1;
      }
    }
  } else {
    // feature maps to single contig
    RawContig *rawContig;
    FeaturePair *newFeature;
    MapperCoordinate *mc;
    MapperRange *mr;

    mr = MapperRangeSet_getRangeAt(mapped,0);

    if (mr->rangeType == MAPPERRANGE_GAP) {
      fprintf(stderr, "Warning: piece of evidence lies on gap\n");
// NIY Not sure what to return
      return emptyVector;
    }

    mc = (MapperCoordinate *)mr;

    // create the new feature
    rawContig = RawContigAdaptor_fetchByDbID(rca, mc->id);
    newFeature = FeaturePair_new();

    FeaturePair_setStart(newFeature,mc->start);
    FeaturePair_setEnd(newFeature,mc->end);
    FeaturePair_setStrand(newFeature,mc->strand);
    FeaturePair_setScore(newFeature,FeaturePair_getScore(fp));
    FeaturePair_setPercId(newFeature,FeaturePair_getPercId(fp));
    FeaturePair_setpValue(newFeature,FeaturePair_getpValue(fp));
    FeaturePair_setHitStart(newFeature,FeaturePair_getHitStart(fp));
    FeaturePair_setHitEnd(newFeature,FeaturePair_getHitEnd(fp));
    FeaturePair_setHitStrand(newFeature,FeaturePair_getHitStrand(fp));
    FeaturePair_setHitSeqName(newFeature,FeaturePair_getHitSeqName(fp));
    FeaturePair_setAnalysis(newFeature,FeaturePair_getAnalysis(fp));
    FeaturePair_setContig(newFeature,rawContig);

    Vector_addElement(out, newFeature);
  }

  return out;
}

void BaseAlignFeature_freePtrs(BaseAlignFeature *baf) {
  if (baf->cigarString) free(baf->cigarString);

  FeaturePair_freePtrs((FeaturePair *)baf);
}


void BaseAlignFeature_freeImpl(BaseAlignFeature *baf) {
  Object_decRefCount(baf);

  if (Object_getRefCount(baf) > 0) {
    return;
  } else if (Object_getRefCount(baf) < 0) {
    fprintf(stderr,"Error: Negative reference count for BaseAlignFeature\n"
                   "       Freeing it anyway\n");
  }

  BaseAlignFeature_freePtrs(baf);

  free(baf);
}

