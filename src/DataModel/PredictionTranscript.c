#define __PREDICTIONTRANSCRIPT_MAIN__
#include "PredictionTranscript.h"
#undef  __PREDICTIONTRANSCRIPT_MAIN__


#include "Exon.h"
#include "StrUtil.h"
#include "SeqUtil.h"

PredictionTranscript *PredictionTranscript_new() {
  PredictionTranscript *transcript;

  if ((transcript = (PredictionTranscript *)calloc(1,sizeof(PredictionTranscript))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for transcripe\n");
    return NULL;
  }

  transcript->objectType = CLASS_PREDICTIONTRANSCRIPT;
  transcript->funcs = &predictionTranscriptFuncs;
  transcript->exons = Vector_new();
  return transcript;
}

char *PredictionTranscript_setType(PredictionTranscript *t, char *type) {
  if ((t->type = (char *)malloc(strlen(type)+1)) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for transcript type\n");
    return NULL;
  }

  strcpy(t->type,type);

  return t->type;
}

void PredictionTranscript_flushExons(PredictionTranscript *trans) {
  PredictionTranscript_setCodingRegionStartIsSet(trans, 0);  
  PredictionTranscript_setCodingRegionEndIsSet(trans, 0);  
  PredictionTranscript_setStartIsSet(trans, 0);  
  PredictionTranscript_setEndIsSet(trans, 0);  
  if (trans->translateableExons) {
    Vector_free(trans->translateableExons,NULL);
    trans->translateableExons = NULL;
  }
// NIY  Transcript_removeAllExons(trans);
// NIY caches 
}

void PredictionTranscript_free(PredictionTranscript *trans) {
  printf("PredictionTranscript_free not implemented\n");
}

int PredictionTranscript_getLength(PredictionTranscript *trans) {
    int length = 0;
    int i;

    for (i=0;i<PredictionTranscript_getExonCount(trans); i++) {
      Exon *ex = PredictionTranscript_getExonAt(trans,i);
 // Check stickies
      if (ex) length += Exon_getLength(ex);
    }
    return length;
}

int PredictionTranscript_setCodingRegionStart(PredictionTranscript *trans, int start) {
  PredictionTranscript_setCodingRegionStartIsSet(trans, 1);  
  trans->codingRegionStart = start;
  return trans->codingRegionStart;
}
  
int PredictionTranscript_getCodingRegionStart(PredictionTranscript *trans) {

  if (!PredictionTranscript_getCodingRegionStartIsSet(trans)) {
    //if the coding start is not defined, use the start of the transcript
    return PredictionTranscript_getStart(trans);
  }

  return trans->codingRegionStart;
}

int PredictionTranscript_setCodingRegionEnd(PredictionTranscript *trans, int end) {
  PredictionTranscript_setCodingRegionEndIsSet(trans, 1);  
  trans->codingRegionEnd = end;
  return trans->codingRegionEnd;
}
  
int PredictionTranscript_getCodingRegionEnd(PredictionTranscript *trans) {

  if (!PredictionTranscript_getCodingRegionEndIsSet(trans)) {
    //if the coding end is not defined, use the end of the transcript
    return PredictionTranscript_getEnd(trans);
  }

  return trans->codingRegionEnd;
}

int PredictionTranscript_setStart(PredictionTranscript *trans, int start) {
  PredictionTranscript_setStartIsSet(trans, 1);  
  trans->start = start;
  return trans->start;
}
  
int PredictionTranscript_getStart(PredictionTranscript *trans) {

  if (!PredictionTranscript_getStartIsSet(trans)) {
    //if the coding start is not defined, use the start of the transcript
    fprintf(stderr, "Error: Trying to call PT getEnd before its set\n");
    exit(1);
  }

  return trans->start;
}

int PredictionTranscript_setEnd(PredictionTranscript *trans, int end) {
  PredictionTranscript_setEndIsSet(trans, 1);  
  trans->end = end;
  return trans->end;
}
  
int PredictionTranscript_getEnd(PredictionTranscript *trans) {
  if (!PredictionTranscript_getEndIsSet(trans)) {
    //if the coding end is not defined, use the end of the transcript
    fprintf(stderr, "Error: Trying to call PT getEnd before its set\n");
    exit(1);
  }

  return trans->end;
}

/*
=head2 add_Exon

  Arg  1    : Bio::EnsEMBL::Exon $exon
  Arg [2]   : int $exon_position
              Use it when you know you dont have exons in the
              beginning. Useful for the Adaptor when retrieving 
              partial PT..
  Function  : Adds given Exon to this prediction transcript. 
              It can be at arbitrary position in the array. Not filled lower
              positions in the exon list are set undef then. Counting starts at 1.
  Returntype: none
  Exceptions: if argument is not Bio::EnsEMBL::Exon
  Caller    : Pipeline runnable Genscan

=cut
*/

void PredictionTranscript_addExon(PredictionTranscript *trans, Exon *exon, int *positionP) {
  if (positionP) {
    Vector_setElementAt(trans->exons, *positionP, exon);
  } else {
    Vector_addElement(trans->exons, exon);
  }

  if (exon && (!PredictionTranscript_getStartIsSet(trans) ||
		Exon_getStart(exon) < PredictionTranscript_getStart(trans))) {
    PredictionTranscript_setStart(trans, Exon_getStart(exon));
  }
  if (exon && (!PredictionTranscript_getEndIsSet(trans) ||
		Exon_getEnd(exon) > PredictionTranscript_getEnd(trans))) {
    PredictionTranscript_setEnd(trans, Exon_getEnd(exon));
  }
}

/*
=head2 get_all_Exons

  Arg [1]   : optional 1 $wish_undefined_exons 
  Function  : Returns all Exons currently in the PredictionTranscript
              in the order 5' to 3'. If this is a partial PredictionTranscript,
              elements of the list will be undef.
  Returntype: listref Bio::EnsEMBL::Exon
  Exceptions: none
  Caller    : self->get_cdna(),Web for display.

=cut
*/

Vector *PredictionTranscript_getAllExons(PredictionTranscript *trans, int wishUndefinedExon) {
   if (wishUndefinedExon) {
     return trans->exons;
   } else {
     return PredictionTranscript_getAllTranslateableExons(trans);
   }
}

/*
=head2 get_all_translateable_Exons

  Arg [1]    : none
  Example    : $exons = $self->get_all_translateable_Exons
  Description: Retreives the same value of get_all_Exons for this prediction
               transcript with the exception that undefined exons (only when
               transcript is in slice coords and exon maps to gap) are not
               returned.  In a prediction transcript there is no UTR and
               thus all exons are entirely translateable.
  Returntype : listref of Bio::EnsEMBL::Exon
  Exceptions : none
  Caller     : general

=cut
*/

Vector *PredictionTranscript_getAllTranslateableExons(PredictionTranscript *trans) {
  int i;

  if (!trans->translateableExons) {
    trans->translateableExons = Vector_new();
    for (i=0; i<PredictionTranscript_getExonCount(trans); i++) {
      Exon *ex = PredictionTranscript_getExonAt(trans,i);
      if (ex) {
        Vector_addElement(trans->translateableExons, ex);
      }
    }
  }

  return trans->translateableExons;
}


/*
=head2 sort

 Function: Sorts the exons by start coordinate
           Sorts forward for forward strand and reverse for reverse strand
           It refills $self->{'exons'} with the sorted exons
 Returns : none
 Args    : none

=cut
*/

void PredictionTranscript_sort(PredictionTranscript *trans) {
  Exon *firstExon = PredictionTranscript_getExonAt(trans,0);
  int strand;
  int i;

  for (i=0;i<PredictionTranscript_getExonCount(trans);i++) {
    if (!PredictionTranscript_getExonAt(trans,i)) {
      return;
    }
  }

  strand = Exon_getStrand(firstExon);

  if (strand == 1) {
    Vector_sort(PredictionTranscript_getExons(trans),Exon_forwardStrandCompFunc);
  } else if (strand == -1) {
    Vector_sort(PredictionTranscript_getExons(trans),Exon_reverseStrandCompFunc);
  }
}

/*
=head2 get_exon_count

  Args      : none
  Function  : How many exons are in this PTranscript. Some might 
              not be in this object, depending on how it was retrieved.
	      (non golden exons missing on Slice->get_predicitonTranscripts()) 
  Returntype: int
  Exceptions: none
  Caller    : general

=cut
*/

int PredictionTranscript_getExonCount(PredictionTranscript *trans) {
  return Vector_getNumElement(trans->exons);
}

/*
=head2 set_exon_count

  Arg 1     : int $number_of_exons
              If the number of exons you put in with add_exon is not the 
              real number of exons in the Transcript, you can set it here.
              Might be necessary in db gets, where you dont get all.
  Function  : sets number of exons, so get_all_Exons returns more undef exons
              at the end. After this, you have to use position argument when you
	      want to insert Exons.
  Returntype: none
  Exceptions: If you set less then already in, you loose what you have :)
  Caller    : $self->adaptor()

=cut
*/

int PredictionTranscript_setExonCount(PredictionTranscript *trans, int count) {

  if (Vector_getNumElement(trans->exons) > count) {
    fprintf(stderr, "Error: Trying to shrink exon vector\n");
  }
  Vector_setNumElement(trans->exons, count);
}

/*
=head2 translate

  Args      : none
  Function  : Give a peptide translation of all exons currently in
              the PT. Gives empty string when none is in.
  Returntype: a Bio::Seq as in transcript->translate()
  Exceptions: if the exons come in two or more groups, with an undef exon
              in the middle, only the first group is translated.
  Caller    : general

=cut
*/

char *PredictionTranscript_translate(PredictionTranscript *trans) {
  char *dna;
  int lenDNA;
  int lengths[6];
  char *frm[6];
  int i;

  dna = PredictionTranscript_getcDNA(trans);

  lenDNA = strlen(dna);

// NIY Case insensitive
  if (lenDNA>=3) {
    char *lastCodon = &(dna[lenDNA-3]);
    if (!strcmp(lastCodon,"TAG") ||
        !strcmp(lastCodon,"TGA") ||
        !strcmp(lastCodon,"TAA")) {
      *lastCodon = '\0';
    }
  }

  for (i=0;i<6;i++) {
    frm[i] = (char *)malloc(lenDNA/3 + 2);
  }

  translate(dna, frm, lengths);

  free(dna);

  // Note don't free frame 0
  for (i=1;i<6;i++) {
    free(frm[i]);
  }

  return frm[0];
}

/*
=head2 get_cdna

  Args      : none
  Function  : Give a concat cdna of all exons currently in
              the PT. Gives empty string when none is in. Pads between not  
              phase matching exons. Builds internal coord translation table.
  Returntype: txt
  Exceptions: if the exons come in two or more groups, with an undef exon
              in the middle, only the first groups cdna is returned.
  Caller    : general, $self->translate()

=cut
*/

char *PredictionTranscript_getcDNA(PredictionTranscript *trans) {
  Vector *exons = PredictionTranscript_getAllExons(trans,0);
  char *cdna = StrUtil_copyString(&cdna, "", 0);
  int lastPhase = 0;
  int i;

  int cdnaStart, cdnaEnd;
  int pepStart, pepEnd;
  int newcdna, pepCount;

  cdnaStart = 1;
  pepStart = 1;

  for (i=0; i<Vector_getNumElement(exons); i++) {
    Exon *exon = Vector_getElementAt(exons, i);
    int phase;
    if (!exon) {
      if (cdna[0] == '\0') {
	continue;
      } else {
	break;
      }
    } 

    phase = 0;

// NIY    if (defined($exon->phase)) {
      phase = Exon_getPhase(exon);
//    }

    if (phase != lastPhase ) {

      if (lastPhase == 1) {
	cdna = StrUtil_appendString(cdna,"NN");
      } else if (lastPhase == 2) {
	cdna = StrUtil_appendString(cdna,"N");
      }

      // startpadding for this exon
      cdna = SeqUtil_addNs(cdna,phase);
    }
    
    cdna = StrUtil_appendString(cdna, Exon_getSeqString(exon));
    lastPhase = Exon_getEndPhase(exon);
  }

// NIY Freeing exons vector?
  return cdna;
}

/*
=head1 pep2genomic

  Arg  1   : integer start - relative to peptide
  Arg  2   : integer end   - relative to peptide

  Function : Provides a list of Bio::EnsEMBL::SeqFeatures which
             is the genomic coordinates of this start/end on the peptide

  Returns  : list of Bio::EnsEMBL::SeqFeature

=cut
*/

MapperRangeSet *PredictionTranscript_pep2Genomic(PredictionTranscript *trans, int start, int end) {
  Exon *firstExon = PredictionTranscript_getExonAt(trans,0);

  // move start end into translate cDNA coordinates now.
  // much easier!
  start = 3*start-2;
  end   = 3*end;

  //
  // Adjust the phase
  //
  if (firstExon) {
    start -= Exon_getPhase(firstExon);
    end   -= Exon_getPhase(firstExon);
  }

  return PredictionTranscript_cDNA2Genomic(trans, start, end);
}

MapperRangeSet *PredictionTranscript_cDNA2Genomic(PredictionTranscript *trans, int start, int end) {
  Mapper *mapper;

  mapper = PredictionTranscript_getcDNACoordMapper(trans);

  return Mapper_mapCoordinates(mapper, (IDType)trans, start, end, 1, CDNA_COORDS );
}

MapperRangeSet *PredictionTranscript_genomic2cDNA(PredictionTranscript *trans, int start, int end, int strand, BaseContig *contig) {
  Mapper *mapper;

  // "ids" in mapper are contigs of exons, so use the same contig that should
  // be attached to all of the exons...
  if (!contig) {
    Vector *translateable = PredictionTranscript_getAllTranslateableExons(trans);
    Exon *firstExon;
    if (!Vector_getNumElement(translateable)) {
      return MapperRangeSet_new();
    }
    firstExon = Vector_getElementAt(translateable, 0);
    contig = Exon_getContig(firstExon);
    Vector_free(translateable, NULL);
  }

  mapper = PredictionTranscript_getcDNACoordMapper(trans);

  return Mapper_mapCoordinates(mapper,(IDType)contig, start, end, strand, GENOMIC_COORDS);
}

Mapper *PredictionTranscript_getcDNACoordMapper(PredictionTranscript *trans) {
  Mapper *mapper;
  int start = 1;
  int i;
  Vector *translateable;

  if (trans->exonCoordMapper) {
    return trans->exonCoordMapper;
  }

  //
  // the mapper is loaded with OBJECTS in place of the IDs !!!!
  //  the objects are the contigs in the exons
  //
  mapper = Mapper_new( CDNA_COORDS, GENOMIC_COORDS );

  translateable = PredictionTranscript_getAllTranslateableExons(trans);
  for (i=0; i<Vector_getNumElement(translateable); i++) {
    Exon *exon = Vector_getElementAt(translateable,i);

    Exon_loadGenomicMapper(exon, mapper, (IDType)trans, start);
    start += Exon_getLength(exon);
  }
  trans->exonCoordMapper = mapper;
  Vector_free(translateable,NULL);
  return mapper;
}
