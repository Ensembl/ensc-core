#include "Transcript.h"
#include "TranscriptAdaptor.h"
#include "TranslationAdaptor.h"
#include "DBEntryAdaptor.h"
#include "StrUtil.h"
#include "SeqUtil.h"
#include "translate.h"

#include <stdlib.h>

Transcript *Transcript_new() {
  Transcript *transcript;

  if ((transcript = (Transcript *)calloc(1,sizeof(Transcript))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for transcripe\n");
    return NULL;
  }

  Transcript_setVersion(transcript,-1);

  transcript->objectType = CLASS_TRANSCRIPT;

  return transcript;
}

Vector *Transcript_getAllDBLinks(Transcript *t) {
  if (!t->dbLinks) {
    TranscriptAdaptor *ta = (TranscriptAdaptor *)Transcript_getAdaptor(t);

    if (ta) {
      DBEntryAdaptor *dbea = DBAdaptor_getDBEntryAdaptor(ta->dba);
      DBEntryAdaptor_fetchAllByTranscript(dbea,t);
    } else {
      t->dbLinks = emptyVector;
    }
  }

  return t->dbLinks;
}

int Transcript_addDBLink(Transcript *t, DBEntry *dbe) {
  if (!t->dbLinks) {
    t->dbLinks = Vector_new();
  }

  Vector_addElement(t->dbLinks, dbe); 
  return 1;
}

char *Transcript_setType(Transcript *t, char *type) {
  if ((t->type = (char *)malloc(strlen(type)+1)) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for transcript type\n");
    return NULL;
  }

  strcpy(t->type,type);

  return t->type;
}

char *Transcript_getStableId(Transcript *transcript) {
  TranscriptAdaptor *ta = (TranscriptAdaptor *)Transcript_getAdaptor(transcript);

  if (StableIdInfo_getStableId(&(transcript->si)) == NULL && ta) {
    TranscriptAdaptor_getStableEntryInfo(ta,transcript);
  }
  return StableIdInfo_getStableId(&(transcript->si));
}

int Transcript_getVersion(Transcript *transcript) {
  TranscriptAdaptor *ta = (TranscriptAdaptor *)Transcript_getAdaptor(transcript);

  if (StableIdInfo_getVersion(&(transcript->si)) == -1 && ta) {
    TranscriptAdaptor_getStableEntryInfo(ta,transcript);
  }
  return StableIdInfo_getVersion(&(transcript->si));
}

Transcript *Transcript_transform(Transcript *trans, IDHash *exonTransforms) {
  Vector *mappedExonVector = Vector_new();
  int i;
  
  for (i=0;i<Transcript_getExonCount(trans);i++) {
    Exon *exon = (Exon *)Transcript_getExonAt(trans,i);
    IDType exonRef = (IDType)exon;

    // the old exon was successfully remapped then store the new exon
/* CHECK */
    if ( IDHash_contains(exonTransforms,exonRef)) {
      Vector_addElement(mappedExonVector,IDHash_getValue(exonTransforms, exonRef));
    }
    // but for the case where the exon was unable to be mapped, as it
    // was outside the bounds of the slice, include the original exon.
    else {
      Vector_addElement(mappedExonVector,exon);
    }
  }

  //Flush the exons and NIY all related internal caches
  Transcript_flushExons(trans);

  // attach the new list of exons to the transcript
  for (i=0; i<Vector_getNumElement(mappedExonVector); i++) {
    Transcript_addExon(trans,(Exon *)Vector_getElementAt(mappedExonVector,i));
  }

  if ( Transcript_getTranslation(trans)) {
    Translation_transform(Transcript_getTranslation(trans), exonTransforms);
  }

  Vector_free(mappedExonVector,NULL);

  return trans;
}

void Transcript_flushExons(Transcript *trans) {
  Transcript_removeAllExons(trans);
// NIY caches 
}

int Transcript_setStart(Transcript *trans, int start) {
  trans->start = start;
  Transcript_setStartIsSet(trans,TRUE);
  return trans->start;
}

int Transcript_getStart(Transcript *trans) {
  if (!Transcript_getStartIsSet(trans)) {
    int start;
    int strand = Exon_getStrand(Transcript_getStartExon(trans));

    if (strand == 1) {
      start = Exon_getStart(Transcript_getStartExon(trans));
    } else {
      start = Exon_getStart(Transcript_getEndExon(trans));
    }
    trans->start = start;
    Transcript_setStartIsSet(trans,TRUE);
  }

  return trans->start;
}


int Transcript_setEnd(Transcript *trans, int end) {
  trans->end = end;
  Transcript_setEndIsSet(trans,TRUE);
  return trans->end;
}

int Transcript_getEnd(Transcript *trans) {
  if (!Transcript_getEndIsSet(trans)) {
    int end;
    int strand = Exon_getStrand(Transcript_getStartExon(trans));

    if (strand == 1) {
      end = Exon_getEnd(Transcript_getEndExon(trans));
    } else {
      end = Exon_getEnd(Transcript_getStartExon(trans));
    }
    trans->end = end;
    Transcript_setEndIsSet(trans,TRUE);
  }

  return trans->end;
}

int Transcript_setCodingRegionEnd(Transcript *trans, int end) {
  trans->codingRegionEnd = end;
  Transcript_setCodingRegionEndIsSet(trans,TRUE);
  return trans->codingRegionEnd;
}

int Transcript_getCodingRegionEnd(Transcript *trans) {
  Translation *translation = Transcript_getTranslation(trans);

  if (!Transcript_getCodingRegionEndIsSet(trans) && translation) {
    int end;
    int strand = Exon_getStrand(Translation_getStartExon(translation));

    if (strand == 1) {
      end = Exon_getStart(Translation_getEndExon(translation));
      end += (Translation_getEnd(translation) - 1);
    } else {
      end = Exon_getEnd(Translation_getStartExon(translation));
      end -= (Translation_getStart(translation) - 1 );
    }
    Transcript_setCodingRegionEnd(trans,end);
  }

  return trans->codingRegionEnd;
}

int Transcript_setCodingRegionStart(Transcript *trans, int start) {
  trans->codingRegionStart = start;
  Transcript_setCodingRegionStartIsSet(trans,TRUE);
  return trans->codingRegionStart;
}

int Transcript_getCodingRegionStart(Transcript *trans) {
  Translation *translation = Transcript_getTranslation(trans);

  if (!Transcript_getCodingRegionStartIsSet(trans) && translation) {
    int start;
    int strand = Exon_getStrand(Translation_getStartExon(translation));

    if (strand == 1) {
      start = Exon_getStart(Translation_getStartExon(translation));
      start += (Translation_getStart(translation) - 1);
    } else {
      start = Exon_getEnd(Translation_getEndExon(translation));
      start -= (Translation_getEnd(translation) - 1 );
    }
    Transcript_setCodingRegionStart(trans,start);
  }

  return trans->codingRegionStart;
}

void Transcript_sort(Transcript *trans) {
  Exon *firstExon = Transcript_getExonAt(trans,0);
  int strand = Exon_getStrand(firstExon);

  if (strand == 1) {
    qsort(Transcript_getExons(trans),Transcript_getExonCount(trans),
          sizeof(void *),Exon_forwardStrandCompFunc);
  } else if (strand == -1) {
    qsort(Transcript_getExons(trans),Transcript_getExonCount(trans),
          sizeof(void *),Exon_reverseStrandCompFunc);
  }
}

Exon *Transcript_getStartExon(Transcript *trans) {
  return Transcript_getExonAt(trans,0); 
}

Exon *Transcript_getEndExon(Transcript *trans) {
  return Transcript_getExonAt(trans,Transcript_getExonCount(trans)-1); 
}

Translation *Transcript_getTranslation(Transcript *trans) {
  if (!trans->translation && trans->translationId ) {
    TranslationAdaptor *ta = DBAdaptor_getTranslationAdaptor(Transcript_getAdaptor(trans)->dba);
    
    trans->translation = TranslationAdaptor_fetchByDbID(ta, trans->translationId, trans);
  }
  return trans->translation;
}


char *Transcript_getSplicedSeq(Transcript *trans) {
  int i;
  char *seqStr;

  
  StrUtil_copyString(&seqStr, "", 0);

  for (i=0; i<Transcript_getExonCount(trans); i++) {
    Exon *ex = Transcript_getExonAt(trans, i);
    StrUtil_appendString(seqStr,Exon_getSeqString(ex));
  }

  return seqStr;
}

char *Transcript_getTranslateableSeq(Transcript *trans) {
  int i;
  int lastPhase = 0;
  int first = 1;
  char *mrna = StrUtil_copyString(&mrna,"",0);
  Vector *translateableExons;

  translateableExons = Transcript_getAllTranslateableExons(trans);
 
  for (i=0; i<Vector_getNumElement(translateableExons); i++) {
    Exon *exon = Vector_getElementAt(translateableExons, i);
    int phase = 0;

// NIY    if (defined($exon->phase)) {
      phase = Exon_getPhase(exon);
//    }

    // startpadding is needed if MONKEY_EXONS are on
    if (first && (!getenv("MONKEY_EXONS"))) {
      mrna = SeqUtil_addNs(mrna,phase);
      first = 0;
    }

    if (phase != lastPhase && getenv("MONKEY_EXONS")) {
      // endpadding for the last exon
      if (lastPhase == 1 ) {
        mrna = StrUtil_appendString(mrna,"NN");
      } else if (lastPhase == 2) {
        mrna = StrUtil_appendString(mrna,"N");
      }
      //startpadding for this exon
      mrna = SeqUtil_addNs(mrna,phase);
    }
    mrna = StrUtil_appendString(mrna, Exon_getSeqString(exon));
    lastPhase = Exon_getEndPhase(exon);
  }

  Vector_free(translateableExons,NULL);

// NIY free adjusted

  return mrna;
}


/*
=head2 cdna_coding_start

  Arg [1]    : (optional) $value
  Example    : $relative_coding_start = $transcript->cdna_coding_start;
  Description: Retrieves the position of the coding start of this transcript
               in cdna coordinates (relative to the start of the 5prime end of
               the transcript, excluding introns, including utrs).
  Returntype : int
  Exceptions : none
  Caller     : five_prime_utr, get_all_snps, general

=cut
*/

int Transcript_setcDNACodingStart(Transcript *trans, int start) {
  trans->cDNACodingStart = start;
  Transcript_setcDNACodingStartIsSet(trans,TRUE);
  return trans->cDNACodingStart;
}

int Transcript_getcDNACodingStart(Transcript *trans) {
  Translation *translation = Transcript_getTranslation(trans);

  if (!Transcript_getcDNACodingStartIsSet(trans) && translation) {
    int start = 0;
    int i;

    for (i=0; i<Transcript_getExonCount(trans); i++) {
      Exon *exon = Transcript_getExonAt(trans,i);

      if (exon == Translation_getStartExon(translation)) {
        //add the utr portion of the start exon
        start += Translation_getStart(translation);
        break;
      } else {
        //add the entire length of this non-coding exon
        start += Exon_getLength(exon);
      }
    }
    Transcript_setcDNACodingStart(trans,start);
  }

  return trans->cDNACodingStart;
}

/*
=head2 cdna_coding_end

  Arg [1]    : (optional) $value
  Example    : $cdna_coding_end = $transcript->cdna_coding_end;
  Description: Retrieves the end of the coding region of this transcript in
               cdna coordinates (relative to the five prime end of the
               transcript, excluding introns, including utrs).
               Note
  Returntype : none
  Exceptions : none
  Caller     : general

=cut
*/

int Transcript_setcDNACodingEnd(Transcript *trans, int end) {
  trans->cDNACodingEnd = end;
  Transcript_setcDNACodingEndIsSet(trans,TRUE);
  return trans->cDNACodingEnd;
}

int Transcript_getcDNACodingEnd(Transcript *trans) {
  Translation *translation = Transcript_getTranslation(trans);

  if (!Transcript_getcDNACodingEndIsSet(trans) && translation) {
    int end = 0;
    int i;

    for (i=0; i<Transcript_getExonCount(trans); i++) {
      Exon *exon = Transcript_getExonAt(trans,i);

      if (exon == Translation_getEndExon(translation)) {
        //add the coding portion of the final coding exon
        end += Translation_getEnd(translation);
        break;
      } else {
        //add the entire length
        end += Exon_getLength(exon);
      }
    }
    Transcript_setcDNACodingEnd(trans,end);
  }

  return trans->cDNACodingEnd;
}

int Transcript_getLength(Transcript *trans) {
  int length = 0;
  int i;

  for (i=0; i<Transcript_getExonCount(trans); i++) {
    Exon *ex = Transcript_getExonAt(trans, i);
    length += Exon_getLength(ex);
  }
  return length;
}


Vector *Transcript_getAllTranslateableExons(Transcript *trans) {
  Vector *translateable;
  Translation *translation;
  Exon *startExon;
  Exon *endExon;
  int tlnStart;
  int tlnEnd;
  int i;
  

  if (!(translation = Transcript_getTranslation(trans))) {
    fprintf(stderr, "Error: No translation attached to transcript object\n");
    exit(1);
  }
  startExon      = Translation_getStartExon(translation);
  endExon        = Translation_getEndExon(translation);
  tlnStart       = Translation_getStart(translation);
  tlnEnd         = Translation_getEnd(translation);

  translateable = Vector_new();

  for (i=0; i<Transcript_getExonCount(trans); i++) {
    Exon *ex = Transcript_getExonAt(trans,i);
    int length;
    int adjustStart = 0;
    int adjustEnd   = 0;

    if (ex != startExon && !Vector_getNumElement(translateable)) {
      continue;   // Not yet in translated region
    }

    length  = Exon_getLength(ex);

    // Adjust to translation start if this is the start exon
    if (ex == startExon ) {
      if (tlnStart < 1 || tlnStart > length) {
        fprintf(stderr,"Error: Translation start %d is outside exon %p length=%d\n", tlnStart,ex,length);
        exit(1);
      }
      adjustStart = tlnStart - 1;
    }

    // Adjust to translation end if this is the end exon
    if (ex == endExon) {
      if (tlnEnd < 1 || tlnEnd > length) {
        fprintf(stderr,"Error: Translation start %d is outside exon %p length=%d\n", tlnEnd,ex,length);
        exit(1);
      }
      adjustEnd = tlnEnd - length;
    }

    // Make a truncated exon if the translation start or
    // end causes the coordinates to be altered.
    if (adjustEnd || adjustStart) {
      Exon *newex = Exon_adjustStartEnd(ex, adjustStart, adjustEnd);

      Vector_addElement(translateable, newex );
    } else {
      Vector_addElement(translateable, ex);
    }

    // Exit the loop when we've found the last exon
    if (ex == endExon) break;
  }
  return translateable;
}


char *Transcript_translate(Transcript *trans) {
  char *mRNA;
  int lenmRNA;
  int lengths[6];
  char *frm[6];
  int i;

  mRNA = Transcript_getTranslateableSeq(trans);

  lenmRNA = strlen(mRNA);

  if (lenmRNA%3 == 0 && lenmRNA) {
    char *lastCodon = &(mRNA[lenmRNA-3]);
    if (!strcmp(lastCodon,"TAG") ||
        !strcmp(lastCodon,"TGA") ||
        !strcmp(lastCodon,"TAA")) {
      *lastCodon = '\0';
    }
  }
  // the above lines will remove the final stop codon from the mRNA
  // sequence produced if it is present, this is so any peptide produced
  // won't have a terminal stop codon
  // if you want to have a terminal stop codon either comment this line out
  // or call translatable seq directly and produce a translation from it

  for (i=0;i<6;i++) {
    frm[i] = (char *)malloc(lenmRNA/3 + 2);
  }

  translate(mRNA, frm, lengths);

  free(mRNA);

  // Note don't free frame 0
  for (i=1;i<6;i++) {
    free(frm[i]);
  }

  return frm[0];
}


#ifdef DONE
char *Transcript_seq(Transcript *trans) {
  int i;
  char *transcriptSeqString;

  for (i=0; i<Transcript_getExonCount(trans); i++) {
    Exon *ex = Transcript_getExonAt(trans,i);
    transcriptSeqString = StrUtil_appendString(transcriptSeqString, Exon_getSeqString(ex));
  }

    my $seq = Bio::Seq->new(
        -DISPLAY_ID => $self->stable_id,
        -MOLTYPE    => 'dna',
        -SEQ        => $transcript_seq_string,
        );

  return seq;
}
#endif



MapperRangeSet *Transcript_cDNA2Genomic(Transcript *trans, int start, int end) {
  Mapper *mapper;

  mapper = Transcript_getcDNACoordMapper(trans);

  return Mapper_mapCoordinates(mapper, (IDType)trans, start, end, 1, CDNA_COORDS );
}

MapperRangeSet *Transcript_genomic2cDNA(Transcript *trans, int start, int end, int strand, BaseContig *contig) {
  Mapper *mapper;

  // "ids" in mapper are contigs of exons, so use the same contig that should
  // be attached to all of the exons...
  Exon *firstExon = Transcript_getExonAt(trans,0);
  if (!contig) {
    contig = Exon_getContig(firstExon);
  }

  mapper = Transcript_getcDNACoordMapper(trans);


  printf("MAPPING %d - %d (%d)\n",start,end,strand);
  printf(" %s = %s\n", BaseContig_getName(contig), 
         BaseContig_getName(Exon_getContig(firstExon)));

  return Mapper_mapCoordinates(mapper,(IDType)contig, start, end, strand, GENOMIC_COORDS);
}


Mapper *Transcript_getcDNACoordMapper(Transcript *trans) {
  Mapper *mapper;
  int start = 1;
  int i;

  if (trans->exonCoordMapper) {
    return trans->exonCoordMapper;
  }

  // 
  // the mapper is loaded with OBJECTS in place of the IDs !!!!
  //  the objects are the contigs in the exons
  // 
  mapper = Mapper_new( CDNA_COORDS, GENOMIC_COORDS );

  for (i=0; i<Transcript_getExonCount(trans); i++) {
    Exon *exon = Transcript_getExonAt(trans,i);

    Exon_loadGenomicMapper(exon, mapper, (IDType)trans, start);
    start += Exon_getLength(exon);
  }
  trans->exonCoordMapper = mapper;
  return mapper;
}


