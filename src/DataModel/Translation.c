#define __TRANSLATION_MAIN__
#include "Translation.h"
#undef __TRANSLATION_MAIN__
#include "TranslationAdaptor.h"
#include "AttributeAdaptor.h"
#include "SeqEdit.h"

#include <strings.h>

Translation *Translation_new() {
  Translation *t;

  if ((t = (Translation *)calloc(1,sizeof(Translation))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for t\n");
    return NULL;
  }

  Translation_setVersion(t,-1);

  t->objectType = CLASS_TRANSLATION;

  t->funcs = &translationFuncs;

  Object_incRefCount(t);

  return t;
}


/*
=head2 length

  Example    : print "Peptide length =", $translation->length();
  Description: Retrieves the length of the peptide sequence (i.e. number of
               amino acids) represented by this Translation object.
  Returntype : int
  Exceptions : none
  Caller     : webcode (protview etc.)
  Status     : Stable

=cut
*/
/* NIY - requires 'seq' method which is fiddly - it requires 'transcript' method
long Translation_getLength(Translation *translation) {
  
  my $seq = $self->seq();

  return ($seq) ? CORE::length($seq) : 0;
}
*/



// NIY:
// Quick hack implementation of this because I got bored getting the gene storing working - revisit
Vector *Translation_getAllDBEntries(Translation *t) {
  if (!t->dbLinks) {
    TranslationAdaptor *tlna = (TranslationAdaptor *)Translation_getAdaptor(t);

    if (tlna) {
      DBEntryAdaptor *dbea = DBAdaptor_getDBEntryAdaptor(tlna->dba);
      DBEntryAdaptor_fetchAllByTranslation(dbea,t);
    } else {
      t->dbLinks = emptyVector;
    }
  }

  return t->dbLinks;
}

char *Translation_getStableId(Translation *translation) {
  TranslationAdaptor *ta = (TranslationAdaptor *)Translation_getAdaptor(translation);

  if (StableIdInfo_getStableId(&(translation->si)) == NULL && ta) {
//    TranslationAdaptor_getStableEntryInfo(ta,translation);
  }
  return StableIdInfo_getStableId(&(translation->si));
}

int Translation_getVersion(Translation *translation) {
  TranslationAdaptor *ta = (TranslationAdaptor *)Translation_getAdaptor(translation);

  if (StableIdInfo_getVersion(&(translation->si)) == -1 && ta) {
//    TranslationAdaptor_getStableEntryInfo(ta,translation);
  }
  return StableIdInfo_getVersion(&(translation->si));
}


/*
=head2 get_all_Attributes

  Arg [1]    : optional string $attrib_code
               The code of the attribute type to retrieve values for.
  Example    : ($sc_attr) = @{$tl->get_all_Attributes('_selenocysteine')};
               @tl_attributes = @{$translation->get_all_Attributes()};
  Description: Gets a list of Attributes of this translation.
               Optionally just get Attrubutes for given code.
               Recognized attribute "_selenocysteine"
  Returntype : listref Bio::EnsEMBL::Attribute
  Exceptions : warning if translation does not have attached adaptor and
               attempts lazy load.
  Caller     : general, modify_translation
  Status     : Stable

=cut
*/
// New

// NIY:
// Because this can filter the results the vector that gets returned must be freeable - so for now
// make a copy of the translation->attributes vector if returning unfiltered so behaviour is 
// consistent. Long term probably want reference count incremented
Vector *Translation_getAllAttributes(Translation *translation, char *attribCode) {
  if (translation->attributes == NULL) {
    TranslationAdaptor *tlna = (TranslationAdaptor *)Translation_getAdaptor(translation);
    if (tlna == NULL) { // No adaptor
// Perl comments out the warning, I'll put it back for now, just in case
      //fprintf(stderr,"Warning: Cannot get attributes without an adaptor.\n");
      return Vector_new();
    }

    AttributeAdaptor *ata = DBAdaptor_getAttributeAdaptor(tlna->dba);
    translation->attributes = AttributeAdaptor_fetchAllByTranslation(ata, translation, NULL);
  }

  if (attribCode != NULL) {
    Vector *results = Vector_new();
    int i;
    for (i=0; i<Vector_getNumElement(translation->attributes); i++) {
      Attribute *attrib = Vector_getElementAt(translation->attributes, i);
      if (!strcasecmp(attrib->code, attribCode)) {
        Vector_addElement(results, attrib);
      }
    }
    return results;
  } else {
// See NIY note above for why I'm making a copy 
    return Vector_copy(translation->attributes);
  }
}

/*
=head2 get_all_SeqEdits

  Example    : my @seqeds = @{$transcript->get_all_SeqEdits()};
  Description: Retrieves all post transcriptional sequence modifications for
               this transcript.
  Returntype : Bio::EnsEMBL::SeqEdit
  Exceptions : none
  Caller     : spliced_seq()
  Status     : Stable

=cut
*/
// New
Vector *Translation_getAllSeqEdits(Translation *translation) {
  char *edits[] = { "initial_met", "_selenocysteine", "amino_acid_sub", NULL };

  Vector *seqEds = Vector_new();

  char **editP = edits;
  while (*editP) {
    char *edit = *editP;

    Vector *attribs = Translation_getAllAttributes(translation, edit);

    // convert attributes to SeqEdit objects
    int i;
    for (i=0; i<Vector_getNumElement(attribs); i++) {
      Attribute *attrib = Vector_getElementAt(attribs, i);
      SeqEdit *seqEd = SeqEdit_newFromAttribute(attrib);

      Vector_addElement(seqEds, seqEd);
    }

    Vector_free(attribs);
    editP++;
  }

  return seqEds;
}


/*
=head2 modify_translation

  Arg [1]    : Bio::Seq $peptide 
  Example    : my $seq = Bio::Seq->new(-SEQ => $dna)->translate();
               $translation->modify_translation($seq);
  Description: Applies sequence edits such as selenocysteines to the Bio::Seq 
               peptide thats passed in
  Returntype : Bio::Seq
  Exceptions : none
  Caller     : Bio::EnsEMBL::Transcript->translate
  Status     : Stable

=cut
*/

char *Translation_modifyTranslation(Translation *translation, char *seq) {
  Vector *seqEds = Translation_getAllSeqEdits(translation);

  // Sort in reverse order to avoid complication of adjusting
  // downstream edits.
  // HACK:   The translation ENSP00000420939 somehow makes the next line
  //         bomb out ($a or $b becomes undef) if the start() method
  //         is used.  I haven't been able to find out why.  It has 10
  //         Selenocysteine seqedits that looks correct.
  //         /Andreas (release 59)
  if (Vector_getNumElement(seqEds)) {
    Vector_sort(seqEds, SeqEdit_reverseStartCompFunc);
  
  //  @seqeds = sort { $b->{'start'} <=> $a->{'start'} } @seqeds;
  
  
    // Apply all edits.
    // Not particularly efficient currently, could improve by precalculating maximum size of new seq prior to applying edits
    int i;
    for (i=0; i<Vector_getNumElement(seqEds); i++) {
      SeqEdit *se = Vector_getElementAt(seqEds, i);
      seq = SeqEdit_applyEdit(se, seq);
    }
  
    //$seq->seq($peptide);
  }

  Vector_free(seqEds);

  return seq;
}

void Translation_transform(Translation *translation, IDHash *exonTransforms) {

  Exon * startExon = Translation_getStartExon(translation);
  Exon * endExon   = Translation_getEndExon(translation);
  IDType startExonRef = (IDType)startExon;
  IDType endExonRef = (IDType)endExon;

/* CHECK */
  if (IDHash_contains(exonTransforms,startExonRef)) {
    Translation_setStartExon(translation,IDHash_getValue(exonTransforms,startExonRef));
  } else {
    // do nothing, the start exon wasnt mapped
  }

  if (IDHash_contains(exonTransforms,endExonRef)) {
    Translation_setEndExon(translation,IDHash_getValue(exonTransforms,endExonRef));
  } else {
    // do nothing, the end exon wasnt mapped
  }
}

void Translation_free(Translation *translation) {
  Object_decRefCount(translation);

  if (Object_getRefCount(translation) > 0) {
    return;
  } else if (Object_getRefCount(translation) < 0) {
    fprintf(stderr,"Error: Negative reference count for Translation\n"
                   "       Freeing it anyway\n");
  }

  
  free(translation);
}

/*
=head2 genomic_start

    Args        : None
    Example     : $translation_genomic_start =
                      $translation->genomic_start();
    Description : Returns the start position of the translation in
                  genomic coordinates on the forward strand.
    Return type : Integer
    Exceptions  : None
    Caller      : General
    Status      : At Risk (Under Development)

=cut
*/
// New (don't bother caching because its pretty quick to calculate
long Translation_getGenomicStart(Translation *translation) {
  long genomicStart;

  if (Exon_getStrand(Translation_getStartExon(translation)) >= 0) {
    genomicStart = Exon_getStart(Translation_getStartExon(translation)) + (Translation_getStart(translation) - 1);
  } else {
    genomicStart = Exon_getEnd(Translation_getEndExon(translation)) - (Translation_getEnd(translation) - 1);
  }

  return genomicStart;
}

/*
=head2 genomic_end

    Args        : None
    Example     : $translation_genomic_end = $translation->genomic_end();
    Description : Returns the end position of the translation in genomic
                  coordinates on the forward strand.
    Return type : Integer
    Exceptions  : None
    Caller      : General
    Status      : At Risk (Under Development)

=cut
*/
// New
long Translation_getGenomicEnd(Translation *translation) {
  long genomicEnd;

  if (Exon_getStrand(Translation_getEndExon(translation)) >= 0) {
    genomicEnd = Exon_getStart(Translation_getEndExon(translation)) + (Translation_getEnd(translation) - 1);
  } else {
    genomicEnd = Exon_getEnd(Translation_getStartExon(translation)) - (Translation_getStart(translation) - 1);
  }

  return genomicEnd;
}

