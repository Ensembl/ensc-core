/*
 * Copyright [1999-2024] EMBL-European Bioinformatics Institute
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *      http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#define __SEQEDIT_MAIN__
#include "SeqEdit.h"
#undef __SEQEDIT_MAIN__
#include "StrUtil.h"
#include <string.h>
#include "EcoString.h"
#include "EnsC.h"
#include "Attribute.h"

/*
=head1 DESCRIPTION

This is a class used to represent post transcriptional
modifications to sequences.  SeqEdit objects are stored as ordinary
Bio::EnsEMBL::Attributes with a parseable value and can be used to
represent RNA editing, selenocysteines etc.

Also see B<Bio::EnsEMBL::Attribute>

*/

/* NOTE: Perl new split into two separate news for different ways it can be initialised
            1) With an attribute (SeqEdit_newFromAttribute())
            2) With values (SeqEdit_new())
=head2 new

  Arg [-ATTRIB] : Bio::EnsEMBL::Attribute
                  Constructs a new SeqEdit from an Attribute.
                  Can only be provided if no other constructor arguments
                  are provided.
  Arg [-START]       : The start position of the edit.
  Arg [-END]         : The end position of the edit.
  Arg [-ALT_SEQ]     : The alternate sequence
  Arg [-CODE]        : A code for this SeqEdit
  Arg [-NAME]        : A name for this SeqEdit
  Arg [-DESCRIPTION] : Arg passed to superclass constructor
  Example    : my $sea = Bio::EnsEMBL::SeqEdit->new(-ATTRIB => $attrib);
               my $sea = Bio::EnsEMBL::SeqEdit->new
                             (-START => 10,
                              -END   => 12,
                              -ALT_SEQ => 'ACG',
                              -CODE    => '_rna_edit',
                              -NAME    => 'RNA Edit',
                              -DESCRIPTION => 'RNA edit');
  Description: Constructs a SeqEdit representing a single edit to a
               sequence, such as an rna modification or a selenocysteine.
  Returntype : Bio::EnsEMBL::SeqEdit
  Exceptions : throws if attribute set and other args aswell
               throws if start and end not set correctly of attribure not set
  Caller     : general
  Status     : Stable

=cut
*/
SeqEdit *SeqEdit_new(long start, long end, char *altSeq, char *name, char *desc, char *code) {
  SeqEdit *seqEdit;
  
  // I'm going to be strict about args
  if (altSeq == NULL || name == NULL || desc == NULL || code == NULL) {
    fprintf(stderr,"Undefined args in SeqEdit_new\n");
    exit(1);
  }

  if ((seqEdit = (SeqEdit *)calloc(1,sizeof(SeqEdit))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for seqEdit\n");
    return NULL;
  }


  if (start > end+1) {
    fprintf(stderr,"SeqEdit start must be less than or equal to end + 1\n");
  }

  if (start < 1) {
    fprintf(stderr,"SeqEdit start must be greater than or equal to 1\n");
  }

  if (end < 0) {
    fprintf(stderr,"SeqEdit end must be greater than or equal to 0\n");
  }

  SeqEdit_setName(seqEdit, name);
  SeqEdit_setCode(seqEdit, code);
  SeqEdit_setDescription(seqEdit, desc);
  SeqEdit_setAltSeq(seqEdit, altSeq);
  SeqEdit_setStart(seqEdit, start);
  SeqEdit_setEnd(seqEdit, end);

  seqEdit->funcs = &seqEditFuncs;

  seqEdit->objectType = CLASS_SEQEDIT;

  return seqEdit;
}


SeqEdit *SeqEdit_newFromAttribute(Attribute *attrib) {
  if (attrib == NULL) {
    fprintf(stderr,"Attrib arg is NULL in SeqEdit_newFromAttribute\n");
    exit(1);
  }

  char **tokens = NULL;
  int nTok;

  StrUtil_tokenize(&tokens, &nTok, Attribute_getValue(attrib));

  long start;
  long end;
  char *altSeq;

  // Note this condition sets start and end as well as checking they are longs!!!!
  if ((nTok < 2 || nTok > 3) || !StrUtil_isLongInteger(&start, tokens[0]) || !StrUtil_isLongInteger(&end, tokens[1])) {
    fprintf(stderr,"Badly formatted SeqEdit value string %s\n", Attribute_getValue(attrib));
    exit(1);
  }

  if (nTok > 2) {
    altSeq = tokens[2];
  } else {
    altSeq = "";
  }

  SeqEdit *seqEdit = SeqEdit_new(start, end, altSeq, Attribute_getName(attrib), Attribute_getDescription(attrib), Attribute_getCode(attrib));

  int i;
  for (i=0;i<nTok;i++) {
    free(tokens[i]);
  }
  free(tokens);

  return seqEdit;
}


/*
=head2 start

  Arg [1]    : (optional) int $start - the new start position
  Example    : $start = $se_attrib->start();
  Description: Getter/Setter for the start position of the region replaced
               by the alt_seq.

               Coordinates are inclusive and one-based, which means that
               inserts are unusually represented by a start 1bp higher than
               the end.

               E.g. start = 1, end = 1 is a replacement of the first base but 
               start = 1, end = 0 is an insert BEFORE the first base.
  Returntype : int
  Exceptions : none
  Caller     : Transcript, Translation
  Status     : Stable

=cut
*/
long SeqEdit_setStart(SeqEdit *seqEd, long start) {
  if (start < 1) {
    fprintf(stderr,"SeqEdit start must be greater than or equal to 1\n");
    exit(1);
  }
  seqEd->start = start;

  return seqEd->start;
}

/*
=head2 end

  Arg [1]    : (optional) int $end - the new end position
  Example    : $end = $se_attrib->end();
  Description: Getter/Setter for the end position of the region replaced
               by the alt_seq.

               Coordinates are inclusive and one-based, which means that
               inserts are unusually represented by a start 1bp higher than
               the end.

               E.g. start = 1, end = 1 is a replacement of the first base but
               start = 1, end = 0 is an insert BEFORE the first base.
  Returntype : int
  Exceptions : throws if end  <= 0
  Caller     : Transcript, Translation
  Status     : Stable

=cut
*/
long SeqEdit_setEnd(SeqEdit *seqEd, long end) {
  if (end < 1) {
    fprintf(stderr,"SeqEdit end must be greater than or equal to 1\n");
    exit(1);
  }
  seqEd->end = end;

  return seqEd->end;
}

/*
=head2 alt_seq

  Arg [1]    : (optional) string $alt_seq
  Example    : my $alt_seq = $se_attrib->alt_seq();
  Description: Getter/Setter for the replacement sequence used by this edit.
               The sequence may either be a string of amino acids or
               nucleotides depending on the context in which this edit is
               used.

               In the case of a deletion the replacement sequence is an empty
               string.
  Returntype : string
  Exceptions : none
  Caller     : Transcript, Translation
  Status     : Stable

=cut
*/
// Get method is in .h
char *SeqEdit_setAltSeq(SeqEdit *seqEdit, char *altSeq) {
  StrUtil_copyString(&(seqEdit->altSeq), altSeq, 0);

  return seqEdit->altSeq;
}


/*
=head2 length_diff

  Arg [1]    : none
  Example    : my $diff = $sea->length_diff();
  Description: Returns the difference in length caused by applying this
               edit to a sequence.  This may be be negative (deletion),
               positive (insertion) or 0 (replacement).

               If either start or end are not defined 0 is returned.
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/
long SeqEdit_getLengthDiff(SeqEdit *seqEd) {

  //return 0 if (!defined($self->{'end'}) || !defined($self->{'start'}));

  return strlen(SeqEdit_getAltSeq(seqEd)) - (SeqEdit_getEnd(seqEd) - SeqEdit_getStart(seqEd) + 1);
}



/*
=head2 name

  Arg [1]    : (optional) string $name
  Example    : my $name = $seqedit->name();
  Description: Getter/Setter for the name of this SeqEdit
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/
// Get method is in .h
ECOSTRING SeqEdit_setName(SeqEdit *seqEdit, char *name) {
  EcoString_copyStr(ecoSTable, &(seqEdit->name),name,0);

  if (seqEdit->name == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for name\n");
    return NULL;
  }

  return seqEdit->name;
}


/*
=head2 code

  Arg [1]    : (optional) string $code
  Example    : my $code = $seqedit->code();
  Description: Getter/Setter for the code of this SeqEdit
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/
// Get method is in .h
ECOSTRING SeqEdit_setCode(SeqEdit *seqEdit, char *code) {
  EcoString_copyStr(ecoSTable, &(seqEdit->code),code,0);

  if (seqEdit->code == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for code\n");
    return NULL;
  }

  return seqEdit->code;
}


/*
=head2 description

  Arg [1]    : (optional) string $desc
  Example    : my $desc = $seqedit->description();
  Description: Getter/Setter for the description of this SeqEdit
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/
// Get method is in .h
ECOSTRING SeqEdit_setDescription(SeqEdit *seqEdit, char *description) {
  EcoString_copyStr(ecoSTable, &(seqEdit->description),description,0);

  if (seqEdit->description == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for description\n");
    return NULL;
  }

  return seqEdit->description;
}

/*
=head2 get_Attribute

  Arg [1]    : none
  Example    : my $attrib = $seqedit->get_Attribute();
               $transcript->add_Attributes($attrib);
  Description: Converts a SeqEdit object into an Attribute object.  This
               allows the SeqEdit to be stored as any other attribute in the
               ensembl database.  The start/end and alt_seq properties
               should be set before calling this method.
  Returntype : Bio::EnsEMBL::Attribute
  Exceptions : warning if start/end or alt_seq properties are not defined
  Caller     : general
  Status     : Stable

=cut
*/
Attribute *SeqEdit_getAttribute(SeqEdit *seqEd) {
  long  start  = SeqEdit_getStart(seqEd);
  long  end    = SeqEdit_getEnd(seqEd);
  char *altSeq = SeqEdit_getAltSeq(seqEd);

  char value[1024];

//  if (defined($start) && defined($end) && defined($alt_seq)) {
  if (altSeq != NULL) {
    sprintf(value,"%ld %ld %s", start, end, altSeq);
  } else {
    // In perl this was a warning but it seems quite serious to me - bail
    fprintf(stderr,"Attribute value cannot be created unless start, end and alt_seq properties are defined\n");
    exit(1);
  }

  Attribute *attrib = Attribute_new();

  Attribute_setCode(attrib, SeqEdit_getCode(seqEd));
  Attribute_setValue(attrib, value);
  Attribute_setName(attrib, SeqEdit_getName(seqEd));
  Attribute_setDescription(attrib, SeqEdit_getDescription(seqEd));

  return attrib;
}


/* Note: C version creates copy of string - it will mostly likely need reallocating
=head2 apply_edit

  Arg [1]    : reference to string $seqref
  Example    : $sequence = 'ACTGAATATTTAAGGCA';
               $seqedit->apply_edit(\$sequence);
               print $sequence, "\n";
  Description: Applies this edit directly to a sequence which is
               passed by reference.  The coordinates of this SeqEdit
               are assumed to be relative to the start of the sequence
               argument.
               If either the start or end of this SeqEdit are not defined
               this function will not do anything to the passed sequence.
  Returntype : reference to the same sequence that was passed in
  Exceptions : none
  Caller     : Transcript, Translation
  Status     : Stable

=cut
*/

char *SeqEdit_applyEdit(SeqEdit *seqEd, char *seq) {
  if (seq == NULL) {
    fprintf(stderr,"Need a sequence to edit\n");
    exit(1);
  }

/* Don't allow undef start or end in C
  if(!defined($self->{'start'}) || !defined($self->{'end'})) {
    return $seqref;
  }
*/

  long len = SeqEdit_getEnd(seqEd) - SeqEdit_getStart(seqEd) + 1;
  long lenDiff = SeqEdit_getLengthDiff(seqEd);

  if (lenDiff) {
    char *newSeq;
  
    if ((newSeq = calloc(strlen(seq) + lenDiff + 1, sizeof(char))) == NULL) {
      fprintf(stderr,"Failed allocating new seq string in SeqEdit_applyEdit\n");
      exit(1);
    }
  
  // Think its always the same
  //    1) First part of seq (upto getStart -1)
  //    2) Next the replacement string
  //    3) Next the rest of seq starting from the end position of the edit
    strncat(newSeq, seq, SeqEdit_getStart(seqEd)-1);
    strcat(newSeq, SeqEdit_getAltSeq(seqEd));
    strcat(newSeq, &seq[SeqEdit_getEnd(seqEd)]);
  
    free(seq);
    seq = newSeq;
  } else { // Replacement is same length as region replaced - no need to make a copy
    fprintf(stderr,"DOING NON REPLACE CHANGE\n");
    memcpy(&seq[SeqEdit_getStart(seqEd)-1], SeqEdit_getAltSeq(seqEd), len);
  }
  
  return seq;
}

int SeqEdit_reverseStartCompFunc(const void *one, const void *two) {
  SeqEdit *se1 = *((SeqEdit**)one);
  SeqEdit *se2 = *((SeqEdit**)two);

  return SeqEdit_getStart(se2) - SeqEdit_getStart(se1);
}

void SeqEdit_free(SeqEdit *seqEdit) {
  fprintf(stderr,"Warning: SeqEdit_free NIY\n");
}
