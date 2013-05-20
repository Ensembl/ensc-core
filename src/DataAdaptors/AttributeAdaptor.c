/*
=head1 NAME

Bio::EnsEMBL::DBSQL::AttributeAdaptor - Provides database interaction for
Bio::EnsEMBL::Attribute objects.
*/


/*
=head1 SYNOPSIS

  # $db is a Bio::EnsEMBL::DBSQL::DBAdaptor object:
  $attribute_adaptor = $db->get_AttributeAdaptor();

  $attributes = $attribute_adaptor->fetch_all_by_MiscFeature($feature);

  $attributes = $attribute_adaptor->fetch_all_by_Slice($slice);

  $attribute_adaptor->store_on_Slice( $slice, \@attributes );

  $attribute_adaptor->store_on_MiscFeature( $misc_feature,
    \@attributes )
*/

/* 
NOTE: In perl this adaptor is implemented using the horrid AUTOLOAD lazyness. Obviously in C 
      I can't do that so I'll actually write the functions!
      Looks like there are 5 object types attributes exist on:
        Gene
        Transcript
        Translation
        Slice (seq_region)
        MiscFeature

      The fetch_all_by_ method seemed to also allow not passing in an object, in which case I 
      think it would return all the attribs of that type ('Gene', 'Slice' etc). I don't think
      that's ever used, so for now I'm not implementing that functionality. If I did I would
      write a 'fetch_all_by_type' method
*/
#include "AttributeAdaptor.h"
#include "Transcript.h"
#include "Gene.h"
#include "Translation.h"
#include "Slice.h"
#include "BaseAdaptor.h"
#include "MysqlUtil.h"

#include "StatementHandle.h"
#include "ResultRow.h"

/*
=head2 new

  Arg [...]  : Superclass args.  See Bio::EnsEMBL::DBSQL::BaseAdaptor
  Description: Instantiates a Bio::EnsEMBL::DBSQL::AttributeAdaptor
  Returntype : Bio::EnsEMBL::AttributeAdaptor
  Exceptions : none
  Caller     : DBAdaptor
  Status     : Stable

=cut
*/
AttributeAdaptor *AttributeAdaptor_new(DBAdaptor *dba) {
  AttributeAdaptor *ata; // Note I don't use the normal 'aa' abbreviation here because I find it difficult to see clearly

  if ((ata = (AttributeAdaptor *)calloc(1,sizeof(AttributeAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for AttributeAdaptor\n");
    exit(1);
  }
  BaseAdaptor_init((BaseAdaptor *)ata, dba, ATTRIBUTE_ADAPTOR);

  return ata;
}


/* Foul and disgusting AUTOLOAD nastiness
use vars '$AUTOLOAD';

sub AUTOLOAD {
  my ($self,@args) = @_;
  my @array_return=();
  my $ref_return = undef;
  $AUTOLOAD =~ /^.*::(\w+_)+(\w+)$/ ;

  my $sub = $1;
  my $type = $2;



#  print STDERR "AUTO".$AUTOLOAD."\n";

#  print STDERR "AUTOLOAD reached with call to $sub of type $type\n";
  if($self->can($sub)){
    return $self->$sub($type,@args);
  }
  else{
    warn("In AttribAdaptor cannot call sub $sub$type\n");
  }
  return undef;
}
*/


/* Truely horrible code
   Will implement stores later (NIY)
sub store_on_ {
  my $self = shift;
  my $type = shift;
  my $object = shift;
  my $attributes = shift;
  my $table;


  my $object_id;
  if($type =~ /[GT][er][na][en]/){
    if (!ref($object)){
      $object_id = $object;
    }
    else {
      $object_id = $object->dbID;
    }
    $table = lc($type);
#    $type = lc($type);
  }
  else{
    if(!ref($object) || !$object->isa('Bio::EnsEMBL::'.$type)) {
      throw("$type argument is required. but you passed $object");
    }
    if($type eq "Slice"){
      $object_id = $object->get_seq_region_id();
      $table = "seq_region"; 
      $type = "seq_region";
    }
    else{
      if($type eq "MiscFeature"){
	$type = "misc_feature";
	$table = "misc"; 
      }
      else{
	$table = lc($type);
      }
      
      $object_id = $object->dbID();
      my $db = $self->db();
      
      if(!$object->is_stored($db)) {
	throw("$type is not stored in this database.");
      }
      
    }
  }
  my $sth = $self->prepare( "INSERT into ".$table."_attrib ".
			    "SET ".$type."_id = ?, attrib_type_id = ?, ".
			    "value = ? " );

  my $undef_circular_cache = 0;
  for my $attrib ( @$attributes ) {
    if(!ref($attrib) && $attrib->isa('Bio::EnsEMBL::Attribute')) {
      throw("Reference to list of Bio::EnsEMBL::Attribute objects " .
            "argument expected.");
    }

#    next if ! $attrib->value;

    my $atid = $self->_store_type( $attrib );
    if ((defined $attrib->code) and ($attrib->code eq 'circular_seq')) {
	$undef_circular_cache = 1;
    }
    $sth->bind_param(1,$object_id,SQL_INTEGER);
    $sth->bind_param(2,$atid,SQL_INTEGER);
    $sth->bind_param(3,$attrib->value,SQL_VARCHAR);
    $sth->execute();
  }

  if($table eq "seq_region") {        
    if ($undef_circular_cache) {
	#the slice is circular
	$object->{'circular'} = 1;
	my $slice_adaptor = $object->adaptor();
        #undefine slice adaptor->is_circular and the circular slice cache
	if (defined $slice_adaptor) {
	    $slice_adaptor->{'is_circular'} = undef;
	    $slice_adaptor->{'circular_sr_id_cache'} = {};
	}
    }
  }

  return;
}

sub remove_from_{
  my $self   = shift;
  my $type   = shift;
  my $object = shift;
  my $code   = shift;
  my $table;

  if(!ref($object) || !$object->isa('Bio::EnsEMBL::'.$type)) {
    throw("$type argument is required or a attrib code. but you passed $object");
  }

  my $object_id;
  if($type eq "Slice"){
    $object_id = $object->get_seq_region_id();
    $table = "seq_region"; 
    $type = "seq_region";
    if ((defined $code) and ($code eq 'circular_seq')) {
	#undefine slice->is_circular, slice adaptor->is_circular and the circular slice cache
	$object->{'circular'} = undef;
	my $slice_adaptor = $object->adaptor();
	if (defined $slice_adaptor) {
	    $slice_adaptor->{'is_circular'} = undef;
	    $slice_adaptor->{'circular_sr_id_cache'} = {};
	}
    }
  }
  else{
    if($type eq "MiscFeature"){
      $type = "misc_feature";
      $table = "misc"; 
    }
    else{
      $table = lc($type);
    }

    $object_id = $object->dbID();
    my $db = $self->db();
    
    if(!$object->is_stored($db)) {
      throw("$type is not stored in this database.");
    }

  }

  if(!defined($object_id)) {
    throw("$type must have dbID.");
  }

  my $sth;
  if(defined($code)){
    $sth = $self->prepare("DELETE a FROM ".$table."_attrib a ,attrib_type at " .
                         "WHERE a.attrib_type_id = at.attrib_type_id AND ".
                         "a.".$type."_id = ? AND ".
			 "at.code like ?");
    $sth->bind_param(1,$object_id,SQL_INTEGER);
    $sth->bind_param(2,$code,SQL_VARCHAR);
  }
  else{
    $sth = $self->prepare("DELETE FROM ".$table."_attrib " .
                         "WHERE ".$type."_id = ?");
    $sth->bind_param(1,$object_id,SQL_INTEGER);
  }
  $sth->execute();

  $sth->finish();

  return;
}
*/

/* Shouldn't need
sub fetch_all {
  throw("Use of method fetch_all not supported for attributes");
}
*/

// Here I explicitly implement the various fetch functions
Vector *AttributeAdaptor_fetchAllByGene(AttributeAdaptor *ata, Gene *gene, char *code) {
  if (gene == NULL) {
    fprintf(stderr,"Error: NULL Gene in AttributeAdaptor_fetchAllByGene\n");
    exit(1);
  }

  char *type  = "gene";
  char *table = "gene";
  IDType id   = Gene_getDbID(gene);

  return AttributeAdaptor_doFetchAllByTypeAndTableAndID(ata, type, table, id, code);
}

Vector *AttributeAdaptor_fetchAllByTranscript(AttributeAdaptor *ata, Transcript *transcript, char *code) {
  if (transcript == NULL) {
    fprintf(stderr,"Error: NULL Transcript in AttributeAdaptor_fetchAllByTranscript\n");
    exit(1);
  }

  char *type  = "transcript";
  char *table = "transcript";
  IDType id   = Transcript_getDbID(transcript);

  return AttributeAdaptor_doFetchAllByTypeAndTableAndID(ata, type, table, id, code);
}

Vector *AttributeAdaptor_fetchAllByTranslation(AttributeAdaptor *ata, Translation *translation, char *code) {
  if (translation == NULL) {
    fprintf(stderr,"Error: NULL Translation in AttributeAdaptor_fetchAllByTranslation\n");
    exit(1);
  }

  char *type  = "translation";
  char *table = "translation";
  IDType id   = Translation_getDbID(translation);

  return AttributeAdaptor_doFetchAllByTypeAndTableAndID(ata, type, table, id, code);
}

Vector *AttributeAdaptor_fetchAllBySlice(AttributeAdaptor *ata, Slice *slice, char *code) {
  if (slice == NULL) {
    fprintf(stderr,"Error: NULL Slice in AttributeAdaptor_fetchAllBySlice\n");
    exit(1);
  }

  char *type  = "seq_region";
  char *table = "seq_region";
  IDType id   = Slice_getSeqRegionId(slice);

  return AttributeAdaptor_doFetchAllByTypeAndTableAndID(ata, type, table, id, code);
}

/* MiscFeature NIY
Vector *AttributeAdaptor_fetchAllByMiscFeature(AttributeAdaptor *ata, MiscFeature *miscFeature, char *code) {
  if (miscFeature == NULL) {
    fprintf(stderr,"Error: NULL MiscFeature in AttributeAdaptor_fetchAllByMiscFeature\n");
    exit(1);
  }

  char *type  = "misc_feature";
  char *table = "misc";
  IDType id   = MiscFeature_getDbID(miscFeature);

  return AttributeAdaptor_doFetchAllByTypeAndTableAndID(ata, type, table, id, code);
}
*/

Vector *AttributeAdaptor_doFetchAllByTypeAndTableAndID(AttributeAdaptor *ata, char *type, char *table, IDType objectId, char *code) {
  char qStr[1024];

  sprintf(qStr, "SELECT at.code, at.name, at.description, t.value "
                  "FROM %s_attrib t, attrib_type at "
                 "WHERE at.attrib_type_id = t.attrib_type_id", table);

  if (code != NULL){
    sprintf(qStr,"%s AND at.code like '%s'", qStr, code);
  }

//  if(defined($object_id)){
  sprintf(qStr,"%s AND t.%s_id = "IDFMTSTR, qStr, type, objectId);
//  }
		   
  StatementHandle *sth = ata->prepare((BaseAdaptor *)ata,qStr,strlen(qStr));
  sth->execute(sth);

  Vector *results = AttributeAdaptor_objectsFromStatementHandle(ata, sth);

  sth->finish(sth);

  return results;
}

/* Equivalent will be AttributeAdaptor_free
sub DESTROY{
}
*/


/* Don't bother - for deprecated use case
#
# _id_check
#
# backwards compatibility check:
# check if $ensID is an object; if so, return $obj->dbID
#

sub _id_check {
  my $self = shift;
  my $ensID = shift;

  if ($ensID =~ /^\d+$/) {
    return $ensID;
  
  } elsif (ref($ensID) eq 'Bio::EnsEMBL::Gene' or
      ref($ensID) eq 'Bio::EnsEMBL::Transcript' or
      ref($ensID) eq 'Bio::EnsEMBL::Translation') {

    warning("You should pass a dbID rather than an ensembl object to store the attribute on");

    if ($ensID->dbID) {
      return $ensID->dbID;
    } else {
      throw("Ensembl object ".$ensID->display_id." doesn't have a dbID, can't store attribute");
    }

  } else {
    throw("Invalid dbID");
  }

}
*/


/*
# _store_type
*/
/* NIY
sub _store_type {
  my $self = shift;
  my $attrib = shift;

  my $sth1 = $self->prepare
    ("INSERT IGNORE INTO attrib_type set code = ?, name = ?, ".
     "description = ?" );


  $sth1->bind_param(1,$attrib->code,SQL_VARCHAR);
  $sth1->bind_param(2,$attrib->name,SQL_VARCHAR);
  $sth1->bind_param(3,$attrib->description,SQL_LONGVARCHAR);

  my $rows_inserted =  $sth1->execute();

  my $atid = $sth1->{'mysql_insertid'};

  if($rows_inserted == 0) {
    # the insert failed because the code is already stored
    my $sth2 = $self->prepare
      ("SELECT attrib_type_id FROM attrib_type " .
       "WHERE code = ?");
    $sth2->bind_param(1,$attrib->code,SQL_VARCHAR);
    $sth2->execute();
    ($atid) = $sth2->fetchrow_array();

    $sth2->finish();

    if(!$atid) {
      throw("Could not store or fetch attrib_type code [".$attrib->code."]\n" .
	    "Wrong database user/permissions?");
    }
  }

  $sth1->finish();

  return $atid;
}
*/


/* 
  In perl note this is misspelled compared to normal _objs_from_sth
  Not sure if that's deliberate to avoid clash or just an error
  Obviously this isn't called through the normal path at all
  (just directly from within the fetch method).
*/
Vector *AttributeAdaptor_objectsFromStatementHandle(AttributeAdaptor *ata, StatementHandle *sth) {
  Vector *results = Vector_new();

  ResultRow *row;
  while (row = sth->fetchRow(sth)) {
    char *code  = row->getStringAt(row, 0);
    char *name  = row->getStringAt(row, 1);
    char *desc  = row->getStringAt(row, 2);
    char *value = row->getStringAt(row, 3);

    Attribute *attr = Attribute_new();
    Attribute_setCode(code);
    Attribute_setName(name);
    Attribute_setDescription(desc);
    Attribute_setValue(value);
    
    Vector_addElement(results, attr);
  }

  return results;
}
