#include "GenomicAlignAdaptor.h"

package Bio::EnsEMBL::Compara::DBSQL::GenomicAlignAdaptor;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Compara::GenomicAlign;
use Bio::EnsEMBL::Compara::DnaFrag;


@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

my $DEFAULT_MAX_ALIGNMENT = 20000;

GenomicAlignAdaptor *GenomicAlignAdaptor_new() {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  my $vals =
    $self->db->get_MetaContainer->list_value_by_key('max_alignment_length');

  if(@$vals) {
    $self->{'max_alignment_length'} = $vals->[0];
  } else {
    $self->warn("Meta table key 'max_alignment_length' not defined\n" .
	       "using default value [$DEFAULT_MAX_ALIGNMENT]");
    $self->{'max_alignment_length'} = $DEFAULT_MAX_ALIGNMENT;
  }

  return $self;
}


=head2 store

  Arg  1     : listref  Bio::EnsEMBL::Compara::GenomicAlign $ga 
               The things you want to store
  Example    : none
  Description: It stores the give GA in the database. Attached
               objects are not stored. Make sure you store them first.
  Returntype : none
  Exceptions : not stored linked dnafrag objects throw.
  Caller     : general

=cut

void GenomicAlignAdaptor_store(GenomicAlignAdaptor *gaa, Vector *genomicAligns) {
  char *qStr;
  StatementHandle *sth;
  char commaStr[2] = {'\0','\0'};
  int i;
  char tmpStr[65556];
  

  StrUtil_copyString(&qStr, "INSERT INTO genomic_align_block"
             " (consensus_dnafrag_id, consensus_start, consensus_end,"
             "  query_dnafrag_id, query_start, query_end, query_strand, method_link_id,"
             "  score, perc_id, cigar_line) VALUES ");
  
  for (i=0; i<Vector_getNumElement(genomicAligns); i++) {
    GenomicAlign *ga = Vector_getElementAt(genomicAligns,i);
    DNAFrag *consDNAFrag  = GenomicAlign_getConsensusDNAFrag(ga);
    DNAFrag *queryDNAFrag = GenomicAlign_getQueryDNAFrag(ga);

    // check that everything has dbIDs
    if (!DNAFrag_getDbID(consDNAFrag) || !DNAFrag_getDbID(queryDNAFrag)) {
      fprintf(stderr, "Error: dna_fragment in GenomicAlign is not in DB\n");
      exit(1);
    }
  }

  // all clear for storing
  
  for (i=0; i<Vector_getNumElement(genomicAligns); i++) {
    GenomicAlign *ga = Vector_getElementAt(genomicAligns,i);
    DNAFrag *consDNAFrag  = GenomicAlign_getConsensusDNAFrag(ga);
    DNAFrag *queryDNAFrag = GenomicAlign_getQueryDNAFrag(ga);

    IDType methodLinkId = GenomicAlignAdaptor_methodLinkIdByAlignmentType(GenomicAlign_getAlignmentType(ga));

    if (!methodLinkId) {
      fprintf(stderr, "Error: There is no method_link with this type [%s] in the DB.\n",
              GenomicAlign_getAlignmentType(ga));
      exit(1);
    }
    
    sprintf(tmpStr," %s(" IDFMTSTR ", %d, %d, " IDFMTSTR ", %d, %d, %d, " IDFMTSTR ", %f, %f, '%s')"   
            commaStr, 
            DNAFrag_getDbID(consDNAFrag),
            GenomicAlign_getConsensusStart(ga),
            GenomicAlign_getConsensusEnd(ga),
            DNAFrag_getDbID(queryDNAFrag),  
            GenomicAlign_getQueryStart(ga),
            GenomicAlign_getQueryEnd(ga),
            GenomicAlign_getQueryStrand(ga),
            methodLinkId,
            GenomicAlign_getScore(ga),
            GenomicAlign_getPercIdScore(ga),
            GenomicAlign_getCigarLine(ga));

    qStr = StrUtil_appendString(qStr, tmpStr);
    commaStr[0] = ','; 
  }
  
  sth = gaa->prepare((BaseAdaptor *)gaa, qStr, strlen(qStr));
  sth->execute(sth);
  sth->finish(sth);

  free(qStr);
}
     
=head2 _fetch_all_by_DnaFrag_GenomeDB_direct

  Arg  1     : Bio::EnsEMBL::Compara::DnaFrag $dnafrag
               All genomic aligns that align to this frag
  Arg [2]    : Bio::EnsEMBL::Compara::GenomeDB $target_genome
               optionally restrict resutls to matches with this
               genome. Has to have a dbID().
  Arg [3]    : int $start
  Arg [4]    : int $end
  Arg [5]    : int $method_link_id
  Example    : none
  Description: Find all GenomicAligns that overlap this dnafrag.
               Return them in a way that this frags are on the
               consensus side of the Alignment.
  Returntype : listref Bio::EnsEMBL::Compara:GenomicAlign
  Exceptions : none
  Caller     : general

=cut


Vector *GenomicAlignAdaptor_fetchAllByDNAFragGenomeDBDirect( GenomicAlignAdaptor *gaa, 
     DNAFrag *dnaFrag, GenomeDB *targetGenome, int *startP, int *endP, methodLinkId) {
  IDType dnaFragId;
  GenomeDB *genomeDB;
  char *qStr;
  char tmpStr[512];
  Vector *results;
  StatementHandle *sth;

  if (!dnaFrag) {
    fprintf(stderr, "Error: Input dnafrag must not be NULL\n");
    exit(1);
  }

  // formatting the dnafrag
  dnaFragId = DNAFrag_getDbID(dnaFrag);

  genomeDB = DNAFrag_getGenomeDB(dnaFrag);

  StrUtil_copyString(&qStr,
     "SELECT gab.consensus_dnafrag_id,"
     "       gab.consensus_start," 
     "       gab.consensus_end,"
     "       gab.query_dnafrag_id," 
     "       gab.query_start," 
     "       gab.query_end,"
     "       gab.query_strand,"
     "       gab.method_link_id,"
     "       gab.score,"
     "       gab.perc_id," 
     "       gab.cigar_line"
     " FROM genomic_align_block gab ";

  if (!target_genome) {
    qStr = StrUtil_appendString(qStr,", dnafrag d");
  }
  sprintf(tmpStr," WHERE gab.method_link_id = " IDFMTSTR, methodLinkId);
  qStr = StrUtil_appendString(qStr,tmpStr);

  results = Vector_new();

  if (!targetGenome ||
      GenomeDB_hasQuery(genomeDB, targetGenome, methodLinkId)) {
    Vector *qres;

    sprintf(tmpStr," AND gab.consensus_dnafrag_id = " IDFMTSTR, dnaFragId);
    qStr = StrUtil_appendString(qStr, tmpStr);

    if (defined $start && defined $end) {
      int lowerBound = *startP - GenomicAlignAdaptor_getMaxAlignmentLength(gaa);
      sprintf(tmpStr,
               " AND gab.consensus_start <= %d"
               " AND gab.consensus_start >= %d"
               " AND gab.consensus_end >= %d", *endP, lowerBound, *startP ) ;
      qStr = StrUtil_appendString(qStr, tmpStr);
    }

    if (targetGenome) {
      sprintf(tmpStr,
              " AND gab.query_dnafrag_id = d.dnafrag_id"
              " AND d.genome_db_id = " IDFMTSTR, GenomeDB_getDbID(targetGenome));
      qStr = StrUtil_appendString(qStr, tmpStr);
    }

    sth = gaa->prepare((BaseAdaptor *)gaa, qStr, strlen(qStr));
    sth->execute(sth);

    qres = GenomicAlignAdaptor_objectsFromStatementHandle(gaa, sth, 0);
    Vector_append(results,qres);
    Vector_free(qres,NULL);

    sth->finish(sth);
  }

  if (!targetGenome ||
      GenomeDB_hasConsensus(genomeDB, targetGenome, methodLinkId)) {
    Vector *cres;

    sprintf(tmpStr," AND gab.query_dnafrag_id = " IDFMTSTR, dnaFragId);
    qStr = StrUtil_appendString(qStr, tmpStr);

    if (startP && endP) {
      int lowerBound = *startP - GenomicAlignAdaptor_getMaxAlignmentLength(gaa);
      sprintf(tmpStr,
               " AND gab.query_start <= %d"
               " AND gab.query_start >= %d"
               " AND gab.query_end >= %d", *endP, lowerBound, *startP ) ;
      qStr = StrUtil_appendString(qStr, tmpStr);
    }
    if (targetGenome) {
      sprintf(tmpStr,
               " AND gab.consensus_dnafrag_id = d.dnafrag_id"
               " AND d.genome_db_id = " IDFMTSTR, GenomeDB_getDbID(targetGenome));
      qStr = StrUtil_appendString(qStr, tmpStr);
    }
    sth = gaa->prepare((BaseAdaptor *)gaa, qStr, strlen(qStr));
    sth->execute(sth);

    cres = GenomicAlignAdaptor_objectsFromStatementHandle(gaa, sth, 1);
    Vector_append(results,cres);
    Vector_free(cres,NULL);

    sth->finish(sth);
  }
  free(qStr);

  return results;
}



=head2 fetch_all_by_DnaFrag_GenomeDB

  Arg  1     : Bio::EnsEMBL::Compara::DnaFrag $dnafrag
  Arg  2     : string $query_species
               The species where the caller wants alignments to
               his dnafrag.
  Arg [3]    : int $start
  Arg [4]    : int $end
  Arg [5]    : string $alignment_type
               The type of alignments to be retrieved
               i.e. WGA or WGA_HCR
  Example    :  ( optional )
  Description: testable description
  Returntype : none, txt, int, float, Bio::EnsEMBL::Example
  Exceptions : none
  Caller     : object::methodname or just methodname

=cut

Vector *GenomicAlignAdaptor_fetchAllbyDNAFragGenomeDB(GenomicAlignAdaptor *gaa,
               DNAFrag *dnaFrag, GenomeDB *targetGenome, int *startP, int *endP, 
               char *alignmentType) {

  GenomeDB *genomeCons;
  IDType methodLinkId;
  GenomeDB *genomeQuery;
  Vector mergedAligns;

  if (!dnaFrag)
    fprintf(stderr, "Error: dnaFrag argument must be non NULL\n");
    exit(1);
  }

  methodLinkId = GenomicAlignAdaptor_methodLinkIdByAlignmentType(gaa, alignmentType);

  genomeCons = DNAFrag_getGenomeDB(dnaFrag);
  genomeQuery = targetGenome;
  
  // direct or indirect ??
  if (GenomeDB_hasConsensus(genomeCons, genomeQuery, methodLinkId) ||
      GenomeDB_hasQuery(genomeCons, genomeQuery, methodLinkId)) {
    return GenomicAlignAdaptor_fetchAllByDnaFragGenomeDBDirect(gaa, 
                   dnaFrag, targetGenome, startP, endP, methodLinkId);
  } else {
    // indirect checks
    Vector *linkedCons  = GenomeDB_linkedGenomesByMethodLinkId(genomeCons, methodLinkId);
    Vector *linkedQuery = GenomeDB_linkedGenomesByMethodLinkId(genomeQuery, methodLinkId);
    
    // there are not many genomes, square effort is cheap
    Vector *linked = Vector_new();
    Vector *set1 = Vector_new();
    int i;

    for (i=0; i<Vector_getNumElement(linkedCons); i++) {
      int j;
      GenomeDB *g1 = Vector_getElementAt(linkedCons, i);

      for (j=0; j<Vector_getNumElement(linkedQuery); j++) {
        GenomeDB *g2 = Vector_getElementAt(linkedQuery, i);
	if (g1 == g2) {
	  Vector_addElement(linked, g1);
	}
      }
    }
    Vector_free(linkedCons,NULL);
    Vector_free(linkedQuery,NULL);

    // collect GenomicAligns from all linked genomes
    for (i=0; i<Vector_getNumElement(linked); i++) {
      GenomeDB *g = Vector_getElementAt(linked, i);

      Vector *gres = GenomicAlignAdaptor_fetchAllByDNAFragGenomeDBDirect(gaa, 
                             dnaFrag, g, startP, endP, methodLinkId);
      Vector_append(set1, gres);

      Vector_free(gres);
    }

    // go from each dnafrag in the result set to target_genome
    // there is room for improvement here: create start end
    // my %frags = map { $_->query_dnafrag->dbID => $_->query_dnafrag } @$set1;
    

    for (i=0; i<Vector_getNumElement(set1); i++) {
      GenomicAlign *alignA = Vector_getElementAt(set1,i);
      DNAFrag *frag = GenomicAlign_getQueryDNAFrag(alignA);
      int qStart = GenomicAlign_getQueryStart(alignA);
      int qEnd   = GenomicAlign_getQueryEnd(alignA);

      Vector dres = GenomicAlignAdaptor_fetchAllByDNAFragGenomeDBDirect(gaa,
                      frag, genomeQuery, &qStart, &qEnd, methodLinkId);
      int j;

      for (j=0; j<Vector_getNumElement(dres); j++) {
        GenomicAlign *alignB = Vector_getElementAt(dres,j);
	GenomicAlignAdaptor_addDerivedAlignments(gaa,  mergedAligns, alignA, alignB);
      } 
      Vector_free(dres,NULL);
    }
// NIY freeing
    return mergedAligns;
  }
}


=head2 _merge_alignsets

  Arg  1     : listref Bio::EnsEMBL::Compara::GenomicAlign $set1
               from consensus to query
  Arg  2     : listref Bio::EnsEMBL::Compara::GenomicAlign $set2
               and over consensus to next species query             
  Example    : none
  Description: set1 contains GAs with consensus species belonging to
               the input dnafragment. Query fragments are the actual reference
               species. In set 2 consensus species is the reference and
               query is the actual target genome. There may be more than
               one reference genome involved.
  Returntype : listref Bio::EnsEMBL::Compara::GenomicAlign
  Exceptions : none
  Caller     : internal

=cut


Vector *GenomicAlignAdaptor_mergeAlignsets(GenomicAlignAdaptor *gaa, Vector *alignSet1, Vector *alignSet2) {

  // sorting of both sets
  // walking through and finding overlapping GAs
  // create GA from overlapping GA
  // return list of those

  // efficiently generating all Aligns that overlap
  // [ key, object, set1 or 2 ]
  // Alignments are twice in big list. They are added to the overlapping
  // set the first time they appear and they are removed the
  // second time they appear. Scanline algorithm

  my @biglist = ();
  for my $align ( @$alignset1 ) {
    push( @biglist, 
          [ $align->query_dnafrag()->dbID(), 
            $align->query_start(), $align, 0 ] );
    push( @biglist, 
          [ $align->query_dnafrag()->dbID(), 
            $align->query_end()+.5, $align, 0 ] );
  }

  for my $align ( @$alignset2 ) {
    push( @biglist, 
          [ $align->consensus_dnafrag()->dbID(), 
            $align->consensus_start(), $align, 1 ] );
    push( @biglist, 
          [ $align->consensus_dnafrag()->dbID(), 
            $align->consensus_end()+.5, $align, 1 ] );
  }
  
  my @sortlist = sort { $a->[0] <=> $b->[0] ||
                        $a->[1] <=> $b->[1] } @biglist;

  // walking from start to end through sortlist and keep track of the 
  // currently overlapping set of Alignments
 
  my @overlapping_sets = ( {}, {} ); 
  my ($align, $setno);
  my $merged_aligns = [];

  for my $aligninfo ( @sortlist ) {
    $align = $aligninfo->[2];
    $setno = $aligninfo->[3];

    if( exists $overlapping_sets[ $setno ]->{ $align } ) {
      // remove from current overlapping set
      delete $overlapping_sets[ $setno ]->{ $align };
    } else {
      // insert into the set and do all the overlap business
      $overlapping_sets[ $setno ]->{ $align } = $align;
      // the other set contains everything this align overlaps with
      for my $align2 ( values %{$overlapping_sets[ 1 - $setno ]} ) {
        if( $setno == 0 ) {
          GenomicAlignAdaptor_addDerivedAlignments(gaa, mergedAligns, align, align2);
        } else {
          GenomicAlignAdaptor_addDerivedAlignments(gaa, mergedAligns, align2, align);
        }
      }
    }
  }

  return mergedAligns;
}



=head2 _add_derived_alignments

  Arg  1     : listref 
    Additional description lines
    list, listref, hashref
  Example    :  ( optional )
  Description: testable description
  Returntype : none, txt, int, float, Bio::EnsEMBL::Example
  Exceptions : none
  Caller     : object::methodname or just methodname

=cut


sub GenomicAlignAdaptor_addDerivedAlignments(GenomicAlignAdaptor *gaa, 
                     Vector *mergedAligns, Vector *alignA, Vector *alignB) {

  // variable name explanation
  // q - query c - consensus s - start e - end l - last
  // o, ov overlap j - jump_in_
  // r - result
  my ( $qs, $qe, $lqs, $lqe, $cs, $ce, $lcs, $lce,
       $ocs, $oce, $oqs, $oqe, $jc, $jq, $ovs, $ove,
       $rcs, $rce, $rqs, $rqe);

  // initialization phase
  

  my @cigA = ( $alignA->cigar_line =~ /(\d*[MDI])/g );
  my @cigB;
  my $line = $alignB->cigar_line();

  if( $alignA->query_strand == -1 ) {
    @cigB = reverse ( $line =~ /(\d*[MDI])/g ); 
  } else {
    @cigB = ( $line =~ /(\d*[MDI])/g ); 
  }

  // need a 'normalized' start for qs, qe, oxs so I dont 
  // have to check strandedness all the time  

  // consensus is strand 1 and is not compared to anything,
  // can keep its original coordinate system
 
  $lce = $alignA->consensus_start() - 1;
  $ce = $lce;
  $cs = $ce + 1;
  
  // alignBs query can be + or - just keep relative coords for now
  $lqe = 0; $lqs = 1;
  $qe = 0; $qs = 1;

  // ocs will be found relative to oce and has to be comparable
  // to oqs. But it could be that we have to move downwards if we
  // are not - strand. thats why coordinates are trnaformed here

  if( $alignA->query_strand == -1 ) {
    // query_end is first basepair of alignment
    if( $alignA->query_end() < $alignB->consensus_end() ) {
      // oqs/e = 0 ocs/e = difference
      $oce = 0; $ocs = 1;
      $oqe = $alignB->consensus_end() - $alignA->query_end();
      $oqs = $oqe + 1;
    } else {
      $oqe = 0; $oqs = 1;
      $oce = $alignA->query_end() - $alignB->consensus_end();
      $ocs = $oce + 1;
    }
  } else {
    // in theory no coordinate magic necessary :-)
    $oqs = $alignA->query_start();
    $oqe = $oqs - 1; 
    $ocs = $alignB->consensus_start();
    $oce = $ocs - 1;
  }

  // initializing result
  $rcs = $rce = $rqs = $rqe = 0;
  my @result_cig= ();

  my $current_match = 0;
  my $new_match;
  

  while( 1 ) {
    // print "ocs $ocs oce $oce oqs $oqs oqe $oqe\n";
    // print "cs $cs ce $ce qs $qs qe $qe\n";
    // print "rcs $rcs rce $rce rqs $rqs rqe $rqe\n";
    // print "\n";

    // exit if you request a new piece of alignment and the cig list is 
    // empty

    if( $oce < $ocs || $oce < $oqs ) {
      // next M area in cigB
      last unless @cigB;
      $self->_next_cig( \@cigB, \$ocs, \$oce, \$qs, \$qe ); 
      next;
    }
    if( $oqe < $oqs || $oqe < $ocs ) {
      // next M area in cigA
      last unless @cigA;
      $self->_next_cig( \@cigA, \$cs, \$ce, \$oqs, \$oqe );
      next;
    }

    // now matching region overlap in reference genome
    $ovs = $ocs < $oqs ? $oqs : $ocs;
    $ove = $oce < $oqe ? $oce : $oqe;
    
    if( $current_match ) {
      $jc = $cs + ( $ovs - $oqs ) - $lce - 1;
      $jq = $qs + ( $ovs - $ocs ) - $lqe - 1;
    } else {
      $jc = $jq = 0;
    }

    $new_match = $ove - $ovs + 1;
    my $new_ga = 0;

    if( $jc == 0 ) {
      if( $jq == 0 ) {
	$current_match += $new_match;
      } else {
        // store current match;
	push( @result_cig, $current_match."M" );
	$jq = "" if ($jq == 1); 
	// jq deletions;
	push( @result_cig, $jq."D" );
	$current_match = $new_match;
      }
    } else {
      if( $jq == 0 ) {
        // store current match;
	push( @result_cig, $current_match."M" );
	// jc insertions;
	$jc = "" if( $jc == 1 );
	push( @result_cig, $jc."I" );
	$current_match = $new_match;
         
      } else {

	push( @result_cig, $current_match."M" );
	// new GA
	my $query_strand = $alignA->query_strand() * $alignB->query_strand();
	my ( $query_start, $query_end );
	if( $query_strand == 1 ) {
	  $query_start = $rqs + $alignB->query_start() - 1;
	  $query_end = $rqe + $alignB->query_start() - 1;
	} else {
	  $query_end = $alignB->query_end() - $rqs + 1;
	  $query_start = $alignB->query_end() - $rqe + 1;
	}
      
	my $score = ( $alignA->score() < $alignB->score()) ? 
	  $alignA->score() : $alignB->score();
	my $perc_id =  int( $alignA->perc_id() * $alignB->perc_id() / 100 );

	my $ga = Bio::EnsEMBL::Compara::GenomicAlign->new
	  ( -consensus_dnafrag => $alignA->consensus_dnafrag,
	    -query_dnafrag => $alignB->query_dnafrag,
	    -cigar_line => join("",@result_cig),
	    -consensus_start => $rcs,
	    -consensus_end => $rce,
	    -query_strand => $query_strand, 
	    -query_start => $query_start,
	    -query_end => $query_end,
	    -adaptor => $self,
	    -perc_id => $perc_id,
	    -score => $score
	  );
	push( @$merged_aligns, $ga );
	$rcs = $rce = $rqs = $rqe = 0;
	@result_cig = ();
	
	$current_match = $new_match;
      }
    }


    
    $rcs = $cs+($ovs-$oqs) unless $rcs;
    $rce = $cs+($ove-$oqs);
    $rqs = $qs+($ovs-$ocs) unless $rqs;
    $rqe = $qs+($ove-$ocs);

    // update the last positions
    $lce = $rce; 
    $lqe = $rqe;

    // next piece on the one that end earlier
    my $cmp = ( $oce <=> $oqe );
 
    if( $cmp <= 0 ) {
      // next M area in cigB
      last unless @cigB;
      $self->_next_cig( \@cigB, \$ocs, \$oce, \$qs, \$qe ); 
    }
    if( $cmp >= 0 ) {
      // next M area in cigA
      last unless @cigA;
      $self->_next_cig( \@cigA, \$cs, \$ce, \$oqs, \$oqe );
    } 
  } // end of while loop

  // if there is a last floating current match
  if( $current_match ) {
    push( @result_cig, $current_match."M" );
    // new GA
    my $query_strand = $alignA->query_strand() * $alignB->query_strand();
    my ( $query_start, $query_end );
    if( $query_strand == 1 ) {
      $query_start = $rqs + $alignB->query_start() - 1;
      $query_end = $rqe + $alignB->query_start() - 1;
    } else {
      $query_end = $alignB->query_end() - $rqs + 1;
      $query_start = $alignB->query_end() - $rqe + 1;
    }
  
    my $score = ( $alignA->score() < $alignB->score()) ? 
      $alignA->score() : $alignB->score();
    my $perc_id =  int( $alignA->perc_id() * $alignB->perc_id() / 100  );
    
    my $ga = Bio::EnsEMBL::Compara::GenomicAlign->new
      ( -consensus_dnafrag => $alignA->consensus_dnafrag,
	-query_dnafrag => $alignB->query_dnafrag,
	-cigar_line => join("",@result_cig),
	-consensus_start => $rcs,
	-consensus_end => $rce,
	-query_strand => $query_strand, 
	-query_start => $query_start,
	-query_end => $query_end,
	-adaptor => $self,
	-perc_id => $perc_id,
	-score => $score
      );
    push( @$merged_aligns, $ga );
  // nothing to return all in merged_aligns
  }
}


sub GenomicAlignAdaptor_nextCig {
  my ( $self, $ciglist, $cs, $ce, $qs, $qe ) = @_;
  
  my ( $cig_elem, $type, $count );
  do {
    $cig_elem = shift( @$ciglist );
    ( $count ) = ($cig_elem =~ /(\d*)/);
    $count || ( $count = 1 );

    ( $type ) = ( $cig_elem =~ /(.)$/ );
    if( $type eq 'D' ) {
      $$qe += $count;
    } elsif( $type eq 'I' ) {
      $$ce += $count;
    } else {
      $$cs = $$ce + 1;
      $$ce = $$cs + $count - 1;
      $$qs = $$qe + 1;
      $$qe = $$qs + $count - 1;
    } 
  } until ( $type eq 'M' || ! ( @$ciglist ));
}

Vector *GenomicAlignAdaptor_objectsFromStatementHandle(GenomicAlignAdaptor *gaa, StatementHandle *sth,
                                                       int reverse) {
  Vector *results = Vector_new();
  ResultRow *row;


  my ( $consensus_dnafrag_id, $consensus_start, $consensus_end, $query_dnafrag_id,
       $query_start, $query_end, $query_strand, $method_link_id, $score, $perc_id, $cigar_string );
  if( $reverse ) {
    $sth->bind_columns
      ( \$query_dnafrag_id, \$query_start, \$query_end,  
	\$consensus_dnafrag_id, \$consensus_start, \$consensus_end, \$query_strand, \$method_link_id,
	\$score, \$perc_id, \$cigar_string );
  } else {
    $sth->bind_columns
      ( \$consensus_dnafrag_id, \$consensus_start, \$consensus_end, 
	\$query_dnafrag_id, \$query_start, \$query_end, \$query_strand, \$method_link_id,
	\$score, \$perc_id, \$cigar_string );
  }

  my $da = $self->db()->get_DnaFragAdaptor();

  while ((row = sth->fetchRow(sth))) {
    GenomicAlign *genomicAlign;
    char *alignment_type = GenomicAlignAdaptor_alignmentTypeByMethodLinkId(gaa, methodLinkId);
    
    if (reverse) {
      $cigar_string =~ tr/DI/ID/;
      if( $query_strand == -1 ) {
	// alignment of the opposite strand

	my @pieces = ( $cigar_string =~ /(\d*[MDI])/g );
	$cigar_string= join( "", reverse( @pieces ));
      }
    }
    
    
    genomicAlign = GenomicAlign_new();
    GenomicAlign_setAdaptor(genomicAlign, gaa);
    GenomicAlign_setConsensusDNAFrag(genomicAlign, );
      (
       -adaptor => $self,
       -consensus_dnafrag => $da->fetch_by_dbID( $consensus_dnafrag_id ),
       -consensus_start => $consensus_start,
       -consensus_end => $consensus_end,
       -query_dnafrag => $da->fetch_by_dbID( $query_dnafrag_id ),
       -query_start => $query_start,
       -query_end => $query_end,
       -query_strand => $query_strand,
       -alignment_type => $alignment_type,
       -score => $score,
       -perc_id => $perc_id,
       -cigar_line => $cigar_string
      );


    Vector_addElement(results, genomicAlign);
  }

  return results;
}


sub GenomicAlignAdaptor_alignmentTypeByMethodLinkId(GenomicAlignAdaptor *gaa, IDType methodLinkId) {
  StatementHandle *sth;
  ResultRow *row;
  char qStr[512];
  char *alignmentType;

  if (!methodLinkId) {
    fprintf(stderr, "Error: methodLinkId has to be defined");
    exit(1);
  } 

  sprintf(qStr,"SELECT type FROM method_link WHERE method_link_id = " IDFMTSTR, methodLinkId);
  sth = gaa->prepare((BaseAdaptor *)gaa, qStr, strlen(qStr));
  sth->execute(sth);

  if ((row = sth->fetchRow(sth))) {
    alignmentType = StrUtil_copyString(&alignmentType, row->getStringAt(row,0), 0); 
  } else {
    fprintf(stderr,"Error: No alignmentType for " IDFMTSTR "\n",methodLinkId);
    exit(1);
  }

  sth->finish(sth);

  return alignmentType;
}

IDType GenomicAlignAdaptor_methodLinkIdByAlignmentType(GenomicAlignAdaptor *gaa, char *alignmentType) {
  StatementHandle *sth;
  ResultRow *row;
  IDType methodLinkId = 0;
  char qStr[512];

  if (!alignmentType) {
    fprintf(stderr, "Error: alignment_type has to be defined\n");
    exit(1);
  }
  
  sprintf(qStr, "SELECT method_link_id FROM method_link WHERE type = '%s'", alignmentType);
  sth = gaa->prepare((BaseAdaptor *)gaa, qStr, strlen(qStr));
  sth->execute(sth);

  if ((row = sth->fetchRow(sth))) {
    methodLinkId = row->getLongLongAt(row,0); 
  } else {
    fprintf(stderr,"Error: No methodLinkId for %s\n",alignmentType);
    exit(1);
  }

  sth->finish(sth);

  return methodLinkId;
}
