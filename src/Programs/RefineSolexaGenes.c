
/*
=head1 DESCRIPTION

This module takes intron spanning dna_align_features and combines them with 
rough transcript models to build refined genes with CDS. The module produces 
transcripts representing all possible combinations of introns and exons which
are then filtered according to the specifications in the config.
The databases containing the various features to combine is defined in
Bio::EnsEMBL::Analysis::Config::Databases and the configuration for the 
module is defined in Bio::EnsEMBL::Analysis::Config::GeneBuild::RefineSolexaGenes
*/


static int  limit = 0;

RefineSolexaGenes *RefineSolexaGenes_new() {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  RefineSolexaGenes *rsg;
  

  if ((rsg = calloc(1,sizeof(RefineSolexaGenes))) == NULL) {
    fprintf(stderr,"Failed allocating RefineSolexaGenes\n");
    exit(1);
  }

  RunnableDB_readAndCheckConfig()
//  $self->read_and_check_config($REFINESOLEXAGENES_CONFIG_BY_LOGIC);

  // Hard limit to the number of possible paths to explore
  RefineSolexaGenes_setRecursiveLimit(rsg, 10000);
  // initialise intron feature cash
  // Doesn't seem to be used
  //my %feature_hash;
  //$self->feature_cash(\%feature_hash);
  return rsg;
}

/*
=head2 fetch_input

    Title        :   fetch_input
    Usage        :   $self->fetch_input
    Returns      :   nothing
    Args         :   none

=cut
*/

void RefineSolexaGenes_fetchInput(RefineSolexaGenes *rsg) {
  char *logic = Analysis_getLogicName(RefineSolexaGenes_getAnalysis(rsg));
  // fetch adaptors and store them for use later
  // explicitly attach the ref db

  RefineSolexaGenes_setDb(rsg, RefineSolexaGenes_getDbAdaptor(rsg, "REFERENCE_DB"));

  DBAdaptor *db = RefineSolexaGenes_getDb(rsg);

  $self->intron_slice_adaptor($self->get_dbadaptor($self->INTRON_DB)->get_SliceAdaptor)
    if $self->INTRON_DB;
  $self->repeat_feature_adaptor($self->db->get_RepeatFeatureAdaptor);
  $self->gene_slice_adaptor($self->get_dbadaptor($self->MODEL_DB)->get_SliceAdaptor);

  // want a slice and a full chromsome to keep everything on in the same coords
  Slice *slice = RefineSolexaGenes_fetchSequence(rsg, RefineSolexaGenes_getInputId());
  Slice *geneSlice =  SliceAdaptor_fetchByRegion(RefineSolexaGenes_getGeneSliceAdaptor(rsg),
                                                 "toplevel",
                                                 Slice_getSeqRegionName(slice),
                                                 Slice_getSeqStart(slice),
                                                 Slice_getSeqEnd(slice),
                                                 1);
  
  Slice *chrSlice = SliceAdaptor_fetchByRegion(DBAdaptor_getSliceAdaptor(db), "toplevel", Slice_getSeqRegionName(slice));

  RefineSolexaGenes_setChrSlice(rsg, chrSlice);

  // we want to fetch and store all the rough transcripts at the beginning - speeds things up
  // also we want to take out tiny introns - merge them into longer structures
  Vector *genes;
  char *modelLogicName = RefineSolexaGenes_getModelLogicName(rsg);
  if (modelLogicName != NULL) {
    genes = Slice_getAllGenes(geneSlice, NULL, modelLogicName, 1, NULL, NULL);
    fprintf(stderr,"Got %d genes with logic name %s\n", Vector_getNumElement(genes), modelLogicName);
  } else {
    genes = Slice_getAllGenes(geneSlice, NULL, NULL, 1, NULL, NULL);
    fprintf(stderr, "Got %d genes\n",Vector_getNumElement(genes)); 
  }

  Vector *prelimGenes = Vector_new();
  int i;
  for (i=0; Vector_getNumElement(genes); i++) {
    Gene *gene = Vector_getElementAt(genes, i);

    // put them on the chromosome
    gene = Gene_transfer(gene, chrSlice);
    // reject genes that are from a different slice that overlap our slice at the start or end
    // say the models has to be > 10% on the slice

    long os = Slice_getStart(slice);
    if (Gene_getStart(gene) > Slice_getStart(slice)) {
      os = Gene_getStart(gene);
    }
    long oe = Slice_getEnd(slice);
    if (Gene_getEnd(gene) < Slice_getEnd(slice)) {
      oe = Gene_getEnd(gene);
    }

    long overlap = oe - os + 1;
    long gc = (overlap / Gene_getLength(gene) * 1000) / 10;
    long sc =  (overlap / Slice_getLength(slice) * 1000) /10;
    printf("Gene has %ld %% overlap with the slice\nSlice has %ld %% overlap with the gene\n", gc, sc);
    if ( gc <= 10 && sc <= 10) {
      printf("Rejecting\n");
      continue;
    } 
    Vector_addElement(prelimGenes, gene);
  }
  fprintf(stderr, "Got %d genes after filtering boundary overlaps\n", Vector_getNumElement(prelimGenes)); 

  // determine strandedness ( including splitting merged genes )
  RefineSolexaGenes_setPrelimGenes(rsg, prelimGenes);

  Vector *intronBamFiles = RefineSolexaGenes_getIntronBamFiles(rsg);

  if (Vector_getNumElement(intronBamFiles) > 0) {
    int i;
    for (i=0; i<Vector_getNumElement(intronBamFiles); i++) {
      // Perl returned a hash but only seemed to use the FILE element
      char *intronFile = Vector_getElementAt(intronBamFiles, i);
    
// NIY IMPLEMENT BAM STUFF HERE
      my $sam = Bio::DB::Sam->new(   -bam => $intron_files->{FILE},
                                     -autoindex => 1,
                                     -expand_flags => 1,
                                 );
      if (sam == NULL) {
        fprintf(stderr, "Bam file " . $intron_files->{FILE} . "  not found \n");
      }
      long long count = 0;

      my $segment = $sam->segment($slice->seq_region_name,$slice->start,$slice->end);

      fprintf(stderr,"Bam file segment not found for slice " .  $slice->name . "\n")
        unless $segment;

      // need to seamlessly merge here with the dna2simplefeatures code
      RefineSolexaGenes_bamToIntronFeatures(rsg,  $self->bam_2_intron_features($segment,$intron_files);
    }
  } else {
    // pre fetch all the intron features
    RefineSolexaGenes_dnaToIntronFeatures(rsg, Slice_getStart(slice), Slice_getEnd(slice));
// NIY ??
//    $self->intron_slice_adaptor->dbc->disconnect_when_inactive(1);
  }
// NIY ??
//  $self->db->disconnect_when_inactive(1);
//  $self->gene_slice_adaptor->dbc->disconnect_when_inactive(1);
}

void  RefineSolexaGenes_run(RefineSolexaGenes *rsg) {
  RefineSolexaGenes_refineGenes(rsg);
}

/*
=head2 refine_genes

    Title        :   refine_genes
    Usage        :   $self->refine_genes
    Returns      :   nothing
    Args         :   none
    Description  :   Combines exons with introns in all possible combinations to 
                 :   Make a series of transcript models

=cut
*/

void RefineSolexaGenes_refineGenes(RefineSolexaGenes *rsg) {

  Vector *prelimGenes = RefineSolexaGenes_getPrelimGenes(rsg);

  int i;
// GENE:
  for (i=0; Vector_getNumElement(prelimGenes)  ;i++) {
    Gene *gene = Vector_getElementAt(prelimGenes, i);
    Transcript *transcript = Gene_getTranscriptAt(gene, 0);
   
    // hack taking out weeny models

    if (Transcript_getLength(gene) < 300) {
      continue;
    }

    Vector *models = Vector_new();

    int singleExon = 0;
    // first run on the fwd strand then on the reverse

//  STRAND: 
    for (strand = -1 ; strand <=1 ; strand+= 2) {
      if ( RefineSolexaGenes_getRecursiveLimit(rsg) > 10000 ) {
        // set recursion to 10000 in case it was raised for a tricky gene
        RefineSolexaGenes_setRecursiveLimit(rsg, 10000);
        fprintf(stderr, "Warning: lowering recursive limit after complex gene\n"); 
      }

      fprintf(stderr,"Running on strand %d\n", strand);

      StringHash *intronCount = StringHash_new(STRINGHASH_SMALL);
      Vector *exonIntron = Vector_new();
      
      StringHash *intronHash = StringHash_new(STRINGHASH_SMALL);
      Vector *exonPrevIntron = Vector_new();
      StringHash *intronExon = StringHash_new(STRINGHASH_SMALL);

      int mostRealIntrons = 0;
      double highestScore = 0;

      fprintf(stderr, "%s : %ld %ld:\n", Gene_getStableId(gene), Gene_getStart(gene), Gene_GetEnd(gene));

      Vector *exons = RefineSolexaGenes_mergeExons(rsg, gene, strand);
      Vector_sort(exons, SeqFeature_posSortFunc);

/*
#      foreach my $exon ( @exons ) {
#        print  "EXTRAEXON: " . 
#          $exon->seq_region_name ." " .
#            ($exon->start +20) ." " .
#              ($exon->end -20)." " .
#                ( $exon->end - $exon->start -40)  ."\n"
#                  if $exon->{"_extra"} ;
#      }
*/
   
      int exonCount = Vector_getNumElement(exons);
      Vector *fakeIntrons = Vector_new();
      StringHash *knownExons = StringHash_new(STRINGHASH_SMALL);

      long offset = 0;

//    EXON:   
      int j;
      for (j=0; j<=exonEount; j++) {
 
        my $exon = clone_Exon($exons[$j]);
        int retainedIntron;
        int leftIntrons = 0;
        int rightIntrons = 0;

// Hack hack hack hack
        $exon->{'left_mask'} = 0;
        $exon->{'right_mask'} = $exon->length;


        printf("%d : %ld %ld\n",j, Exon_getStart(exon), Exon_getEnd(exon));
        // make intron features by collapsing the dna_align_features

        Vector *introns;

        ($introns,$offset) = $self->fetch_intron_features($exon->seq_region_start,$exon->seq_region_end,$offset);

        my @left_c_introns;
        my @right_c_introns;
        my @left_nc_introns;
        my @right_nc_introns;
        my @filtered_introns;
        my $intron_overlap;
        my @retained_introns;

//      INTRON: 
        int k;
        for (j=0; i<Vector_getNumElement(introns); k++) {
          DNAAlignFeature *intron = Vector_getElementAt(introns, k);

          if (DNAAlignFeature_getStrand(intron) != strand) {
            continue;
          } 
          if (DNAAlignFeature_getLength(intron) > RefineSolexaGenes_getMinIntronSize(rsg)) {
            continue;
          } 
          if (DNAAlignFeature_getLength(intron) <= RefineSolexaGenes_getMaxIntronSize(rsg)) {
            continue;
          } 

          // discard introns that splice over our exon
          if (DNAAlignFeature_getStart(intron) < Exon_getStart(exon) && 
              DNAAlignFeature_getEnd(intron) > Exon_getEnd(exon) ) {
            intronOverlap++;
            continue;
          }

          // check to see if this exon contains a retained intron
          if (DNAAlignFeature_getStart(intron) > Exon_getStart(exon) && 
              DNAAlignFeature_getEnd(intron) < Exon_getEnd(exon) ) {
            retainedIntron = 1;
// Hack hack hack
            $exon->{'retained'} =1;

            Vector_addElement(retainedIntrons, intron);
          } else {
            // we need to know how many consensus introns we have to the 
            // left and  right in order to determine whether to put in 
            // a non consensus intron
            if (DNAAlignFeature_getEnd(intron) <= Exon_getEnd(exon)) {
              if (strstr(DNAAlignFeature_getHitSeqName(intron), "non canonical")) {
                if (Intron_getScore(intron) > 1) {
                  Vector_addElement(leftNonConsIntrons, intron);
                }
              } else {
                Vector_addElement(leftConsIntrons, intron);
              }
            }
            if (DNAAlignFeature_getStart(intron) >= Exon_getStart(exon)) {
              if (strstr(DNAAlignFeature_getHitSeqName(intron), "non canonical")) {
                if (Intron_getScore(intron) > 1) {
                  Vector_addElement(rightNonConsIntrons, intron);
                }
              } else {
                Vector_addElement(rightConsIntrons, intron);
              }
            }
          }
        }
        
        // Restrict internal exons splice sites to most common
        // that way our alt splices will all share the same boundaries
        // but have different combinations of exons
        if ( $self->STRICT_INTERNAL_SPLICE_SITES && 
             // either we apply it equeally to all exons
             ( $self->STRICT_INTERNAL_END_EXON_SPLICE_SITES or
               // only apply to internal exons, leave out end exons
               ( !$self->STRICT_INTERNAL_END_EXON_SPLICE_SITES && 
                 ( scalar(@left_c_introns)  + scalar(@left_nc_introns) ) > 0 &&
                 ( scalar(@right_c_introns) + scalar(@right_nc_introns)) > 0 ))){
          // pick best left splice
          my $best_left_splice;
          my $best_left_score = 0;
          my @all_left_introns =  @left_c_introns;
          push @all_left_introns, @left_nc_introns;

          foreach my $intron ( @all_left_introns ) {
            if ( $best_left_score < $intron->score ) {
              $best_left_score = $intron->score;
              $best_left_splice = $intron->end;
            }
          }
          
          // pick best right  splice
          my $best_right_splice;
          my $best_right_score = 0;
          my @all_right_introns =  @right_c_introns;
          push @all_right_introns, @right_nc_introns;
          foreach my $intron ( @all_right_introns ) {
            if ( $best_right_score < $intron->score ) {
              $best_right_score = $intron->score;
              $best_right_splice = $intron->start;
            }
          }
          // filter out introns that pick other splice sites
          foreach my $intron ( @all_left_introns ) {
            push @filtered_introns, $intron if $intron->end == $best_left_splice;
          }

          foreach my $intron ( @all_right_introns ) {
            push @filtered_introns, $intron if $intron->start == $best_right_splice;
          }
          
        } else {
          
          // add non consensus introns only where there are no consensus introns
          Vector_append(filteredIntrons, leftConsIntrons);
          Vector_append(filteredIntrons, rightConsIntrons);

          if (Vector_getNumElement(leftConsIntrons) == 0) {
            Vector_append(filteredIntrons, leftNonConsIntrons);
          }
          if (Vector_getNumElement(rightConsIntrons) == 0) {
            Vector_append(filteredIntrons, rightNonConsIntrons);
          }
        }

/* Does nothing
        if ( scalar(@left_c_introns)  == 0 && scalar(@left_nc_introns)  > 0) {
         # print STDERR "using " . scalar(@left_nc_introns) . " NC left \n";
        } 
        if ( scalar(@right_c_introns)  == 0 && scalar(@right_nc_introns)  > 0 ) {
         # print STDERR "using " . scalar(@right_nc_introns) . " NC right \n";
        }
*/
        
        // single exon models are a special case
        if ( Vector_getNumElement(exons) == 1 && 
             Vector_getNumElement(filteredIntrons) == 0 &&  
             Vector_getNumElement(retainedIntrons) == 0 ) {
          //# at least on this strand this model looks like a single exon
          singleExon += 1;
        }
        
        // we dont want to allow left and right introns to overlap - 
        // it leads to -ve length exons
        
        // we put all the retained introns in at the end we want to do all the 
        // entrances and exits to each exon before we look at whether its 
        // retained or not
        Vector_sort(retainedIntrons, SeqFeature_posSortFunc);

        // push @filtered_introns, @retained_introns;
//      INTRON:  
        int k;
        for (k=0; k<Vector_getNumElement(filteredIntrons); k++) {
          DNAAlignFeature *intron = Vector_getElementAt(filteredIntrons, k);
          //# print STDERR "\t" . $intron->start . " " . $intron->end . " " . $intron->strand . " " . $intron->hseqname . " " . $intron->score . "\n";
          // becasue we make a new exons where we have a reatained intron to 
          // stop circular references we need to allow the final 
          // intron splicing out of the exon to be used more than once
          // by each new exon in fact

          StringHash_add(intronCount, DNAAlignFeature_getHitSeqName(intron),
          $intron_count{$intron->hseqname}++ unless $retained_intron;

          StringHash_add(intronHash, DNAAlignFeature_getHitSeqName(intron), intron);

          // only use each intron twice once at the end and once at the start of
          // an exon
          // exon_intron links exons to the intron on their right ignoring strand
          if (DNAAlignFeature_getEnd(intron) > Exon_getEnd(exon)) {
            push @{ $exon_intron[$j]}  , $intron
          }

          // intron exon links introns to exons on their right ignoring strand
          if (DNAAlignFeature_getStart(intron) < Exon_getStart(exon)) {
// NOT SURE WHAT THIS LINE SHOULD BE!!!!!
            push @{$intron_exon{$intron->hseqname}} , $intron;
            // exon_prev_intron links exons to introns on their left ignoring strand
            push @{ $exon_prev_intron[$j]}  , $intron ;
          }
        }
        if (Vector_getNumElement(retainedIntrons) > 0) {
          //#print STDERR "Dealing with " . scalar( @retained_introns ) . " retained introns \n";
          my @new_exons;
          push @new_exons,  $exon;
          // sort first by start then by end where start is the same
          @retained_introns =  sort {$a->start <=> $b->start } @retained_introns;
          for ( my $i = 0; $i < $#retained_introns ; $i++ ) {
            if ( $retained_introns[$i]->start ==  $retained_introns[$i+1]->start &&
                 $retained_introns[$i]->end >  $retained_introns[$i+1]->end ) {
              // reverse the order
              my $temp =  $retained_introns[$i];
              $retained_introns[$i] = $retained_introns[$i+1];
              $retained_introns[$i+1] = $temp;
            }
          }
          // now lets deal with any retained introns we have
//  RETAINED: 
          foreach my $intron ( @retained_introns ) {
            // we dont need to make all new exons for each alternate splice
            // check the intron is still retained given the new exons
            my $retained = 1;
            foreach my $new_exon ( @new_exons ) {
              if  (  $intron->start > $new_exon->start && $intron->end < $new_exon->end ) {
              } else {
                $retained = 0;
              }
              next RETAINED unless $retained;
            }
            my $reject_score = 0;
            // intron is within the exon - this is not a true exon but a retained intron
            if (  $intron->start > $exon->start && $intron->end < $exon->end && $intron->length > $self->MIN_INTRON_SIZE ) {
              // we are going to make a new exon and chop it up
              // add intron penalty
              //#print STDERR "RETAINED INTRON PENALTY for " . $intron->display_id ." before " . $intron->score . " ";
              $reject_score = $intron->score - $self->RETAINED_INTRON_PENALTY;
              // intron penalty is doubled for nc introns 
              if ( $intron->hseqname  =~ /non canonical/ ) {
                $reject_score = $reject_score - $self->RETAINED_INTRON_PENALTY;
              }
              //#print STDERR " after " . $reject_score ."\n";
              if ( $reject_score < 1 ) {
                // treat as single exon
                if ( scalar(@exons) == 1 ) {
                  // at least on this strand this model looks like a single exon
                  $single_exon += 1;
                }
                next;
              }
              //#print STDERR  "Exon " . $exon->start ."\t". $exon->end . " has retained intron:\n     " . $intron->start ."\t" .  $intron->end ." "."\n";
              // dont have circular references to exons or the paths
              // will be infinite so clone this exon instead
              // I guess we also want to keep the original exon too?
              my $new_exon1 = clone_Exon( $exon );
              my $new_exon2 = clone_Exon( $exon );
              // chop it up a bit so it no longer overlaps the other introns
              fprintf(stderr, "TRIMMING EXON\n");
              my $length = $intron->end - $intron->start;
              $new_exon1->end( $intron->start + int( $length / 2 ) - 2 );
              $new_exon2->start( $intron->end - int( $length / 2 ) + 2 );
              
              push @new_exons,$new_exon1 unless $known_exons{$new_exon1->start."-".$new_exon1->end."-".$new_exon1->strand} ;
              push @new_exons,$new_exon2 unless $known_exons{$new_exon2->start."-".$new_exon2->end."-".$new_exon2->strand};
              
              $known_exons{$new_exon1->start."-".$new_exon1->end."-".$new_exon1->strand} = 1;
              $known_exons{$new_exon2->start."-".$new_exon2->end."-".$new_exon2->strand} = 1;
            }
          }
          if ( scalar (@new_exons > 1 ) ) {
            // we want to split the score equally across the new exons
            foreach my $e ( @new_exons ) {
              foreach my $d ( @{$e->get_all_supporting_features} ) {
                $d->score($d->score / scalar(@new_exons));
              }
            }
            
            splice( @exons,$i,1,@new_exons);
            for ( my $i = 0 ; $i<= $#exons ; $i++ ) {
              my $e = $exons[$i];
            }
            print "ADDED " . scalar( @new_exons) . " new exons\n";
            $exon_count+= $#new_exons;
            // make sure they are all stil sorted
            @exons = sort { $a->start <=> $b->start }  @exons;
          }
        }
      }
      
      next unless @exon_intron;
      // Loop around the path generation, 
      // if there are too many paths to process return undef
      // then re-run the path processing but with increasing strictness
      // where strictness = elimianating alternate low scoring introns
      my $paths;
      my $strict = 0;
      while ( !$paths ) {
        $paths = $self->process_paths( \@exons, \@exon_intron, \%intron_exon, $strict );
        next GENE if $paths && $paths eq 'Give up';
        $strict++;
      }  
      print STDERR "STRAND $strand BEFORE COLLAPSING  PATHS  = " . scalar( keys %$paths ) . "\n";
      // lets collapse redundant paths
      foreach my $path ( sort keys %$paths ) {
        //  print "PATHS $path\n";
        my @array = split ( /\./,$path);
        my ($start,$end,$middle);
        for ( my $j = 0 ; $j < scalar( @array ) ; $j++ )  {
          $start .= $array[$j] ."."  unless $j < 2;
          $middle .= $array[$j] ."." unless $j < 2  or $j >= $#array-1 ;
          $end .= $array[$j] . "." unless $j >= $#array-1;
        }
        // remove redunancy from the array
        delete $paths->{$start} if $start && $paths->{$start};
        delete $paths->{$end} if $end && $paths->{$end};
        delete $paths->{$middle} if $middle && $paths->{$middle};
      }
      
      print STDERR "AFTER COLLAPSING PATHS  = " . scalar( keys %$paths ) . "\n";
      push @models, @{$self->make_models($paths,$strand, \@exons,$gene,\%intron_hash )};
      print STDERR "Now have " . scalar ( @models ) ." models \n";
    }

    // recursively recluster the models to identify 'other' models 
    // with no overlap to the 'best' model
    my $model_count = 0;
    my  ($clustered_models,$new_clusters) = $self->recluster_models(\@models);
    my @clean_clusters;
    if ( $new_clusters ) {
      while (  scalar(@{$new_clusters}) ){
        push @clean_clusters,@{$clustered_models};
        ($clustered_models,$new_clusters) = $self->recluster_models($new_clusters);
//        print "Now have " .  scalar(@{$new_clusters}) ." new clusters after reclustering\n";
      }
    }
    push @clean_clusters,@{$clustered_models} if $clustered_models;
    
    // filter to identify 'best', 'other' and 'bad' models
    $self->filter_models(\@clean_clusters);

    // process single exon models
    // if it has no introns on either strand
    if ( $self->SINGLE_EXON_MODEL && $single_exon == 2 ) {
      my $exon = $self->merge_exons($gene,1)->[0];
      my $single_exon_model;
      //print STDERR " Single exon = $single_exon\n";
      next unless $exon->length+40 >= $self->MIN_SINGLE_EXON;
      //print STDERR "Passed length filter " . $exon->length ."\n";
      // trim padding 
      $exon->start($exon->start + 20);
      $exon->end  ($exon->end   - 20);
      // trim away strings of Ns from the start and  end
      // check start
      my $ex_seq = $exon->seq->seq;
      if ( $ex_seq =~ /^(N+)(\S+)$/ ) {
        $exon->start($exon->start + length($1));
      }
      // check end
      $ex_seq = reverse $ex_seq;
      if ( $ex_seq =~ /^(N+)(\S+)$/ ) {
        $exon->end($exon->end - length($1));
      }
     
      // get the cds
      my $fwd_exon =  clone_Exon($exon);
      $fwd_exon->strand(1);
      my $rev_exon = clone_Exon($exon);
      my $fwd_t =  new Bio::EnsEMBL::Transcript(-EXONS => [$fwd_exon]);
      my $fwd_tran = compute_translation(clone_Transcript($fwd_t));
      my $fwd_t_len;
      $fwd_t_len = $fwd_tran->translation->genomic_end - $fwd_tran->translation->genomic_start 
        if $fwd_tran->translateable_seq;
      //print STDERR "FWD t length $fwd_t_len\n";
      my $rev_t =  new Bio::EnsEMBL::Transcript(-EXONS => [$rev_exon]);
      my $rev_tran = compute_translation(clone_Transcript($rev_t));
      my $rev_t_len;
      $rev_t_len = $rev_tran->translation->genomic_end - $rev_tran->translation->genomic_start
        if $rev_tran->translateable_seq;;
      //print STDERR "REV t length $rev_t_len\n";
      if ( $fwd_tran->translateable_seq &&  
           ( $fwd_t_len / $fwd_tran->length )* 100 >= $self->SINGLE_EXON_CDS &&
           $fwd_t_len >  $rev_t_len ) {
        // keep this one
        $single_exon_model =  $fwd_tran;
      }
      if ( $rev_tran->translateable_seq &&  
           ( $rev_t_len / $rev_tran->length )* 100 >= $self->SINGLE_EXON_CDS &&
           $rev_t_len >  $fwd_t_len ) {
        // keep this one
        $single_exon_model = $rev_tran;
      }
      if ( $single_exon_model ) {
        $single_exon_model->analysis($self->analysis);
        $single_exon_model->version(1);
        my ( $new_gene ) = @{convert_to_genes(($single_exon_model),$gene->analysis)};
        $new_gene->biotype($self->SINGLE_EXON_MODEL); 
        // score comes from exon supporting feature;
        my $score =  $exon->get_all_supporting_features->[0]->score;
        $exon->flush_supporting_features;
        $new_gene->stable_id($gene->stable_id . "-v1-" . int($score) );
        push @{$self->output} , $new_gene;
      }
    }
  }
}

Analysis *RefineSolexaGenes_createAnalysisObject(RefineSolexaGenes *rsg, char *logicName) {
  DBAdaptor *outDb = RefineSolexaGenes_getDbAdaptor(rsg, RefineSolexaGenes_getOutputDb(rsg));

  AnalysisAdaptor *aa = DBAdaptor_getAnalysisAdaptor(outDb);
  Analysis *analysis = AnalysisAdaptor_fetchByLogicName(aa, logicName);

  if (analysis != NULL) {
     return analysis ;
  }
  // need to create analysis object first
  analysis = Analysis_new();
  Analysis_setLogicName(logicName);

  return analysis;
}

/*
=head2 recluster_models

    Title        :   recluster_models
    Usage        :   $self->recluster_models(\@models);
    Returns      :   nothing
    Args         :   Array ref of array references of Bio::EnsEMBL::Transcript
    Description  :   reclusters 'other' models that have no overlap with 'best' models

=cut
*/
// called with 'models' which was immediately changed to clusters, so switched method arg to clusters
Vector *RefineSolexaGenes_reclusterModels(RefineSolexaGenes *rsg, Vector *clusters) {
  Vector *newClusters   = Vector_new();
  Vector *finalClusters = Vector_new();

  int strand;
  for  (strand=1; strand >= -1; strand-=2) {
   // print "Running on strand $strand\n";
    
    int i;
    for (i=0; i<Vector_getNumElements(clusters); i++) {
      GeneCluster *cluster = Vector_getElementAt(clusters, i);

      if (cluster->finalModels == NULL) {
        continue;
      }

      Vector *strandedGenes = Vector_new();

      int j;
      for (j=0; j<Vector_getNumElement(cluster->finalModels); j++) {
        Gene *gene = Vector_getElementAt(cluster->finalModels, j);
        if (Gene_getStrand(gene) == strand) {
        //#print "GENE STRAND " . $gene->strand ." vs $strand \n";
          Vector_addElement(strandedGenes, gene);
        }
      }

      if (Vector_getNumElement(strandedGenes) == 0) {
        Vector_free(strandedGenes);
        continue;
      }

      GeneCluster *strandedCluster = RefineSolexaGenes_recalculateCluster(rsg, strandedGenes);

      Gene *best = NULL;
      Vector *genes;
      Vector *otherGenes;

      for (j=0; j<Vector_getNumElement(strandedCluster->finalModels); j++) {
        Gene *gene = Vector_getElementAt(strandedCluster->finalModels);

        if (Gene_getStrand(gene) == strand) {
          if (!strcmp(Gene_getBioType(gene), RefineSolexaGenes_getBestScoreType(rsg))) {
            //$best =  $gene if $gene->biotype eq $self->BEST_SCORE;
            best = gene;
            Vector_addElement(genes, gene);
          }
      }

      if (best == NULL) {
        fprintf(stderr,"No best model found\n");
        exit(1);
      }

      // now recluster 
    //OTHERGENE: 
      for (j=0; j<Vector_getNumElement(strandedCluster->finalModels); j++) {
        Gene *gene = Vector_getElementAt(strandedCluster->finalModels);

        if (!strcmp(Gene_getBioType(gene), RefineSolexaGenes_getBestScoreType(rsg))) {
          continue;
        }

        my @best_exons =  sort { $a->start <=> $b->start } @{$best->get_all_Transcripts->[0]->get_all_Exons};
        my @other_exons = sort { $a->start <=> $b->start } @{$gene->get_all_Transcripts->[0]->get_all_Exons};
        // exon overlap with a best model

      //BESTEXON: 
        int k;
        for ( k=0; k < Vector_getNumElement(bestExons); k++ ) {
          Exon *be = Exon_getElementAt(bestExons, k);

        //OTHEREXON: 
          int m;
          for (m=0 ; m < Vector_getNumElement(otherExons); m++) {
            Exon *oe = Exon_getElementAt(otherExons, m);

            if (Exon_getEnd(oe) < Exon_getStart(be)) {
              next OTHEREXON;
            }
// Big jump here - think about how to do
            if (Exon_getStart(oe) > Exon_getEnd(be)) {
              next BESTEXON;
            }

            // does it have exon overlap with the best model
            if (Exon_getStart(be) <= Exon_getEnd(oe) && 
                Exon_getEnd(be)   >= Exon_getStart(oe)) {
              // yes - store it and move on 
              Vector_addElement(genes, gene);
              
           //   print "Overlap " . $gene->start . " " , $gene->end . " " . $gene->strand ."\n";
// Big jump here - think about how to do
              next OTHERGENE;
            }
          }
        }
        // other model has no exon overlap with best model it needs to be in a new cluster
        Vector_addElement(otherGenes, gene);
        //#print "No overlap " . $gene->start . " " , $gene->end . " " . $gene->strand ."\n";
      }

      // now we need to fix the clusters
      if (Vector_getNumElement(otherGenes) > 0) {
        Vector_addElement(finalClusters, RefineSolexaGenes_recalculateCluster(rsg, genes),
        Vector_addElement(newClusters,   RefineSolexaGenes_recalculateCluster(rsg, otherGenes),
// NIY: Free strandedCluster???
      } else {
        // keep it as it was
        Vector_addElement(finalClusters, strandedCluster);
      }
// NIY: Free stuff???

    }
  }

  // print "CLUSTERS " . scalar(@final_clusters) ." " .  scalar(@new_clusters) ."\n";
  return (\@final_clusters,\@new_clusters);
}

GeneCluster *RefineSolexaGenes_recalculateCluster(RefineSolexaGenes *rsg, Vector *genes) {

  GeneCluster *cluster = GeneCluster_new();
  Gene *firstGene = Vector_getElementAt(genes, 0);

  long start   = Gene_getStart(firstGene);
  long end     = Gene_getEnd(firstGene);
  long strand  = Gene_getStrand(firstGene);
  double score = 0;

  int i;
  for (i=0; i<Vector_getNumElement(genes); i++) {
    Gene *g = Vector_getElementAt(genes, i);
    if (Gene_getStart(g) < start) {
      start = Gene_getStart(g);
    }
    if (Gene_getEnd(g) > end) {
      end = Gene_getEnd(g);
    }

    Transcript *t = Gene_getTranscriptAt(g, i);
    if (Transcript_getScore(t) > score) {
      score = Transcript_getScore(t);
    }
  }

  cluster->start  = start;
  cluster->end    = end;
  cluster->strand = strand;

  for (i=0; i<Vector_getNumElement(genes); i++) {
    Gene *g = Vector_getElementAt(genes, i);
    Transcript *t = Gene_getTranscriptAt(g, i);

    if (Transcript_getScore(t) == score ) {
      Gene_setBiotype(g, RefineSolexaGenes_getBestScoreType(rsg));
      // print "BEST SCORE " . $g->start ." " . $g->end . " " . $g->strand . "\n";
    } else {
      Gene_setBiotype(g, RefineSolexaGenes_getOtherIsoformsType(rsg));
    }
  }

  cluster->finalModels = genes;

  return cluster;
}


/*
=head2 filter_models

    Title        :   filter_models
    Usage        :   $self->filter_models(\@models);
    Returns      :   nothing
    Args         :   Array ref of array references of Bio::EnsEMBL::Transcript
    Description  :   Labels or removes models overlapping better scoring models on the 
                     opposite strand

=cut
*/
// NOTE: Modifies output vector of RunnableDB
// Perl had models instead of clusters but then immediately assigned clusters = models, but then uses models for something else, so just call
// it clusters from start
void RefineSolexaGenes_filterModels(RefineSolexaGenes *rsg, Vector *clusters) {

  Vector *fwd = Vector_new();
  Vector *rev = Vector_new();

  Vector *models = Vector_new();

  int i;
  for (i=0; i<Vector_getNumElement(clusters); i++) {
    GeneCluster *cluster = Vector_getElementAt(clusters, i);

    if (cluster->finalModels == NULL) {
      fprintf(stderr,"Had a cluster with no final models\n");
      continue;
    }

    if (cluster->strand == 1) {
      Vector_addElement(fwd, cluster);
    }
    if (cluster->strand == -1) {
      Vector_addElement(rev, cluster);
    }
    Vector_addElement(models, cluster);
  }

  // overlaps
  for (i=0; i<Vector_getNumElement(fwd); i++) {
    GeneCluster *fc = Vector_getElementAt(fwd, i);
    int j;
    for (j=0; j<Vector_getNumElement(rev); j++) {
      GeneCluster *rc = Vector_getElementAt(fwd, i);

      // one is within the other  or they are the same
      // they proably need to be rejected on the basis of coding overlap
      if (( fc->start >= rc->start && fc->end <= rc->end) || 
          ( rc->start >= fc->start && rc->end <= fc->end))  {
        
        // do they have coding overlap?
        Vector *fgs = fc->finalModels;
        Vector *rgs = rc->finalModels;

        // do they have coding overlap?
//      FG: 
        int k;
        for (k=0; k<Vector_getNumElement(fgs); k++) {
          Gene *fg = Vector_getElementAt(fgs, k);
          Transcript *ft = Gene_getTranscriptAt(fg, 0);
          
// Next what?????
          next unless $ft->translateable_seq;
          if (  $ft->translation->length <=  100 ) {
            Gene_setBioType(fg, "bad");
            // Do an else instead next FG;
          } else {

            int m;
            Vector *ftTranslateableExons = Transcript_getAllTranslateableExons(ft);
            for (m=0; m<Vector_getNumElement(ftTranslateableExons); m++) {
              Exon *fe = Vector_getElementAt(ftTranslateableExons,  m);

//            RG: 
              int n;
              for (n=0; n<Vector_getNumElement(rgs); n++) {
                Gene *rg = Vector_getElementAt(rgs, n);
                Transcript *rt = Gene_getTranscriptAt(rg, 0);

// Next what?????
                next unless $rt->translateable_seq;
                if (  $rt->translation->length <=  100 ) {
                  Gene_setBioType(rg, "bad");
                  // Do an else instead next RG;
                } else {
                  int p;
                  Vector *rtTranslateableExons = Transcript_getAllTranslateableExons(rt);
                  for (p=0; p<Vector_getNumElement(rtTranslateableExons); p++) {
                    Exon *re = Vector_getElementAt(rtTranslateableExons, p);

                    if ( Exon_getStart(fe) <= Exon_getEnd(re) && 
                         Exon_getEnd(fe)  >=  Exon_getStart(re)) {
                      // coding overlap        

// Hack hack hack
                      if (Transcript_getScore(ft) < Transcript_getScore(rt)) {
                        //# get rid of / label the reverse genes 
                        Gene_setBioType(fg, "bad");
                      } else {
                        Gene_setBioType(rg, "bad");
                      }
// Big jump - think about how to do
                      next FG;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  for (i=0; i<Vector_getNumElement(models); i++) {
    GeneCluster *cluster = Vector_getElementAt(models, i);
    my %exon_use_hash;
    my %exon_starts;
    my %exon_ends;
    my %exon_pattern;

    int count = 0;
    long translationStart = 2000000000;
    long translationEnd   = 0;

    int j;
    for (j=0; j<Vector_getNumElement(cluster->finalModels); j++) {
      Gene *gene = Vector_getElementAt(cluster->finalModels, j);
      Transcript *transcript = Gene_getTranscriptAt(gene, 0);

     if ( $transcript->translateable_seq ) {
        if (  $transcript->coding_region_start < $translation_start ) {
          $translation_start =  $transcript->coding_region_start;
        }
        if (  $transcript->coding_region_end > $translation_end ) {
          $translation_end =  $transcript->coding_region_end;
        }
      }

      my $exon_use = $transcript->{'_exon_use'};
      if ( $exon_use_hash{$exon_use}  ) {
        Gene_setBiotype(gene, "bad");
      }

      $exon_use_hash{$exon_use} = 1 ;
      //print "TRANSCRIPT " . $gene->get_all_Transcripts->[0]->{'_depth'} .           " Exon use $exon_use Biotype " . $gene->biotype ."\n";
      int es = 0;
      int ee = 0;
      my $pattern;
      foreach my $exon ( @{$gene->get_all_Exons} ) {
        $pattern .= $exon->start.":".$exon->end.":";
        $es++ if $exon_starts{$exon->start};
        $ee++ if $exon_ends{$exon->end};
        $exon_starts{$exon->start} = 1;
        $exon_ends{$exon->end} = 1;
      }
      if ( $ee == scalar(  @{$gene->get_all_Exons} ) &&
           $es == scalar(  @{$gene->get_all_Exons} ) ) {
        // seen it before - or something very much like it
        Gene_setBiotype(gene, "bad");
        //        print "CALLING it bad\n";
      }
      if ( $exon_pattern{$pattern} ) {
        // seen it before - or something very much like it
        Gene_setBiotype(gene, "duplicate");
        //print "CALLING it bad\n";
      }
      $exon_pattern{$pattern} = 1;
    } 

    // promote "bad" models that have a cds as long as the best cds to 
    // alt isoforms
    my @final_models = @{$cluster->{'final_models'}};
    my $best_cds = 0;
    for (  my $g = 0; $g < scalar(@final_models) ; $g++ ) {
      my $gene = $final_models[$g];
      my $transcript =  $gene->get_all_Transcripts->[0];
      fprintf(stdout, "%d - %f tran length %ld\n", 
              g, Transcript_getScore(transcript), Transcript_getcDNACodingEnd(transcript) - Transcript_getcDNACodingStart(transcript));

      if (g == 0) {
        // best scoring model 
        if  ( $transcript->translateable_seq ) {
          $best_cds =  $transcript->cdna_coding_end - $transcript->cdna_coding_start;
        }
      }
      if ( $transcript->translateable_seq ) {
        if ( $gene->biotype eq 'bad' && 
             $transcript->coding_region_start == $translation_start &&
             $transcript->coding_region_end == $translation_end ) {
          Gene_setBiotype(gene, RefineSolexaGenes_getOtherIsoformsType(rsg));
        }
        if (!strcmp(Gene_getBiotype(gene), "bad") && 
            $transcript->cdna_coding_end - $transcript->cdna_coding_start  > $best_cds ) {
          Gene_setBiotype(gene, RefineSolexaGenes_getOtherIsoformsType(rsg));
        }
      } 
      if (!strcmp(Gene_getBiotype(gene), "bad")) {
        // change type to  a bad model if it is bad 
        // store it on output array if the bad type is defined
        if ( $self->BAD_MODELS ) {
          $gene->biotype( $self->BAD_MODELS ) ;
          push @{$self->output} , $gene if $count <= $self->OTHER_NUM ;
        }
      } else {
        // if NOT duplicate 
        if (strcmp(Gene_getBiotype(gene), "duplicate")) {
          if ( $gene->biotype eq $self->BEST_SCORE ) {
            // trim the UTR
            $self->prune_UTR($gene);
            push @{$self->output} , $gene ;
          } else {
            if ( $self->OTHER_NUM  && $self->OTHER_ISOFORMS && $count <= $self->OTHER_NUM ) {
              // trim the UTR
              $self->prune_UTR($gene);
              push @{$self->output} , $gene;
            }
          }
        }
      }
      $count++ if $gene->biotype eq $self->OTHER_ISOFORMS ;
      $count++ if $gene->biotype eq $self->BEST_SCORE ;
    }
  }
}


/*
=head2 make_models

    Title        :   
    Usage        :   $self->make_models($paths, $strand ,$exons,$gene, $intron_hash);
    Returns      :   Array ref of array references of Bio::EnsEMBL::Transcript
    Description  :   Turns abstract paths into Bio::EnsEMBL::Gene models. Paths are
                     clustered and sorted by score - only the top X models for 
                     each cluster of paths get built ( X is defined in config )

=cut
*/

sub RefineSolexaGenes_makeModels(RefineSolexaGenes *rsg, StringHash *paths, int strand, Vector *exons, Gene *gene, StringHash *intronHash) {
  // paths are stored as text - turn them into arrays of features "models"
  my @clusters;
  my @models;
  my @genes;

  foreach my $path ( keys %$paths ) {
    my $exon_use;
    my @model;
    double exonScore = 0;
    double intronScore = 0;

    foreach my $feature ( split(/\./,$path ) ) {
      if ( $feature =~ /canonical/ ) {
        push @model, $intron_hash->{$feature};
        $intron_score+= $intron_hash->{$feature}->score;
      } else {
        $exon_use.= "$feature,";
        push @model, $exons->[$feature];
        foreach my $daf ( @{$exons->[$feature]->get_all_supporting_features} ) {
          $exon_score += $daf->score;
        }
      }
    }

    double totalScore = int($exon_score)/100 + intronScore;
    // last elements are the strand and score
    push @model, $exon_use;
    push @model, $total_score;
    push @models,\@model;
  }

  print "Starting model_cluster\n";
  // now lets cluster the models so that they are non overlapping
  // and return the clusters arranged by score
  my @model_clusters = @{$self->model_cluster(\@models, $strand)};
  
  print "Starting gene cycle\n";
  // Now we cycle through all the models and turn them into genes
  // we start with the highest scoring modes and work backwards
  // until we have enough 
  int clusterCount = 0;
  foreach my $cluster ( @model_clusters ) {
    clusterCount++;
    my @trans;
    int version = 0;

    my $strand = $cluster->{'strand'};
    // we want the array in the reverse order, highest scoring first
    // the score is the last array element 
    my @models_by_score =sort {$b->[-1] <=> $a->[-1]}  @{$cluster->{'models'}};
    
    // all the models with a particualar score highest first
//  MODEL:   
    foreach my $model (@models_by_score) {
      my @ises;
      // the score is the last array element 
      my $s = pop(@{$model});
      // list of the rough exons used in the model
      my $exon_use = pop(@{$model});

      my @introns;

      int intronCount = 0;
      double intronScore = 0;
      double exonScore = 0;
      int nonConIntrons = 0;

      my @new_exons;

      // make an array containing cloned exons
      for ( my $i = 0 ; $i < scalar(@$model) ; $i++ ) {
        unless ( $model->[$i]->isa("Bio::EnsEMBL::DnaDnaAlignFeature") ) {        
          my $new_exon = clone_Exon($model->[$i]);
          // add in exon coverage scores from supporting features
          foreach my $daf ( @{$model->[$i]->get_all_supporting_features} ) {
            $exon_score += $daf->score;
          }
          $new_exon->strand($strand);
          push @new_exons,$new_exon;
        } else {
          push @new_exons, $model->[$i];
        }
      }

      // trim the exons using the intron features to get the splice sites correct
      for ( my $i = 0 ; $i < scalar(@$model) ; $i++ ) {
        if ( $model->[$i]->isa("Bio::EnsEMBL::DnaDnaAlignFeature") ) {
          my $intron = $model->[$i];
          next unless $intron->strand == $strand;
          next unless $new_exons[$i-1] && $new_exons[$i+1];
          push @introns,$intron;
          // its an intron trim the exons accordingly
          $new_exons[$i-1]->end( $intron->start );
          $new_exons[$i+1]->start( $intron->end );
          if ( $new_exons[$i-1]->start >=  $new_exons[$i-1]->end ) {
            next MODEL;
          }
          $intron_count++;
          $intron_score+= $intron->score;
          $non_con_introns++ if $intron->hseqname =~ /non canonical/;
          // use the new intron feature code to store introns
          // provided we have not seen them before
          my $if;
          if ( $strand == 1 ){
            $if =  Bio::EnsEMBL::Intron->new( $new_exons[$i-1] , $new_exons[$i+1] );
          } else {
            $if =  Bio::EnsEMBL::Intron->new( $new_exons[$i+1] , $new_exons[$i-1] ); 
          }
          my $ise = Bio::EnsEMBL::IntronSupportingEvidence->new(
                                                                  -ANALYSIS => $intron->analysis,
                                                                -INTRON   => $if, 
                                                                -HIT_NAME => $intron->hseqname,
                                                                -SCORE    => $intron->score,
                                                                -SCORE_TYPE  => 'DEPTH',
                                                               );
          if ( $intron->hseqname =~ /non canonical/ ) {
            $ise->is_splice_canonical(0);
          }
          push @ises, $ise if $ise;
        }
      }
      next MODEL unless $intron_count;
      
      // trim padding from the start and end exons
      $new_exons[0]->start($new_exons[0]->start + 20) ;
      $new_exons[-1]->end ($new_exons[-1]->end  - 20) ;

      // get rid of impossibly small exons
      foreach my $e ( @new_exons){
        if ( $e->end - $e->start <= 0 ) {
          next MODEL;
        }
      }
      
      // trim away strings of Ns from the 1st and last exons
      // use same regex for 1st and last exon and reverse the
      // sequence accordingly depending on the strand
      my $ex_seq = $new_exons[0]->seq->seq;
      $ex_seq = reverse $ex_seq if $strand == - 1;
      if ( $ex_seq =~ /^(N+)(\S+)$/ ) {
        $new_exons[0]->start($new_exons[0]->start + length($1));
      }
      $ex_seq = $new_exons[-1]->seq->seq;
      $ex_seq = reverse $ex_seq if $strand ==  1;
      if ( $ex_seq =~ /^(N+)(\S+)$/ ) {
        $new_exons[-1]->end($new_exons[-1]->end - length($1));
      }
     
      // get rid of impossibly small exons again after the N trimming
      foreach my $e ( @new_exons){
        if ( $e->end - $e->start <= 0 ) {
          next MODEL;
        }
      }

      // make it into a gene
      my @modified_exons;
      foreach my $exon ( @new_exons ) {
        next if $exon->isa("Bio::EnsEMBL::DnaDnaAlignFeature");
        push @modified_exons, clone_Exon($exon);
      }
      if ( $strand == 1 ) {
        @modified_exons = sort { $a->start <=> $b->start } @modified_exons;
      } else {
        @modified_exons = sort { $b->start <=> $a->start } @modified_exons;
      }

      // make it into a gene
      my $t =  new Bio::EnsEMBL::Transcript(-EXONS => \@modified_exons);
      foreach my $ise ( @ises ) {
        $t->add_IntronSupportingEvidence($ise);
      }

//# check for dna
//#      my $check = $t->seq->seq ;
//#      my $Ns =  $check =~  s/N//g;
//#      if( length($t->seq->seq) == $Ns ){
//#        $self->throw("There does not appear to be ay DNA in the database, transcript seq is all N's\n");
//#      }
      // add a translation 
      my $tran = compute_translation(clone_Transcript($t));
      // stop spam coming from the Exon module
      $tran->dbID(0) ;        
      // store the introns along with the transcript so we can use them later for UTR trimming 
      $tran->{'_introns'} = \@introns;
      // keep track of the scores for this transcript
      $tran->analysis($self->analysis);
      $tran->version(1);

      // favor longer cds by adding doubling the score for coding exons
      // only use exons that are completely coding otherwise you also 
      // end up adding in score which is really UTR for long terminal exons
      // that have a bit of coding in them
      double codingBonus = 0;
      int codingExons = 0;
      if ( Transcript_getTranslation(tran) != NULL ) {
        foreach my $ce ( @{$tran->get_all_Exons} ) {
          unless ( $ce->phase == -1 or $ce->end_phase == -1 ) {
            $coding_bonus += $ce->get_all_supporting_features->[0]->score;
            codingExons++;
          }
        }
      }
      
      // remove any supporting features from the exons
      foreach my $e ( @{$tran->get_all_Exons} ) {
        $e->flush_supporting_features;
      }

// There's a bunch of hackery here adding random stuff to the transcript

      // print "Coding Bonus of $coding_bonus from $coding_exons completely coding exons \n";
      $tran->{'_score'} =  ( (int ( $intron_score + $exon_score ) / 10 ) + $coding_bonus  );
      // print "Final score = $intron_score + int( $exon_score / 100 ) + $coding_bonus = " . $tran->{'_score'} ;
      $tran->{'_depth'} =  ( $intron_score + $exon_score );
      // print " for tran " .$tran->{'_depth'} . "\n";
      $tran->{'_NC_introns'} =  $non_con_introns ;
      $tran->{'_exon_use'} = $exon_use;
      // print STDERR " EXON count $exon_count\n";
      $tran->{'_intron_count'} = $intron_count;

      push @trans, $tran;

      // we want X number of models
      if ( $self->BEST_SCORE &&  $self->MAX_NUM  ) {
        last MODEL if scalar(@trans)  >= ( $self->MAX_NUM +1 )  ;
      }
    }

    // re-sort the transcripts to take account of the revised scores
    Vector_sort(trans, SeqFeature_revScoreSortFunc);
    //@trans = sort { $b->{'_score'} <=> $a->{'_score'} } @trans;

    Transcript *best;
    int j;
    for (j=0; j<Vector_getNumElement(trans); j++) {
      Transcript *tran = Vector_getElementAt(trans, j);

      my ( $new_gene ) = @{convert_to_genes(($tran),$gene->analysis)};

      version++;
      Gene_setBiotype(newGene, RefineSolexaGenes_getOtherIsoformsType(rsg));

      if (version == 1) {
        Gene_setBiotype(newGene, RefineSolexaGenes_getBestScoreType(rsg));
        best = tran;
      }

      char stableId[2048];
      sprintf(stableId, "%s-v%d.%d-%d-%d-%d-NC-%d-%d",
              Gene_getStableId(gene),
              clusterCount,
              version,
              (int)(Transcript_getScore(trans)),
              (int)(Transcript_getDepth(trans)),
// NIY - not sure how to do these
              $tran->{'_intron_count'},
              $tran->{'_NC_introns'}, 
              Transcript_getStrand(tran));
      Gene_setStableId(newGene, stableId);
      Vector_addElement(cluster->finalModels, newGene);
    }
  }

  printf("Done gene cycle\n");
  return modelClusters;
}

Gene *RefineSolexaGenes_pruneUTR(RefineSolexaGenes *rsg, Gene *gene) {

  if (!RefineSolexaGenes_trimUTR(rsg)) {
    return gene;
  }

  my $transcript = $gene->get_all_Transcripts->[0];
  unless ( $transcript->translateable_seq ) {
    return $gene;
  }
 
  // fetch introns 
  my $introns = $transcript->{'_introns'};

  // otherwise trim the UTR according to the values set out in the config
  my %intron_hash;

  foreach my $intron ( @{$introns} ) {
    my $key = $intron->start  .":". $intron->end .":". $intron->strand;
    $intron_hash{$key} = $intron;
  }
  
  my @new_fivep;
  my @new_threep;
  my @new_exons;
  my @features;
  my @exons = sort {$a->start <=> $b->start }  @{$transcript->get_all_Exons};
  
  // put everything into the features array
  push @features, @exons;
  for ( my $i =0 ; $i < $#exons  ; $i++ ) {
    my $key = ($exons[$i]->end) .":". ($exons[$i+1]->start ) . ":" . $exons[$i]->strand;    
    if ( my $intron = $intron_hash{$key}  ) {
      push @features, $intron;
    }
  }
  @features = sort { $a->start <=> $b->start } @features;
  // so now we should have an array of alternating introns and exons
  print "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
Trimming UTR
Transcript " .  $transcript->seq_region_name ." " . 
    $transcript->start ." " . 
      $transcript->end ." " . 
        $transcript->strand ." " . 
          scalar(@{$transcript->get_all_Exons}) ."

";
  throw("Something is wrong we are missing " . scalar(@{$introns}) ." introns " . scalar(@exons) . "  exons " . scalar(@features) . " exons and introns\n")
    unless scalar(@features) == (scalar(@exons) * 2) -1 ;
  my $average_intron = 0;
  my $intron_count = 0;
  // leave single exon genes alone for now
  if ( scalar(@features) == 1 or scalar(@{$transcript->get_all_translateable_Exons}) == 1 )  {
    // lets strip the UTR
    my $trimmed_tran =  $self->modify_transcript($transcript,$transcript->get_all_translateable_Exons);
    // The naughty bit!
    $gene->{_transcript_array} = [];
    $gene->add_Transcript($trimmed_tran);
    return   $gene;
  }
  // first calculate the average
  foreach my $f ( @features ) {
    if ( $f->isa("Bio::EnsEMBL::DnaDnaAlignFeature")) {
      $average_intron += $f->score;
      $intron_count++;
    }
  }
  $average_intron /= $intron_count;

  foreach my $f ( @features ) {
   //# print $f->start . "\t" . $f->end ."\t$f\t";
    if ( $f->isa("Bio::EnsEMBL::DnaDnaAlignFeature")) {
    //#  print $f->score;
      if ( $average_intron && ($f->score /  $average_intron) * 100 <= $self->REJECT_INTRON_CUTOFF ) {
        print " Potentially bad ";
      }
    }
  //#  print "\n";
  }
  throw("Something is wrong we are missing introns " . scalar(@exons) . "  exons  and $intron_count introns\n")
    unless $intron_count == scalar(@exons) -1 ;
  print  "Average intron depth = $average_intron \n";
  
  
  my @fivep;
  my @threep;
  my $coding =0 ;
  // need to account for strand
  @features = sort { $b->start <=> $a->start } @features if $transcript->strand == -1;

  for ( my $i = 0 ; $i <=  $#features ; $i += 2  ) {
    my $e = $features[$i];
    throw("Got a DNA align feature where I should have an exon\n")
      unless $e->isa("Bio::EnsEMBL::Exon");
    if ( $e->coding_region_start($transcript) ) {
      // first coding exon
      for ( my $j = 0 ; $j <= $i ; $j++ ) {
        push @fivep,$features[$j];
      }
      last;
    }
  }
  for ( my $i = $#features ; $i >= 0 ;  $i -= 2  ) {
    my $e = $features[$i];
    throw("Got a DNA align feature where I should have an exon\n")
      unless $e->isa("Bio::EnsEMBL::Exon");
    if ( $e->coding_region_end($transcript) ) {
          # first coding exon
      for ( my $j = $i ; $j <= $#features ; $j++ ) {
        push @threep,$features[$j];
      }
      last;
    }
  }
  
  // want to start at last coding exon and work outwards so....
  @fivep = reverse @fivep;
  // now we should be good
  print "FIVE P \n";
  @new_exons = @{$transcript->get_all_translateable_Exons};
  my $fivep_cds = shift(@new_exons);
  my $threep_cds = pop(@new_exons);
  my $fiveplen;
  my $threeplen;
  my $fivepc = 0 ;
  my $threepc = 0 ;
  my $nmd;
  
  // FIVE PRIME RULES
  
 FIVEP: for ( my $i = 0 ; $i <= $#fivep ; $i++ ) {
    my $f =  $fivep[$i];
    print $f->start . "\t" . $f->end ."\t$f\t";
    if ( $i == 0 ) {
      throw("First feature is not an exon \n")
        unless $f->isa("Bio::EnsEMBL::Exon");
      // UTR starts in this exon - how long is it?
      my $cds_start = $f->coding_region_start($transcript);
      $cds_start = $f->coding_region_end($transcript)  if $transcript->strand == -1;
      throw("First coding exon has no CDS \n") unless $cds_start;
      print "CDS START $cds_start\t";
      $fiveplen = $cds_start - $f->start +1 if $transcript->strand == 1;
      $fiveplen = $f->end - $cds_start   +1 if $transcript->strand == -1;
      // is the coding exon too long
      if ( $fiveplen > $self->MAX_5PRIME_LENGTH ) {
        // replace it with the cds
        @new_fivep = ($fivep_cds);
        print " 5p too long $fiveplen \n";
        last FIVEP;
      }
      push @new_fivep,$f;
      $fivepc++;
    } else {
      if (  $f->isa("Bio::EnsEMBL::Exon") ) {
        $fivepc++;
        $fiveplen+= $f->end - $f->start +1;
        // does it make the UTR too long?
        if ( $fiveplen > $self->MAX_5PRIME_LENGTH ) {
          // dont add it
          print " 5p too long $fiveplen \n";
          last FIVEP;
        }
        // is it too many exons?
        if ( $fivepc > $self->MAX_5PRIME_EXONS ) {
          // dont add it
          print " too many 5p  $fivepc cut them all as we are not sure \n";
          @new_fivep = ($fivep_cds);
          last FIVEP;
        }
        push @new_fivep,$f;
      }
    }
    // Does the intron score well enough to include the exon
    // apply rules and add successful exons into the mix
    if ( $f->isa("Bio::EnsEMBL::DnaDnaAlignFeature")) {
      print $f->score;
      if ( $average_intron && ($f->score /  $average_intron) *100 <= $self->REJECT_INTRON_CUTOFF ) {
        print " Rejecting on score cutoff " . $f->score ." vs  $average_intron\n";
        // dont add any more 
        last FIVEP;
      }
    }
    print "\n";
  }
  
  // three P
  print "THREE P \n";
 THREEP:   for ( my $i = 0 ; $i <= $#threep ; $i++ ) {
    my $f = $threep[$i];
    print $f->start . "\t" . $f->end ."\t$f\t";
    if ( $i == 0 ) {
      throw("First feature is not an exon \n")
        unless $f->isa("Bio::EnsEMBL::Exon");
      // UTR starts in this exon - how long is it?
      my $cds_end = $f->coding_region_end($transcript);
      $cds_end = $f->coding_region_start($transcript)  if $transcript->strand == -1;
      throw("last coding exon has no CDS \n") unless $cds_end;
      print "CDS END $cds_end\t";
      $threeplen = $cds_end - $f->start +1 if $transcript->strand == -1;
      $threeplen = $f->end - $cds_end   +1 if $transcript->strand == 1;
      # is the coding exon too long
      if ( $threeplen > $self->MAX_3PRIME_LENGTH ) {
        // replace it with the cds
        @new_threep = ($threep_cds);
        print " 3p too long $threeplen \n";
        last THREEP;
      }
      push @new_threep,$f;
      $nmd = $threeplen ;
      $threepc++;
    } else {
      if (  $f->isa("Bio::EnsEMBL::Exon") ) {
        // does it break the NMD rule?
        if ( $nmd > 55 ) {
          print " splice is after $nmd bp from stop codon - rejected on NMD rule of maximum 55 bp \n";
          @new_threep = ($threep_cds);
          last THREEP;
        }
        $threepc++;
        $threeplen+= $f->end - $f->start +1;
        // does it make the UTR too long?
        if ( $threeplen > $self->MAX_3PRIME_LENGTH ) {
          // dont add it
          print " 3p too long $threeplen \n";
          last THREEP;
        }
        // is it too many exons?
        if ( $threepc > $self->MAX_3PRIME_EXONS ) {
          // dont add it
          print " too many 3p  $threepc cut them all as we are not sure \n";
          @new_threep = ($threep_cds);
          last THREEP;
        }
        push @new_threep,$f;
      }
    }
    // Does the intron score well enough to include the exon
    // apply rules and add successful exons into the mix
    if ( $f->isa("Bio::EnsEMBL::DnaDnaAlignFeature")) {
      print $f->score;
      if ( $average_intron && ($f->score /  $average_intron) * 100 <= $self->REJECT_INTRON_CUTOFF ) {
        print " Rejecting on score cutoff " . $f->score ." vs  $average_intron\n";
        // dont add any more 
        last THREEP;
      }
    }
    print "\n";
  }
  
  push @new_exons, @new_fivep;
  push @new_exons, @new_threep;
  print " New transript has " . scalar(@new_exons) , " exons\n";
  my @clones;
  foreach my $e ( @new_exons ) {
    throw("Not is not an exon " . $e->start ." " . $e->end . " $e\n")
      unless $e->isa("Bio::EnsEMBL::Exon");
    push @clones, clone_Exon($e);
  }
  @clones = sort { $a->start <=> $b->start } @clones;
  @clones =  reverse(@clones) if $transcript->strand == -1;
  my $trimmed_tran =  $self->modify_transcript($transcript,\@clones);
  // The naughty bit!
  $gene->{_transcript_array} = [];
  $gene->add_Transcript($trimmed_tran);
  return $gene;
}


Transcript *RefineSolexaGenes_modifyTranscript(RefineSolexaGenes *rsg, Transcript *tran, Vector *exons) {
  long cdsStart;
  long cdsEnd;

  if (Transcript_getStrand(tran) == 1) {
    cdsStart = Transcript_getCodingRegionStart(tran);
    cdsEnd   = Transcript_getCodingRegionEnd(tran);
  } else {
    cdsEnd   = Transcript_getCodingRegionStart(tran);
    cdsStart = Transcript_getCodingRegionEnd(tran);
  }
  
  printf("CDS START END %ld %ld\n", cdsStart,  cdsEnd);
  Translation *tln = Transcript_getTranslation(tran); 
  printf("PHASE %ld %ld\n", Tranlation_getStart(tln), Translation_getEnd(tln));

// t created here and again below - not sure if need both
  Transcript *t = Transcript_new();
  Transcript_setExons(exons);

  Exon *se;
  Exon *ee;
  int i;
  for (i=0; i<Transcript_getExonCount(t); i++) {
    Exon *e = Vector_getExonAt(t, i);
    if (Exon_getStart(e) <= cdsEnd && Exon_getEnd(e) >= cdsEnd) {
      ee = e;
    }
    if (Exon_getStart(e) <= cdsStart && Exon_getEnd(e) >= cdsStart) {
      se = e;
    }
  }

  Transcript *ts;
  Transcript *te;
  if ( $tran->strand == -1 ) {
    ts =  Exon_getEnd(se) - cdsStart + 1;
    te =  Exon_getEnd(ee) - cdsEnd   + 1;
  } else {
    ts =   cdsStart - Exon_getStart(se) + 1;
    te =   cdsEnd   - Exon_getStart(ee) + 1;
  }

// Why created again?????
/*
  Transcript *t = Transcript_new();
  Transcript_setExons(exons);
*/

  // transfer the intron supporting evidence
  // except for where we have trimmed the intron
  Vector *intronSupport = Transcript_getAllIntronSupportingEvidence(tran);
  if (intronSupport!=NULL) {
    for (i=0; i<Vector_getNumElement(intronSupport, i); i++) {
      Vector *ise = Vector_getElementAt(intronSupport, i);
      if (IntronSupportingEvidence_getSeqRegionStart(ise) > Transcript_getStart(t) &&  
          IntronSupportingEvidence_getSeqRegionEnd(ise) < Transcript_getEnd(t)) {
        Transcript_addIntronSupportingEvidence(t, ise);
      }
    }
  }
  //#  my $start_phase = $se->phase;
  Translation *translation = Translation_new();
  Translation_setStartExon(se); 
  Translation_setEndExon(se); 
  Translation_setSeqStart(ts); 
  Translation_setSeqEnd(te); 

  printf("S-E %ld %ld\n",ts, te); //#START PHASE $start_phase\n";
  printf("GS %ld %ld\n", Translation_getGenomicStart(translation), Translation_getGenomicEnd(translation));
  Transcript_setTranslation(t, translation);
  //# calculate_exon_phases($t,$start_phase);

  unless (  $tran->translation->seq eq $t->translation->seq ) {
    $self->throw("Translations do not match: Before " . $tran->translation->seq ."\nAfter  " .           $t->translation->seq ."\n");
  }
  return t;
}

void RefineSolexaGenes_writeOutput(RefineSolexaGenes *rsg) {

  DBAdaptor *outdb = RefineSolexaGenes_getDbAdaptor(rsg, RefineSolexaGenes_getOutputDb(rsg));

  GeneAdaptor *geneAdaptor = DBAdaptor_getGeneAdaptor(outdb);
// NIY:
//  $outdb->dbc->disconnect_when_inactive(0);

  Vector *output = RefineSolexaGenes_getOutput(rsg);
  
  int fails = 0;
  int total = 0;
//  GENE: 
  int i;
  for (i=0; i<Vector_getNumElement(output); i++) {
    Vector *gene = Vector_getElementAt(output, i);
    Analysis *anal = RefineSolexaGenes_getAnalysis(rsg);

    Gene_setAnalysis(gene, anal);
    Gene_setSource(gene, Analysis_getLogicName(anal));

    int j;
    for (j=0; j<Gene_getTranscriptCount(gene); j++) {
      Transcript *tran = Gene_getTranscriptAt(gene, j);
      Transcript_setAnalysis(tran, anal);
    }

    // filter single exon genes that may have been made through UTR trimming
    my @exons = @{$gene->get_all_Exons};
    if ( scalar(@exons == 1 )) {
      if ( $self->SINGLE_EXON_MODEL ) {
        next GENE unless $exons[0]->length >= $self->MIN_SINGLE_EXON;
        Gene_setBiotype(gene, RefineSolexaGenes_getSingleExonModelType(rsg));
      } else {
        // dont store it
        continue;
      }
    }
    
// Was 'eval'ed
    GeneAdaptor_store(geneAdaptor, gene);
/*
    if ($@){
      warning("Unable to store gene!!\n$@");
      $fails++;
    }
*/
    total++;
  }

// Can't be any currently, but maybe should be possible to trap errors
  if (fails > 0) {
    fprintf(stderr, "Not all genes could be written successfully (%d fails out of %d)\n", fails, total);
  }

  DNAAlignFeatureAdaptor *intronAdaptor = DBAdaptor_getDNAALignFeatureAdaptor(outdb);
  fails = 0;
  total = 0;
 
  Vector *intronFeatures = RefineSolexaGenes_getIntronFeatures(rsg);

  for (i=0; i<Vector_getNumElement(intronFeatures); i++) {
    DNAAlignFeature *intron = Vector_getElementAt(intronFeatures, i);

# SMJS For now leave as exon end and exon start so that it edge matches in apollo
#    $intron->start($intron->start+1);
#    $intron->end($intron->end-1);
//    eval {
        DNAAlignFeatureAdaptor_store(intronAdaptor, intron);
//    };
/*
    if ($@){
      warning("Unable to store DnaAlignFeature!!\n$@");
      $fails++;
    }
*/
    total++;
  }
  if (fails > 0) {
    fprintf(stderr, "Not all introns could be written successfully (%d fails out of %d)\n", fails, total);
  }
}


/*
=head2 ProcessTree

    Title        :   ProcessTree
    Usage        :   $self->ProcessTree
    Returns      :   String containing paths through the gene
    Args         :   A hash reference contianing the possible intron exons
                 :   Integer key for the hashref
                 :   String containing keys used up to this point
                 :   String containing the paths under construction
    Description  :   Recursive method that creates paths that explore all possible
                 :   routes through a hashref, uses a configurable recursion limit 
                 :   to prevent it running out of memory if the paths cannot be solved
                 :   or are too large to be practical

=cut
*/

char *RefineSolexaGenes_processTree(RefineSolexaGenes *rsg, hashref, index, sofar, paths) {
  my ($self,$hashref,$index,$sofar,$paths) = @_;
  // dont let it go on for ever eating memory
  if ($limit > $self->recursive_limit){
    print STDERR  "Too many recursive possibilities\n";
    return "ERROR";
  }
  my @node =  keys %{$hashref->{$index}} ;
  $sofar.= "$index.";
  foreach my $child (@node){
    $limit++;
     my $result =  $self->ProcessTree($hashref,$child,$sofar,$paths);
    if ( $result eq "ERROR" ) {
      $limit = 0;
      return "ERROR"; 
    }
   //# $result->{$sofar} = 1;
  }
  if ( scalar(@node == 0) ) {
    //#print "$sofar\n";
    $paths->{$sofar} = 1;
  }
  return $paths;
}

/*
=head2 ProcessTree

    Title        :   ProcessTree
    Usage        :   $self->process_paths
    Returns      :   String containing paths through the gene
    Args         :   A hash reference contianing the possible intron exons
                 :   Integer key for the hashref
                 :   String containing keys used up to this point
                 :   Integer flag indicating filtering should take place
    Description  :   Filters paths to remove the lowest scoring intron
                 :   for a given pair of exons where more than one intron
                 :   is possible. Filters progressivley if the paths cannot be
                 :   made for the model until the paths can be created or the 
                 :   model cannot be filtered any more, in this case the number
                 :   of recursions can be raised and the process repeated
                 :   untill the max_recursions limit is reached

=cut
*/

char *RefineSolexaGenes_processPaths(RefineSolexaGenes *rsg, Vector *exons, ??? exonIntron, ??? intronExon, int strict) {

  my $variants;
  my $removed;

  //# now lets make a hash of hashes holding which exons connect to which
  int i;
  for (i=0; i<Vector_getNumElement(exons); i++) {
    Exon *exon = Vector_getElementAt(exons, i);
    if ( $exon_intron->[$i] ) {
      if (strict) {
        //# Throw out exons that have retained introns for a start
        next if  $exons->[$i]->{'retained'};
        if () {
          continue;
        }
      }

      if (strict > 1) {
        // group all the introns by exon pairs
        my %intron_groups;
        foreach my $intron ( @{$exon_intron->[$i]} ) {
          next if  $intron->hseqname =~ /REMOVED/  ;
          push @{$intron_groups{$i}}, $intron;
        }

        //# now lets sort these groups by score
        foreach my $group ( keys %intron_groups ) {
          @{$intron_groups{$group}}  = sort {$b->score <=> $a->score} @{$intron_groups{$group}};
        }

        //# now lets see what they look like
        fprintf(stderr,"EXON %d: %ld - %ld - %d\n",
                i, Exon_getStart(exon), Exon_getEnd(exon), Exon_getStrand(exon));

        foreach my $group ( keys %intron_groups ) {
          foreach my $intron ( @{$intron_groups{$group}} ) {
            print STDERR "$group " . $intron->hseqname . " " . $intron->score ."\n";
          }
          if ( scalar( @{$intron_groups{$group}} ) > 1 ) {
            //  remove the lowest scoring one
            my $intron = pop( @{$intron_groups{$group}} ) ;
            print STDERR "Eliminating " . $intron->hseqname . " " .
              $intron->score . "\n";
            unless (  $intron->hseqname =~ /REMOVED/ ) {
              $intron->hseqname($intron->hseqname."-REMOVED") ;
              $removed++;
            }
          }        
        }
      }
      
      foreach my $intron ( @{$exon_intron->[$i]} ) {
        //# only allow each intron to connect to 1 exon
        foreach my $exon ( @{$intron_exon->{$intron->hseqname}} ) {
          if ( $intron->end > $exons->[$exon]->end ) {
            next ;
          }
          //# store the possible paths as a hash (splice)variants
          $variants->{$i}->{$intron->hseqname} = 1;
          $variants->{$intron->hseqname}->{$exon} = 1;
          //# check 
          if ( $intron->end > $exons->[$exon]->end ) {
            throw(" exon $i start end " . $intron->end ." - " . $exons->[$exon]->start );
          }
        }
      }
    }
  }
  
  if (strict && !removed) {
    Exon *firstExon = Vector_getElementAt(exons, 0);

// ??? why &! rather than && ! if ($strict &! $removed ) {
    fprintf(stderr,"Cannot simplify this gene any more EXON 0: %ld - %ld - %d\n",
            Exon_getStart(firstExon), Exon_getEnd(firstExon), Exon_getStrand(firstExon));

    if (RefineSolexaGenes_getRecursiveLimit(rsg) < RefineSolexaGenes_getMaxRecursions(rsg)) {
      RefineSolexaGenes_setRecursiveLimit(rsg, RefineSolexaGenes_getRecursiveLimit(rsg) * 10);
      fprintf(stderr, "Upping recursive limit to %d to see if it helps\n", RefineSolexaGenes_getRecursiveLimit(rsg));
    } else {
      fprintf(stderr,"Giving up on EXON 0: %ld - %ld - %d\n",
              Exon_getStart(firstExon), Exon_getEnd(firstExon), Exon_getStrand(firstExon));
//!!!!!!!!!!! NIY What to return
      return "Give up";
    }
  }
  
  // work out all the possible paths given the features we have
  char *result;
  my %paths;
  for (i=0 ; i<Vector_getNumElement(exons); i++) {
    limit = 0;
    result = RefineSolexaGenes_processTree(rsg, variants, i, NULL, ??? paths);

    if (!strcmp(result,"ERROR")){
      fprintf(stderr, "Could not process cluster %d trying again with simpler cluster\n", i);
      return NULL;
    }
  }

  return result;
}

/*
=head2 model_cluster

    Title        :   model_cluster
    Usage        :   $self->model_cluster($models,$strand);
    Returns      :   2D array ref of exons and intron features
    Args         :   Array ref of  exons and intron features
                 :   Integer indicating strand
    Description  :   Clusters the initial models by start end
                 :   orders the models in each cluster by score

=cut
*/

typedef struct GeneModelStruct {
  Exon *firstExon;
  Exon *lastExon;
  int exonIUse;
  double score;
} GeneModel;

Vector *RefineSolexaGenes_modelCluster(RefineSolexaGenes *rsg, Vector *models, int strand) {
  Vector *clusters = Vector_new();

  // sort them by the start of the fist exon ( the first array element )
  Vector_sort(models, GeneCluster_posSortFunc);

  //my @models = sort { $a->[0]->start <=> $b->[0]->start }  @$models ;

  //# $model->[0]  = 1st exon 
  //# $model->[-3] = last exon 
  //# $model->[-2] = exon iuse
  //# $model->[-1] = score
  
  printf("Have %d to cluster\n", Vector_getNumElement(models));

  int startInd = 0;

  int i;
  for (i=0; i<Vector_getNumElement(models); i++) {
    GeneModel *model = Vector_getElementAt(models, i);
    int clustered = 0;

    int nClust = Vector_getNumElement(clusters);

    int j;
    for (j=startInd; j<nClust; j++) {
//#    foreach my $cluster ( @clusters ) {
      // do they overlap?
      GeneCluster *cluster = Vector_getElementAt(clusters, j);

      if ( Exon_getStart(model->firstExon) <= cluster->end &&  Exon_getEnd(model->lastExon) >= cluster->start) {
        // Expand the cluster
        if (Exon_getStart(model->firstExon) < cluster->start) cluster->start = Exon_getStart(model->firstExon);
        if (Exon_getEnd(model->lastExon)    > cluster->end)   cluster->end   = Exon_getEnd(model->lastExon);

        Vector_addElement(cluster->models, model);
        clustered = 1;
      } else if (Exon_getStart(model->firstExon) > cluster->start) {
        startInd++;
      }
    }

    if (!clustered) {
      GeneCluster *cluster = GeneCluster_new();

// NIY: Do we need to add
      cluster->models = Vector_new();

      Vector_addElement(cluster->models, model);

      cluster->start  = Exon_getStart(model->firstExon);
      cluster->end    = Exon_getEnd(model->lastExon);
      cluster->strand = strand;

      Vector_addElement(clusters, cluster);
    }
  }

  printf("Have %d after clustering them\n", Vector_getNumElement(clusters));
  return clusters;
}


/*
=head2 merge_exons

    Title        :   merge_exons
    Usage        :   $self->merge_exons($gene)
    Returns      :   Array ref of Bio::EnsEMBL::Exon
    Args         :   Bio::EnsEMBL::Gene
    Description  :   Merges adjacent exons where the intron is covered by repeats or
                 :   is very small

=cut
*/
// lets us merge exons with tiny  introns between them  unless they contain an intron

Vector *RefineSolexaGenes_mergeExons(RefineSolexaGenes *rsg, Gene *gene, int strand) {

  my @exons;
  next unless $gene->get_all_Transcripts->[0];
  foreach my $exon ( @{$gene->get_all_Transcripts->[0]->get_all_Exons} ) {
    push @exons, clone_Exon($exon);
  }

  my $ec = scalar(@exons) ;
  my $extra_exons = $self->extra_exons;

  // the extra exon is a list of start end coords of the spliced intron sections
  // ie: end:start:end:start where the 1st and last coords are anchors to tie it 
  // into our rough model both must match before we can try and add any potentialy 
  // novel exons in
  my @sortedexons =  sort { $a->start <=> $b->start } @exons;
  foreach my $key ( keys %$extra_exons ) {
    my @coords = split(/:/,$key);
    my $start_anchor = shift(@coords);
    my $end_anchor = pop(@coords);
    #print "START AND END $start_anchor  $end_anchor \n";
    if ($start_anchor > $end_anchor) {
      fprintf(stderr, "START AND END %d %d\n", startAnchor,  endAnchor);
      exit(1);
    }
    //# do the anchors lie within the model?
    //# SMJS Did some optimisation here
//#   foreach my $exon ( @exons ) {
//#      if ( $start_anchor <= $exon->end &&
//#           $start_anchor >= $exon->start ) {
//#        $start_anchor = -1;
//#      }
//#      if ( $end_anchor <= $exon->end &&
//#           $end_anchor >= $exon->start ) {
//#        $end_anchor = -1;
//#      }
//#    }

    if ($self->bin_search_for_overlap(\@sortedexons, $start_anchor)) {
      $start_anchor = -1;
    }
    if ($self->bin_search_for_overlap(\@sortedexons, $end_anchor)) {
      $end_anchor = -1;
    }
       
//#    my $i;
//#    for ($i=0;$i<scalar(@sortedexons); $i++) {
//#      my $exon = $sortedexons[$i];
//#      if ( $start_anchor <= $exon->end && 
//#           $start_anchor >= $exon->start ) {
//#        $start_anchor = -1;
//#        last;
//#      } elsif ($exon->start > $start_anchor) {
//#        last;
//#      }
//#    }
//#
//#    for ( ;$i<scalar(@sortedexons); $i++) {
//#      my $exon = $sortedexons[$i];
//#      if ( $end_anchor <= $exon->end && 
//#           $end_anchor >= $exon->start ) {
//#        $end_anchor = -1;
//#        last;
//#      } elsif ($exon->start > $end_anchor) {
//#        last;
//#      }
//#    }
    if ( $start_anchor == -1 && $end_anchor == -1 ) {
      // now to make new the exon(s)
      for ( my $i = 0 ; $i <= $#coords ; $i += 2 ) {
        my $left = $coords[$i];
        my $right = $coords[$i+1];
        my $extra = $self->make_exon(undef,$left,$right,$extra_exons->{$key},$key );
        $extra->{"_extra"} = 1;
        push @exons,$extra;
      }
    }
  }
//  print "After extras add have " . scalar(@exons) . "\n";
  
//  print "Merging exons - done extras\n";

  
//#  @exons =  sort { $a->start <=> $b->start } @exons;
//#  # want to get rid of any overlapping exons
//#  while  ( $ec != scalar(@exons) ) {
//#    $ec = scalar(@exons);
//#    for ( my $i = 1 ; $i < scalar(@exons) ; $i++ ) {
//#      my $left_exon = $exons[$i-1];
//#      my $right_exon = $exons[$i];
//#      # do they overlap 
//#      if ( $left_exon->start <= $right_exon->end && 
//#           $left_exon->end >= $right_exon->start ) {
//#        # merge them 
//#        if (   $right_exon->end >= $left_exon->end &&
//#               $right_exon->start <= $left_exon->start ){
//#          print "HERE\n";
//#          $left_exon->{"_extra"} = 0;
//#        }
//#        $left_exon->start($right_exon->start) 
//#          if $right_exon->start < $left_exon->start;
//#        $left_exon->end($right_exon->end) 
//#          if $right_exon->end > $left_exon->end;
//#        # get rid of right exon
//#        splice(@exons,$i,1);
//#        $i-- ;
//#        @exons =  sort { $a->start <=> $b->start } @exons;
//#      }
//#    }
//#  }

  @exons =  sort { $a->start <=> $b->start } @exons;
  for ( my $i = 1 ; $i < scalar(@exons) ; $i++ ) {
     my $left_exon = $exons[$i-1];
     my $right_exon = $exons[$i];
     // do they overlap 
     if ( $left_exon->start <= $right_exon->end && 
          $left_exon->end >= $right_exon->start ) {
      // merge them 
      if ( $right_exon->end >= $left_exon->end &&
           $right_exon->start <= $left_exon->start ){
        $left_exon->{"_extra"} = 0;
      }
//# Don't see how this is possible - they are sorted on start!
//#      $left_exon->start($right_exon->start) 
//#        if $right_exon->start < $left_exon->start;
      
      if ($right_exon->start < $left_exon->start) {
        die "HOW IS THIS POSSIBLE\n";
      }

      $left_exon->end($right_exon->end) 
        if $right_exon->end > $left_exon->end;
      // get rid of right exon
      splice(@exons,$i,1);
      $i-- ;
//#      @exons =  sort { $a->start <=> $b->start } @exons;
    }
  }
  
  printf("Merging exons - done overlap filter\n");
  
  long offset = 0;
  for (i = 1; i<Vector_getNumElement(exons); i++) {
    Exon *exon     = Vector_getElementAt(exons, i);
    Exon *prevExon = Vector_getElementAt(exons, i-1);

    int intronCount = 0;

    Vector *introns = RefineSolexaGenes_fetchIntronFeatures(rsg, Exon_getEnd(prevExon), Exon_getStart(exon), &offset);
    
    // is the intron tiny?    
    // we know it lies across the boundary does it lie within the 2 exons?
    int j;
    for (j=0; j<Vector_getNumElement(introns); j++) {
      DNAAlignFeature *intron = Vector_getElementAt(introns, j);
 //#     print "INTRON " . $intron->start . " " . $intron->end . " " , $intron->strand ." " , $intron->score ."\n";

      //# ignore non consensus introns at this point
      if (!strstr(DNAAlignFeature_getHitSeqName(intron), "non canonical")) continue;

      if ($intron->start > $prev_exon->start &&
          $intron->end <  $exon->end && 
          $intron->strand == $strand ){
        intronCount++;
      }
    }

    // remove very small introns if there is no direct evidence for them
    if (Exon_getStart(exon) - Exon_getEnd(prevExon) <= 20  && intronCount == 0)   {
      Exon_setStart(exon, Exon_getStart(prevExon));

      //splice(@exons,$i-1,1);
      Vector_removeElementAt(exons, i-1);

// NIY: Free prevExon

      i--;
// Is this necessary???      next;
    }
  }
  printf("Merging exons done\n");
  
  return exons;
}

Exon *RefineSolexaGenes_binSearchForOverlap(RefineSolexaGenes *rsg, Vector *exons, int pos) {
  int iMin = 0;
  int iMax = Vector_getNumElement(exons) - 1;

   while (iMax >= iMin) {
    int iMid = (iMax+iMin) / 2;

    Exon *e = Vector_getElementAt(exons, iMid);

    if (pos > Exon_getEnd(e)) {
      iMin = iMid + 1;
    } else if (pos < Exon_getStart(e)) {
      iMax = iMid - 1;
    } else {
      return e;
    }
  }

  return NULL;
}

/*
=head2 bam_2_intron_features
    Title        :   bam_2_intron_features
    Usage        :   $self->bam_2_intron_features($segment)
    Returns      :   None
    Args         :   Bam file segement
    Description  :   Fetches all alignments from the bam file segment, collapses them down into a 
                 :   non redundant set and builds a Bio::EnsEMBL::DnaDnaAlignFeature to 
                 :   represent it, then stores it in $self->intron_features
                 :   analyses splice sites for consensus and non consensus splices as this data is 
                 :   not stored in the BAM.
                 :   Also checks for small exons defined by a single read splicing at least twice
                 :   stores any additional exons found this way in $self->extra_exons

=cut
*/

sub bam_2_intron_features {
  my ($self,$segment,$intron_files) = @_;
  my $slice_adaptor = $self->gene_slice_adaptor;
  my @ifs;
  my $extra_exons = $self->extra_exons;
  my %id_list;
  my %read_groups;
  if (  $intron_files->{GROUPNAME} && scalar(@{$intron_files->{GROUPNAME}} > 0 ) ) {
    my @groups = @{$intron_files->{GROUPNAME}};
    print "Limiting to read groups ";
    foreach my $group ( @groups ) {
      print " $group";
      $read_groups{$group} = 1;
    }
    print "\n";
  }


  my $iterator = $segment->features(-iterator=>1);
 READ:  while (my $read = $iterator->next_seq) {
    my $spliced;
    // ignore unspliced reads if the bam file is a mixture of spliced and 
    // unspliced reads
    if ( $intron_files->{MIXED_BAM} ) {
      $spliced = $read->get_tag_values('XS');
      next READ unless $spliced; 
    }
    // filter by read group if needed
    
    // need to recreate the ungapped features code as the
    // auto splitting code does not seem to work with > 2 features
    if ( $intron_files->{GROUPNAME}  && scalar(@{$intron_files->{GROUPNAME}} > 0 )) {
      next unless ($read_groups{$read->get_tag_values('RG')}) ;
    }
    my @mates = sort { $a->[2] <=> $b->[2] } @{$self->ungapped_features($read)};

    // if mates > 2 then we have a possibility of adding in some extra exons into our rough models
    // as the read has spliced into and out of an exon
    // lets make them unique
    if ( scalar(@mates) > 2 ) {
      my $string;
      for ( my $i = 0 ; $i <= $#mates  ; $i++ ) {

        my $start  = $mates[$i]->[2];
        my $end    = $mates[$i]->[3];
        my $hstrand = $read->strand;
        $string .= $start .":" if $i > 0 ;
        $string .= $end .":" if $i < $#mates ;
      }
      $extra_exons->{$string} ++;
      //# print "Not doing extra_exon stuff for now\n";
    }
    my $strand = $read->target->strand;
    if   ($intron_files->{MIXED_BAM} ) {
      $strand = 1 if $spliced eq '+';
      $strand = -1 if $spliced eq '-'; 
    } 
    //# print "\nREAD " . $read->cigar_str;
    my $offset;
    for ( my $i = 0 ; $i <= $#mates  ; $i++ ) {
      //#   print "\n";
      // intron reads should be split according to the CIGAR line
      // the default split function seems to ad
      // we want the ungapped features to make our introns
      my $name   = $read->seq_id;
      // we dont allow . in the seq region name as we use them to delimit our paths
// NOT REALLY COMMENTED OUT - JUST TO STOP SYNTAX HIGHLIGHTER GOING NUTS      $name =~ s/\./*/g;   
      my $start  = $mates[$i]->[2];
      my $end    = $mates[$i]->[3];
      my $cigar  = $mates[$i]->[4];
      my $hstrand = $read->strand;
   //#   print "$name $start $end $strand $hstrand $cigar\t";
   //#   print $read->get_tag_values('FIRST_MATE') ."\t";
      next if $i == $#mates;
      my $unique_id = $name . ":" . 
        $mates[$i]->[3] . ":" .
          $mates[$i+1]->[2] . ":" . 
            $strand ;
      $id_list{$unique_id} ++;
  //#    print "$unique_id";
    }
 //#   print "\n";
  }

  // For param testing store introns with different anal to results
  my $conslim = $ENV{CONSLIM};
  my $nonconslim = $ENV{NONCONSLIM};
  my $intron_anal = $self->create_analysis_object("intron_c" . $conslim . "_nc" . $nonconslim);

  //# collapse them down and make them into simple features
  foreach my $key ( keys %id_list ) {
    // filter on score if appropriate
    if ( $intron_files->{DEPTH} ) {
      if ( $intron_files->{DEPTH} > $id_list{$key} ) {
        //#print "Rejecting on score " . $id_list{$key} ."\n";
        next;
      }
    }
    my @data = split(/:/,$key) ;
    my $length =  $data[2] - $data[1] -1;
    next unless $length > 0 ;
    my $name = $data[0]. ":" . ($data[1]+1) . ":" . ($data[2] -1) .":" . $data[3] .":";

    my $if = Bio::EnsEMBL::DnaDnaAlignFeature->new
      (
       -start => $data[1],
       -end => $data[2],
       -strand => $data[3],
       -hstart => 1,
       -hend => $length,
       -hstrand => 1,
       -slice => $self->chr_slice,
       -analysis => $intron_anal,#$self->analysis,
       -score =>  $id_list{$key},
       -hseqname => "$name",
       -cigar_string => $length ."M",
      );
    my $canonical = 1;
    // figure out if its cannonical or not
    my $left_splice = $slice_adaptor->fetch_by_region('toplevel',
                                                      $if->seq_region_name,
                                                      $if->start+1,
                                                      $if->start+2,
                                                      $if->strand
                                                     );
    my $right_splice = $slice_adaptor->fetch_by_region('toplevel',
                                                       $if->seq_region_name,
                                                       $if->end-2,
                                                       $if->end-1,
                                                       $if->strand
                                                      );
    //#  print "KEY $key " . $if->score ."\n";;
    //#  print "LEFT  " . $left_splice->start ." " . $left_splice->end  ." " . $left_splice->strand ." " . $left_splice->seq . "\n";
    //#  print "RIGHT " . $right_splice->start ." " . $right_splice->end  ." " . $right_splice->strand  ." " . $right_splice->seq ."\n\n";
    
    
    if ( $left_splice->seq eq 'NN' && $right_splice->seq eq 'NN' ) {
      warn("Cannot find dna sequence for " . $key .
           " this is used in detecting non cannonical splices\n");
    } else {
      //# is it cannonical
      if ( $if->strand  == 1 ) {
        #        print "Splice type " . $left_splice->seq ."-".  $right_splice->seq ." ";
        # is it GTAG?
        unless ( $left_splice->seq eq 'GT' && $right_splice->seq eq 'AG' ) {
          $canonical = 0;
        }
      } else {
        //#        print "Splice type " . $right_splice->seq ."-".  $left_splice->seq ." ";
        //# is it GTAG?
        unless ( $right_splice->seq eq 'GT' && $left_splice->seq eq 'AG' ) {
          $canonical = 0;
        }
      }
    }
    if ( $canonical ) {
       $if->hseqname($if->hseqname."canonical");
    } else {
      $if->hseqname($if->hseqname."non canonical");
    }
    //#print "INTRONS ".  $if->hseqname ."\n";
    push @ifs , $if;
  }
  // sort them
  @ifs = sort {$a->start <=> $b->start} @ifs;
  if ($self->FILTER_ON_OVERLAP) {
      my @tmp_array;
      my $threshold = $self->FILTER_ON_OVERLAP;
      my $array_length = scalar(@ifs); 
      if ($array_length > 1) {
          for (my $j = 0; $j < $array_length-1; $j++) {
              my $k = 0;
              my $count = 1;
              my $overlapped_support = 0;
              while () {
                  ++$k;
                  if ($count > $threshold) {
                      if ($overlapped_support < $ifs[$j]->score) {
//#                          print STDERR "\t",$ifs[$j+$k]->hseqname, ': ', $ifs[$j+$k]->start, ':', $ifs[$j+$k]->end, "\n";
                          push (@tmp_array, $ifs[$j]);
                      }
                      else {
//#                          print STDERR 'THROWING: ', $ifs[$j]->hseqname, ': ', $ifs[$j]->start, ':', $ifs[$j]->end, "\n";
                      }
                      last;
                  }
                  $overlapped_support += $ifs[$j+$k]->score;
                  if (($ifs[$j]->end < $ifs[$j+$k]->start) or (($j+$k) == $array_length-1)) {
//#                      print STDERR "\t",$ifs[$j+$k]->hseqname, ': ', $ifs[$j+$k]->start, ':', $ifs[$j+$k]->end, "\n";
                      push (@tmp_array, $ifs[$j]);
                      last;
                  }
//#                      print STDERR $ifs[$j+$k]->hseqname, "\n";
                  next unless ($ifs[$j]->strand == $ifs[$j+$k]->strand);
                  ++$count;
              }
          }
          @ifs = @tmp_array;
      }
  }
//#  print STDERR 'RES: ', scalar(@ifs), "\n";
//# SMJS Filter here 

  if (!defined($conslim) || !defined($nonconslim)) {
    die "Env vars for CONSLIM and NONCONSLIM not set\n";
  }

  print "Filter parameters:  Consensus splice coverage $conslim    Non consensus splice coverage $nonconslim\n";

  my @tmp;
  foreach my $f (@ifs) {
    if($f->hseqname =~ /non canonical/) {
      if ($f->score > $nonconslim && $f->end-$f->start < 50000) {
        push @tmp,$f;
      } else {
        print STDERR "Rejected non_canonical feature with score " . $f->score . "\n";
      }
    } else {
      if ($f->score > $conslim && $f->end-$f->start < 150000) {
        push @tmp,$f;
      } else {
        print STDERR "Rejected feature with score " . $f->score . "\n";
      }
    }
  }
  @ifs = @tmp;

  $self->intron_features(\@ifs);
  $self->extra_exons($extra_exons);
  print STDERR "Got " . scalar(@ifs)  . " unique introns  " ;
  print STDERR " and " . scalar(keys %$extra_exons) . " potential novel exons from " . $intron_files->{FILE} . "\n";
  return;
}

sub ungapped_features {
  my ($self,$read) = @_;
  my @ugfs;
  my @tmp_ugfs;
  my $string = $read->cigar_str;
  my $start = $read->start;
  my $end = $read->end;

  //#rint "THINGS $start $end $string\n";
  my @pieces = ( $string =~ /(\d*[MDN])/g );

  for my $piece ( @pieces ) {
    my ($length) = ( $piece =~ /^(\d*)/ );
    if( $length eq "" ) { $length = 1 }
    if( $piece =~ /M$/ ) {
      //
      // MATCH
      //
      my ( $qstart, $qend);
      $qstart = $start;
      $qend = $start + $length - 1;
      $start = $qend + 1;
      
      my $ugf;
      $ugf->[0] = $read->query->name;
      $ugf->[1] = $read->seq_id;
      $ugf->[2] = $qstart;
      $ugf->[3] = $qend;
      $ugf->[4] = $length."M";
      push @tmp_ugfs, $ugf;
      //#print "UNGAPPED " .$ugf->[2] .
        //#" " . $ugf->[3] . " " . $ugf->[4] ."\n";
    } elsif( $piece =~ /N$/ ) {
      //
      // INSERT
      //
      $start += $length;
      push @tmp_ugfs,"intron";
    } elsif( $piece =~ /D$/ ) {
      //
      // DELETION
      //
      $start += $length;
      push @tmp_ugfs,"deletion";
    } else {
      throw( "Illegal cigar line $string!" );
    }
  }
  // only return the UGFS either side of splices
  my %used_pieces;
  foreach ( my $i = 0 ; $i < scalar(@pieces); $i++ )  {
    my $piece = $pieces[$i];
    if ( $piece =~ /\d*N/) {
      // it's a splice push the Matches either side of it
      for ( my $j = $i-1 ; $j >= 0 ; $j-- ) {
        if ( $tmp_ugfs[$j] && $pieces[$j] =~ /\d*M/ )  {
          unless ( $used_pieces{$j} ) {
            my $ugf =  $tmp_ugfs[$j];
            $self->throw("Cannot find ugf $j\n") unless $ugf;
            push @ugfs, $ugf;
            $used_pieces{$j} =1;
            last ;
          }
        }
      }
      for ( my $j = $i+1 ; $j < scalar(@pieces)  ; $j++ ) {
        if ( $tmp_ugfs[$j] && $pieces[$j] =~ /\d*M/ )  {
          unless ( $used_pieces{$j} ) {
            my $ugf =  $tmp_ugfs[$j];
            $self->throw("Cannot find ugf $j\n") unless $ugf;
            push @ugfs, $ugf;
            $used_pieces{$j} =1;
            last ;
          }
        }
      }
    }
  }
  return \@ugfs;
}

/*
=head2 dna_2_extra_exons
    Title        :   dna_2_extra_exons
    Usage        :   $self->dna_2_extra_exons($start,$end)
    Returns      :   None
    Args         :   Int start
                 :   Int end
    Description  :   Fetches all dna_align_features from the small intron db that lie within
                 :   the range determined by start and end, collapses them down into a 
                 :   non redundant set and builds extra exons to account for short exons
                 :   missed by genomic alignemnt of reads

=cut
*/

Vector *RefineSolexaGenes_dnaToExtraExons(RefineSolexaGenes *rsg, long start, long end) {

  Vector *extraIntrons = Vector_getIntronFeatures(rsg);

  SliceAdaptor *smallIntronSliceAdaptor = RefineSolexaGenes_getSmallIntronSliceAdaptor(rsg);
  Slice *smallIntronSlice = SliceAdaptor_fetchByRegion(smallIntronSliceAdaptor, "toplevel", Slice_getSeqRegionName(chrSlice), start, end); 

  my @ifs;
  push @ifs, @$extra_introns if $extra_introns;
  my $rough_genes = $self->prelim_genes;
  my %exon_list;
  my %id_list;
  
  // fetch all the dna_align_features for this slice by logic name
  my @extra_reads;
  // look for extra introns from ESLA

  // featch all the dna_align_features for this slice by logic name
  Vector *logicNames = RefineSolexaGenes_getLogicNames(rsg);
  Vector *reads;

  if (Vector_getNumElement(logicNames)) {
    reads = Vector_new();
    fprintf(stderr,"Fetching reads with logic names: ");

    int i;
    for (i; i<Vector_getNumElement(logicNames); i++) {
      char *logicName = Vector_getElementAt(logicNames, i);

      fprintf(stderr,"%s ", logicName);
      Vector *feats = Slice_getAllDNAAlignFeatures(smallIntronSlice, logicName);
      Vector_append(reads, feats); 
      Vector_free(feats);
    }
    fprintf(stderr,"\n");
  } else {
  //# fetch them all if no logic name is Supplied
    reads =  Slice_getAllDNAAlignFeatures(smallIntronSlice);
  }
  fprintf(stderr, "Got %d extra reads\n", Vector_getNumElement(reads));

  // process extra reads and assign them to rough models
  while ( scalar @extra_reads > 0 ) {
    my $read = pop(@extra_reads);
    my $type = 'canonical';
    my $roughid;
    $type = 'non canonical' if ( $read->hseqname =~ /\:NC$/ ) ;
    $read = $read->transfer($self->chr_slice);
    // assign a rough model to the reads
    my @ugfs = $read->ungapped_features;
    @ugfs = sort { $a->start <=> $b->start } @ugfs;
    next if ( scalar(@ugfs) == 0 ) ;
    // first and or last ugf should overlap an exon
  ROUGH:  foreach my $rough ( @$rough_genes ) {
      foreach my $exon ( @{$rough->get_all_Exons} ) {
        if( ( $exon->start <= $ugfs[0]->end &&  $exon->end >= $ugfs[0]->start ) or 
            ( $exon->start <= $ugfs[-1]->end &&  $exon->end >= $ugfs[-1]->start ) ) {
          //ONLY USE THIS READ WITH THIS MODEL
          $roughid = $rough->stable_id;
          last ROUGH;
        }
      }
    }
    unless ( $roughid ) {
      print STDERR "Ignoring read " . $read->hseqname . " cannot pair it with a rough model\n";
    }
    for ( my $i = 0 ; $i < scalar(@ugfs) - 1 ; $i++ ) {
      //# one read can span several exons so make all the features 
      // cache them by internal boundaries
      // we use . to deliminate entries in our paths so dont allow them in the seq_region_name or it wont work
      my $name = $read->seq_region_name;
// NOT REALLY COMMENTED OUT, JUST TO STOP SYNTAX HIGHLIGHTING GOING NUTS      $name =~ s/\./*/g;
      my $unique_id = $name . ":" . 
        $ugfs[$i]->end . ":" .
          $ugfs[$i+1]->start . ":" . 
            $read->strand .":$type";
      $id_list{$unique_id} ++;
      if  ( $i > 0 && $i < scalar(@ugfs) - 1 ) {
        // if the read splices into and out of an exon, store that exon in case it is not 
        // found in the rough model but is an extra short exon determined by ExonerateSolexaLocalAlignment
        // check that it is introns on both sides and not inserts, need to be bigger than the min intron length
        next unless $ugfs[$i+1]->start - $ugfs[$i]->end > $self->MIN_INTRON_SIZE;
        next unless $ugfs[$i]->start - $ugfs[$i-1]->end > $self->MIN_INTRON_SIZE;
        my $name = $read->seq_region_name;
// NOT REALLY COMMENTED OUT, JUST TO STOP SYNTAX HIGHLIGHTING GOING NUTS      $name =~ s/\./*/g;
        my $unique_id = $name . ":" . 
          $ugfs[$i]->start . ":" .
            $ugfs[$i]->end . ":" . 
              $read->strand .":$type";
        $exon_list{$roughid}{$unique_id} = $ugfs[$i];        
      }
    }
  }
  
  //collapse down the intron feaures and make them into simple features
  foreach my $key ( keys %id_list ) {
   // print "KEY $key\n";
    my @data = split(/:/,$key) ;
    my $length =  $data[2] - $data[1] -1;
    next unless $length > 0 ;
    my $name = $self->chr_slice->seq_region_name . ":" . ($data[1]+1) . ":" . ($data[2] -1) .":" . $data[3] .":".$data[4];
    my $if = Bio::EnsEMBL::DnaDnaAlignFeature->new
      (
       -start => $data[1],
       -end => $data[2],
       -strand => $data[3],
       -hstart => 1,
       -hend => $length,
       -hstrand => 1,
       -slice => $self->chr_slice,
       -analysis => $self->analysis,
       -score =>  $id_list{$key},
       -hseqname => $name,
       -cigar_string => $length ."M",
      );
    push @ifs , $if;
  }

  // sort them
  @ifs = sort {$a->start <=> $b->start} @ifs;
  $self->intron_features(\@ifs);
  print STDERR "Got " . scalar(@ifs) . " intron features\n";
  // make the exon features
  my @extra_exons;
  foreach my $rough ( keys %exon_list ) {
    foreach my $key ( keys %{$exon_list{$rough}} ) {
      my $extra_exon = $self->make_exon($exon_list{$rough}{$key});
      $extra_exon->{"_extra"} = 1;
      $extra_exon->{"_model"} = $rough;
      push @extra_exons, $extra_exon;
    }
  }
  print STDERR "Got " . scalar(@extra_exons) . " extra exons\n";
  // make the exon features
  $self->extra_exons(\@extra_exons);
  return;
}


/*
=head2 dna_2_intron_features
    Title        :   dna_2_intron_features
    Usage        :   $self->dna_2_intron_features($start,$end)
    Returns      :   None
    Args         :   Int start
                 :   Int end
    Description  :   Fetches all dna_align_features from the intron db that lie within
                 :   the range determined by start and end, collapses them down into a 
                 :   non redundant set and builds a Bio::EnsEMBL::DnaAlignFeature to 
                 :   represent it, then stores it in $self->intron_features
                 :   also checks for small exons defined by a single read splicing at least twice
                 :   stores any additional exons found this way in $self->extra_exons

=cut
*/
void RefineSolexaGenes_dnaToIntronFeatures(RefineSolexaGenes *rsg, long start, long end) {
// Unused  my $rough_genes = $self->prelim_genes;
  SliceAdaptor *intronSliceAdaptor = RefineSolexaGenes_getIntronSliceAdaptor(rsg);
  Slice *chrSlice = RefineSolexaGenes_getChrSlice(rsg);

  Slice *intronSlice = SliceAdaptor_fetchByRegion(intronSliceAdaptor, "toplevel", Slice_getSeqRegionName(chrSlice), start, end); 

  // featch all the dna_align_features for this slice by logic name
  Vector *logicNames = RefineSolexaGenes_getLogicNames(rsg);
  Vector *reads;

  if (Vector_getNumElement(logicNames)) {
    reads = Vector_new();
    fprintf(stderr,"Fetching reads with logic names: ");

    int i;
    for (i; i<Vector_getNumElement(logicNames); i++) {
      char *logicName = Vector_getElementAt(logicNames, i);

      fprintf(stderr,"%s ", logicName);
      Vector *feats = Slice_getAllDNAAlignFeatures(intronSlice, logicName);
      Vector_append(reads, feats); 
      Vector_free(feats);
    }
    fprintf(stderr,"\n");
  } else {
  //# fetch them all if no logic name is Supplied
    reads =  Slice_getAllDNAAlignFeatures(intronSlice);
  }
  fprintf(stderr, "Got %d reads\n", Vector_getNumElement(reads));

  StringHash *idList = StringHash_new(STRINGHASH_MEDIUM);

  int i;
  for (i=0; i<Vector_getNumElement(reads); i++) {
    DNAAlignFeature *read = Vector_getElementAt(reads, i);

    char *type = "canonical";
 // Was  $type = 'non canonical' if ( $read->hseqname =~ /\:NC$/ ) ;
    if (!strstr(DNAAlignFeature_getHitSeqName(read), ":NC")) {
      type = "non canonical";
    }
    
    // Do we actually need to tranfer this???
    DNAAlignFeature *origRead = read;
    read = DNAAlignFeature_transfer(read, chrSlice);
// NIY: Need to free

    Vector *ugfs = DNAAlignFeature_getUngappedFeatures(read);

    if (Vector_getNumElement(ugfs) > 0) {
      Vector_sort(ugfs, SeqFeature_posSortFunc);
     
      //# one read can span several exons so make all the features 
      // cache them by internal boundaries
      // we use . to deliminate entries in our paths so dont allow them in the seq_region_name or it wont work
      char name[1024];
      strcpy(name, DNAAlignFeature_getSeqRegionName(read));
      StrUtil_strReplChr(name, '.', '*');
      
      char uniqueId[2048];

      for (j=0 ; j<Vector_getNumElement(ugfs)-1; j++) {
        DNAAlignFeature *ugf   = Vector_getElementAt(ugfs, j);
        DNAAlignFeature *ugfP1 = Vector_getElementAt(ugfs, j+1);
        
        sprintf(uniqueId, "%s:%ld:%ld:%d:%s", name, DNAAlignFrature_getEnd(ugf), 
                DNAAlignFeature_getStart(ugfP1), DNAAlignFeature_getStrand(read), type)

        if (!StringHash_contains(idList, uniqueId)) {
          StringHash_add(idList, uniqueId, long_new());
        }
        *(StringHash_getValue(idList, uniqueId)) ++;
      }
    }
  }

  fprintf(stderr, "Got %d collapsed introns\n", StringHash_getNumValues(idList));

  //# collapse them down and make them into simple features

  Vector *intFeats = Vector_new();
  Analysis *analysis = RefineSolexaGenes_getAnalysis(rsg);

  foreach my $key ( keys %id_list ) {
    my @data = split(/:/,$key) ;
    my $length =  $data[2] - $data[1] -1;

    if (length > 0) {
      my $name = $self->chr_slice->seq_region_name . ":" . ($data[1]+1) . ":" . ($data[2] -1) .":" . $data[3] .":".$data[4];

      DNAAlignFeature *intFeat = DNAAlignFeature_new();

      DNAAlignFeature_setStart      (intFeat, );
      DNAAlignFeature_setEnd        (intFeat, );
      DNAAlignFeature_setStrand     (intFeat, );
      DNAAlignFeature_setHitStart   (intFeat, 1);
      DNAAlignFeature_setHitEnd     (intFeat, length);
      DNAAlignFeature_setHitStrand  (intFeat, 1);
      DNAAlignFeature_setSlice      (intFeat, chrSlice);
      DNAAlignFeature_setAnalysis   (intFeat, analysis);
      DNAAlignFeature_setScore      (intFeat, );
      DNAAlignFeature_setHitSeqName (intFeat, name);
      char cigStr[256]
      sprintf(cigStr, "%dM", length);
      DNAAlignFeature_setCigarString(intFeat, cigStr);

      my $if = Bio::EnsEMBL::DnaDnaAlignFeature->new
        (
         -start => $data[1],
         -end => $data[2],
         -strand => $data[3],
         -hstart => 1,
         -hend => $length,
         -hstrand => 1,
         -slice => $self->chr_slice,
         -analysis => $self->analysis,
         -score =>  $id_list{$key},
         -hseqname => $name,
         -cigar_string => $length ."M",
        );

      Vector_addElement(intFeats, intFeat);
    }
  }

  // sort them
  Vector_sort(intFeats, SeqFeature_posSortFunc);

  RefineSolexaGenes_setIntronFeatures(intFeats);

  fprintf(stderr, "Got %d intron features\n", Vector_getNumElement(intFeats));
  return;
}

/*
=head2 fetch_intron_features

    Title        :   fetch_intron_features
    Usage        :   $self->fetch_intron_features($start,$end)
    Returns      :   Array ref of Bio::EnsEMBL::DnaAlignFeature
    Args         :   Int start
                 :   Int end
    Description  :   Accesses the pre computed simple features representing introns
                 :   Filters out non consensus models that overlap consensus models

=cut
*/

Vector *RefineSolexaGenes_fetchIntronFeatures(RefineSolexaGenes *rsg, long start, long end, long offset) {
  Vector *sfs = RefineSolexaGenes_getIntronFeatures(rsg);

  long intronStart = 0;
  if (offset) intronStart = offset;
  int index = -1;

  Vector *chosenSf = Vector_new();

  // sfs is a sorted array
  int i;
  for (i=intronStart; i<Vector_getNumElement(sfs); i++) {
    DNAAlignFeature *intron = Vector_getElementAt(sfs, i);

    if (DNAAlignFeature_getStart(intron) > end) break;

    if (DNAAlignFeature_getStart(intron) <= end && DNAAlignFeature_getEnd(intron) >= start) {
      Vector_addElement(chosenSf, intron);

      // remember the position of the 1st intron to overlap
      // this exon - we will start counting from here next time
      if (index == -1) {
        index = i;
      }
    }
  }

  Vector *filteredIntrons = Vector_new();
// INTRON: 
  for (i=0; i<Vector_getNumElement(chosenSf); i++) {
    DNAAlignFeature *intron = Vector_getElementAt(chosenSf, i);

    if (strstr(DNAAlignFeature_getHitSeqName(intron), "non canonical")) {
      // check it has no overlap with any consensus introns
      // unless it out scores a consensus intron
      int j;
      for (j=0; j<Vector_getNumElement(chosenSf); j++) {
        DNAAlignFeature *compIntron = Vector_getElementAt(chosenSf, j);
        if (strstr(DNAAlignFeature_getHitSeqName(compIntron), "non canonical") == NULL) {
          if (DNAAlignFeature_getEnd(intron)    >  DNAAlignFeature_getStart(compIntron) && 
              DNAAlignFeature_getStart(intron)  <  DNAAlignFeature_getEnd(compIntron) &&
              DNAAlignFeature_getStrand(intron) == DNAAlignFeature_getStrand(compIntron)) {
            next INTRON if $intron->score <= $i->score;
          }
        }
      }
      Vector_addElement(filteredIntrons, intron);
    } else {
      Vector_addElement(filteredIntrons, intron);
    }
  }

// NIY Need to return index
  return (\@filtered_introns,$index);
}


/*

=head2 make_exon
    Title        :   pad_exons
    Usage        :   $self->($ungapped_feature)
    Returns      :   Bio::EnsEMBL::Exon
    Args         :   Bio::EnsEMBL::FeaturePair 
    Description  :   Takes an ungapped feature, pads it and builds a 
                 :   Exon from it 
=cut
*/

Exon *RefineSolexaGenes_makeExon(RefineSolexaGenes *rsg, UngappedFeature *ugf, long start, long end, double score, char *diplayId) {
  Slice *chrSlice = RefineSolexaGenes_getChrSlice(rsg);

  if (ugf) {
    start     = ugf->start;
    end       = ugf->end;
    displayId = ugf->displayId;
    score     = ugf->score;
  }

  Exon *paddedExon = creatExon(start-20, end+20, -1, -1, -1, RefineSolexaGenes_getAnalysis(rsg), NULL, NULL, RefineSolexaGenes_getChrSlice(rsg));

  // dont let it fall of the slice because of padding
  if (Exon_getStart(paddedExon) <= 0) {
    Exon_setStart(paddedExon, 1);
  }

  if (Exon_getEnd(paddedExon) >= Slice_getLength(chrSlice)) {
    Exon_setEnd(paddedExon, Slice_getLength(chrSlice)-1);
  }
  
  DNAAlignFeature *feat = DNAAlignFeature_new();
  DNAAlignFeature_setSlice     (feat, RefineSolexaGenes_getChrSlice(rsg));
  DNAAlignFeature_setStart     (feat, Exon_getStart(paddedExon));
  DNAAlignFeature_setEnd       (feat, Exon_getEnd(paddedExon));
  DNAAlignFeature_setStrand    (feat, -1); // Why -1??????????????
  DNAAlignFeature_setHitSeqName(feat, displayId);
  DNAAlignFeature_setHitStart  (feat, 1);
  DNAAlignFeature_setHitStrand (feat, 1);
  DNAAlignFeature_setHitEnd    (feat, Exon_getLength(paddedExon));
  DNAAlignFeature_setAnalysis  (feat, RefineSolexaGenes_getAnalysis(rsg));
  DNAAlignFeature_setScore     (feat, score);

  char tmpStr[256];
  sprintf(tmpStr, "%ldM", Exon_getLength(paddedExon));

  my @feats;
  push @feats,$feat;
  $padded_exon->add_supporting_features(@feats);

  return paddedExon;
}

 
//##################################################################
//# Containers

RefineSolexaGenes_setRecursiveLimit(RefineSolexaGenes *rsg, int limit) {
  rsg->recursiveLimit = limit;
}

RefineSolexaGenes_getRecursiveLimit(RefineSolexaGenes *rsg) {
  return rsg->recursiveLimit;
}

/* Not used
sub repeat_feature_adaptor {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_rfa} = $val;
  }

  return $self->{_rfa};
}
*/

sub gene_slice_adaptor {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_gsa} = $val;
  }

  return $self->{_gsa};
}

sub intron_slice_adaptor {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_isa} = $val;
  }

  return $self->{_isa};
}

sub small_intron_slice_adaptor {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_sisa} = $val;
  }

  return $self->{_sisa};
}

/* Not used
sub feature_cash {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_feature_cash} = $val;
  }

  return $self->{_feature_cash};
}
*/

sub prelim_genes {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_prelim_genes} = $val;
  }

  return $self->{_prelim_genes};
}


sub chr_slice {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_chr_slice} = $val;
  }

  return $self->{_chr_slice};
}

sub intron_features {
  my ($self, $val) = @_;

  if (defined $val) {
    push @{$self->{_introns}}, @$val;
    # make sure it is still sorted
    @{$self->{_introns}} = sort { $a->start <=> $b->start } @{$self->{_introns}};
  }
  return $self->{_introns};
}


sub extra_exons {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_extra_exons} = $val;
  }

  return $self->{_extra_exons};
}

//####################################
//# config variable holders
//####################################

sub INTRON_DB {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_INTRON_DB'} = $value;
  }
  
  if (exists($self->{'_CONFIG_INTRON_DB'})) {
    return $self->{'_CONFIG_INTRON_DB'};
  } else {
    return undef;
  }
}


sub OUTPUT_DB {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_OUTPUT_DB'} = $value;
  }
  
  if (exists($self->{'_CONFIG_OUTPUT_DB'})) {
    return $self->{'_CONFIG_OUTPUT_DB'};
  } else {
    return undef;
  }
}


sub MODEL_DB {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_MODEL_DB'} = $value;
  }
  
  if (exists($self->{'_CONFIG_MODEL_DB'})) {
    return $self->{'_CONFIG_MODEL_DB'};
  } else {
    return undef;
  }
}

sub LOGICNAME {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_LOGICNAME'} = $value;
  }
  
  if (exists($self->{'_CONFIG_LOGICNAME'})) {
    return $self->{'_CONFIG_LOGICNAME'};
  } else {
    return undef;
  }
}

sub RETAINED_INTRON_PENALTY {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_RETAINED_INTRON_PENALTY'} = $value;
  }
  
  if (exists($self->{'_CONFIG_RETAINED_INTRON_PENALTY'})) {
    return $self->{'_CONFIG_RETAINED_INTRON_PENALTY'};
  } else {
    return undef;
  }
}


sub MIN_INTRON_SIZE {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_MIN_INTRON_SIZE'} = $value;
  }
  
  if (exists($self->{'_CONFIG_MIN_INTRON_SIZE'})) {
    return $self->{'_CONFIG_MIN_INTRON_SIZE'};
  } else {
    return undef;
  }
}


sub MAX_INTRON_SIZE {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_MAX_INTRON_SIZE'} = $value;
  }
  
  if (exists($self->{'_CONFIG_MAX_INTRON_SIZE'})) {
    return $self->{'_CONFIG_MAX_INTRON_SIZE'};
  } else {
    return undef;
  }
}


sub BEST_SCORE {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_BEST_SCORE'} = $value;
  }
  
  if (exists($self->{'_CONFIG_BEST_SCORE'})) {
#    return $self->{'_CONFIG_BEST_SCORE'};
    my $conslim = $ENV{CONSLIM};
    my $nonconslim = $ENV{NONCONSLIM};
    return $self->{'_CONFIG_BEST_SCORE'} . "_c" . $conslim . "_nc" . $nonconslim;
  } else {
    return undef;
  }
}


sub OTHER_NUM {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_OTHER_NUM'} = $value;
  }
  
  if (exists($self->{'_CONFIG_OTHER_NUM'})) {
    return $self->{'_CONFIG_OTHER_NUM'};
  } else {
    return undef;
  }
}

sub OTHER_ISOFORMS {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_OTHER_ISOFORMS'} = $value;
  }
  
  if (exists($self->{'_CONFIG_OTHER_ISOFORMS'})) {
    return $self->{'_CONFIG_OTHER_ISOFORMS'};
  } else {
    return undef;
  }
}

sub MODEL_LN {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_MODEL_LN'} = $value;
  }
  
  if (exists($self->{'_CONFIG_MODEL_LN'})) {
    return $self->{'_CONFIG_MODEL_LN'};
  } else {
    return undef;
  }
}

sub BAD_MODELS {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_BAD_MODELS'} = $value;
  }
  
  if (exists($self->{'_CONFIG_BAD_MODELS'})) {
    return $self->{'_CONFIG_BAD_MODELS'};
  } else {
    return undef;
  }
}

sub MAX_NUM {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_MAX_NUM'} = $value;
  }
  
  if (exists($self->{'_CONFIG_MAX_NUM'})) {
    return $self->{'_CONFIG_MAX_NUM'};
  } else {
    return undef;
  }
}

sub MAX_RECURSIONS {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_MAX_RECURSIONS'} = $value;
  }
  
  if (exists($self->{'_CONFIG_MAX_RECURSIONS'})) {
    return $self->{'_CONFIG_MAX_RECURSIONS'};
  } else {
    return undef;
  }
}

sub MIN_SINGLE_EXON {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_MIN_SINGLE_EXON'} = $value;
  }
  
  if (exists($self->{'_CONFIG_MIN_SINGLE_EXON'})) {
    return $self->{'_CONFIG_MIN_SINGLE_EXON'};
  } else {
    return undef;
  }
}

sub SINGLE_EXON_CDS {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_SINGLE_EXON_CDS'} = $value;
  }
  
  if (exists($self->{'_CONFIG_SINGLE_EXON_CDS'})) {
    return $self->{'_CONFIG_SINGLE_EXON_CDS'};
  } else {
    return undef;
  }
}

sub SINGLE_EXON_MODEL {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_SINGLE_EXON_MODEL'} = $value;
  }
  
  if (exists($self->{'_CONFIG_SINGLE_EXON_MODEL'})) {
    return $self->{'_CONFIG_SINGLE_EXON_MODEL'};
  } else {
    return undef;
  }
}

sub STRICT_INTERNAL_SPLICE_SITES{
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_STRICT_INTERNAL_SPLICE_SITES'} = $value;
  }
  
  if (exists($self->{'_CONFIG_STRICT_INTERNAL_SPLICE_SITES'})) {
    return $self->{'_CONFIG_STRICT_INTERNAL_SPLICE_SITES'};
  } else {
    return undef;
  }
}

sub STRICT_INTERNAL_END_EXON_SPLICE_SITES {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_STRICT_INTERNAL_END_EXON_SPLICE_SITES'} = $value;
  }
  
  if (exists($self->{'_CONFIG_STRICT_INTERNAL_END_EXON_SPLICE_SITES'})) {
    return $self->{'_CONFIG_STRICT_INTERNAL_END_EXON_SPLICE_SITES'};
  } else {
    return undef;
  }
}

sub INTRON_BAM_FILES {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_INTRON_BAM_FILE'} = $value;
  }
  
  if (exists($self->{'_CONFIG_INTRON_BAM_FILE'})) {
    return $self->{'_CONFIG_INTRON_BAM_FILE'};
  } else {
    return undef;
  }
}


sub WRITE_INTRONS {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_WRITE_INTRONS'} = $value;
  }
  
  if (exists($self->{'_WRITE_INTRONS'})) {
    return $self->{'_WRITE_INTRONS'};
  } else {
    return undef;
  }
}

sub TRIM_UTR {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_TRIM_UTR'} = $value;
  }
  
  if (exists($self->{'_CONFIG_TRIM_UTR'})) {
    return $self->{'_CONFIG_TRIM_UTR'};
  } else {
    return undef;
  }
}


sub MAX_3PRIME_EXONS {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_MAX_3PRIME_EXONS'} = $value;
  }
  
  if (exists($self->{'_CONFIG_MAX_3PRIME_EXONS'})) {
    return $self->{'_CONFIG_MAX_3PRIME_EXONS'};
  } else {
    return undef;
  }
}


sub MAX_3PRIME_LENGTH {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_MAX_3PRIME_LENGTH'} = $value;
  }
  
  if (exists($self->{'_CONFIG_MAX_3PRIME_LENGTH'})) {
    return $self->{'_CONFIG_MAX_3PRIME_LENGTH'};
  } else {
    return undef;
  }
}


sub MAX_5PRIME_EXONS {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_MAX_5PRIME_EXONS'} = $value;
  }
  
  if (exists($self->{'_CONFIG_MAX_5PRIME_EXONS'})) {
    return $self->{'_CONFIG_MAX_5PRIME_EXONS'};
  } else {
    return undef;
  }
}


sub FILTER_ON_OVERLAP {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_FILTER_ON_OVERLAP'} = $value;
  }
  
  if (exists($self->{'_CONFIG_FILTER_ON_OVERLAP'})) {
    return $self->{'_CONFIG_FILTER_ON_OVERLAP'};
  } else {
    return undef;
  }
}

sub MAX_5PRIME_LENGTH {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_MAX_5PRIME_LENGTH'} = $value;
  }
  
  if (exists($self->{'_CONFIG_MAX_5PRIME_LENGTH'})) {
    return $self->{'_CONFIG_MAX_5PRIME_LENGTH'};
  } else {
    return undef;
  }
}


sub REJECT_INTRON_CUTOFF {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_REJECT_INTRON_CUTOFF'} = $value;
  }
  
  if (exists($self->{'_CONFIG_REJECT_INTRON_CUTOFF'})) {
    return $self->{'_CONFIG_REJECT_INTRON_CUTOFF'};
  } else {
    return undef;
  }
}


