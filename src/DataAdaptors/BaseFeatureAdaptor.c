#include "BaseFeatureAdaptor.h"

#include "AnalysisAdaptor.h"
#include "AssemblyMapperAdaptor.h"
#include "StrUtil.h"
#include "Error.h"
#include "SliceAdaptor.h"
#include "ProjectionSegment.h"
#include "MetaContainer.h"
#include "MetaCoordContainer.h"

// For testing
//#include "DNAAlignFeature.h"


/*
=head1 NAME

Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor - An Abstract Base class for all
FeatureAdaptors

=head1 SYNOPSIS

Abstract class - should not be instantiated.  Implementation of
abstract methods must be performed by subclasses.

=head1 DESCRIPTION

This is a base adaptor for feature adaptors. This base class is simply a way
of eliminating code duplication through the implementation of methods
common to all feature adaptors.
*/


static int SLICE_FEATURE_CACHE_SIZE    = 4;
static int MAX_SPLIT_QUERY_SEQ_REGIONS = 3;
static int SILENCE_CACHE_WARNINGS      = 0;


/*
=head2 new

  Arg [1]    : list of args @args
               Superclass constructor arguments
  Example    : none
  Description: Constructor which warns if caching has been switched off
  Returntype : Bio::EnsEMBL::BaseFeatureAdaptor
  Exceptions : none
  Caller     : implementing subclass constructors
  Status     : Stable

=cut
*/

void BaseFeatureAdaptor_init(BaseFeatureAdaptor *bfa, DBAdaptor *dba, int adaptorType) {
  BaseAdaptor_init((BaseAdaptor *)bfa,dba,adaptorType);

  bfa->sliceFeatureCache = Cache_new(SLICE_FEATURE_CACHE_SIZE);

  if ( DBAdaptor_noCache(bfa->dba) && ! SILENCE_CACHE_WARNINGS) {
    fprintf(stderr, "You are using the API without caching most recent features. "
                     "Performance might be affected.\n");
  }
/*
  bfa->objectsFromStatementHandle = BaseFeatureAdaptor_objectsFromStatementHandle;
  bfa->getTables                  = BaseFeatureAdaptor_getTables;
  bfa->getColumns                 = BaseFeatureAdaptor_getColumns;
  bfa->finalClause                = BaseFeatureAdaptor_finalClause;
  bfa->leftJoin                   = BaseFeatureAdaptor_leftJoin;
  bfa->defaultWhereClause         = BaseFeatureAdaptor_defaultWhereClause;
  bfa->store                      = BaseFeatureAdaptor_store;
*/

  return;
}
/*
sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  if ( defined $self->db->no_cache() && $self->db->no_cache() && ! $SILENCE_CACHE_WARNINGS) {
    warning(  "You are using the API without caching most recent features. "
            . "Performance might be affected." );
  }
  return $self;
}
*/

/*
=head2 start_equals_end

  Arg [1]    : (optional) boolean $newval
  Example    : $bfa->start_equals_end(1);
  Description: Getter/Setter for the start_equals_end flag.  If set
               to true sub _slice_fetch will use a simplified sql to retrieve 1bp slices.
  Returntype : boolean
  Exceptions : none
  Caller     : EnsemblGenomes variation DB build
  Status     : Stable
  
=cut
*/

// Do in header
/*
sub start_equals_end {
  my ( $self, $value ) = @_;

  if ( defined($value) ) {
    $self->{'start_equals_end'} = $value;
  }
  return $self->{'start_equals_end'};
}
*/


/*
=head2 clear_cache

  Args      : None
  Example   : my $sa =
                $registry->get_adaptor( 'Mus musculus', 'Core',
                                        'Slice' );
              my $ga =
                $registry->get_adaptor( 'Mus musculus', 'Core',
                                        'Gene' );

              my $slice =
                $sa->fetch_by_region( 'Chromosome', '1', 1e8,
                                      1.05e8 );

              my $genes = $ga->fetch_all_by_Slice($slice);

              $ga->clear_cache();

  Description   : Empties the feature cache associated with this
                  feature adaptor.
  Return type   : None
  Exceptions    : None
  Caller        : General
  Status        : At risk (under development)

=cut
*/

void BaseFeatureAdaptor_clearCache(BaseFeatureAdaptor *bfa) {
  
  BaseFeatureAdaptor_clearSliceFeatureCache(bfa);

// Not implemented ID caching in C because didn't seem to be used much in perl
//  if(!$self->_no_id_cache()) {
//    $self->_id_cache()->clear_cache();
//  }
  return;
}

void BaseFeatureAdaptor_clearSliceFeatureCache(BaseFeatureAdaptor *bfa) {
  Cache_empty(bfa->sliceFeatureCache);

  return;
}

/*
=head2 _slice_feature_cache
 
  Description  : Returns the feature cache if we are allowed to cache and
                will build it if we need to. We will never return a reference
                to the hash to avoid unintentional auto-vivfying caching
  Returntype   : Bio::EnsEMBL::Utils::Cache
  Exceptions   : None
  Caller       : Internal

=cut
*/

/* Created in BFA_init in C
sub _slice_feature_cache {
  my ($self) = @_;
  return if $self->db()->no_cache();
  if(! exists $self->{_slice_feature_cache}) {
    tie my %cache, 'Bio::EnsEMBL::Utils::Cache', $SLICE_FEATURE_CACHE_SIZE;
    $self->{_slice_feature_cache} = \%cache;
  }
  return $self->{_slice_feature_cache};
}
*/

/*
=head2 fetch_all_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               the slice from which to obtain features
  Arg [2]    : (optional) string $logic_name
               the logic name of the type of features to obtain
  Example    : $fts = $a->fetch_all_by_Slice($slice, 'Swall');
  Description: Returns a listref of features created from the database 
               which are on the Slice defined by $slice. If $logic_name is 
               defined only features with an analysis of type $logic_name 
               will be returned. 
               NOTE: only features that are entirely on the slice's seq_region
               will be returned (i.e. if they hang off the start/end of a
               seq_region they will be discarded). Features can extend over the
               slice boundaries though (in cases where you have a slice that
               doesn't span the whole seq_region).
  Returntype : listref of Bio::EnsEMBL::SeqFeatures in Slice coordinates
  Exceptions : none
  Caller     : Bio::EnsEMBL::Slice
  Status     : Stable

=cut
*/

Vector *BaseFeatureAdaptor_fetchAllBySlice(BaseFeatureAdaptor *bfa, Slice *slice,
                                           char *logicName) {
  // fetch by constraint with empty constraint
  return BaseFeatureAdaptor_fetchAllBySliceConstraint(bfa, slice, "", logicName);
}


/*
=head2 fetch_Iterator_by_Slice_method

  Arg [1]    : CODE ref of Slice fetch method
  Arg [2]    : ARRAY ref of parameters for Slice fetch method
  Arg [3]    : Optional int: Slice index in parameters array
  Arg [4]    : Optional int: Slice chunk size. Default=500000
  Example    : my $slice_iter = $feature_adaptor->fetch_Iterator_by_Slice_method
                                       ($feature_adaptor->can('fetch_all_by_Slice_Arrays'),
                                     \@fetch_method_params,
                                     0,#Slice idx
                                    );

               while(my $feature = $slice_iter->next && defined $feature){
                 #Do something here
               }

  Description: Creates an Iterator which chunks the query Slice to facilitate
               large Slice queries which would have previously run out of memory
  Returntype : Bio::EnsEMBL::Utils::Iterator
  Exceptions : Throws if mandatory params not valid
  Caller     : general
  Status     : at risk

=cut

#Does not support Collections. See Funcgen ResultFeatureAdaptor::fetch_collection_Iterator_by_Slice_method

sub fetch_Iterator_by_Slice_method{
  my ($self, $slice_method_ref, $params_ref, $slice_idx, $chunk_size) = @_;

  if(! ( defined $slice_method_ref &&
     ref($slice_method_ref) eq 'CODE')
  ){
  throw('Must pass a valid Slice fetch method CODE ref');
  }

  if (! ($params_ref && 
     ref($params_ref) eq 'ARRAY')) {
  #Don't need to check size here so long as we have valid Slice
  throw('You must pass a method params ARRAYREF');
  }
  
  $slice_idx    = 0 if(! defined $slice_idx);
  my $slice     = $params_ref->[$slice_idx];
  $chunk_size ||= 1000000;
    
  my @feat_cache;
  my $finished     = 0;
  my $start        = 1;  #local coord for sub slice
  my $end          = $slice->length;
  my $num_overlaps = 0;
  
  my $coderef = 
  sub {
    
    while (scalar(@feat_cache) == 0 &&
       ! $finished) {
    
    my $new_end = ($start + $chunk_size - 1);
    
    if ($new_end >= $end) {
      # this is our last chunk
      $new_end = $end;
      $finished = 1;  
    }
    
    #Chunk by sub slicing
    my $sub_slice             = $slice->sub_Slice($start, $new_end);
    $params_ref->[$slice_idx] = $sub_slice;
    @feat_cache = @{ $slice_method_ref->($self, @$params_ref)};
    
    #Remove & count overlapping features
    splice(@feat_cache, 0, $num_overlaps) if($num_overlaps);
    my $i;
    
    if (scalar(@feat_cache) > 0) {
      
      my $feat_end  = $feat_cache[$#feat_cache]->seq_region_end;
      my $slice_end = $sub_slice->end;
      $num_overlaps = 0;
      
      for ($i = $#feat_cache; $i >=0; $i--) {
      
      if ($feat_end > $slice_end) {
        $feat_end  = $feat_cache[$i]->end;
        $num_overlaps ++;
      } else {
        last;
      }
      
      }
    }
    
    # update the start coordinate
    $start = $new_end + 1;
    }
    
    #this maybe returning from an undef cache
    #Need to sub this out even more?
    return shift @feat_cache;
  };

  return Bio::EnsEMBL::Utils::Iterator->new($coderef);
}


=head2 fetch_Iterator_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice
  Arg [2]    : Optional string: logic name of analysis
  Arg [3]    : Optional int: Chunk size to iterate over. Default is 500000
  Example    : my $slice_iter = $feature_adaptor->fetch_Iterator_by_Slice($slice);

               while(my $feature = $slice_iter->next && defined $feature){
                 #Do something here
               }

  Description: Creates an Iterator which chunks the query Slice to facilitate
               large Slice queries which would have previously run out of memory
  Returntype : Bio::EnsEMBL::Utils::Iterator
  Exceptions : None
  Caller     : general
  Status     : at risk

=cut

sub fetch_Iterator_by_Slice{
  my ($self, $slice, $logic_name, $chunk_size) = @_;

  my $method_ref = $self->can('fetch_all_by_Slice');

  return $self->fetch_Iterator_by_Slice_method($method_ref, [$slice, $logic_name], 0, $chunk_size);
}
*/


/*
=head2 fetch_all_by_Slice_and_score

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               the slice from which to obtain features
  Arg [2]    : (optional) float $score
               lower bound of the the score of the features retrieved
  Arg [3]    : (optional) string $logic_name
               the logic name of the type of features to obtain
  Example    : $fts = $a->fetch_all_by_Slice_and_score($slice,90,'Swall');
  Description: Returns a list of features created from the database which are 
               are on the Slice defined by $slice and which have a score 
               greater than $score. If $logic_name is defined, 
               only features with an analysis of type $logic_name will be 
               returned. 
  Returntype : listref of Bio::EnsEMBL::SeqFeatures in Slice coordinates
  Exceptions : none
  Caller     : Bio::EnsEMBL::Slice
  Status     : Stable

=cut
*/

Vector *BaseFeatureAdaptor_fetchAllBySliceAndScore(BaseFeatureAdaptor *bfa, Slice *slice,
                                                   double *scoreP, char *logicName) {
  char constraintStr[256];
  NameTableType *tables = bfa->getTables();

  constraintStr[0] = '\0';
  // Perl does a defined check on score
  if (scoreP) {
    sprintf(constraintStr,"%s.score > %f",(*tables)[0][SYN], *scoreP);
  }

  return BaseFeatureAdaptor_fetchAllBySliceConstraint(bfa, slice, constraintStr, logicName);
}

/*
=head2 fetch_all_by_Slice_constraint

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               the slice from which to obtain features
  Arg [2]    : (optional) string $constraint
               An SQL query constraint (i.e. part of the WHERE clause)
  Arg [3]    : (optional) string $logic_name
               the logic name of the type of features to obtain
  Example    : $fs = $a->fetch_all_by_Slice_constraint($slc, 'perc_ident > 5');
  Description: Returns a listref of features created from the database which 
               are on the Slice defined by $slice and fulfill the SQL 
               constraint defined by $constraint. If logic name is defined, 
               only features with an analysis of type $logic_name will be 
               returned. 
  Returntype : listref of Bio::EnsEMBL::SeqFeatures in Slice coordinates
  Exceptions : thrown if $slice is not defined
  Caller     : Bio::EnsEMBL::Slice
  Status     : Stable

=cut
*/

Vector *BaseFeatureAdaptor_fetchAllBySliceConstraint(BaseFeatureAdaptor *bfa, Slice *slice, char *constraint, char *logicName) {
  Vector *result = Vector_new();
  char allConstraint[655500];
  char key[655500];
  allConstraint[0] = '\0';
  key[0] = '\0';

  if (constraint != NULL) {
    strcpy(allConstraint, constraint);
  }

  if ( ! BaseFeatureAdaptor_logicNameToConstraint(bfa, allConstraint, logicName)) {
  // If the logic name was invalid, undef was returned
    return result;
  }
 
  // Will only use feature_cache if hasn't got no_cache attribute set
  if ( !DBAdaptor_noCache(bfa->dba)) {

/*
    // strain test and add to constraint if so to stop caching.
    if ( $slice->isa('Bio::EnsEMBL::StrainSlice') ) {
      my $string =
        $self->dbc()->db_handle()->quote( $slice->strain_name() );

      if ( $constraint ne "" ) {
        $constraint .= " AND $string = $string ";
      } else {
        $constraint .= " $string = $string ";
      }
    }
*/

    // Check the cache and return the cached results if we have already
    // done this query.  The cache key is the made up from the slice
    // name, the constraint, and the bound parameters (if there are any).
    sprintf(key, "%s:%s", Slice_getName(slice), constraint);

    StrUtil_strupr(key);
    
 /* In C I don't have bound params, I put them into the constraint, so there should be no need for this bit of the key
    if ( defined($bind_params) ) {
      $key .= ':'
        . join( ':', map { $_->[0] . '/' . $_->[1] } @{$bind_params} );
    }
 */

    if (Cache_contains(bfa->sliceFeatureCache, key)) {
      return Cache_findElem(bfa->sliceFeatureCache, key);
    }
  }

  Vector *projVec = BaseFeatureAdaptor_getAndFilterSliceProjections(bfa, slice);
  int nBound = 0;
  long *bounds = BaseFeatureAdaptor_generateFeatureBounds(bfa, slice, &nBound); 

  // fetch features for the primary slice AND all symlinked slices
  int i;
  for (i=0; i<Vector_getNumElement(projVec); i++) {
    ProjectionSegment *seg = Vector_getElementAt(projVec, i);

    long offset     = ProjectionSegment_getFromStart(seg);
    Slice *segSlice = ProjectionSegment_getToSlice(seg);

    Vector *features = BaseFeatureAdaptor_sliceFetch(bfa, segSlice, allConstraint);

    // If this was a symlinked slice offset the feature coordinates as
    // needed.
    if ( EcoString_strcmp(Slice_getName(segSlice), Slice_getName(slice))) {
      int j;
      for (j=0; j<Vector_getNumElement(features); j++) {
        SeqFeature *f = Vector_getElementAt(features, j);
        if ( offset != 1 ) {
          SeqFeature_setStart(f, (SeqFeature_getStart(f) + (offset-1)));
          SeqFeature_setEnd(f, (SeqFeature_getEnd(f) + (offset-1)));
        }

        // discard boundary crossing features from symlinked regions
        int k;
        int skipFlag = 0;
        for (k=0; k<nBound && !skipFlag; k++) {
          long bound = bounds[k];
          if ( SeqFeature_getStart(f) < bound && SeqFeature_getEnd(f) >= bound ) {
            skipFlag = 1;
          }
        }

        // NIY: Do I need to free slice f was on????
        if (!skipFlag) {
          SeqFeature_setSlice(f, slice);
          Vector_addElement(result, f);
        } else {
          // NIY: Free feature if it was out of bounds
        }
      }
      Vector_free(features);
    } else {
      Vector_append(result, features);
      Vector_free(features);
    }
  }
  // NIY: I need to free ProjectionSegments

  // Will only use feature_cache when set attribute no_cache in DBAdaptor
  // Condition looks slightly odd, but key will have only been set to something if
  // the code entered the noCache controlled condition above
  if (key[0]) {
    Cache_addElement(bfa->sliceFeatureCache, key, result, NULL);
  }

  return result;
}


/*
=head2 fetch_all_by_logic_name

  Arg [1]    : string $logic_name
               the logic name of the type of features to obtain
  Example    : $fs = $a->fetch_all_by_logic_name('foobar');
  Description: Returns a listref of features created from the database.
               only features with an analysis of type $logic_name will
               be returned.  If the logic name is invalid (not in the
               analysis table), a reference to an empty list will be
               returned.
  Returntype : listref of Bio::EnsEMBL::SeqFeatures
  Exceptions : thrown if no $logic_name
  Caller     : General
  Status     : Stable

=cut
*/

Vector *BaseFeatureAdaptor_fetchAllByLogicName(BaseFeatureAdaptor *bfa, char *logicName) {
  if ( logicName == NULL ) {
    fprintf(stderr,"Need a logic_name\n");
    exit(1);
  }
  char constraint[1024];

  if (!BaseFeatureAdaptor_logicNameToConstraint(bfa, constraint, logicName )) {
    fprintf(stderr, "Invalid logic name: %s\n", logicName);
    return Vector_new();
  }

  return BaseAdaptor_genericFetch((BaseAdaptor *)bfa, constraint, NULL, NULL);
}

/*
=head2 fetch_all_by_stable_id_list

  Arg [1]    : string $logic_name
               the logic name of the type of features to obtain
  Arg [2]    : Bio::EnsEMBL::Slice $slice
               the slice from which to obtain features
  Example    : $fs = $a->fetch_all_by_stable_id_list(["ENSG00001","ENSG00002", ...]);
  Description: Returns a listref of features identified by their stable IDs.
               This method only fetches features of the same type as the calling
               adaptor. 
               Results are constrained to a slice if the slice is provided.
  Returntype : listref of Bio::EnsEMBL::Feature
  Exceptions : thrown if no stable ID list is provided.
  Caller     : General
  Status     : Stable

=cut
*/

// Adapted from BaseAdaptor->uncached_fetch_all_by_dbID_list
Vector *BaseFeatureAdaptor_fetchAllByStableIdList(BaseFeatureAdaptor *bfa, Vector *ids, Slice *slice) {

  fprintf(stderr," BaseFeatureAdaptor_fetchAllByStableIdList not implemented until I decide how to do it\n");
  exit(1);
  //return BaseFeatureAdaptor_uncachedFetchAllByIdList(bfa, ids, slice, "stable_id");
}

/*
# Method that creates an object.  Called by the _objs_from_sth() method
# in the sub-classes (the various feature adaptors).  Overridden by the
# feature collection classes.
*/

/* Try to do it without these
sub _create_feature {
  my ( $self, $feature_type, $args ) = @_;
  return $feature_type->new( %{$args} );
}
*/

/*
# This is the same as the above, but calls the new_fast() constructor of
# the feature type.
*/

/* Try to do it without these
sub _create_feature_fast {
  my ( $self, $feature_type, $args ) = @_;
  return $feature_type->new_fast($args);
}
*/

/*
=head2 count_by_Slice_constraint

    Arg [1]     : Bio::EnsEMBL::Slice
    Arg [2]     : String Custom SQL constraint
    Arg [3]     : String Logic name to search by
    Description : Finds all features with at least partial overlap to the given
                  slice and sums them up
    Returntype  : Integer
=cut
*/

int BaseFeatureAdaptor_countBySliceConstraint(BaseFeatureAdaptor *bfa, Slice *slice, char *constraint, char *logicName) {
  int count = 0;
  
  //Table synonym
  NameTableType *tables = bfa->getTables();
  char **primTab = (*tables)[0];
  char *tableName = primTab[NAME];
  char *tableSynonym = primTab[SYN];
  
  //Constraints
  char allConstraint[655500];
  allConstraint[0] = '\0';

  if (constraint != NULL) {
    strcpy(allConstraint, constraint);
  }

  if ( ! BaseFeatureAdaptor_logicNameToConstraint(bfa, allConstraint, logicName)) {
  // If the logic name was invalid, undef was returned
    return 0;
  }

  //Query logic
  // Doesn't seem to be used in perl SliceAdaptor *sa = Slice_getAdaptor(slice);
  Vector *projVec = BaseFeatureAdaptor_getAndFilterSliceProjections(bfa, slice);

  int projVecLen = Vector_getNumElement(projVec);
  //Manual loop to support look-ahead/behind
  int i;
  for (i = 0; i < projVecLen; i++) {
    ProjectionSegment *seg = Vector_getElementAt(projVec, i);
    Slice *segSlice = ProjectionSegment_getToSlice(seg);
    
    // Because we cannot filter boundary crossing features in code we need to
    // do it in SQL. So we detect when we are not on the original query Slice
    // we do manual filtering by the seg_Slice's start and end. If we are on
    // the *1st* section we only limit by the *end* and when we are on the *last*
    // we filter by the *start*
    // 
    // This is the same as fetch_all_by_Slice_constraint()'s in-memory filtering
    // except we need to alter projected features in that code with an offset
    if ( EcoString_strcmp(Slice_getName(segSlice), Slice_getName(slice))) {
      //limit both by default
      int limitStart = 1;
      int limitEnd   = 1; 

      // Special cases at ends
      if (i == 0) {
        limitStart = 0; // don't check start as we are on the first projection
      } else if (i == (projVecLen-1)) {
        limitEnd   = 0; // don't check end as we are on the final projection
      }

      if (allConstraint[0]) {
        strcat(allConstraint, " AND ");
      }

      char tmpStr[1024];
      //Do not cross the start boundary so our feature must be less than slice end on all counts
      if (limitStart) {
        sprintf(tmpStr,"%s.seq_region_start <= %ld AND %s.seq_region_end <= %ld", 
                       tableSynonym, Slice_getEnd(segSlice), tableSynonym, Slice_getEnd(segSlice));
        strcat(allConstraint, tmpStr);
      }
      //Do not cross the start boundary so our feature must be larger than slice start on all counts
      if (limitEnd) {
        if (limitStart) strcat(allConstraint, " AND ");

        sprintf(tmpStr, "%s.seq_region_start >= %ld AND %s.seq_region_end >= %ld",
                       tableSynonym, Slice_getStart(segSlice), tableSynonym, Slice_getStart(segSlice));
        strcat(allConstraint, tmpStr);
      }
    }
    
    Vector *countVec = BaseFeatureAdaptor_getBySlice(bfa, segSlice, allConstraint, "count");

    // Data comes out as an array
    int j;
    for (j=0; j<Vector_getNumElement(countVec); j++) {
      int *cnt = Vector_getElementAt(countVec, j);
      count += *cnt;
    }
  }
  
  return count;
}

/*
=head2 _get_and_filter_Slice_projections


    Arg [1]     : Bio::EnsEMBL::Slice
    Description : Delegates onto SliceAdaptor::fetch_normalized_slice_projection() 
                  with filtering on
    Returntype  : ArrayRef Bio::EnsEMBL::ProjectionSegment; Returns an array
                  of projected segments
=cut
*/
Vector *BaseFeatureAdaptor_getAndFilterSliceProjections(BaseFeatureAdaptor *bfa, Slice *slice) {
  SliceAdaptor *sa = (SliceAdaptor *)Slice_getAdaptor(slice);
  int filterProjections = 1;
  return SliceAdaptor_fetchNormalizedSliceProjection(sa, slice, filterProjections);
}

/*
=head2 _generate_feature_bounds

    Arg [1]     : Bio::EnsEMBL::Slice
    Description : Performs a projection of Slice and records the bounds
                  of that projection. This can be used later on to restrict
                  Features which overlap into unwanted areas such as
                  regions which exist on another HAP/PAR region.
                  
                  Bounds are defined as projection_start - slice_start + 1.
    Example     : my $bounds = $self->_generate_feature_bounds($slice);
    Returntype  : ArrayRef Integer; Returns the location of the bounds.
=cut
*/

long *BaseFeatureAdaptor_generateFeatureBounds(BaseFeatureAdaptor *bfa, Slice *slice, int *nBound) {
  SliceAdaptor *sa = (SliceAdaptor *)Slice_getAdaptor(slice);

  // construct list of Hap/PAR boundaries for entire seq region
  IDType srId = Slice_getSeqRegionId(slice);

  Slice *entSlice = SliceAdaptor_fetchBySeqRegionId(sa, srId, POS_UNDEF, POS_UNDEF, STRAND_UNDEF);

  if (Slice_getStrand(slice) == -1 ) {
    Slice *toDel = entSlice;
    entSlice = Slice_invert(entSlice);
    Slice_free(toDel);
  }

  Vector *entProj = SliceAdaptor_fetchNormalizedSliceProjection(sa, entSlice, 0);

  long *bounds;

  if ((bounds = calloc(Vector_getNumElement(entProj), sizeof(long))) == NULL) {
    fprintf(stderr, "Failed allocating bounds array for %d bounds\n", Vector_getNumElement(entProj));
    exit(1);
  }

  *nBound = Vector_getNumElement(entProj)-1;

  int i;

  // skip first; 1st does not have bounds normally; may change if we ever have a patch a pos 1 
  for (i=1; i<Vector_getNumElement(entProj); i++) {
    ProjectionSegment *seg = Vector_getElementAt(entProj, i);
    bounds[i-1] = ProjectionSegment_getFromStart(seg) - Slice_getStart(slice) + 1;

    // Take the opportunity to free seg
    ProjectionSegment_free(seg);
  }
  Vector_free(entProj);
  Slice_free(entSlice);

  return bounds;
}

/*
=head2 _get_by_Slice
    Arg [0]    : Bio::EnsEMBL::Slice to find all the features within
    Arg [1]    : SQL constraint string
    Arg [2]    : Type of query to run. Default behaviour is to select, but 
                 'count' is also valid
    Description: Abstracted logic from _slice_fetch
    Returntype : Listref of Bio::EnsEMBL::Feature, or integers for counting mode
=cut
*/

typedef struct queryAccumDataStruct {
  char *constraint;
  AssemblyMapper *mapper;
  Slice *slice;
} QueryAccumData;

QueryAccumData *QueryAccumData_new(char *constraint, AssemblyMapper *mapper, Slice *slice) {
  QueryAccumData *qad;
  if ((qad = calloc(1,sizeof(QueryAccumData))) == NULL) {
    fprintf(stderr, "Failed allocating space for qad\n");
    exit(1);
  }
  qad->constraint = constraint;
  qad->mapper = mapper;
  qad->slice  = slice;

  return qad;
}

void QueryAccumData_free(QueryAccumData *qad) {
  free(qad->constraint); 
  // Don't free mapper or slice (we don't own those, well certainly not the mapper)
  free(qad);
}

Vector *BaseFeatureAdaptor_getBySlice(BaseFeatureAdaptor *bfa, Slice *slice, char *origConstraint, char *queryType) {
  // features can be scattered across multiple coordinate systems
  NameTableType *tables = bfa->getTables();
  char **primTab = (*tables)[0];
  char *tableName = primTab[NAME];
  char *tableSynonym = primTab[SYN];

  AssemblyMapper *mapper;

  MetaContainer *metaContainer = DBAdaptor_getMetaContainer(bfa->dba);
  MetaCoordContainer *metaCoordContainer = DBAdaptor_getMetaCoordContainer(bfa->dba);

  Vector *featureCoordSystems;

  char tmpStr[1024];
  sprintf(tmpStr, "%sbuild.level", tableName);  
  Vector *metaValues = MetaContainer_listValueByKey(metaContainer, tmpStr);

  if (Vector_getNumElement(metaValues) && Slice_isTopLevel(slice)) {
    featureCoordSystems = Vector_new();
    Vector_addElement(featureCoordSystems, Slice_getCoordSystem(slice));
  } else {
    featureCoordSystems = MetaCoordContainer_fetchAllCoordSystemsByFeatureType(metaCoordContainer, tableName);
  }
  
  AssemblyMapperAdaptor *ama = DBAdaptor_getAssemblyMapperAdaptor(bfa->dba);

  Vector *panCoordFeatures = Vector_new();
        
  int i;
  for (i=0; i<Vector_getNumElement(featureCoordSystems); i++) {
    int doneCoordSystem = 0;

    CoordSystem *coordSystem = Vector_getElementAt(featureCoordSystems, i);

    Vector *queryAccumulator = Vector_new();

    // Build up a combination of query constraints that will quickly establish the result set
    char *constraint;
    // Note allocated because stored in struct for later execution
    if ((constraint = calloc(655500, sizeof(char))) == NULL) {
      fprintf(stderr,"Failed allocating constraint\n");
      exit(1);
    }

    constraint[0] = '\0';
    if (origConstraint) {
      strcpy(constraint, origConstraint);
    }

    if ( ! CoordSystem_compare(coordSystem, Slice_getCoordSystem(slice))) {
      long maxLen = BaseFeatureAdaptor_getMaxFeatureLength(bfa);
      if (!maxLen) {
        maxLen = MetaCoordContainer_fetchMaxLengthByCoordSystemFeatureType(metaCoordContainer,  coordSystem, tableName);
      }
                       
      IDType seqRegionId;

      if (Slice_getAdaptor(slice)) {
        seqRegionId = SliceAdaptor_getSeqRegionId((SliceAdaptor *)Slice_getAdaptor(slice), slice);
      } else {
        SliceAdaptor *sa = DBAdaptor_getSliceAdaptor(bfa->dba);

        seqRegionId = SliceAdaptor_getSeqRegionId(sa, slice);
      }
            
/* Not bothering about the external seq region ids for now
      my @seq_region_ids = ($seq_region_id);
      while (1) {
        my $ext_seq_region_id = $self->get_seq_region_id_external($seq_region_id);
        
        if ( $ext_seq_region_id == $seq_region_id ) { last }
        
        push( @seq_region_ids, $ext_seq_region_id );
        $seq_region_id = $ext_seq_region_id;
      }
*/
      if (constraint[0]) {
        strcat(constraint," AND ");
      }

      sprintf(tmpStr, "%s.seq_region_id = "IDFMTSTR" AND ", tableSynonym, seqRegionId);
      strcat(constraint, tmpStr);
            
      //faster query for 1bp slices where SNP data is not compressed
      if (BaseFeatureAdaptor_getStartEqualsEnd(bfa) && Slice_getStart(slice) == Slice_getEnd(slice)) {
        sprintf(tmpStr, " AND %s.seq_region_start = %ld AND %s.seq_region_end = %ld",
                tableSynonym, Slice_getEnd(slice), tableSynonym, Slice_getStart(slice));;
        strcat(constraint, tmpStr);
      } else {
        //if ( !$slice->is_circular() ) 
        if (1) {
          // Deal with the default case of a non-circular chromosome.
          sprintf(tmpStr,"%s.seq_region_start <= %ld AND %s.seq_region_end >= %ld", 
                  tableSynonym, Slice_getEnd(slice), tableSynonym, Slice_getStart(slice));
          strcat(constraint, tmpStr);
            
          if (maxLen) {
            long minStart = Slice_getStart(slice) - maxLen;
            sprintf(tmpStr," AND %s.seq_region_start >= %ld", tableSynonym, minStart);
            strcat(constraint, tmpStr);
          }
        } else {
/* Don't deal with circular - not supported in C implementation
          // Deal with the case of a circular chromosome.
          if ( $slice->start > $slice->end ) {
            $constraint .= " ( ".$table_synonym.".seq_region_start >= ".$slice->start
                            . " OR ".$table_synonym.".seq_region_start <= ".$slice->end
                            . " OR ".$table_synonym.".seq_region_end >= ".$slice->start
                            . " OR ".$table_synonym.".seq_region_end <= ".$slice->end
                            . " OR ".$table_synonym.".seq_region_start > ".$table_synonym.".seq_region_end";
          } else {
            $constraint .= " ((".$table_synonym.".seq_region_start <= ".$slice->end
                            . " AND ".$table_synonym.".seq_region_end >= ".$slice->start.") "
                            . "OR (".$table_synonym.".seq_region_start > ".$table_synonym.".seq_region_end"
                            . " AND (".$table_synonym.".seq_region_start <= ".$slice->end
                            . " OR ".$table_synonym.".seq_region_end >= ".$slice->start.")))";
          }
*/
        }
      }

      Vector_addElement(queryAccumulator, QueryAccumData_new(constraint, NULL, slice));
    } else { 
      //coordinate systems do not match
      mapper = AssemblyMapperAdaptor_fetchByCoordSystems(ama, Slice_getCoordSystem(slice), coordSystem);

      if (mapper == NULL) continue;

      // Get list of coordinates and corresponding internal ids for
      // regions the slice spans
      MapperRangeSet *coords = AssemblyMapper_map(mapper, Slice_getSeqRegionName(slice), Slice_getStart(slice), Slice_getEnd(slice), 
                                                  Slice_getStrand(slice), Slice_getCoordSystem(slice), 0, NULL);

      MapperRangeSet_removeGaps(coords);
//      @coords = grep { !$_->isa('Bio::EnsEMBL::Mapper::Gap') } @coords;

// This is a very big jump - I'm going to have to think about how best to do this - there's tidying to do
      //next CORD_SYSTEM if ( !@coords );
      if ( ! MapperRangeSet_getNumRange(coords)) {
        // Do the 'next COORD_SYSTEM' with this flag
        doneCoordSystem = 1;

      } else {
        if ( MapperRangeSet_getNumRange(coords) > MAX_SPLIT_QUERY_SEQ_REGIONS //&& ! $slice->isa('Bio::EnsEMBL::LRGSlice') 
                      && strcmp(Slice_getCoordSystemName(slice), "lrg")) {

// Think this is already done above         $constraint = $orig_constraint;
          if (constraint[0]) {
            strcat(constraint, " AND ");
          }

          sprintf(tmpStr, "%s.seq_region_id IN (", tableSynonym);
          int length = MapperRangeSet_getNumRange(coords);
          int j;
          for (j = 0; j < length; j++ ) {
            MapperCoordinate *coord = (MapperCoordinate *)MapperRangeSet_getRangeAt(coords, j);
            if (j!=0) {
              strcat(tmpStr,", ");
            }
            sprintf(tmpStr,"%s"IDFMTSTR, tmpStr, coord->id);
          }
          strcat(tmpStr, ")");
          strcat(constraint, tmpStr);
                  
          Vector_addElement(queryAccumulator, QueryAccumData_new(constraint, mapper, slice));

        } else if (strcmp(Slice_getCoordSystemName(slice), "lrg")) {
          long maxLen = BaseFeatureAdaptor_getMaxFeatureLength(bfa);
          if (!maxLen) {
            maxLen = MetaCoordContainer_fetchMaxLengthByCoordSystemFeatureType(metaCoordContainer,  coordSystem, tableName);
          }
  
          // Free constraint because its going to be reallocated each time through loop
          free(constraint);

          int length = MapperRangeSet_getNumRange(coords);
          int j;
          for (j = 0; j < length; j++ ) {
            MapperCoordinate *coord = (MapperCoordinate *)MapperRangeSet_getRangeAt(coords, j);
  
            // Note allocated because stored in struct for later execution
            if ((constraint = calloc(655500, sizeof(char))) == NULL) {
              fprintf(stderr,"Failed allocating constraint\n");
              exit(1);
            }
        
            constraint[0] = '\0';
            if (origConstraint) {
              strcpy(constraint, origConstraint);
            }

            if (constraint[0]) {
              strcat(constraint, " AND ");
            }
            sprintf(tmpStr," %s.seq_region_id = "IDFMTSTR" AND %s.seq_region_start <= %ld AND %s.seq_region_end >= %ld",
                    tableSynonym, coord->id, tableSynonym, coord->end, tableSynonym, coord->start);
            strcat(constraint, tmpStr);
  
            if (maxLen) {
              long minStart = coord->start - maxLen;
              sprintf(tmpStr, " AND %s.seq_region_start >= %ld", tableSynonym, minStart);
              strcat(constraint, tmpStr);
            }
                      
            Vector_addElement(queryAccumulator, QueryAccumData_new(constraint, mapper, slice));
          }
        } else { // LRG - ignore this stuff
          fprintf(stderr,"Ignoring lrg coord system\n");
        }
      }
      // NIY: Free mapper range set
    }
     
    int j;
    for (j=0; j<Vector_getNumElement(queryAccumulator) && !doneCoordSystem; j++) {
      QueryAccumData *qad = Vector_getElementAt(queryAccumulator, j);

      if (!strcmp(queryType, "count")) {
        int count = BaseAdaptor_genericCount((BaseAdaptor *)bfa, qad->constraint);

        int *countP;
        if ((countP = calloc(1,sizeof(int))) == NULL) {
          fprintf(stderr, "Failed allocating space for a count\n");
          exit(1);
        }
    
        *countP = count;
        Vector_addElement(panCoordFeatures, countP);

      } else {
        Vector *features = BaseAdaptor_genericFetch((BaseAdaptor *)bfa, qad->constraint, qad->mapper,  qad->slice);
        fprintf(stderr,"Here!!!!!!!!!!!!!!!!!! with %d features and %d pan coord features\n", Vector_getNumElement(features),  Vector_getNumElement(panCoordFeatures));
        Vector *remappedFeatures = BaseFeatureAdaptor_remap(bfa, features, qad->mapper, qad->slice);

        Vector_append(panCoordFeatures, remappedFeatures);
        fprintf(stderr,"Here 22!!!!!!!!!!!!!!!!!! with %d features and %d remapped features and %d pan coord features\n", Vector_getNumElement(features), Vector_getNumElement(remappedFeatures), Vector_getNumElement(panCoordFeatures));

        // NIY: Free stuff (features vector etc)
        // NIY: Note hack hack hack to unset any freeFunc that has been set on the features vector - to stop the features we've transferred to another vector being freed
        Vector_setFreeFunc(features, NULL);
        Vector_free(features);
        Vector_free(remappedFeatures);
      }

      QueryAccumData_free(qad);
    }
    // NIY: Free query accumulator stuff
    Vector_free(queryAccumulator);

    mapper = NULL;
  }

  return panCoordFeatures;
}
/*
#
# helper function used by fetch_all_by_Slice_constraint method
#
*/
// NIY: Figure out if this method is needed at all - it seems to just pass straight through to getBySlice
Vector *BaseFeatureAdaptor_sliceFetch(BaseFeatureAdaptor *bfa, Slice *slice, char *origConstraint) {
  Vector *features = BaseFeatureAdaptor_getBySlice(bfa, slice, origConstraint, "fetch");

  return features;
}


/*
#for a given seq_region_id, gets the one used in an external database, if present, otherwise, returns the internal one
sub get_seq_region_id_external {
  my ( $self, $sr_id ) = @_;
  my $cs_a = $self->db()->get_CoordSystemAdaptor();
  return ( exists( $cs_a->{'_internal_seq_region_mapping'}->{$sr_id} )
           ? $cs_a->{'_internal_seq_region_mapping'}->{$sr_id}
           : $sr_id );
}

*/
/*
#for a given seq_region_id and coord_system, gets the one used in the internal (core) database
sub get_seq_region_id_internal{
  my ( $self, $sr_id ) = @_;
  my $cs_a = $self->db()->get_CoordSystemAdaptor();
  return (  exists $cs_a->{'_external_seq_region_mapping'}->{$sr_id} 
            ? $cs_a->{'_external_seq_region_mapping'}->{$sr_id} 
            : $sr_id);
}
*/

/*
#
# Helper function containing some common feature storing functionality
#
# Given a Feature this will return a copy (or the same feature if no changes 
# to the feature are needed) of the feature which is relative to the start
# of the seq_region it is on. The seq_region_id of the seq_region it is on
# is also returned.
#
# This method will also ensure that the database knows which coordinate
# systems that this feature is stored in.
#

sub _pre_store {
  my $self    = shift;
  my $feature = shift;

  if(!ref($feature) || !$feature->isa('Bio::EnsEMBL::Feature')) {
    throw('Expected Feature argument.');
  }
  my $slice = $feature->slice();

  $self->_check_start_end_strand($feature->start(),$feature->end(),
                                 $feature->strand(), $slice);


  my $db = $self->db();

  my $slice_adaptor = $db->get_SliceAdaptor();

  if(!ref($slice) || !($slice->isa('Bio::EnsEMBL::Slice') or $slice->isa('Bio::EnsEMBL::LRGSlice'))  ) {
    throw('Feature must be attached to Slice to be stored.');
  }

  # make sure feature coords are relative to start of entire seq_region

  if($slice->start != 1 || $slice->strand != 1) {
    #move feature onto a slice of the entire seq_region
    $slice = $slice_adaptor->fetch_by_region($slice->coord_system->name(),
                                             $slice->seq_region_name(),
                                             undef, #start
                                             undef, #end
                                             undef, #strand
                                             $slice->coord_system->version());

    $feature = $feature->transfer($slice);

    if(!$feature) {
      throw('Could not transfer Feature to slice of ' .
            'entire seq_region prior to storing');
    }
  }

  # Ensure this type of feature is known to be stored in this coord system.
  my $cs = $slice->coord_system;

  my ($tab) = $self->_tables();
  my $tabname = $tab->[0];

  my $mcc = $db->get_MetaCoordContainer();

  $mcc->add_feature_type($cs, $tabname, $feature->length);

  my $seq_region_id = $slice_adaptor->get_seq_region_id($slice);

  if(!$seq_region_id) {
    throw('Feature is associated with seq_region which is not in this DB.');
  }

  return ($feature, $seq_region_id);
}


# The same function as _pre_store
# This one is used to store user uploaded features in XXX_userdata db

sub _pre_store_userdata {
  my $self    = shift;
  my $feature = shift;

  if(!ref($feature) || !$feature->isa('Bio::EnsEMBL::Feature')) {
    throw('Expected Feature argument.');
  }

  my $slice = $feature->slice();
  my $slice_adaptor = $slice->adaptor;
  
  $self->_check_start_end_strand($feature->start(),$feature->end(),
                                 $feature->strand(), $slice);


  if(!ref($slice) || !($slice->isa('Bio::EnsEMBL::Slice') or $slice->isa('Bio::EnsEMBL::LRGSlice')) ) {
    throw('Feature must be attached to Slice to be stored.');
  }

  # make sure feature coords are relative to start of entire seq_region

  if($slice->start != 1 || $slice->strand != 1) {
    #move feature onto a slice of the entire seq_region
    $slice = $slice_adaptor->fetch_by_region($slice->coord_system->name(),
                                             $slice->seq_region_name(),
                                             undef, #start
                                             undef, #end
                                             undef, #strand
                                             $slice->coord_system->version());

    $feature = $feature->transfer($slice);

    if(!$feature) {
      throw('Could not transfer Feature to slice of ' .
            'entire seq_region prior to storing');
    }
  }

  # Ensure this type of feature is known to be stored in this coord system.
  my $cs = $slice->coord_system;

  my ($tab) = $self->_tables();
  my $tabname = $tab->[0];

  my $db = $self->db;
  my $mcc = $db->get_MetaCoordContainer();

  $mcc->add_feature_type($cs, $tabname, $feature->length);

  my $seq_region_id = $slice_adaptor->get_seq_region_id($slice);

  if(!$seq_region_id) {
    throw('Feature is associated with seq_region which is not in this DB.');
  }

  return ($feature, $seq_region_id);
}
*/


/*
#
# helper function used to validate start/end/strand and 
# hstart/hend/hstrand etc.
#
*/
/*
sub _check_start_end_strand {
  my $self = shift;
  my $start = shift;
  my $end   = shift;
  my $strand = shift;
  my $slice = shift;

  // 
  // Make sure that the start, end, strand are valid
  //
  if(int($start) != $start) {
    throw("Invalid Feature start [$start].  Must be integer.");
  }
  if(int($end) != $end) {
    throw("Invalid Feature end [$end]. Must be integer.");
  }
  if(int($strand) != $strand || $strand < -1 || $strand > 1) {
    throw("Invalid Feature strand [$strand]. Must be -1, 0 or 1.");
  }
  if($end < $start && !$slice->is_circular()) {
    throw("Invalid Feature start/end [$start/$end]. Start must be less " .
          "than or equal to end.");
  }

  return 1;
}
*/


/*
#
# Given a list of features checks if they are in the correct coord system
# by looking at the first features slice.  If they are not then they are
# converted and placed on the slice.
#
*/
Vector *BaseFeatureAdaptor_remap(BaseFeatureAdaptor *bfa, Vector *features, AssemblyMapper *mapper, Slice *slice) {
  //check if any remapping is actually needed
  if (Vector_getNumElement(features) > 0) {
    SeqFeature *f = Vector_getElementAt(features, 0);
// NIY: Check this should actually be == rather than slice equivalence
    if (SeqFeature_getSlice(f) == slice) {
      return features;
    }
  }

  // remapping has not been done, we have to do our own conversion from
  // to slice coords
  Vector *out = Vector_new();

  long sliceStart      = Slice_getStart(slice);
  long sliceEnd        = Slice_getEnd(slice);
  int sliceStrand      = Slice_getStrand(slice);
  CoordSystem *sliceCs = Slice_getCoordSystem(slice);

  //char *seqRegion;
  IDType seqRegionId;
  long start;
  long end;
  int strand;

  IDType sliceSeqRegionId = Slice_getSeqRegionId(slice);
  char *sliceSeqRegion = Slice_getSeqRegionName(slice);

  int nFeat = Vector_getNumElement(features);

  int i;
  for (i=0; i<nFeat; i++) {
    SeqFeature *f = Vector_getElementAt(features, i);

    // since feats were obtained in contig coords, attached seq is a contig
    Slice *fSlice = SeqFeature_getSlice(f);

    if (fSlice == NULL) {
      fprintf(stderr, "Feature does not have attached slice.\n");
      exit(1);
    }

    char *fSeqRegion = Slice_getSeqRegionName(fSlice);
    IDType fSeqRegionId = Slice_getSeqRegionId(fSlice);
    CoordSystem *fCs = Slice_getCoordSystem(fSlice);

    if (CoordSystem_compare(sliceCs, fCs)) {
      // slice of feature in different coord system, mapping required

      // Note my implementation of fastMap is different to perl, in that it returns a MapperRangeSet rather than a list specifying a single location
      MapperRangeSet *mrs = AssemblyMapper_fastMap(mapper, fSeqRegion, SeqFeature_getStart(f), SeqFeature_getEnd(f), SeqFeature_getStrand(f), fCs, NULL);

      // empty set means gap
      if (MapperRangeSet_getNumRange(mrs) == 0) {
        MapperRangeSet_free(mrs);
        continue;
      }

      MapperCoordinate *mc = (MapperCoordinate *)MapperRangeSet_getRangeAt(mrs, 0);
    
      //seqRegion = mc->
      seqRegionId = mc->id;
      start       = mc->start;
      end         = mc->end;
      strand      = mc->strand;

      MapperRangeSet_free(mrs);
    } else {
      start       = SeqFeature_getStart(f);
      end         = SeqFeature_getEnd(f);
      strand      = SeqFeature_getStrand(f);
      seqRegionId = Slice_getSeqRegionId(SeqFeature_getSlice(f));
      //seqRegion  = Slice_getSeqRegionName(SeqFeature_getSlice(f));
    }
    
    // maps to region outside desired area
    if ((start > sliceEnd) || (end < sliceStart) || seqRegionId != sliceSeqRegionId) { // was strcmp(sliceSeqRegion, seqRegion))
      // NIY: Anything need freeing???
      continue;
    }

    // shift the feature start, end and strand in one call
    if (sliceStrand == -1) {
      SeqFeature_setStart (f, sliceEnd - end + 1);
      SeqFeature_setEnd   (f, sliceEnd - start + 1);
      SeqFeature_setStrand(f, strand * -1 );
    } else {
      SeqFeature_setStart (f, start - sliceStart + 1);
      SeqFeature_setEnd   (f, end - sliceStart + 1);
      SeqFeature_setStrand(f, strand);
    }
      
    SeqFeature_setSlice(f, slice);

    //fprintf(stderr, "Cigar after remap = %s\n", DNAAlignFeature_getCigarString((DNAAlignFeature *)f));

    Vector_addElement(out, f);
  }

  return out;
}


/*
#
# Given a logic name and an existing constraint this will
# add an analysis table constraint to the feature.  Note that if no
# analysis_id exists in the columns of the primary table then no
# constraint is added at all
#
*/
char *BaseFeatureAdaptor_logicNameToConstraint(BaseFeatureAdaptor *bfa, char *constraint, char *logicName) {

  if (logicName == NULL) {
    return constraint;
  }

  // make sure that an analysis_id exists in the primary table
  NameTableType *tables = bfa->getTables();
  char **primTab = (*tables)[0];
  char *primSynonym = primTab[SYN];

  int foundAnalysis=0;
  
  char **columns = bfa->getColumns();
  int i=0;
  while (columns[i] != NULL) {
    char syn[1024];
    char colName[1024];

    char *dotP = strchr(columns[i], '.');
    if (dotP == NULL) {
      fprintf(stderr, "No '.' in columns string %s\n", columns[i]);
      exit(1);
    }
    int synLen = dotP-columns[i];
    strncpy(syn, columns[i], synLen);
    syn[synLen] = '\0';
    strcpy(colName, dotP+1);
    
    if (!strcmp(syn, primSynonym)) {
      if (!strcmp(colName, "analysis_id")) {
        foundAnalysis = 1;
        break;
      }
    }
    i++;
  }

  if (!foundAnalysis) {
    fprintf(stderr,"This feature is not associated with an analysis.\n"
                   "Ignoring logic_name argument = [%s].\n", logicName);
    return constraint;
  }

  AnalysisAdaptor *aa = DBAdaptor_getAnalysisAdaptor(bfa->dba);
  Analysis *an = AnalysisAdaptor_fetchByLogicName(aa, logicName);

  if (an == NULL) {
  // Perl returned undef, I'm going to throw for now
  //  return undef;
    fprintf(stderr, "Unknown analysis %s in BFA_logicNameToConstraint - exiting\n", logicName);
    exit(1);
  }

  IDType anId = Analysis_getDbID(an);

  if (constraint[0]) {
    strcat(constraint, " AND");
  }
  char tmpStr[1024];
  sprintf(tmpStr, " %s.analysis_id = "IDFMTSTR, primSynonym, anId);
  strcat(constraint, tmpStr);

  return constraint;
}


/*
=head2 store

  Arg [1]    : list of Bio::EnsEMBL::SeqFeature
  Example    : $adaptor->store(@feats);
  Description: ABSTRACT  Subclasses are responsible for implementing this 
               method.  It should take a list of features and store them in 
               the database.
  Returntype : none
  Exceptions : thrown method is not implemented by subclass
  Caller     : general
  Status     : At Risk
             : throws if called.

=cut
*/

/* Moved to BaseAdaptor
void BaseFeatureAdaptor_store(BaseFeatureAdaptor *bfa) {

  fprintf(stderr, "Abstract method store not defined by implementing subclass\n");
  exit(1);
}
*/


/*
=head2 remove

  Arg [1]    : A feature $feature 
  Example    : $feature_adaptor->remove($feature);
  Description: This removes a feature from the database.  The table the
               feature is removed from is defined by the abstract method
               _tablename, and the primary key of the table is assumed
               to be _tablename() . '_id'.  The feature argument must 
               be an object implementing the dbID method, and for the
               feature to be removed from the database a dbID value must
               be returned.
  Returntype : none
  Exceptions : thrown if $feature arg does not implement dbID(), or if
               $feature->dbID is not a true value
  Caller     : general
  Status     : Stable

=cut


sub remove {
  my ($self, $feature) = @_;

  if(!$feature || !ref($feature) || !$feature->isa('Bio::EnsEMBL::Feature')) {
    throw('Feature argument is required');
  }

  if(!$feature->is_stored($self->db)) {
    throw("This feature is not stored in this database");
  }

  my @tabs = $self->_tables;
  my ($table) = @{$tabs[0]};

  my $sth = $self->prepare("DELETE FROM $table WHERE ${table}_id = ?");
  $sth->bind_param(1,$feature->dbID,SQL_INTEGER);
  $sth->execute();

  #unset the feature dbID ad adaptor
  $feature->dbID(undef);
  $feature->adaptor(undef);

  return;
}


=head2 remove_by_Slice

  Arg [1]    : Bio::Ensembl::Slice $slice
  Example    : $feature_adaptor->remove_by_Slice($slice);
  Description: This removes features from the database which lie on a region
               represented by the passed in slice.  Only features which are
               fully contained by the slice are deleted; features which overlap
               the edge of the slice are not removed.
               The table the features are removed from is defined by
               the abstract method_tablename.
  Returntype : none
  Exceptions : thrown if no slice is supplied
  Caller     : general
  Status     : Stable

=cut

sub remove_by_Slice {
  my ($self, $slice) = @_;

  if(!$slice || !ref($slice) || !($slice->isa('Bio::EnsEMBL::Slice') or $slice->isa('Bio::EnsEMBL::LRGSlice')) ) {
    throw("Slice argument is required");
  }

  my @tabs = $self->_tables;
  my ($table_name) = @{$tabs[0]};

  my $seq_region_id = $self->db->get_SliceAdaptor->get_seq_region_id($slice);
  my $start = $slice->start();
  my $end   = $slice->end();

  #
  # Delete only features fully on the slice, not overlapping ones
  #
  my $sth = $self->prepare("DELETE FROM $table_name " .
                           "WHERE seq_region_id = ? " .
                           "AND   seq_region_start >= ? " .
                           "AND   seq_region_end <= ?");

  $sth->bind_param(1,$seq_region_id,SQL_INTEGER);
  $sth->bind_param(2,$start,SQL_INTEGER);
  $sth->bind_param(3,$end,SQL_INTEGER);
  $sth->execute();
  $sth->finish();
}
*/


/*
#
# Internal function. Allows the max feature length which is normally
# retrieved from the meta_coord table to be overridden.  This allows
# for some significant optimizations to be put in when it is known
# that requested features will not be over a certain size.
#
*/
// Do in header
/*
sub _max_feature_length {
  my $self = shift;
  return $self->{'_max_feature_length'} = shift if(@_);
  return $self->{'_max_feature_length'};
}
*/


/*
#
# Lists all seq_region_ids that a particular feature type is found on.
# Useful e.g. for finding out which seq_regions have genes.
# Returns a listref of seq_region_ids.
#
*/
Vector *BaseFeatureAdaptor_listSeqRegionIds(BaseFeatureAdaptor *bfa, char *table) {
  Vector *out = Vector_new();
  
  char qStr[1024];
  sprintf(qStr,"SELECT DISTINCT sr.seq_region_id "
               "FROM      seq_region sr, "
                         "%s a, "
                         "coord_system cs "
               "WHERE     sr.seq_region_id = a.seq_region_id "
                 "AND     sr.coord_system_id = cs.coord_system_id "
                 "AND     cs.species_id = "IDFMTSTR, table, BaseFeatureAdaptor_getSpeciesId(bfa));

  StatementHandle *sth = bfa->prepare((BaseAdaptor *)bfa,qStr,strlen(qStr));

  sth->execute(sth);

  ResultRow *row;
  while (row = sth->fetchRow(sth)) {
    IDType id = row->getLongLongAt(row, 0);
    IDType *idP;

    if ((idP = calloc(1,sizeof(IDType))) == NULL) {
      fprintf(stderr, "Failed allocating space for a id\n");
      exit(1);
    }

    *idP = id;
    Vector_addElement(out, idP);
  }

  sth->finish(sth);

  return out;
}


/*
=head1 DEPRECATED METHODS

=cut


=head2 fetch_all_by_RawContig_constraint

  Description: DEPRECATED use fetch_all_by_RawContig_constraint instead

=cut

sub fetch_all_by_RawContig_constraint {
  my $self = shift;
  deprecate('Use fetch_all_by_Slice_constraint() instead.');
  return $self->fetch_all_by_slice_constraint(@_);
}

=head2 fetch_all_by_RawContig

  Description: DEPRECATED use fetch_all_by_Slice instead

=cut

sub fetch_all_by_RawContig {
  my $self = shift;
  deprecate('Use fetch_all_by_Slice() instead.');
  return $self->fetch_all_by_Slice(@_);
}

=head2 fetch_all_by_RawContig_and_score

  Description: DEPRECATED use fetch_all_by_Slice_and_score instead

=cut

sub fetch_all_by_RawContig_and_score{
  my $self = shift;
  deprecate('Use fetch_all_by_Slice_and_score() instead.');
  return $self->fetch_all_by_Slice_and_score(@_);
}

=head2 remove_by_RawContig

  Description: DEPRECATED use remove_by_Slice instead

=cut

sub remove_by_RawContig {
  my $self = shift;
  deprecate("Use remove_by_Slice instead");
  return $self->remove_by_Slice(@_);
}


sub remove_by_analysis_id {
  my ($self, $analysis_id) = @_;

  $analysis_id or throw("Must call with analysis id");

  my @tabs = $self->_tables;
  my ($tablename) = @{$tabs[0]};

  my $sql = "DELETE FROM $tablename WHERE analysis_id = $analysis_id";
#  warn "SQL : $sql";
      
  my $sth = $self->prepare($sql);
  $sth->execute();
  $sth->finish();
}

sub remove_by_feature_id {
  my ($self, $features_list) = @_;

  my @feats = @$features_list or throw("Must call store with features");

  my @tabs = $self->_tables;
  my ($tablename) = @{$tabs[0]};

  my $sql = sprintf "DELETE FROM $tablename WHERE ${tablename}_id IN (%s)", join ', ', @feats;
#  warn "SQL : $sql";
      
  my $sth = $self->prepare($sql);
  $sth->execute();
  $sth->finish();
}

*/

