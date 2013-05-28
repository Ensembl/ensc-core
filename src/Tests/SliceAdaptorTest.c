#include <stdio.h>

#include "SliceAdaptor.h"
#include "ProjectionSegment.h"
#include "DBAdaptor.h"
#include "EnsC.h"
#include "SeqUtil.h"
#include "StrUtil.h"

#include "BaseTest.h"

void printSlices(Vector *slices);
void testTopLevelLocation(SliceAdaptor *sa, char *location, char *csName, char *seqRegionName, long start, long end, int strand);
void testLocationSlice(char *location, Slice *incomingSlice, char *csName, char *seqRegionName, long start, long end, int strand);

int nTest = 1;

void printSlices(Vector *slices) {
  int i;
  for (i=0;i<Vector_getNumElement(slices); i++) {
    Slice *slice = Vector_getElementAt(slices, i);
    fprintf(stderr,"  %s\n", slice ? Slice_getName(slice) : "UNDEF");
  } 
  fprintf(stderr,"Got %d slices\n", Vector_getNumElement(slices));
}


void testTopLevelLocation(SliceAdaptor *sa, char *location, char *csName, char *seqRegionName, long start, long end, int strand) {
  Slice *incomingSlice = SliceAdaptor_fetchByTopLevelLocation(sa, location, 1, 0);
  fprintf(stderr,"Location test with %s\n", location);
  testLocationSlice(location, incomingSlice, csName, seqRegionName, start, end, strand);
  return;
}

void testLocationSlice(char *location, Slice *incomingSlice, char *csName, char *seqRegionName, long start, long end, int strand) {
  if (strand == STRAND_UNDEF) strand = 1;

  int okFlag = 0;
  okFlag =  ( Slice_getStart(incomingSlice) == start &&
              Slice_getEnd(incomingSlice) == end &&
              Slice_getStrand(incomingSlice) == strand &&
              !strcmp(Slice_getSeqRegionName(incomingSlice), seqRegionName) &&
              !strcmp(Slice_getCoordSystemName(incomingSlice), csName));
  ok(nTest++, okFlag);

  if (!okFlag) {
    fprintf(stderr, " slice name is %s\n", Slice_getName(incomingSlice)); 
  }
  return;
}

static char *CHR = "20";
static long START = 30252000;
static long END   = 31252001;
static long FLANKING = 1000;


int main(int argc, char *argv[]) {
  DBAdaptor *dba;
  Slice *slice;
  SliceAdaptor *sa;

  initEnsC(argc, argv);

  //dba = DBAdaptor_new("ensembldb.ensembl.org","anonymous",NULL,"homo_sapiens_core_70_37",5306,NULL);
  dba = DBAdaptor_new("ens-livemirror.internal.sanger.ac.uk","ensro",NULL,"homo_sapiens_core_70_37",3306,NULL);
  //dba = DBAdaptor_new("genebuild2.internal.sanger.ac.uk","ensadmin","ensembl","steve_hs_testdb",3306,NULL);

  ok(nTest++,!strcmp("NCBI34",DBAdaptor_getAssemblyType(dba)));

  sa = DBAdaptor_getSliceAdaptor(dba);

  ok(nTest++, sa!=NULL);

  slice = SliceAdaptor_fetchByRegion(sa,"chromosome","1",1,250000000,1,NULL,0);

  ok(nTest++, slice!=NULL);

  printf("slice name = %s\n",Slice_getName(slice));
  ok(nTest++, !strcmp(Slice_getName(slice),"chromosome:NCBI33:1:1:250000000:1"));

//  slice = SliceAdaptor_fetchByRegion(sa,"chromosome","Y",1,3000000,1,NULL,0);
  slice = SliceAdaptor_fetchByRegion(sa,"chromosome","1",1,250000000,1,NULL,0);
  char *seq = Slice_getSeq(slice);
  //printf("slice seq = %s\n", seq);
 
  SeqUtil_writeFasta(stdout, Slice_getName(slice), seq, 60);

  fflush(stdout);
  exit(1);

// Tests from perl

//
// fetch_by_region
//
  slice = SliceAdaptor_fetchByRegion(sa, "chromosome",CHR, START, END, STRAND_UNDEF, NULL, 0);
  ok(nTest++, !strcmp(Slice_getSeqRegionName(slice), CHR));
  ok(nTest++, Slice_getStart(slice) == START);
  ok(nTest++, Slice_getEnd(slice) == END);
  ok(nTest++, Slice_getSeqRegionLength(slice) == 62842997);
  fprintf(stderr,"slice seq_region length = %ld\n", Slice_getSeqRegionLength(slice));


//
// fetch_by_contig_name
//

  Vector *projection = Slice_project(slice,"seqlevel", NULL);

//it is important to get a contig not cut off by slice start or end
  if (Vector_getNumElement(projection) <= 2) {
    fprintf(stderr,"Warning: There aren't enough tiles in this path for this test to work\n");
  }

  ProjectionSegment *p1 = Vector_getElementAt(projection, 1);
  
  long chrStart = ProjectionSegment_getFromStart(p1);
  long chrEnd   = ProjectionSegment_getFromEnd(p1);
  Slice *contig = ProjectionSegment_getToSlice(p1);

  ok(nTest++, Slice_getLength(contig) == (chrEnd - chrStart + 1));

//  char *seq1 = Slice_getSubSeq(slice, chrStart, chrEnd, 1);
//  char *seq2 = Slice_getSeq(contig);

//  ok(10, !strcmp(seq1,seq2));


#ifdef DONE
#
# 12-13 fetch_by_fpc_name
#
#my $fpc_name = 'NT_011387';
#$slice = $slice_adaptor->fetch_by_supercontig_name($fpc_name);
#ok($new_slice->chr_start);
#ok($new_slice->chr_end);



#
# 14 - 15 fetch_by_clone_accession
#
#my $clone_acc = 'AL031658';
#$slice = $slice_adaptor->fetch_by_clone_accession($clone_acc);
#$new_slice = $slice_adaptor->fetch_by_clone_accession($clone_acc, $FLANKING);
#ok($new_slice->chr_start == $slice->chr_start - $FLANKING);
#ok($new_slice->chr_end   == $slice->chr_end   + $FLANKING);


#
# 16-17 fetch by transcript_stable_id
#
my $t_stable_id = 'ENST00000217315';
$slice = $slice_adaptor->fetch_by_transcript_stable_id($t_stable_id);
my $new_slice = $slice_adaptor->fetch_by_transcript_stable_id($t_stable_id,
                                                           $FLANKING);

ok($new_slice->start == $slice->start - $FLANKING);
ok($new_slice->end   == $slice->end   + $FLANKING);


#
# 18-19 fetch by transcript_id
#
my $transcript = $db->get_TranscriptAdaptor->fetch_by_stable_id($t_stable_id);
my $tid = $transcript->dbID;
$slice = $slice_adaptor->fetch_by_transcript_id($tid);
$new_slice = $slice_adaptor->fetch_by_transcript_id($tid, $FLANKING);
ok($new_slice->start == $slice->start - $FLANKING);
ok($new_slice->end   == $slice->end   + $FLANKING);
ok($slice->seq_region_length == 62842997);
debug("new slice seq_region length = " . $new_slice->seq_region_length());

#
# 20-23 fetch_by_gene_stable_id
#
my $g_stable_id = 'ENSG00000125964';
$slice = $slice_adaptor->fetch_by_gene_stable_id($g_stable_id);
$new_slice = $slice_adaptor->fetch_by_gene_stable_id($g_stable_id, $FLANKING);
ok($new_slice->start == $slice->start - $FLANKING);
ok($new_slice->end   == $slice->end   + $FLANKING);

#verify we can retrieve the gene from this slice
my $gene_found = 0;
foreach my $g (@{$slice->get_all_Genes}) {
  if($g->stable_id eq $g->stable_id) {
    $gene_found = 1;
    last;
  }
}
ok($gene_found);

# same test for flanking slice
$gene_found = 0;
foreach my $g (@{$new_slice->get_all_Genes}) {
  if($g->stable_id eq $g->stable_id) {
    $gene_found = 1;
    last;
  }
}
ok($gene_found);

#endif


//
//  fetch_by_region (entire region)
//
  slice = SliceAdaptor_fetchByRegion(sa, "chromosome", CHR, POS_UNDEF, POS_UNDEF, STRAND_UNDEF, NULL, 0);
  ok(nTest++, !strcmp(Slice_getSeqRegionName(slice), CHR));
  ok(nTest++, Slice_getStart(slice) == 1);

#ifdef DONE

#
# fetch_by_misc_feature_attribute
#
my $flanking= 1000;
$slice = $slice_adaptor->fetch_by_misc_feature_attribute('superctg',
                                                         'NT_030871',
                                                         $flanking);

ok($slice->seq_region_name eq '20');
ok($slice->start == 59707812 - $flanking);
ok($slice->end   == 60855021 + $flanking);

#endif

//
// normalized projected slice
//

//
// a slice with a PAR region
// 24,25
//
  slice = SliceAdaptor_fetchByRegion(sa, "chromosome", "Y", 9000000, 11000000, 1, NULL, 0);

  Vector *results = SliceAdaptor_fetchNormalizedSliceProjection(sa, slice, 0 );

  fprintf(stderr,"Pseudo autosomal region results\n");
  int i;
  for (i=0; i<Vector_getNumElement(results); i++) {
    ProjectionSegment *projection = Vector_getElementAt(results, i);
    fprintf(stderr, "Start: %ld\n",  ProjectionSegment_getFromStart(projection) );
    fprintf(stderr, "End: %ld\n",    ProjectionSegment_getFromEnd(projection) );
    fprintf(stderr, "To Slice: %s\n",Slice_getName(ProjectionSegment_getToSlice(projection)) );
    fprintf(stderr, "-----------\n" );
  }

  ok(nTest++,  Vector_getNumElement(results) == 3 );
  
  if (Vector_getNumElement(results) > 1) {
    ProjectionSegment *par1 = Vector_getElementAt(results, 1);
    ok(nTest++,  !strcmp(Slice_getSeqRegionName(ProjectionSegment_getToSlice(par1)), "20") );
  } else {
    ok(nTest++,  0);
  }

//
// a slice with a haplotype 
// 26,27
//

  slice = SliceAdaptor_fetchByRegion(sa, "chromosome", "20_HAP1", 30000000, 31000000, 1, NULL, 0);

  results = SliceAdaptor_fetchNormalizedSliceProjection(sa, slice, 0 );

  fprintf(stderr,"Haplotype projection results\n");
  for (i=0; i<Vector_getNumElement(results); i++) {
    ProjectionSegment *projection = Vector_getElementAt(results, i);
    fprintf(stderr, "Start: %ld\n",  ProjectionSegment_getFromStart(projection) );
    fprintf(stderr, "End: %ld\n",    ProjectionSegment_getFromEnd(projection) );
    fprintf(stderr, "To Slice: %s\n",Slice_getName(ProjectionSegment_getToSlice(projection)) );
    fprintf(stderr, "-----------\n" );
  }

  ok(nTest++,  Vector_getNumElement(results) == 3 );
  char *correctHapNames[] = {"20","20_HAP1","20"};
  for (i=0; i<Vector_getNumElement(results); i++) {
    ProjectionSegment *projection = Vector_getElementAt(results, i);
    ok(nTest++,  !strcmp(Slice_getSeqRegionName(ProjectionSegment_getToSlice(projection)), correctHapNames[i]) );
  }



//try a projection from chromosome 20 to supercontigs
  slice = SliceAdaptor_fetchByRegion(sa, "chromosome", "20", 29252000, 31252001, 1, NULL, 0);

  fprintf(stderr,"Projection from chromosome 20 to supercontig\n");
  results = Slice_project(slice, "supercontig", NULL);
  ok(nTest++, Vector_getNumElement(results) == 1);
  ProjectionSegment *p0 = Vector_getElementAt(results, 0);
  ok(nTest++,  !strcmp(Slice_getSeqRegionName(ProjectionSegment_getToSlice(p0)), "NT_028392") );

  for (i=0; i<Vector_getNumElement(results); i++) {
    ProjectionSegment *projection = Vector_getElementAt(results, i);
    fprintf(stderr, "Start: %ld\n",  ProjectionSegment_getFromStart(projection) );
    fprintf(stderr, "End: %ld\n",    ProjectionSegment_getFromEnd(projection) );
    fprintf(stderr, "To Slice: %s\n",Slice_getName(ProjectionSegment_getToSlice(projection)) );
    fprintf(stderr, "-----------\n" );
  }

//try a projection from clone to supercontig
  slice = SliceAdaptor_fetchByRegion(sa, "clone", "AL121583.25", POS_UNDEF,POS_UNDEF,STRAND_UNDEF, NULL, 0);

  fprintf(stderr, "Projection from clone AL121583.25 to supercontig\n");

  results = Slice_project(slice,"supercontig", NULL);
  ok(nTest++, Vector_getNumElement(results) == 1);
  p0 = Vector_getElementAt(results, 0);
  ok(nTest++,  !strcmp(Slice_getSeqRegionName(ProjectionSegment_getToSlice(p0)), "NT_028392") );
  for (i=0; i<Vector_getNumElement(results); i++) {
    ProjectionSegment *projection = Vector_getElementAt(results, i);
    fprintf(stderr, "Start: %ld\n",  ProjectionSegment_getFromStart(projection) );
    fprintf(stderr, "End: %ld\n",    ProjectionSegment_getFromEnd(projection) );
    fprintf(stderr, "To Slice: %s\n",Slice_getName(ProjectionSegment_getToSlice(projection)) );
    fprintf(stderr, "-----------\n" );
  }

#ifdef DONE
#
# test storing a couple of different slices
#
my $csa = $db->get_CoordSystemAdaptor();
my $ctg_cs  = $csa->fetch_by_name('contig');

$multi->save('core', 'seq_region', 'dna', 'dnac');

my $ctg_len = 50;
my $name = 'testregion';

#
# Store a slice with sequence
#

my $ctg_slice = Bio::EnsEMBL::Slice->new(-COORD_SYSTEM    => $ctg_cs,
                                         -SEQ_REGION_NAME => $name,
                                         -SEQ_REGION_LENGTH => $ctg_len,
                                         -START           => 1,
                                         -END             => $ctg_len,
                                         -STRAND          => 1); 

my $seq   = 'A' x $ctg_len;



$slice_adaptor->store($ctg_slice, \$seq);

$ctg_slice = $slice_adaptor->fetch_by_region('contig', $name);

ok($ctg_slice->length == $ctg_len);
ok($ctg_slice->seq eq $seq);
ok($ctg_slice->seq_region_name eq $name);

#
# Store a slice without sequence
#

my $chr_cs  = $csa->fetch_by_name('chromosome');

my $chr_len = 50e6;
$name = 'testregion2';
my $chr_slice = Bio::EnsEMBL::Slice->new(-COORD_SYSTEM    => $chr_cs,
                                         -SEQ_REGION_NAME => $name,
                                         -SEQ_REGION_LENGTH => $chr_len,
                                         -START           => 1,
                                         -END             => $chr_len,
                                         -STRAND          => 1); 

$slice_adaptor->store($chr_slice);

$chr_slice = $slice_adaptor->fetch_by_region('chromosome', $name);
ok($chr_slice->length() == $chr_len);
ok($chr_slice->seq_region_length() == $chr_len);
ok($chr_slice->seq_region_name eq $name);

#
# Store an assembly between the slices
#
my $asm_start = 9999;
my $asm_slice = $chr_slice->sub_Slice( $asm_start, $asm_start + $ctg_len - 1 );
my $str = $slice_adaptor->store_assembly( $asm_slice, $ctg_slice );

ok( $str eq "chromosome:NCBI33:testregion2:9999:10048:1<>".
            "contig::testregion:1:50:1" );

my $ctg_map = $chr_slice->project( $ctg_cs->name, $ctg_cs->version );
# Test currently fails as assembly cached somewhere.
#ok( @$ctg_map == 1 and
#    $ctg_map->[0]->[0] == $asm_slice->start and
#    $ctg_map->[0]->[1] == $asm_slice->end and
#    $ctg_map->[0]->[2]->name eq $ctg_slice->name );

my $chr_map = $ctg_slice->project( $chr_cs->name, $chr_cs->version );
# Test currently fails as assembly cached somewhere.
#ok( @$chr_map == 1 and
#    $chr_map->[0]->[0] == $ctg_slice->start and
#    $chr_map->[0]->[1] == $ctg_slice->end and
#    $chr_map->[0]->[2]->name eq $chr_slice->name );


$multi->restore('core', 'seq_region', 'dna', 'dnac');
#endif


//
// There was a bug such that features were not being retrieved
// from slices that had a start < 1.  This is a test for that case.
//
  slice = SliceAdaptor_fetchByRegion(sa, "chromosome", "20", 1, 35000000, STRAND_UNDEF, NULL, 0);
  fprintf(stderr,"debug: slice start = %ld\n", Slice_getStart(slice));
  fprintf(stderr,"debug: slice end = %ld\n", Slice_getEnd(slice));

#ifdef DONE
my $sfs1 = $slice->get_all_SimpleFeatures();
print_features($sfs1);
#endif

  slice = SliceAdaptor_fetchByRegion(sa, "chromosome", "20", -10, 35000000, STRAND_UNDEF, NULL, 0);
  fprintf(stderr,"debug: slice start = %ld\n", Slice_getStart(slice));
  fprintf(stderr,"debug: slice end = %ld\n", Slice_getEnd(slice));

#ifdef DONE
my $sfs2 = $slice->get_all_SimpleFeatures();
print_features($sfs2);

ok(@$sfs1 == @$sfs2);

#endif

//
// test fetch_by_name
//
  slice = SliceAdaptor_fetchByName(sa, Slice_getName(slice));

  ok(nTest++, !strcmp(CoordSystem_getName(Slice_getCoordSystem(slice)), "chromosome"));
  ok(nTest++, !strcmp(Slice_getSeqRegionName(slice), "20"));
  ok(nTest++, Slice_getStart(slice) == -10);
  ok(nTest++, Slice_getStrand(slice) == 1);
  ok(nTest++, Slice_getEnd(slice) == 35000000);

  slice = SliceAdaptor_fetchByName(sa, "clone::AL121583.25:1:10000:-1");

  ok(nTest++, !strcmp(CoordSystem_getName(Slice_getCoordSystem(slice)), "clone"));
  ok(nTest++, !strcmp(Slice_getSeqRegionName(slice), "AL121583.25"));
  ok(nTest++, Slice_getStart(slice) == 1);
  ok(nTest++, Slice_getStrand(slice) == -1);
  ok(nTest++, Slice_getEnd(slice) == 10000);

  slice = SliceAdaptor_fetchByName(sa, "clone::AL121583.25:1::");

  fprintf(stderr,"debug: slice name %s\n", Slice_getName(slice));
  ok(nTest++, !strcmp(CoordSystem_getName(Slice_getCoordSystem(slice)), "clone"));
  ok(nTest++, !strcmp(Slice_getSeqRegionName(slice), "AL121583.25"));
  ok(nTest++, Slice_getStart(slice) == 1);
  ok(nTest++, Slice_getStrand(slice) == 1);
  ok(nTest++, Slice_getEnd(slice) == 84710);


//
// test fetch_all
//

//default no duplicates and reference only
  Vector *slices = SliceAdaptor_fetchAll(sa, "chromosome", NULL, 0);
  printSlices(slices);

  ok(nTest++, Vector_getNumElement(slices) == 63); //, 'References slices for coord system chromosome');

// include duplicates
  slices = SliceAdaptor_fetchAll(sa, "chromosome", NULL, SA_INCLUDE_DUPLICATES);
  printSlices(slices);

  ok(nTest++, Vector_getNumElement(slices) == 62); //is(@$slices, 62, 'References slices for coord system chromosome when including duplicates (Y should become 1 region not 2)'); 

  slices = SliceAdaptor_fetchAll(sa, "contig", NULL, 0);

  ok(nTest++, Vector_getNumElement(slices) == 13);
  printSlices(slices);



  slices = SliceAdaptor_fetchAll(sa, "toplevel", NULL, 0);

  slice = Vector_getElementAt(slices,0);
  ok(nTest++, Vector_getNumElement(slices) == 1 && !strcmp(Slice_getSeqRegionName(slice), "20"));
  printSlices(slices);

//
// test the fuzzy matching of clone accessions
//
  char *cloneName = "AL031658";
  slice = SliceAdaptor_fetchByRegion(sa, "clone", cloneName, POS_UNDEF,POS_UNDEF,STRAND_UNDEF, NULL, 0);

  char *regionName = Slice_getSeqRegionName(slice);
  int lenCloneName = strlen(cloneName);

  fprintf(stderr,"debug: Fuzzy matched clone name %s Got %s\n",cloneName, Slice_getSeqRegionName(slice));

  int okFlag = 0;
  long ver;
  if (strstr(regionName, cloneName) == regionName && strlen(regionName) > strlen(cloneName)+2 && StrUtil_isLongInteger(&ver,&(regionName[lenCloneName+2]))) {
    okFlag = 1;
  }
  ok(nTest++, okFlag);

//make sure that it does not fuzzy match too much
  slice = SliceAdaptor_fetchByRegion(sa, "contig", cloneName, POS_UNDEF,POS_UNDEF,STRAND_UNDEF, NULL, 0);
  ok(nTest++, slice == NULL);
//print_slices([$slice]);

//make sure that you can fetch a seq_region without knowing its version
  slice = SliceAdaptor_fetchByRegion(sa, NULL, "20", POS_UNDEF,POS_UNDEF,STRAND_UNDEF, NULL, 0);
  ok(nTest++, slice!=NULL && !strcmp(Slice_getSeqRegionName(slice), "20"));

  slice = SliceAdaptor_fetchByRegion(sa, "toplevel", "20", POS_UNDEF,POS_UNDEF,STRAND_UNDEF, NULL, 0);
  ok(nTest++, slice!=NULL && !strcmp(Slice_getSeqRegionName(slice), "20"));

  slice = SliceAdaptor_fetchByRegion(sa, "toplevel", "20", 10, 20, STRAND_UNDEF, NULL, 0);
  ok(nTest++, slice!=NULL && Slice_getStart(slice) == 10 && Slice_getEnd(slice) == 20);

  slice = SliceAdaptor_fetchByRegion(sa, NULL, "20", 10, 20, 1,"NCBI33", 0);
  ok(nTest++, slice!=NULL && !strcmp(Slice_getSeqRegionName(slice), "20"));

  slice = SliceAdaptor_fetchByRegion(sa, NULL, "20", 10, 20, 1,"bogus", 0);
  ok(nTest++, slice == NULL);

  slice = SliceAdaptor_fetchByRegion(sa, "toplevel", "20", 10, 20, 1, "bogus", 0);
  ok(nTest++, slice!=NULL && !strcmp(Slice_getSeqRegionName(slice), "20"));

  // try fuzzy matching in conjunction with coord system guessing
  slice = SliceAdaptor_fetchByRegion(sa, NULL, cloneName, POS_UNDEF,POS_UNDEF,STRAND_UNDEF, NULL, 0);
  regionName = Slice_getSeqRegionName(slice);
  okFlag = 0;
  if (strstr(regionName, cloneName) == regionName && strlen(regionName) > strlen(cloneName)+2 && StrUtil_isLongInteger(&ver,&(regionName[lenCloneName+2]))) {
    okFlag = 1;
  }
  ok(nTest++, okFlag);

#ifdef DONE

# Testing synonym fetching
{
  my $syn_slice = $slice_adaptor->fetch_by_region(undef, 'anoth_20');
  is($syn_slice->seq_region_name(), '20', 'Ensuring slice is Chr20 as expected');
  my $chr_syn_slice = $slice_adaptor->fetch_by_region('chromosome', 'anoth_20');
  is($chr_syn_slice->seq_region_name(), '20', 'Ensuring slice is Chr20 as expected');
}

#{
#  my @slices = @{$slice_adaptor->fetch_all_by_synonym('anoth_20', 'UniGene')};
#  is(scalar(@slices), 0, 'Checking querying with a bad external name means no Slices');
#  
#  @slices = @{$slice_adaptor->fetch_all_by_synonym('anoth_20', 'RFAM')}; #Yeah ... RFAM
#  is(scalar(@slices), 1, 'Checking querying with a good external name means Slices');
#  is($slices[0]->seq_region_name(), '20', 'Ensuring slice is Chr20 as expected');
#}

# test that with multiple sequence regions with the same name, the
# highest (lowest-numbered) ranked comes out first
$multi->hide('core', 'seq_region');

my $sth = $db->dbc->prepare(qq{INSERT INTO seq_region (coord_system_id, name,
                                                  length)
                SELECT cs.coord_system_id, 'TESTREGION', 1000000
                FROM coord_system cs
                WHERE cs.name in ('supercontig', 'chromosome')});

$sth->execute();
$sth->finish();

$slice = $slice_adaptor->fetch_by_region('toplevel', 'TESTREGION');

ok($slice->seq_region_name() eq 'TESTREGION');
ok($slice->coord_system()->name() eq 'chromosome');


$multi->restore('core', 'seq_region');
#endif

//###### FETCH BY LOCATION
  testTopLevelLocation(sa, "1:1-1000", "chromosome",  "1", 1, 1000, 1);
// Silly  testTopLevelLocation(sa, "1:1-", "chromosome",  "1", 1, 246874334, 1);
// Silly  testTopLevelLocation(sa, "1:-10", "chromosome",  "1", 1, 10, 1);
  testTopLevelLocation(sa, "1:100", "chromosome",  "1", 100, 246874334, 1);
// Silly  testTopLevelLocation(sa, "1:", "chromosome",  "1", 1, 246874334, 1);
  testTopLevelLocation(sa, "1", "chromosome",  "1", 1, 246874334, 1);
  
  testTopLevelLocation(sa, "1:1..1000", "chromosome",  "1", 1, 1000, 1);
// Silly  testTopLevelLocation(sa, "1:1..", "chromosome",  "1", 1, 246874334, 1);
// Silly  testTopLevelLocation(sa, "1:..10", "chromosome",  "1", 1, 10, 1);
  testTopLevelLocation(sa, "1:100", "chromosome",  "1", 100, 246874334, 1);
// Silly  testTopLevelLocation(sa, "1:", "chromosome",  "1", 1, 246874334, 1);
  testTopLevelLocation(sa, "1", "chromosome",  "1", 1, 246874334, 1);
  
  testTopLevelLocation(sa, "1: 1-1,000", "chromosome",  "1", 1, 1000, 1);
  testTopLevelLocation(sa, "1: 1-1,000,000", "chromosome",  "1", 1, 1000000, 1);
  testTopLevelLocation(sa, "1: 1-1 000 000", "chromosome",  "1", 1, 1000000, 1);
  testTopLevelLocation(sa, "1: 1", "chromosome",  "1", 1, 246874334, 1);
// Silly  testTopLevelLocation(sa, "1: -10", "chromosome",  "1", 1, 10, 1);
  testTopLevelLocation(sa, "1: 100", "chromosome",  "1", 100, 246874334, 1);
  testTopLevelLocation(sa, "1:100..2_000_000_000", "chromosome",  "1", 100, 246874334, 1);
  testTopLevelLocation(sa, "1:100..2E9", "chromosome",  "1", 100, 246874334, 1);
  
  //#Try strands
  testTopLevelLocation(sa, "1:1-1000:1", "chromosome",  "1", 1, 1000, 1);
  testTopLevelLocation(sa, "1:1-1000:-1", "chromosome",  "1", 1, 1000, -1);
  testTopLevelLocation(sa, "1:1-1000:+", "chromosome",  "1", 1, 1000, 1);
  testTopLevelLocation(sa, "1:1-1000:-", "chromosome",  "1", 1, 1000, -1);
// Silly  testTopLevelLocation(sa, "1:1-1000..1", "chromosome",  "1", 1, 1000, 1);
// Silly  testTopLevelLocation(sa, "1:1-1000--1", "chromosome",  "1", 1, 1000, -1);

#ifdef DONE
dies_ok { $slice_adaptor->fetch_by_toplevel_location(); } 'Checking calling without a location fails';
dies_ok { $slice_adaptor->fetch_by_toplevel_location('', 1); } 'Checking calling with a blank location fails';
dies_ok { $slice_adaptor->fetch_by_toplevel_location('1:1_000_000_000..100', 1); } 'Checking calling with an excessive start throws an error';
ok(!defined $slice_adaptor->fetch_by_toplevel_location('wibble', 1), 'Checking with a bogus region returns undef');
ok(!defined $slice_adaptor->fetch_by_toplevel_location('1:-100--50', 1), 'Checking with a bogus region with negative coords returns undef');

# Try without toplevel_location
{
  my $location = 'AL359765.6.1.13780:2-100';
  
  note "Testing $location by asking for seqlevel";
  my $seqlevel_slice = $slice_adaptor->fetch_by_location($location, 'seqlevel');
  test_slice($location, $seqlevel_slice, 'contig', 'AL359765.6.1.13780', 2, 100, 1);
  
  note "Testing $location by asking for contig";
  my $contig_slice = $slice_adaptor->fetch_by_location($location, 'contig');
  test_slice($location, $contig_slice, 'contig', 'AL359765.6.1.13780', 2, 100, 1);
}

{
  #Non standard name check
  my ($name, $start, $end, $strand) = $slice_adaptor->parse_location_to_values('GL21446.1');
  is($name, 'GL21446.1', 'Name parses out');
  ok(!defined $start, 'Start is undefined');
  ok(!defined $end, 'End is undefined');
  ok(!defined $strand, 'Strand is undefined');
}

############# METHODS BELOW HERE 
#endif


#ifdef DONE
sub print_features {
  my $fs = shift;
  foreach my $f (@$fs) {
    my $start  = $f->start();
    my $end    = $f->end();
    my $strand = $f->strand();
    debug("  $start-$end($strand)");
  }
}

done_testing();

#endif


  return 0;
}
