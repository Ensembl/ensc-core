#include "DNAAlignFeatureAdaptor.h"

package Bio::EnsEMBL::Compara::DBSQL::DnaAlignFeatureAdaptor;
use strict;
use vars qw(@ISA);

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Cache; #CPAN LRU cache
use Bio::EnsEMBL::DnaDnaAlignFeature;

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

my $CACHE_SIZE = 4;

=head2 new

  Arg [1]    : list of args to super class constructor
  Example    : $dafa = new Bio::EnsEMBL::Compara::Genomi
  Description: Creates a new DnaAlignFeatureAdaptor.  The superclass 
               constructor is extended to initialise an internal cache.  This
               class should be instantiated through the get method on the 
               DBAdaptor rather than calling this method directory.
  Returntype : none
  Exceptions : none
  Caller     : Bio::EnsEMBL::DBSQL::DBConnection

=cut

sub DNAAlignFeatureAdaptor_new {
  my ($class, @args) = @_;

  my $self = $class->SUPER::new(@args);

  #initialize internal LRU cache
  tie(%{$self->{'_cache'}}, 'Bio::EnsEMBL::Utils::Cache', $CACHE_SIZE);

  return $self;
}




=head2 fetch_all_by_species_region

 Arg [1]    : string $cs_species
              e.g. "Homo sapiens"
 Arg [2]    : string $cs_assembly
              e.g. "NCBI_31"
 Arg [3]    : string $qy_species
              e.g. "Mus musculus"
 Arg [4]    : string $qy_assembly
              e.g. "MGSC_3"
 Arg [5]    : string $chr_name
              the name of the chromosome to retrieve alignments from (e.g. 'X')
 Arg [6]    : int start
 Arg [7]    : int end
 Arg [8]    : string $alignment_type
              The type of alignments to be retrieved
              e.g. WGA or WGA_HCR
 Example    : $gaa->fetch_all_by_species_region("Homo sapiens", "NCBI_31",
						"Mus musculus", "MGSC_3",
                                                "X", 250_000, 750_000,"WGA");
 Description: Retrieves alignments between the consensus and query species
              from a specified region of the consensus genome.
 Returntype : an array reference of Bio::EnsEMBL::DnaDnaAlignFeature objects
 Exceptions : none
 Caller     : general

=cut

sub DNAAlignFeatureAdaptor_fetchAllBySpeciesRegion(DNAAlignFeatureAdaptor *dafa,
                                                   char *csSpecies, char *csAssembly,
                                                   char *qySpecies, char *qyAssembly,
                                                   char *chrName, int start, int end,
                                                   char *alignmentType) {

  char *dnaFragType = "Chromosome";
  Vector *out;

  //get the genome database for each species
  GenomeDBAdaptor *gdba = $self->db->get_GenomeDBAdaptor;
  GenomeDB *csGdb = GenomeDBAdaptor_fetchByNameAssembly(gdba, csSpecies, csAssembly);
  GenomeDB *qyGdb = GenomeDBAdaptor_fetchByNameAssembly(gdba, qySpecies, qyAssembly);

  //retrieve dna fragments from the subjects species region of interest
  DNAFragAdaptor *dfa = $self->db->get_DnaFragAdaptor;
  Vector dnaFrags = DNAFragAdaptor_fetchAllByGenomeDBRegion(dfa, 
                                                    csGdb,
						    dnaFragType,
						    chrName,
						    start,
						    end);

  GenomicAlignAdaptor *gaa = $self->db->get_GenomicAlignAdaptor;
  int i;

  out = Vector_new();

  for (i=0; i<Vector_getNumElement(dnaFrags); i++) {
    // caclulate coords relative to start of dnafrag
    Vector *genomicAligns;
    DNAFrag df = Vector_getElementAt(dnaFrags,i);
    int dfStart = start - DNAFrag_getStart(df) + 1;
    int dfEnd   = end   - DNAFrag_getStart(df) + 1;

    // constrain coordinates so they are completely within the dna frag
    int len = DNAFrag_getEnd(df) - DNAFrag_getStart(df) + 1;
    int j;

    dfStart = (dfStart < 1)  ? 1 : dfStart;
    dfEnd   = (dfEnd > len) ? len : dfEnd;

    // fetch all alignments in the region we are interested in
    genomicAligns = GenomicAlignAdaptor_fetchAllByDNAFragGenomeDB(gaa,
                                                             df,
							     qyGdb,
							     dfStart,
							     dfEnd,
							     alignmentType);

    // convert genomic aligns to dna align features
    for (j=0; j<Vector_getNumElement(genomicAligns); j++) {
      GenomicAlign *ga = Vector_getElementAt(genomicAligns, j);
      DNAAlignFeature *f = DNAAlignFeature_new();
      DNAAlignFeature_setCigarString(GenomicAlign_getCigarString(ga));

      DNAFrag *qdf = GenomicAlign_getQueryDNAFrag(ga);

      // calculate chromosomal coords
      int cStart = DNAFrag_getStart(df) + GenomicAlign_getConsensusStart(ga) - 1;
      int cEnd   = DNAFrag_getStart(df) + GenomicAlign_getConsensusEnd(ga) - 1;

      // skip features which do not overlap the requested region
      // next if ($cstart > $end || $cend < $start); 

      DNAAlignFeature_setSeqName(f, ;
      $f->seqname($df->contig->chr_name);
      DNAAlignFeature_setStart(f, cStart);
      DNAAlignFeature_setEnd(f, cEnd);
      DNAAlignFeature_setStrand(f, 1);
      DNAAlignFeature_setSpecies(f, csSpecies);
      DNAAlignFeature_setScore(f, GenomicAlign_getScore(ga));
      DNAAlignFeature_setPercId(f, GenomicAlign_getPercId(ga));

      DNAAlignFeature_setHitStart(f, DNAFrag_getStart(qdf) + GenomicAlign_getQueryStart(ga) - 1);
      DNAAlignFeature_setHitEnd(f, DNAFrag_getStart(qdf) + GenomicAlign_getQueryEnd(ga) - 1);
      DNAAlignFeature_setHitStrand(f,  GenomicAlign_getQueryStrand(ga));
      DNAAlignFeature_setHitSeqName(f, ;
      $f->hseqname($qdf->contig->chr_name);
      DNAAlignFeature_setHitSpecies(f, qySpecies);

      Vector_addElement(out, f);
    }
  }

  return out;
}




=head2 fetch_all_by_Slice

 Arg [1]    : Bio::EnsEMBL::Slice
 Arg [2]    : string $qy_species
              The query species to retrieve alignments against
 Arg [3]    : string $qy_assembly
 Arg [4]    : string $$alignment_type
              The type of alignments to be retrieved
              e.g. WGA or WGA_HCR
 Example    : $gaa->fetch_all_by_Slice($slice, "Mus musculus","WGA");
 Description: find matches of query_species in the region of a slice of a 
              subject species
 Returntype : an array reference of Bio::EnsEMBL::DnaDnaAlignFeature objects
 Exceptions : none
 Caller     : general

=cut

sub DNAAlignFeatureAdaptor_fetchAllBySlice(DNAAlignFeatureAdaptor *dafa,
              Slice *slice, char *qySpecies, char *qyAssembly, char *assemblyType) {
  int sliceStart;
  int sliceEnd;
  int sliceStrand;
  char *csAssembly;
  Vector *features;

  if (!slice || slice->objectType==CLASS_SLICE) {
    fprintf(stderr, "Error: Invalid slice argument\n");
    exit(1);
  }

  if (!qySpecies || !qyAssembly) {
    fprintf(stderr, "Error: Query species argument is required\n");
    exit(1);
  }

  csSpecies =
      $slice->adaptor->db->get_MetaContainer->get_Species->binomial;
  csAssembly = Slice_getAssemblyType(slice);

  my $key = uc(join(':', "SLICE", $slice->name,
		 $cs_species,$cs_assembly,
		 $qy_species, $qy_assembly,$assembly_type));

  if (exists $self->{'_cache'}->{$key}) {
    return $self->{'_cache'}->{$key};
  }

  sliceStart  = Slice_getChrStart(slice);
  sliceEnd    = Slice_getChrEnd(slice);
  sliceStrand = Slice_getStrand(slice);

  features = DNAAlignFeatureAdaptor_fetchAllBySpeciesRegion(dafa, 
                                                    csSpecies, csAssembly,
						    qySpecies,qyAssembly,
						    Slice_getChrName(slice),
						    sliceStart, sliceEnd, assemblyType);

  if (sliceStrand == 1) {
    int i;
    for (i=0; i<Vector_getNumElement(features); i++) {
      DNAAlignFeature *f = Vector_getElementAt(features,i);
      int start  = DNAAlignFeature_getStart(f) - sliceStart + 1;
      int end    = DNAAlignFeature_getEnd(f)   - sliceStart + 1;

      DNAAlignFeature_setStart(f, start);
      DNAAlignFeature_setEnd(f, end);
      DNAAlignFeature_setContig(f, slice);
    }
  } else {
    int i;
    for (i=0; i<Vector_getNumElement(features); i++) {
      DNAAlignFeature *f = Vector_getElementAt(features,i);
      int start  = sliceEnd - DNAAlignFeature_getStart(f) + 1;
      int end    = sliceEnd - DNAAlignFeature_getEnd(f)   + 1;
      int strand = DNAAlignFeature_getStrand(f) * -1;

      DNAAlignFeature_setStart(f, start);
      DNAAlignFeature_setEnd(f, end);
      DNAAlignFeature_setStrand(f, strand);
      DNAAlignFeature_setContig(f, slice);
    }
  }

  // update the cache
  $self->{'_cache'}->{$key} = $features;

  return features;
}
