/*
 * Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
 * Copyright [2016-2017] EMBL-European Bioinformatics Institute
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

#include "ComparaDNAAlignFeatureAdaptor.h"
#include "GenomeDBAdaptor.h"
#include "GenomicAlignAdaptor.h"
#include "DNAFragAdaptor.h"
#include "DBAdaptor.h"
#include "Slice.h"
#include "MetaContainer.h"

int COMPARA_DAFA_CACHE_SIZE = 4;

ComparaDNAAlignFeatureAdaptor *ComparaDNAAlignFeatureAdaptor_new(ComparaDBAdaptor *dba) {
  ComparaDNAAlignFeatureAdaptor *cdafa = NULL;

  if ((cdafa = (ComparaDNAAlignFeatureAdaptor *)
            calloc(1,sizeof(ComparaDNAAlignFeatureAdaptor))) == NULL) {
    fprintf(stderr,"Error: Failed allocating cdafa\n");
    cdafa = NULL;
  } else {
    BaseComparaAdaptor_init((BaseComparaAdaptor *)cdafa, dba, COMPARADNAALIGNFEATURE_ADAPTOR);

    cdafa->regionCache = Cache_new(COMPARA_DAFA_CACHE_SIZE);
  }

  return cdafa;
}

Vector *ComparaDNAAlignFeatureAdaptor_fetchAllBySpeciesRegion(ComparaDNAAlignFeatureAdaptor *dafa,
                                                   char *csSpecies, char *csAssembly,
                                                   char *qySpecies, char *qyAssembly,
                                                   char *chrName, int start, int end,
                                                   char *alignmentType) {
  char *dnaFragType = "Chromosome";
  Vector *out;

  //get the genome database for each species
  GenomeDBAdaptor *gdba = ComparaDBAdaptor_getGenomeDBAdaptor(dafa->dba);
  GenomeDB *csGdb = GenomeDBAdaptor_fetchByNameAssembly(gdba, csSpecies, csAssembly);
  GenomeDB *qyGdb = GenomeDBAdaptor_fetchByNameAssembly(gdba, qySpecies, qyAssembly);

  //retrieve dna fragments from the subjects species region of interest
  DNAFragAdaptor *dfa = ComparaDBAdaptor_getDNAFragAdaptor(dafa->dba);
  Vector *dnaFrags = DNAFragAdaptor_fetchAllByGenomeDBRegion(dfa, 
                                                    csGdb,
						    dnaFragType,
						    chrName,
						    &start,
						    &end);

  GenomicAlignAdaptor *gaa = ComparaDBAdaptor_getGenomicAlignAdaptor(dafa->dba);
  int i;

  out = Vector_new();

  for (i=0; i<Vector_getNumElement(dnaFrags); i++) {
    // caclulate coords relative to start of dnafrag
    Vector *genomicAligns;
    DNAFrag *df = Vector_getElementAt(dnaFrags,i);
    Slice *slice = (Slice *)DNAFrag_getContig(df);
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
							     &dfStart,
							     &dfEnd,
							     alignmentType);

    // convert genomic aligns to dna align features
    for (j=0; j<Vector_getNumElement(genomicAligns); j++) {
      GenomicAlign *ga = Vector_getElementAt(genomicAligns, j);
      DNAAlignFeature *f = DNAAlignFeature_new();
      DNAAlignFeature_setCigarString(f, GenomicAlign_getCigarString(ga));

      DNAFrag *qdf = GenomicAlign_getQueryDNAFrag(ga);
      Slice *qSlice = (Slice *)DNAFrag_getContig(qdf);

      // calculate chromosomal coords
      int cStart = DNAFrag_getStart(df) + GenomicAlign_getConsensusStart(ga) - 1;
      int cEnd   = DNAFrag_getStart(df) + GenomicAlign_getConsensusEnd(ga) - 1;

      // skip features which do not overlap the requested region
      // next if ($cstart > $end || $cend < $start); 

      DNAAlignFeature_setSeqName((SeqFeature*)f, Slice_getChrName(slice));
      DNAAlignFeature_setStart(f, cStart);
      DNAAlignFeature_setEnd(f, cEnd);
      DNAAlignFeature_setStrand(f, 1);
      DNAAlignFeature_setSpecies((FeaturePair*)f, csSpecies);
      DNAAlignFeature_setScore(f, GenomicAlign_getScore(ga));
      DNAAlignFeature_setPercId(f, GenomicAlign_getPercentId(ga));

      DNAAlignFeature_setHitStart(f, DNAFrag_getStart(qdf) + GenomicAlign_getQueryStart(ga) - 1);
      DNAAlignFeature_setHitEnd(f, DNAFrag_getStart(qdf) + GenomicAlign_getQueryEnd(ga) - 1);
      DNAAlignFeature_setHitStrand(f,  GenomicAlign_getQueryStrand(ga));
      DNAAlignFeature_setHitSeqName(f, Slice_getChrName(qSlice));
      DNAAlignFeature_setHitSpecies((FeaturePair*)f, qySpecies);

      Vector_addElement(out, f);
    }
  }

// NIY Freeing genomic aligns
  return out;
}

Vector *ComparaDNAAlignFeatureAdaptor_fetchAllBySlice(ComparaDNAAlignFeatureAdaptor *dafa,
              Slice *slice, char *qySpecies, char *qyAssembly, char *assemblyType) {
  int sliceStart;
  int sliceEnd;
  int sliceStrand;
  char *csAssembly;
  char *csSpecies;
  Vector *features;
  char cacheKey[1024];
  void *val;
  MetaContainer *mc;

  if (!slice || slice->objectType!=CLASS_SLICE) {
    fprintf(stderr, "Error: Invalid slice argument\n");
    return NULL;
  }

  if (!qySpecies || !qyAssembly) {
    fprintf(stderr, "Error: Query species argument is required\n");
    return NULL;
  }

  mc = DBAdaptor_getMetaContainer(Slice_getAdaptor(slice)->dba);

  csSpecies = Species_getBinomialName(MetaContainer_getSpecies(mc));
  csAssembly = Slice_getAssemblyType(slice);

  sprintf(cacheKey,"%s:%s:%s:%s:%s:%s", 
            Slice_getName(slice),
            csSpecies,
            csAssembly,
            qySpecies,
            qyAssembly,
            assemblyType);
           
  if ((val = Cache_findElem(dafa->regionCache, cacheKey)) != NULL) {
    return (Vector *)val;
  }

  sliceStart  = Slice_getChrStart(slice);
  sliceEnd    = Slice_getChrEnd(slice);
  sliceStrand = Slice_getStrand(slice);

  features = ComparaDNAAlignFeatureAdaptor_fetchAllBySpeciesRegion(dafa, 
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
// Was setContig
      DNAAlignFeature_setSlice(f, slice);
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
// Was setContig
      DNAAlignFeature_setSlice(f, slice);
    }
  }

  // update the cache
  Cache_addElement(dafa->regionCache, cacheKey, features, NULL);

  return features;
}
