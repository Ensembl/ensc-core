/*
 * Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

#include "SeqFeatureFactory.h"

#include "BaseAlignFeature.h"
#include "DNAAlignFeature.h"
#include "DNAPepAlignFeature.h"
#include "FeaturePair.h"
#include "PredictionTranscript.h"
#include "SeqFeature.h"
#include "Exon.h"
#include "PredictionExon.h"
#include "Transcript.h"
#include "SimpleFeature.h"
#include "Class.h"

SeqFeature *SeqFeatureFactory_newFeature(ClassType type) {
  void *feature;

  switch (type) {
    case CLASS_BASEALIGNFEATURE:
      feature = BaseAlignFeature_new(); 
      break;
    case CLASS_DNADNAALIGNFEATURE:
      feature = DNAAlignFeature_new(); 
      break;
    case CLASS_DNAPEPALIGNFEATURE:
      feature = DNAPepAlignFeature_new(); 
      break;
    case CLASS_FEATUREPAIR:
      feature = FeaturePair_new(); 
      break;
    case CLASS_PREDICTIONTRANSCRIPT:
      feature = PredictionTranscript_new(); 
      break;
    case CLASS_SIMPLEFEATURE:
      feature = SimpleFeature_new(); 
      break;
    case CLASS_SEQFEATURE:
      feature = SeqFeature_new(); 
      break;
    case CLASS_EXON:
      feature = Exon_new(); 
      break;
    case CLASS_TRANSCRIPT:
      feature = Transcript_new(); 
      break;
    case CLASS_PREDICTIONEXON:
      feature = PredictionExon_new(); 
      break;
    case CLASS_INTRONSUPPORTINGEVIDENCE:
      fprintf(stderr,"!!!!!! ise new in seqfeaturefactory\n");
      feature = IntronSupportingEvidence_new(); 
      break;
    default:
      fprintf(stderr,"Error: Unknown feature type %s SeqFeatureFactory_newFeature\n",Class_findByType(type)->name);
      exit(1);
  }

  return (SeqFeature *)feature;
}

SeqFeature *SeqFeatureFactory_newFeatureFromFeature(SeqFeature *sf) {
  void *feature;

  switch (sf->objectType) {
    case CLASS_BASEALIGNFEATURE:
      feature = BaseAlignFeature_new(); 
      break;
    case CLASS_DNADNAALIGNFEATURE:
      feature = DNAAlignFeature_shallowCopy((DNAAlignFeature *)sf); 
      break;
    case CLASS_DNAPEPALIGNFEATURE:
      feature = DNAPepAlignFeature_shallowCopy((DNAPepAlignFeature *)sf); 
      break;
    case CLASS_FEATUREPAIR:
      feature = FeaturePair_new(); 
      break;
    case CLASS_PREDICTIONTRANSCRIPT:
      feature = PredictionTranscript_new(); 
      break;
    case CLASS_SIMPLEFEATURE:
      feature = SimpleFeature_shallowCopy((SimpleFeature *)sf); 
      break;
    case CLASS_SEQFEATURE:
      feature = SeqFeature_new(); 
      break;
    case CLASS_EXON:
      feature = Exon_shallowCopy(sf); 
      break;
    case CLASS_TRANSCRIPT:
      feature = Transcript_shallowCopy(sf); 
      break;
    case CLASS_PREDICTIONEXON:
      feature = PredictionExon_shallowCopy(sf); 
      break;
    case CLASS_GENE:
      feature = Gene_shallowCopy(sf); 
      break;
    case CLASS_INTRONSUPPORTINGEVIDENCE:
      feature = IntronSupportingEvidence_shallowCopy(sf);
      break;
    default:
      fprintf(stderr,"Error: Unknown feature type %s in SeqFeatureFactory_newFeatureFromFeature\n",Class_findByType(sf->objectType)->name);
      exit(1);
  }

  return (SeqFeature *)feature;
}
