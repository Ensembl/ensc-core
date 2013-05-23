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
    default:
      fprintf(stderr,"Error: Unknown feature type %d\n",type);
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
    default:
      fprintf(stderr,"Error: Unknown feature type %d\n",sf->objectType);
      exit(1);
  }

  return (SeqFeature *)feature;
}
