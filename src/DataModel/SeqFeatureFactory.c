#include "SeqFeatureFactory.h"

#include "BaseAlignFeature.h"
#include "DNAAlignFeature.h"
#include "DNAPepAlignFeature.h"
#include "FeaturePair.h"
#include "PredictionTranscript.h"
#include "SeqFeature.h"

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
    case CLASS_SEQFEATURE:
      feature = SeqFeature_new(); 
      break;
    default:
      fprintf(stderr,"Error: Unknown feature type %d\n",type);
      exit(1);
  }

  return (SeqFeature *)feature;
}
