#ifndef __SEQFEATUREFACTORY_H__
#define __SEQFEATUREFACTORY_H__

#include "SeqFeature.h"

SeqFeature *SeqFeatureFactory_newFeature(ClassType type);
SeqFeature *SeqFeatureFactory_newFeatureFromFeature(SeqFeature *sf);

#endif
