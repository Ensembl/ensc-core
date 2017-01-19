/*
 * Copyright [1999-2017] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

#ifndef __INTRONSUPPORTINGEVIDENCE_H__
#define __INTRONSUPPORTINGEVIDENCE_H__

#include "DataModelTypes.h"
#include "SeqFeature.h"

#define INTRONSUPPORTINGEVIDENCEFUNCS_DATA(CLASSTYPE) \
  SEQFEATUREFUNCS_DATA(CLASSTYPE)

SEQFEATUREFUNC_TYPES(IntronSupportingEvidence)

typedef struct IntronSupportingEvidenceFuncsStruct {
  INTRONSUPPORTINGEVIDENCEFUNCS_DATA(IntronSupportingEvidence)
} IntronSupportingEvidenceFuncs;
  

#define INTRONSUPPORTINGEVIDENCE_DATA \
  SEQFEATURE_DATA \
  char      isSpliceCanonical; \
  char *    hitName; \
  ECOSTRING scoreType;


#define FUNCSTRUCTTYPE IntronSupportingEvidenceFuncs
struct IntronSupportingEvidenceStruct {
  INTRONSUPPORTINGEVIDENCE_DATA
};
#undef FUNCSTRUCTTYPE

IntronSupportingEvidence *IntronSupportingEvidence_new(void);
Intron *IntronSupportingEvidence_getIntron(IntronSupportingEvidence *ise, Transcript *transcript);
void IntronSupportingEvidence_setValuesFromIntron(IntronSupportingEvidence *ise, Intron *intron);
int IntronSupportingEvidence_hasLinkedTranscripts(IntronSupportingEvidence *ise);
Exon *IntronSupportingEvidence_findPreviousExon(IntronSupportingEvidence *ise, Transcript *transcript);
Exon *IntronSupportingEvidence_findNextExon(IntronSupportingEvidence *ise, Transcript *transcript);

#define IntronSupportingEvidence_isStored(ise, db) Storable_isStored(&((ise)->st), (db))

ECOSTRING IntronSupportingEvidence_setScoreType(IntronSupportingEvidence *ise, char *scoreType);
#define IntronSupportingEvidence_getScoreType(ise)  (ise)->scoreType

char *IntronSupportingEvidence_setHitName(IntronSupportingEvidence *ise, char *str);
#define IntronSupportingEvidence_getHitName(ise)  (ise)->hitName

#define IntronSupportingEvidence_setIsSpliceCanonical(ise,flag) (ise)->isSpliceCanonical = (flag)
#define IntronSupportingEvidence_getIsSpliceCanonical(ise)  (ise)->isSpliceCanonical

#define IntronSupportingEvidence_setStart(ise,start) SeqFeature_setStart((ise),(start))
#define IntronSupportingEvidence_getStart(ise) SeqFeature_getStart((ise))

#define IntronSupportingEvidence_setEnd(ise,end) SeqFeature_setEnd((ise),(end))
#define IntronSupportingEvidence_getEnd(ise) SeqFeature_getEnd((ise))

#define IntronSupportingEvidence_setStrand(ise,strand) SeqFeature_setStrand((ise),(strand))
#define IntronSupportingEvidence_getStrand(ise) SeqFeature_getStrand((ise))

#define IntronSupportingEvidence_setDbID(ise,dbID) SeqFeature_setDbID((ise),(dbID))
#define IntronSupportingEvidence_getDbID(ise) SeqFeature_getDbID((ise))

#define IntronSupportingEvidence_setAnalysis(ise,anal) SeqFeature_setAnalysis((ise),(anal))
#define IntronSupportingEvidence_getAnalysis(ise) SeqFeature_getAnalysis((ise))

#define IntronSupportingEvidence_setAdaptor(ise,adaptor) SeqFeature_setAdaptor((ise),(adaptor))
#define IntronSupportingEvidence_getAdaptor(ise) SeqFeature_getAdaptor((ise))

#define IntronSupportingEvidence_setSlice(ise,contig) SeqFeature_setSlice((ise),(contig))
#define IntronSupportingEvidence_getSlice(ise) SeqFeature_getSlice((ise))

#define IntronSupportingEvidence_setScore(ise,score) SeqFeature_setScore((ise),(score))
#define IntronSupportingEvidence_getScore(ise) SeqFeature_getScore((ise))

#define IntronSupportingEvidence_free(ise) SeqFeature_free((ise))

void IntronSupportingEvidence_freeImpl(IntronSupportingEvidence *ise);

#define IntronSupportingEvidence_getSeqRegionStart(ise) SeqFeature_getSeqRegionStart((SeqFeature *)(ise))
#define IntronSupportingEvidence_getSeqRegionEnd(ise) SeqFeature_getSeqRegionEnd((SeqFeature *)(ise))
#define IntronSupportingEvidence_getSeqRegionStrand(ise) SeqFeature_getSeqRegionStrand((SeqFeature *)(ise))

IntronSupportingEvidence *IntronSupportingEvidence_shallowCopyImpl(IntronSupportingEvidence *ise);
#define IntronSupportingEvidence_shallowCopy(ise) SeqFeature_shallowCopy((ise))

#ifdef __INTRONSUPPORTINGEVIDENCE_MAIN__
  IntronSupportingEvidenceFuncs 
    intronSupportingEvidenceFuncs = {
                        IntronSupportingEvidence_freeImpl, // free
                        IntronSupportingEvidence_shallowCopyImpl, // shallowCopy
                        NULL, // deepCopy
                        NULL, // getStart
                        NULL, // setStart
                        NULL, // getEnd
                        NULL  // setEnd
                       };
#else 
  extern IntronSupportingEvidenceFuncs intronSupportingEvidenceFuncs;
#endif

#endif
