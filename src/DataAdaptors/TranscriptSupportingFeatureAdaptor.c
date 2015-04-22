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

/*
=head1 NAME

Bio::EnsEMBL::DBSQL::TranscriptSupportingFeatureAdaptor - Retrieves
supporting features from the database.
*/

#include "TranscriptSupportingFeatureAdaptor.h"

#include "DNAAlignFeatureAdaptor.h"
#include "ProteinAlignFeatureAdaptor.h"
#include "BaseAlignFeature.h"
#include "DNAAlignFeature.h"
#include "DNAPepAlignFeature.h"
#include "DBAdaptor.h"
#include "BaseFeatureAdaptor.h"
#include "SliceAdaptor.h"

#include <string.h>

TranscriptSupportingFeatureAdaptor *TranscriptSupportingFeatureAdaptor_new(DBAdaptor *dba) {
  TranscriptSupportingFeatureAdaptor *tsfa;

  if ((tsfa = (TranscriptSupportingFeatureAdaptor *)calloc(1,sizeof(TranscriptSupportingFeatureAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for TranscriptSupportingFeatureAdaptor\n");
    return NULL;
  }
  BaseAdaptor_init((BaseAdaptor *)tsfa, dba, TRANSCRIPTSUPPORTINGFEATURE_ADAPTOR);

  return tsfa;
}

/*
=head2 fetch_all_by_Transcript

  Arg [1]    : Bio::EnsEMBL::Transcript $transcript 
               The transcript to fetch supporting features for
  Example    : @sfs = @{$supporting_feat_adaptor->fetch_all_by_Transcript($transcript)};
  Description: Retrieves supporting features (evidence) for a given transcript. 
  Returntype : list of Bio::EnsEMBL::BaseAlignFeatures in the same coordinate
               system as the $transcript argument
  Exceptions : warning if $transcript is not in the database (i.e. dbID not defined)
               throw if a retrieved supporting feature is of unknown type 
  Caller     : Bio::EnsEMBL::Transcript
  Status     : Stable

=cut
*/
Vector *TranscriptSupportingFeatureAdaptor_fetchAllByTranscript(TranscriptSupportingFeatureAdaptor *tsfa, Transcript *transcript) {
  StatementHandle *sth;
  char qStr[512];
  ResultRow *row;
  DNAAlignFeatureAdaptor *dafa;
  ProteinAlignFeatureAdaptor *pafa;


  if (!Transcript_getDbID(transcript)) {
    fprintf(stderr,"WARNING: transcript has no dbID can't fetch evidence from db "
                   "no relationship exists\n");
    return Vector_new();
  }

  Vector *out = Vector_new();
  sprintf(qStr,"SELECT tsf.feature_type, tsf.feature_id "
               "FROM   transcript_supporting_feature tsf "
               "WHERE  transcript_id = " IDFMTSTR, Transcript_getDbID(transcript));

  sth = tsfa->prepare((BaseAdaptor *)tsfa, qStr, strlen(qStr));

  sth->execute(sth);

  pafa = DBAdaptor_getProteinAlignFeatureAdaptor(tsfa->dba);
  dafa = DBAdaptor_getDNAAlignFeatureAdaptor(tsfa->dba);

  while ((row = sth->fetchRow(sth))) {
    SeqFeature *sf = NULL;
    char *type = row->getStringAt(row,0);
    IDType dbId = row->getLongLongAt(row,1);

// sf is HACK HACK HACK
    if (!strcmp(type,"protein_align_feature")) {
      sf = (SeqFeature *)ProteinAlignFeatureAdaptor_fetchByDbID(pafa, dbId);
    } else if (!strcmp(type,"dna_align_feature")) {
      sf = (SeqFeature *)DNAAlignFeatureAdaptor_fetchByDbID(dafa, dbId);
    } else {
      fprintf(stderr,"Error: Unknown feature type [%s]\n",type);
      exit(1);
    }

    if (sf == NULL) {
      fprintf(stderr,"Warning: Transcript supporting feature %s "IDFMTSTR" does not exist in DB\n", type, dbId);
    } else {
      SeqFeature *newSf = SeqFeature_transfer(sf, (Slice *)Transcript_getSlice(transcript));
      // NIY: Free pretranferred one??
      if (newSf) {
        Vector_addElement(out, newSf);
      } else {
        fprintf(stderr,"Warning: Failed to transfer transcript supporting feature %s "IDFMTSTR" onto transcript slice\n", type, dbId);
      }
    }
  }

  sth->finish(sth);
  return out;
}


/*
=head2 store
  Arg [2]    : Int $transID
               The dbID of an EnsEMBL transcript to associate with supporting
               features
  Arg [1]    : Ref to array of Bio::EnsEMBL::BaseAlignFeature (the support)
  Example    : $dbea->store($transcript_id, \@features);
  Description: Stores a set of alignment features and associates an EnsEMBL transcript
               with them
  Returntype : none
  Exceptions : thrown when invalid dbID is passed to this method
  Caller     : TranscriptAdaptor
  Status     : Stable

=cut
*/
void TranscriptSupportingFeatureAdaptor_store(TranscriptSupportingFeatureAdaptor *tsfa, IDType transcriptDbID, Vector *alnObjs) {
  char pepCheckSql[1024];
  char dnaCheckSql[1024];
  char assocCheckSql[1024];
  char assocWriteSql[1024];

  if (alnObjs == NULL) return;

  sprintf(pepCheckSql,
      "SELECT protein_align_feature_id "
      "FROM protein_align_feature "
      "WHERE seq_region_id = %"IDFMTSTR
      " AND   seq_region_start = %%ld"
      " AND   seq_region_end   = %%ld"
      " AND   seq_region_strand = %%d"
      " AND   hit_name = '%%s'"
      " AND   hit_start = %%d"
      " AND   hit_end   = %%d"
      " AND   analysis_id = %"IDFMTSTR
      " AND   hcoverage = %%f"
      " AND   cigar_line = '%%s'");

  sprintf(dnaCheckSql,
      "SELECT dna_align_feature_id "
      "FROM dna_align_feature "
      "WHERE seq_region_id = %"IDFMTSTR
      " AND   seq_region_start = %%ld"
      " AND   seq_region_end   = %%ld"
      " AND   seq_region_strand = %%d"
      " AND   hit_name = '%%s'"
      " AND   hit_start = %%d"
      " AND   hit_end   = %%d"
      " AND   analysis_id = %"IDFMTSTR
      " AND   cigar_line = '%%s'"
      " AND   hcoverage = %%f"
      " AND   hit_strand = %%d");

  sprintf(assocCheckSql, 
      "SELECT * "
      "FROM  transcript_supporting_feature "
      "WHERE  transcript_id = "IDFMTSTR
      " AND   feature_type = '%%s'"
      " AND   feature_id   = %"IDFMTSTR, transcriptDbID);

  sprintf(assocWriteSql, 
          "INSERT into transcript_supporting_feature (transcript_id, feature_id, feature_type) values(%"IDFMTSTR", %"IDFMTSTR", '%%s')");

  StatementHandle *pepCheckSth   = tsfa->prepare((BaseAdaptor *)tsfa, pepCheckSql, strlen(pepCheckSql));
  StatementHandle *dnaCheckSth   = tsfa->prepare((BaseAdaptor *)tsfa, dnaCheckSql, strlen(dnaCheckSql));
  StatementHandle *assocCheckSth = tsfa->prepare((BaseAdaptor *)tsfa, assocCheckSql, strlen(assocCheckSql));
  StatementHandle *sfSth         = tsfa->prepare((BaseAdaptor *)tsfa, assocWriteSql, strlen(assocWriteSql));

  DNAAlignFeatureAdaptor *dnaAdaptor     = DBAdaptor_getDNAAlignFeatureAdaptor(tsfa->dba);
  ProteinAlignFeatureAdaptor *pepAdaptor = DBAdaptor_getProteinAlignFeatureAdaptor(tsfa->dba);
  SliceAdaptor *sliceAdaptor             = DBAdaptor_getSliceAdaptor(tsfa->dba);

  int i;
  for (i=0; i < Vector_getNumElement(alnObjs); i++) {
    BaseAlignFeature *f = Vector_getElementAt(alnObjs, i);

/* I really wish I understood why this is necessary - until I do, I'm not doing these transfers
    // check that the feature is in toplevel coords
    if($f->slice->start != 1 || $f->slice->strand != 1) {
    #move feature onto a slice of the entire seq_region
      my $tls = $self->db->get_sliceAdaptor->fetch_by_region($f->slice->coord_system->name(),
                                                             $f->slice->seq_region_name(),
                                                             undef, #start
                                                             undef, #end
                                                             undef, #strand
                                                             $f->slice->coord_system->version());
      $f = $f->transfer($tls);

      if(!$f) {
        throw('Could not transfer Feature to slice of ' .
              'entire seq_region prior to storing');
      }
    }
*/

    Class_assertType(CLASS_BASEALIGNFEATURE, f->objectType);

/*
    if(!$f->isa("Bio::EnsEMBL::BaseAlignFeature")){
      throw("$f must be an align feature otherwise" .
            "it can't be stored");
    }
*/
    
    IDType sfDbID;
    char *type;
    BaseFeatureAdaptor *adap;
    StatementHandle *checkSth;

    IDType seqRegionId = SliceAdaptor_getSeqRegionId(sliceAdaptor, BaseAlignFeature_getSlice(f));
    
// Note - moved the checkSth execute into the condition because I can't do the variable args
    if (f->objectType == CLASS_DNADNAALIGNFEATURE) {
      adap     = (BaseFeatureAdaptor*)dnaAdaptor;      
      type     = "dna_align_feature";

      checkSth = dnaCheckSth;
      checkSth->execute(checkSth, seqRegionId, 
                        BaseAlignFeature_getSeqRegionStart((SeqFeature*)f), 
                        BaseAlignFeature_getSeqRegionEnd((SeqFeature*)f), 
                        BaseAlignFeature_getSeqRegionStrand((SeqFeature*)f), 
                                  BaseAlignFeature_getHitSeqName(f), 
                                  BaseAlignFeature_getHitStart(f), 
                                  BaseAlignFeature_getHitEnd(f), 
                                  Analysis_getDbID(BaseAlignFeature_getAnalysis(f)), 
                                  BaseAlignFeature_getCigarString(f), 
                                  BaseAlignFeature_gethCoverage(f),
                                  DNAAlignFeature_getHitStrand((DNAAlignFeature *)f));

    } else if (f->objectType == CLASS_DNAPEPALIGNFEATURE) {
      adap     = (BaseFeatureAdaptor*)pepAdaptor;
      type     = "protein_align_feature";

      checkSth = pepCheckSth;
      checkSth->execute(checkSth, seqRegionId, 
                        BaseAlignFeature_getSeqRegionStart((SeqFeature*)f), 
                        BaseAlignFeature_getSeqRegionEnd((SeqFeature*)f), 
                        BaseAlignFeature_getSeqRegionStrand((SeqFeature*)f), 
                                  BaseAlignFeature_getHitSeqName(f), 
                                  BaseAlignFeature_getHitStart(f), 
                                  BaseAlignFeature_getHitEnd(f), 
                                  Analysis_getDbID(BaseAlignFeature_getAnalysis(f)), 
                                  BaseAlignFeature_getCigarString(f),
                                  BaseAlignFeature_gethCoverage(f));

    } else {
      fprintf(stderr, "Warning: Supporting feature of unknown type. Skipping\n");
      continue;
    }

/// HOW?? - moved into conditions above
//    checkSth->execute(@check_args);

    if (checkSth->numRows(checkSth) > 0) {
      ResultRow *row = checkSth->fetchRow(checkSth);
      sfDbID = row->getLongLongAt(row, 0);

    } else {
      Vector *vec = Vector_new();
      Vector_addElement(vec, f);
      //BaseFeatureAdaptor_store(adap, vec);
      adap->store((BaseAdaptor*)adap, vec);
      Vector_free(vec);
     
      sfDbID = BaseAlignFeature_getDbID(f);
    }

    // now check association
    assocCheckSth->execute(assocCheckSth, type, sfDbID);

    if (checkSth->numRows(assocCheckSth) == 0) {
      sfSth->execute(sfSth, transcriptDbID, sfDbID, type);
    }
  }

  dnaCheckSth->finish(dnaCheckSth);
  pepCheckSth->finish(pepCheckSth);
  assocCheckSth->finish(assocCheckSth);
  sfSth->finish(sfSth);
  
  return;
}
