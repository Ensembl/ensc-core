/*
 * Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

#include "GenomicAlignAdaptor.h"
#include "StrUtil.h"
#include "CigarStrUtil.h"
#include "DNAFragAdaptor.h"
#include "IDHash.h"

void GenomicAlignAdaptor_nextCig(GenomicAlignAdaptor *gaa,
    Vector *cigList, int *cigListPos, int *cs, int *ce, int *qs, int *qe);

#define DEFAULT_MAX_ALIGNMENT 20000

GenomicAlignAdaptor *GenomicAlignAdaptor_new(ComparaDBAdaptor *dba) {
  GenomicAlignAdaptor *gaa;
  Vector *vals;
  int maxAlLen;
  MetaContainer *mc;

  if ((gaa = (GenomicAlignAdaptor *)calloc(1,sizeof(GenomicAlignAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for GenomicAlignAdaptor\n");
    return NULL;
  }
  BaseComparaAdaptor_init((BaseComparaAdaptor *)gaa, dba, GENOMICALIGN_ADAPTOR);


  mc = ComparaDBAdaptor_getMetaContainer(gaa->dba);

  if (MetaContainer_getIntValueByKey(mc, "max_alignment_length", &maxAlLen)) {
    gaa->maxAlignmentLength = maxAlLen;
  } else {
    fprintf(stderr, "Warning: Meta table key 'max_alignment_length' not defined\n" 
	       "using default value [%d]", DEFAULT_MAX_ALIGNMENT);
    gaa->maxAlignmentLength = DEFAULT_MAX_ALIGNMENT;
  }

  return gaa;
}

void GenomicAlignAdaptor_store(GenomicAlignAdaptor *gaa, Vector *genomicAligns) {
  char *qStr;
  StatementHandle *sth;
  char commaStr[2] = {'\0','\0'};
  int i;
  char tmpStr[65556];
  

  StrUtil_copyString(&qStr, "INSERT INTO genomic_align_block"
             " (consensus_dnafrag_id, consensus_start, consensus_end,"
             "  query_dnafrag_id, query_start, query_end, query_strand, method_link_id,"
             "  score, perc_id, cigar_line) VALUES ",0);
  
  for (i=0; i<Vector_getNumElement(genomicAligns); i++) {
    GenomicAlign *ga = Vector_getElementAt(genomicAligns,i);
    DNAFrag *consDNAFrag  = GenomicAlign_getConsensusDNAFrag(ga);
    DNAFrag *queryDNAFrag = GenomicAlign_getQueryDNAFrag(ga);

    // check that everything has dbIDs
    if (!DNAFrag_getDbID(consDNAFrag) || !DNAFrag_getDbID(queryDNAFrag)) {
      fprintf(stderr, "Error: dna_fragment in GenomicAlign is not in DB\n");
      exit(1);
    }
  }

  // all clear for storing
  
  for (i=0; i<Vector_getNumElement(genomicAligns); i++) {
    GenomicAlign *ga = Vector_getElementAt(genomicAligns,i);
    DNAFrag *consDNAFrag  = GenomicAlign_getConsensusDNAFrag(ga);
    DNAFrag *queryDNAFrag = GenomicAlign_getQueryDNAFrag(ga);

    IDType methodLinkId = GenomicAlignAdaptor_methodLinkIdByAlignmentType(gaa, GenomicAlign_getAlignmentType(ga));

    if (!methodLinkId) {
      fprintf(stderr, "Error: There is no method_link with this type [%s] in the DB.\n",
              GenomicAlign_getAlignmentType(ga));
      exit(1);
    }
    
    sprintf(tmpStr," %s(" IDFMTSTR ", %d, %d, " IDFMTSTR ", %d, %d, %d, " IDFMTSTR ", %f, %f, '%s')", 
            commaStr, 
            DNAFrag_getDbID(consDNAFrag),
            GenomicAlign_getConsensusStart(ga),
            GenomicAlign_getConsensusEnd(ga),
            DNAFrag_getDbID(queryDNAFrag),  
            GenomicAlign_getQueryStart(ga),
            GenomicAlign_getQueryEnd(ga),
            GenomicAlign_getQueryStrand(ga),
            methodLinkId,
            GenomicAlign_getScore(ga),
            GenomicAlign_getPercentId(ga),
            GenomicAlign_getCigarString(ga));

    qStr = StrUtil_appendString(qStr, tmpStr);
    commaStr[0] = ','; 
  }
  
  sth = gaa->prepare((BaseAdaptor *)gaa, qStr, strlen(qStr));
  sth->execute(sth);
  sth->finish(sth);

  free(qStr);
}
     
/*
=head2 _fetch_all_by_DnaFrag_GenomeDB_direct

  Arg  1     : Bio::EnsEMBL::Compara::DnaFrag $dnafrag
               All genomic aligns that align to this frag
  Arg [2]    : Bio::EnsEMBL::Compara::GenomeDB $target_genome
               optionally restrict resutls to matches with this
               genome. Has to have a dbID().
  Arg [3]    : int $start
  Arg [4]    : int $end
  Arg [5]    : int $method_link_id
  Example    : none
  Description: Find all GenomicAligns that overlap this dnafrag.
               Return them in a way that this frags are on the
               consensus side of the Alignment.
  Returntype : listref Bio::EnsEMBL::Compara:GenomicAlign
  Exceptions : none
  Caller     : general

=cut
*/


Vector *GenomicAlignAdaptor_fetchAllByDNAFragGenomeDBDirect( GenomicAlignAdaptor *gaa, 
     DNAFrag *dnaFrag, GenomeDB *targetGenome, int *startP, int *endP, IDType methodLinkId) {
  IDType dnaFragId;
  GenomeDB *genomeDB;
  char *qStr;
  char tmpStr[512];
  Vector *results;
  StatementHandle *sth;

  if (!dnaFrag) {
    fprintf(stderr, "Error: Input dnafrag must not be NULL\n");
    exit(1);
  }

  // formatting the dnafrag
  dnaFragId = DNAFrag_getDbID(dnaFrag);

  genomeDB = DNAFrag_getGenomeDB(dnaFrag);

  StrUtil_copyString(&qStr,
     "SELECT gab.consensus_dnafrag_id,"
     "       gab.consensus_start," 
     "       gab.consensus_end,"
     "       gab.query_dnafrag_id," 
     "       gab.query_start," 
     "       gab.query_end,"
     "       gab.query_strand,"
     "       gab.method_link_id,"
     "       gab.score,"
     "       gab.perc_id," 
     "       gab.cigar_line"
     " FROM genomic_align_block gab ",0);

  if (targetGenome) {
    qStr = StrUtil_appendString(qStr,", dnafrag d");
  }
  sprintf(tmpStr," WHERE gab.method_link_id = " IDFMTSTR, methodLinkId);
  qStr = StrUtil_appendString(qStr,tmpStr);

  results = Vector_new();

  if (!targetGenome ||
      GenomeDB_hasQuery(genomeDB, targetGenome, methodLinkId)) {
    Vector *qres;

    sprintf(tmpStr," AND gab.consensus_dnafrag_id = " IDFMTSTR, dnaFragId);
    qStr = StrUtil_appendString(qStr, tmpStr);

    if (startP && endP) {
      int lowerBound = *startP - gaa->maxAlignmentLength;
      sprintf(tmpStr,
               " AND gab.consensus_start <= %d"
               " AND gab.consensus_start >= %d"
               " AND gab.consensus_end >= %d", *endP, lowerBound, *startP ) ;
      qStr = StrUtil_appendString(qStr, tmpStr);
    }

    if (targetGenome) {
      sprintf(tmpStr,
              " AND gab.query_dnafrag_id = d.dnafrag_id"
              " AND d.genome_db_id = " IDFMTSTR, GenomeDB_getDbID(targetGenome));
      qStr = StrUtil_appendString(qStr, tmpStr);
    }

    sth = gaa->prepare((BaseAdaptor *)gaa, qStr, strlen(qStr));
    sth->execute(sth);

    qres = GenomicAlignAdaptor_objectsFromStatementHandle(gaa, sth, 0);
    Vector_append(results,qres);
    Vector_free(qres);

    sth->finish(sth);
  }

  if (!targetGenome ||
      GenomeDB_hasConsensus(genomeDB, targetGenome, methodLinkId)) {
    Vector *cres;

    sprintf(tmpStr," AND gab.query_dnafrag_id = " IDFMTSTR, dnaFragId);
    qStr = StrUtil_appendString(qStr, tmpStr);

    if (startP && endP) {
      int lowerBound = *startP - gaa->maxAlignmentLength;
      sprintf(tmpStr,
               " AND gab.query_start <= %d"
               " AND gab.query_start >= %d"
               " AND gab.query_end >= %d", *endP, lowerBound, *startP ) ;
      qStr = StrUtil_appendString(qStr, tmpStr);
    }
    if (targetGenome) {
      sprintf(tmpStr,
               " AND gab.consensus_dnafrag_id = d.dnafrag_id"
               " AND d.genome_db_id = " IDFMTSTR, GenomeDB_getDbID(targetGenome));
      qStr = StrUtil_appendString(qStr, tmpStr);
    }
    sth = gaa->prepare((BaseAdaptor *)gaa, qStr, strlen(qStr));
    sth->execute(sth);

    cres = GenomicAlignAdaptor_objectsFromStatementHandle(gaa, sth, 1);
    Vector_append(results,cres);
    Vector_free(cres);

    sth->finish(sth);
  }
  free(qStr);

  return results;
}

/*

=head2 fetch_all_by_DnaFrag_GenomeDB

  Arg  1     : Bio::EnsEMBL::Compara::DnaFrag $dnafrag
  Arg  2     : string $query_species
               The species where the caller wants alignments to
               his dnafrag.
  Arg [3]    : int $start
  Arg [4]    : int $end
  Arg [5]    : string $alignment_type
               The type of alignments to be retrieved
               i.e. WGA or WGA_HCR
  Example    :  ( optional )
  Description: testable description
  Returntype : none, txt, int, float, Bio::EnsEMBL::Example
  Exceptions : none
  Caller     : object::methodname or just methodname

=cut
*/

Vector *GenomicAlignAdaptor_fetchAllByDNAFragGenomeDB(GenomicAlignAdaptor *gaa,
               DNAFrag *dnaFrag, GenomeDB *targetGenome, int *startP, int *endP, 
               char *alignmentType) {

  GenomeDB *genomeCons;
  IDType methodLinkId;
  GenomeDB *genomeQuery;
  Vector *mergedAligns;

  if (!dnaFrag) {
    fprintf(stderr, "Error: dnaFrag argument must be non NULL\n");
    exit(1);
  }

  methodLinkId = GenomicAlignAdaptor_methodLinkIdByAlignmentType(gaa, alignmentType);

  genomeCons = DNAFrag_getGenomeDB(dnaFrag);
  genomeQuery = targetGenome;
  
  // direct or indirect ??
  if (GenomeDB_hasConsensus(genomeCons, genomeQuery, methodLinkId) ||
      GenomeDB_hasQuery(genomeCons, genomeQuery, methodLinkId)) {
    return GenomicAlignAdaptor_fetchAllByDNAFragGenomeDBDirect(gaa, 
                   dnaFrag, targetGenome, startP, endP, methodLinkId);
  } else {
    // indirect checks
    Vector *linkedCons  = GenomeDB_linkedGenomesByMethodLinkId(genomeCons, methodLinkId);
    Vector *linkedQuery = GenomeDB_linkedGenomesByMethodLinkId(genomeQuery, methodLinkId);
    
    // there are not many genomes, square effort is cheap
    Vector *linked = Vector_new();
    Vector *set1 = Vector_new();
    mergedAligns = Vector_new();
    int i;

    for (i=0; i<Vector_getNumElement(linkedCons); i++) {
      int j;
      GenomeDB *g1 = Vector_getElementAt(linkedCons, i);

      for (j=0; j<Vector_getNumElement(linkedQuery); j++) {
        GenomeDB *g2 = Vector_getElementAt(linkedQuery, i);
	if (g1 == g2) {
	  Vector_addElement(linked, g1);
	}
      }
    }
    Vector_free(linkedCons);
    Vector_free(linkedQuery);

    // collect GenomicAligns from all linked genomes
    for (i=0; i<Vector_getNumElement(linked); i++) {
      GenomeDB *g = Vector_getElementAt(linked, i);

      Vector *gres = GenomicAlignAdaptor_fetchAllByDNAFragGenomeDBDirect(gaa, 
                             dnaFrag, g, startP, endP, methodLinkId);
      Vector_append(set1, gres);

      Vector_free(gres);
    }

    // go from each dnafrag in the result set to target_genome
    // there is room for improvement here: create start end
    // my %frags = map { $_->query_dnafrag->dbID => $_->query_dnafrag } @$set1;
    

    for (i=0; i<Vector_getNumElement(set1); i++) {
      GenomicAlign *alignA = Vector_getElementAt(set1,i);
      DNAFrag *frag = GenomicAlign_getQueryDNAFrag(alignA);
      int qStart = GenomicAlign_getQueryStart(alignA);
      int qEnd   = GenomicAlign_getQueryEnd(alignA);

      Vector *dres = GenomicAlignAdaptor_fetchAllByDNAFragGenomeDBDirect(gaa,
                      frag, genomeQuery, &qStart, &qEnd, methodLinkId);
      int j;

      for (j=0; j<Vector_getNumElement(dres); j++) {
        GenomicAlign *alignB = Vector_getElementAt(dres,j);
	GenomicAlignAdaptor_addDerivedAlignments(gaa,  mergedAligns, alignA, alignB);
      } 
      Vector_free(dres);
    }
// NIY freeing
    return mergedAligns;
  }
}

/*
=head2 _merge_alignsets

  Arg  1     : listref Bio::EnsEMBL::Compara::GenomicAlign $set1
               from consensus to query
  Arg  2     : listref Bio::EnsEMBL::Compara::GenomicAlign $set2
               and over consensus to next species query             
  Example    : none
  Description: set1 contains GAs with consensus species belonging to
               the input dnafragment. Query fragments are the actual reference
               species. In set 2 consensus species is the reference and
               query is the actual target genome. There may be more than
               one reference genome involved.
  Returntype : listref Bio::EnsEMBL::Compara::GenomicAlign
  Exceptions : none
  Caller     : internal

=cut
*/

typedef struct GenomicAlignListElemStruct {
  IDType queryDbID;
  double queryStart;
  GenomicAlign *align;
  int setNum;
} GenomicAlignListElem;

GenomicAlignListElem *GenomicAlignListElem_new(IDType id, int start, GenomicAlign *ga, int setNum) {
  GenomicAlignListElem *gale;

  if ((gale = (GenomicAlignListElem *)calloc(1,sizeof(GenomicAlignListElem))) == NULL) {
    fprintf(stderr, "Error: Failed allocating gale\n");
    exit(1);
  }
  gale->queryDbID = id;
  gale->queryStart = start;
  gale->align = ga;
  gale->setNum = setNum;
 
  return gale;
}

int GenomicAlignListElem_compFunc(const void *a, const void *b) {
  GenomicAlignListElem **one = (GenomicAlignListElem **)a; 
  GenomicAlignListElem **two = (GenomicAlignListElem **)b; 

  if ((*one)->queryDbID == (*two)->queryDbID) {
    if ((*one)->queryStart < (*two)->queryStart) {
      return -1;
    } else if ((*one)->queryStart > (*two)->queryStart) {
      return 1;
    } else {
      return 0;
    }
  } else if ((*one)->queryDbID < (*two)->queryDbID) {
    return -1;
  } else {
    return 1;
  }
}
  // sorting of both sets
  // walking through and finding overlapping GAs
  // create GA from overlapping GA
  // return list of those

  // efficiently generating all Aligns that overlap
  // [ key, object, set1 or 2 ]
  // Alignments are twice in big list. They are added to the overlapping
  // set the first time they appear and they are removed the
  // second time they appear. Scanline algorithm

Vector *GenomicAlignAdaptor_mergeAlignsets(GenomicAlignAdaptor *gaa, Vector *alignSet1, Vector *alignSet2) {
  int i;
  Vector *bigList = Vector_new();
  IDHash *overlappingSets[2];
  Vector *mergedAligns;


  for (i=0;i<Vector_getNumElement(alignSet1); i++) {
    GenomicAlign *align = Vector_getElementAt(alignSet1, i);
    Vector_addElement(bigList, GenomicAlignListElem_new(DNAFrag_getDbID(GenomicAlign_getQueryDNAFrag(align)),
                                                        GenomicAlign_getQueryStart(align), align, 0));
    Vector_addElement(bigList, GenomicAlignListElem_new(DNAFrag_getDbID(GenomicAlign_getQueryDNAFrag(align)),
                                                        GenomicAlign_getQueryEnd(align)+0.5, align, 0));
  }

  for (i=0;i<Vector_getNumElement(alignSet2); i++) {
    GenomicAlign *align = Vector_getElementAt(alignSet2, i);
    Vector_addElement(bigList, GenomicAlignListElem_new(DNAFrag_getDbID(GenomicAlign_getConsensusDNAFrag(align)),
                                                        GenomicAlign_getConsensusStart(align), align, 1));
    Vector_addElement(bigList, GenomicAlignListElem_new(DNAFrag_getDbID(GenomicAlign_getConsensusDNAFrag(align)),
                                                        GenomicAlign_getConsensusEnd(align)+0.5, align, 1));
  }
  
  Vector_sort(bigList, GenomicAlignListElem_compFunc);

  // walking from start to end through sortlist and keep track of the 
  // currently overlapping set of Alignments
 
  overlappingSets[0] = IDHash_new(IDHASH_SMALL);
  overlappingSets[1] = IDHash_new(IDHASH_SMALL);

  mergedAligns = Vector_new();

  for (i=0; i<Vector_getNumElement(bigList); i++) {
    GenomicAlignListElem *gale  = Vector_getElementAt(bigList,i);

    GenomicAlign *align = gale->align;
    int setNo           = gale->setNum;

    if (IDHash_contains(overlappingSets[setNo], align)) {
      // remove from current overlapping set
      IDHash_remove(overlappingSets[setNo], align, NULL);
    } else {
      int j;
      void **values = IDHash_getValues(overlappingSets[1-setNo]);

      // insert into the set and do all the overlap business
      IDHash_add(overlappingSets[setNo], align, align);

      // the other set contains everything this align overlaps with
      for (j=0; j<IDHash_getNumValues(overlappingSets[1-setNo]); j++) {
        GenomicAlign *align2 = values[j];
        if (setNo == 0) {
          GenomicAlignAdaptor_addDerivedAlignments(gaa, mergedAligns, align, align2);
        } else {
          GenomicAlignAdaptor_addDerivedAlignments(gaa, mergedAligns, align2, align);
        }
      }
      free(values);
    }
  }

// NIY Free gale

  return mergedAligns;
}

/*
=head2 _add_derived_alignments

  Arg  1     : listref 
    Additional description lines
    list, listref, hashref
  Example    :  ( optional )
  Description: testable description
  Returntype : none, txt, int, float, Bio::EnsEMBL::Example
  Exceptions : none
  Caller     : object::methodname or just methodname

=cut
*/

void GenomicAlignAdaptor_addDerivedAlignments(GenomicAlignAdaptor *gaa, 
                     Vector *mergedAligns, GenomicAlign *alignA, GenomicAlign *alignB) {

  // variable name explanation
  // q - query c - consensus s - start e - end l - last
  // o, ov overlap j - jump_in_
  // r - result

  int  qs, qe, lqs, lqe, cs, ce, lcs, lce,
       ocs, oce, oqs, oqe, jc, jq, ovs, ove,
       rcs, rce, rqs, rqe;
  int currentMatch = 0;
  int newMatch;
  int cigAPos = 0, cigBPos = 0;
  char *resultCig;
  char tmpStr[128];

  // initialization phase
  Vector *cigA = CigarStrUtil_getPieces(GenomicAlign_getCigarString(alignA));
  Vector *cigB = CigarStrUtil_getPieces(GenomicAlign_getCigarString(alignB));

  if (GenomicAlign_getQueryStrand(alignA) == -1 ) {
    Vector_reverse(cigB);
  }

  // need a 'normalized' start for qs, qe, oxs so I dont 
  // have to check strandedness all the time  

  // consensus is strand 1 and is not compared to anything,
  // can keep its original coordinate system
 
  lce = GenomicAlign_getConsensusStart(alignA) - 1;
  ce = lce;
  cs = ce + 1;
  
  // alignBs query can be + or - just keep relative coords for now
  lqe = 0; lqs = 1;
  qe = 0; qs = 1;

  // ocs will be found relative to oce and has to be comparable
  // to oqs. But it could be that we have to move downwards if we
  // are not - strand. thats why coordinates are transformed here

  if (GenomicAlign_getQueryStrand(alignA) == -1 ) {
    // query_end is first basepair of alignment
    if (GenomicAlign_getQueryEnd(alignA) < GenomicAlign_getConsensusEnd(alignB)) {
      oce = 0; ocs = 1;
      oqe = GenomicAlign_getConsensusEnd(alignB) - GenomicAlign_getQueryEnd(alignA);
      oqs = oqe + 1;
    } else {
      oqe = 0; oqs = 1;
      oce = GenomicAlign_getQueryEnd(alignA) - GenomicAlign_getConsensusEnd(alignB);
      ocs = oce + 1;
    }
  } else {
    // in theory no coordinate magic necessary :-)
    oqs = GenomicAlign_getQueryStart(alignA);
    oqe = oqs - 1; 
    ocs = GenomicAlign_getConsensusStart(alignB);
    oce = ocs - 1;
  }

  // initializing result
  rcs = rce = rqs = rqe = 0;
  resultCig= StrUtil_copyString(&resultCig,"",0);

  while (1) {
    int newGa;
    // exit if you request a new piece of alignment and the cig list is 
    // empty

    if (oce < ocs || oce < oqs) {
      // next M area in cigB
      if (cigBPos == Vector_getNumElement(cigB)) break;
      GenomicAlignAdaptor_nextCig(gaa, cigB, &cigBPos, &ocs, &oce, &qs, &qe ); 
      continue;
    }
    if (oqe < oqs || oqe < ocs) {
      // next M area in cigA
      if (cigAPos == Vector_getNumElement(cigA)) break;
      GenomicAlignAdaptor_nextCig(gaa, cigA, &cigAPos, &cs, &ce, &oqs, &oqe );
      continue;
    }

    // now matching region overlap in reference genome
    ovs = ocs < oqs ? oqs : ocs;
    ove = oce < oqe ? oce : oqe;
    
    if (currentMatch) {
      jc = cs + (ovs - oqs) - lce - 1;
      jq = qs + (ovs - ocs) - lqe - 1;
    } else {
      jc = jq = 0;
    }

    newMatch = ove - ovs + 1;
    newGa = 0;

    if (jc==0) {
      if (jq==0) {
	currentMatch += newMatch;
      } else {
        // store current match;
        sprintf(tmpStr,"%dM",currentMatch);
        resultCig = StrUtil_appendString(resultCig,tmpStr);

	// jq deletions;
	if (jq == 1) {
          resultCig = StrUtil_appendString(resultCig,"D");
        } else {
          sprintf(tmpStr,"%dD",jq);
          resultCig = StrUtil_appendString(resultCig,tmpStr);
        }
	currentMatch = newMatch;
      }
    } else {
      if (jq==0) {
        // store current match;
        sprintf(tmpStr,"%dM",currentMatch);
        resultCig = StrUtil_appendString(resultCig,tmpStr);

	// jc insertions;
	if (jc==1) {
          resultCig = StrUtil_appendString(resultCig,"I");
        } else {
          sprintf(tmpStr,"%dI",jc);
          resultCig = StrUtil_appendString(resultCig,tmpStr);
        }
	currentMatch = newMatch;
         
      } else {
        double percId;
        double score;
        GenomicAlign *ga;

        sprintf(tmpStr,"%dM",currentMatch);
        resultCig = StrUtil_appendString(resultCig,tmpStr);

	// new GA
	int queryStrand = GenomicAlign_getQueryStrand(alignA) * GenomicAlign_getQueryStrand(alignB);
	int queryStart, queryEnd;
	if (queryStrand == 1) {
	  queryStart = rqs + GenomicAlign_getQueryStart(alignB) - 1;
	  queryEnd = rqe + GenomicAlign_getQueryStart(alignB) - 1;
	} else {
	  queryEnd = GenomicAlign_getQueryEnd(alignB) - rqs + 1;
	  queryStart = GenomicAlign_getQueryEnd(alignB) - rqe + 1;
	}
      
        score = (GenomicAlign_getScore(alignA) < GenomicAlign_getScore(alignB)) ? 
          GenomicAlign_getScore(alignA) : GenomicAlign_getScore(alignB);
        percId =  (int)(GenomicAlign_getPercentId(alignA)*GenomicAlign_getPercentId(alignB)/100.0);
        
        ga = GenomicAlign_new();
    
        GenomicAlign_setConsensusDNAFrag(ga, GenomicAlign_getConsensusDNAFrag(alignA));
        GenomicAlign_setQueryDNAFrag(ga, GenomicAlign_getQueryDNAFrag(alignB));
        GenomicAlign_setCigarString(ga, resultCig);
        GenomicAlign_setConsensusStart(ga, rcs);
        GenomicAlign_setConsensusEnd(ga, rce);
        GenomicAlign_setQueryStrand(ga, queryStrand);
        GenomicAlign_setQueryStart(ga, queryStart);
        GenomicAlign_setQueryEnd(ga, queryEnd);
        GenomicAlign_setAdaptor(ga, (BaseAdaptor *)gaa);
        GenomicAlign_setPercentId(ga, percId);
        GenomicAlign_setScore(ga, score);

	Vector_addElement(mergedAligns, ga);

        rcs = rce = rqs = rqe = 0;
	resultCig[0] = '\0';
	
	currentMatch = newMatch;
      }
    }


    
    if (!rcs) rcs = cs+(ovs-oqs);
    rce = cs+(ove-oqs);
    if (!rqs) rqs = qs+(ovs-ocs);
    rqe = qs+(ove-ocs);

    // update the last positions
    lce = rce; 
    lqe = rqe;

    // next piece on the one that end earlier
 
    if (oce <= oqe) {
      // next M area in cigB
      if (cigBPos == Vector_getNumElement(cigB)) break;
      GenomicAlignAdaptor_nextCig(gaa, cigB, &cigBPos, &ocs, &oce, &qs, &qe ); 
    }
    if (oce >= oqe) {
      // next M area in cigA
      if (cigAPos == Vector_getNumElement(cigA)) break;
      GenomicAlignAdaptor_nextCig(gaa, cigA, &cigAPos, &cs, &ce, &oqs, &oqe );
    } 
  } // end of while loop

  // if there is a last floating current match
  if (currentMatch) {
    
    // new GA
    int queryStrand = GenomicAlign_getQueryStrand(alignA) * GenomicAlign_getQueryStrand(alignB);
    int queryStart, queryEnd;
    double percId;
    double score;
    GenomicAlign *ga;

    sprintf(tmpStr,"%dM",currentMatch);
    resultCig = StrUtil_appendString(resultCig, tmpStr);

    if (queryStrand == 1) {
      queryStart = rqs + GenomicAlign_getQueryStart(alignB) - 1;
      queryEnd = rqe + GenomicAlign_getQueryStart(alignB) - 1;
    } else {
      queryEnd = GenomicAlign_getQueryEnd(alignB) - rqs + 1;
      queryStart = GenomicAlign_getQueryEnd(alignB) - rqe + 1;
    }
  
    score = (GenomicAlign_getScore(alignA) < GenomicAlign_getScore(alignB)) ? 
      GenomicAlign_getScore(alignA) : GenomicAlign_getScore(alignB);
    percId =  (int)(GenomicAlign_getPercentId(alignA)*GenomicAlign_getPercentId(alignB)/100.0);
    
    ga = GenomicAlign_new();

    GenomicAlign_setConsensusDNAFrag(ga, GenomicAlign_getConsensusDNAFrag(alignA));
    GenomicAlign_setQueryDNAFrag(ga, GenomicAlign_getQueryDNAFrag(alignB));
    GenomicAlign_setCigarString(ga, resultCig);
    GenomicAlign_setConsensusStart(ga, rcs);
    GenomicAlign_setConsensusEnd(ga, rce);
    GenomicAlign_setQueryStrand(ga, queryStrand);
    GenomicAlign_setQueryStart(ga, queryStart);
    GenomicAlign_setQueryEnd(ga, queryEnd);
    GenomicAlign_setAdaptor(ga, (BaseAdaptor *)gaa);
    GenomicAlign_setPercentId(ga, percId);
    GenomicAlign_setScore(ga, score);

    Vector_addElement(mergedAligns, ga);
  }

  free(resultCig);

  Vector_free(cigA);
  Vector_free(cigB);

  // nothing to return all in merged_aligns
}


void GenomicAlignAdaptor_nextCig(GenomicAlignAdaptor *gaa,
    Vector *cigList, int *cigListPos, int *cs, int *ce, int *qs, int *qe)  {
  int count;
  char type;
  char *cigElem;
  int lenElem;
  
  do {
    cigElem = Vector_getElementAt(cigList, *cigListPos);
    (*cigListPos)++;
    lenElem = strlen(cigElem);
    type = cigElem[lenElem-1];

    
    if (type!='M' && type!='I' && type!='D') {
      fprintf(stderr,"Error: Cigar string format error for %s\n",cigElem);
      exit(1);
    }
  
    if (lenElem > 1) {
      cigElem[lenElem-1] = '\0';
      count = atol(cigElem);
    } else {
      count = 1;
    }


    switch (type) {
      case 'D':
        *qe += count;
        break;
      case 'I':
        *ce += count;
        break;
      case 'M':
        *cs = *ce + 1;
        *ce = *cs + count - 1;
        *qs = *qe + 1;
        *qe = *qs + count - 1;
    } 
  } while (type != 'M' && *cigListPos!=Vector_getNumElement(cigList));
}

Vector *GenomicAlignAdaptor_objectsFromStatementHandle(GenomicAlignAdaptor *gaa, StatementHandle *sth,
                                                       int reverse) {
  Vector *results = Vector_new();
  ResultRow *row;
  DNAFragAdaptor *dfa;

  IDType consensusDNAFragId;
  IDType queryDNAFragId;
  int consensusStart;
  int consensusEnd;
  int queryStart;
  int queryEnd;
  int queryStrand;
  IDType methodLinkId;
  double score;
  double percId;
  char *cigarString;

  dfa = ComparaDBAdaptor_getDNAFragAdaptor(gaa->dba);

  while ((row = sth->fetchRow(sth))) {
    GenomicAlign *genomicAlign;
    char *alignmentType;
    
    if (reverse) {
      queryDNAFragId     = row->getLongLongAt(row,0);
      queryStart         = row->getIntAt(row,1);
      queryEnd           = row->getIntAt(row,2);
      consensusDNAFragId = row->getLongLongAt(row,3);
      consensusStart     = row->getIntAt(row,4);
      consensusEnd       = row->getIntAt(row,5);
    } else {
      consensusDNAFragId = row->getLongLongAt(row,0);
      consensusStart     = row->getIntAt(row,1);
      consensusEnd       = row->getIntAt(row,2);
      queryDNAFragId     = row->getLongLongAt(row,3);
      queryStart         = row->getIntAt(row,4);
      queryEnd           = row->getIntAt(row,5);
    }
    queryStrand  = row->getIntAt(row,6);
    methodLinkId = row->getLongLongAt(row,7);
    score        = row->getDoubleAt(row,8);
    percId       = row->getDoubleAt(row,9);
    cigarString  = row->getStringAt(row,10);

    alignmentType = GenomicAlignAdaptor_alignmentTypeByMethodLinkId(gaa, methodLinkId);

    if (reverse) {
      StrUtil_strReplChrs(cigarString,"DI","ID");

      // alignment of the opposite strand
      if (queryStrand == -1) {
        cigarString = CigarStrUtil_reverse(cigarString, strlen(cigarString));
      }
    }
    
    
    genomicAlign = GenomicAlign_new();
    GenomicAlign_setAdaptor(genomicAlign, (BaseAdaptor *)gaa);
    GenomicAlign_setConsensusDNAFrag(genomicAlign, DNAFragAdaptor_fetchByDbID(dfa,consensusDNAFragId));
    GenomicAlign_setConsensusStart(genomicAlign, consensusStart);
    GenomicAlign_setConsensusEnd(genomicAlign, consensusEnd);
    GenomicAlign_setQueryDNAFrag(genomicAlign, DNAFragAdaptor_fetchByDbID(dfa,queryDNAFragId));
    GenomicAlign_setQueryStart(genomicAlign, queryStart);
    GenomicAlign_setQueryEnd(genomicAlign, queryEnd);
    GenomicAlign_setQueryStrand(genomicAlign, queryStrand);
    GenomicAlign_setAlignmentType(genomicAlign, alignmentType);
    GenomicAlign_setScore(genomicAlign, score);
    GenomicAlign_setPercentId(genomicAlign, percId);
    GenomicAlign_setCigarString(genomicAlign, cigarString);

    Vector_addElement(results, genomicAlign);
  }

  return results;
}


char *GenomicAlignAdaptor_alignmentTypeByMethodLinkId(GenomicAlignAdaptor *gaa, IDType methodLinkId) {
  StatementHandle *sth;
  ResultRow *row;
  char qStr[512];
  char *alignmentType;

  if (!methodLinkId) {
    fprintf(stderr, "Error: methodLinkId has to be defined");
    exit(1);
  } 

  sprintf(qStr,"SELECT type FROM method_link WHERE method_link_id = " IDFMTSTR, methodLinkId);
  sth = gaa->prepare((BaseAdaptor *)gaa, qStr, strlen(qStr));
  sth->execute(sth);

  if ((row = sth->fetchRow(sth))) {
    alignmentType = StrUtil_copyString(&alignmentType, row->getStringAt(row,0), 0); 
  } else {
    fprintf(stderr,"Error: No alignmentType for " IDFMTSTR "\n",methodLinkId);
    exit(1);
  }

  sth->finish(sth);

// NIY switch to using passed in string
  return alignmentType;
}

IDType GenomicAlignAdaptor_methodLinkIdByAlignmentType(GenomicAlignAdaptor *gaa, char *alignmentType) {
  StatementHandle *sth;
  ResultRow *row;
  IDType methodLinkId = 0;
  char qStr[512];

  if (!alignmentType) {
    fprintf(stderr, "Error: alignment_type has to be defined\n");
    exit(1);
  }
  
  sprintf(qStr, "SELECT method_link_id FROM method_link WHERE type = '%s'", alignmentType);
  sth = gaa->prepare((BaseAdaptor *)gaa, qStr, strlen(qStr));
  sth->execute(sth);

  if ((row = sth->fetchRow(sth))) {
    methodLinkId = row->getLongLongAt(row,0); 
  } else {
    fprintf(stderr,"Error: No methodLinkId for %s\n",alignmentType);
    exit(1);
  }

  sth->finish(sth);

  return methodLinkId;
}
