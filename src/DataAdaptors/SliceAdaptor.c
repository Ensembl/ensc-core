#include "SliceAdaptor.h"
#include "BaseAdaptor.h"
#include "ChromosomeAdaptor.h"
#include "MysqlUtil.h"
#include "GenomicRange.h"

/* private methods */
GenomicRange *SliceAdaptor_getChrStartEndOfTranscript(SliceAdaptor *sa, char *fromClause, char *whereClause);
GenomicRange *SliceAdaptor_getChrStartEndOfGene(SliceAdaptor *sa, char *geneStableId);
GenomicRange *SliceAdaptor_getChrStartEndOfContig(SliceAdaptor *sa, char *contigName);

SliceAdaptor *SliceAdaptor_new(DBAdaptor *dba) {
  SliceAdaptor *sa;

  if ((sa = (SliceAdaptor *)calloc(1,sizeof(SliceAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for SliceAdaptor\n");
    return NULL;
  }
  BaseAdaptor_init((BaseAdaptor *)sa, dba, SLICE_ADAPTOR);

  return sa;
}

Slice *SliceAdaptor_fetchByChrStartEnd(SliceAdaptor *sa, char *chr, int start, int end) {
  Slice *slice;
  char *assemblyType;

  if (!chr) {
    fprintf(stderr,"ERROR: chromosome name argument must be defined and not ''\n");
    exit(1);
  }

  if(start > end) {
    fprintf(stderr,"ERROR: start must be less than end: parameters %s:%d:%d\n",chr,start,end);
    exit(1);
  }

  assemblyType = DBAdaptor_getAssemblyType(sa->dba);

  slice = Slice_new(chr, start, end, 1, assemblyType, sa, 0, FALSE);

  return slice;
}

Slice *SliceAdaptor_fetchByChrName(SliceAdaptor *sa, char *chr) {
  Slice *slice;
  int chrStart = 1;
  int chrEnd;
  ChromosomeAdaptor *ca;
  Chromosome *chromosome;
  char *assemblyType;

  // set the end of the slice to the end of the chromosome

  ca = DBAdaptor_getChromosomeAdaptor(sa->dba);
  chromosome = ChromosomeAdaptor_fetchByChrName(ca,chr); 

  if (!chromosome) {
    fprintf(stderr, "ERROR: Unknown chromosome %s\n",chr);
    exit(1);
  }

  chrEnd = Chromosome_getLength(chromosome);

  assemblyType = DBAdaptor_getAssemblyType(sa->dba);

  slice = Slice_new(chr, chrStart, chrEnd, 1, assemblyType, sa, 0, FALSE);

  return slice;
}

Slice *SliceAdaptor_fetchByContigName(SliceAdaptor *sa, char *name, int *sizeP) {
  GenomicRange *gr;
  Slice *slice;
  
  gr = SliceAdaptor_getChrStartEndOfContig(sa, name);
  
  if (sizeP) GenomicRange_expand(gr,*sizeP);
  
  slice = SliceAdaptor_fetchByChrStartEnd(sa,
                                          GenomicRange_getChrName(gr),
                                          GenomicRange_getChrStart(gr),
                                          GenomicRange_getChrEnd(gr));
  GenomicRange_free(gr);
  return slice;
}


/*
=head2 fetch_by_supercontig_name

  Arg [1]    : string $supercontig_name
  Example    : $slice = $slice_adaptor->fetch_by_supercontig_name('NT_004321');
  Description: Creates a Slice on the region of the assembly where 
               the specified super contig lies.  Note that this slice will
               have the same orientation as the supercontig. If the supercontig
               has a negative assembly orientation, the slice will also have
               a negative orientation relative to the assembly.
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : none
  Caller     : general

=cut

*/

Slice *SliceAdaptor_fetchBySupercontigName(SliceAdaptor *sa, char *supercontigName) {
  char *assemblyType;
  StatementHandle *sth;
  ResultRow *row;
  char qStr[1024];
  Slice *slice = NULL;
  
  assemblyType = DBAdaptor_getAssemblyType(sa->dba);

  if (!assemblyType) {
    fprintf(stderr,"Error: No assembly type defined\n");
    exit(1);
  }
  
  sprintf(qStr,
       "SELECT chr.name, a.superctg_ori, MIN(a.chr_start), MAX(a.chr_end)"
       " FROM assembly a, chromosome chr"
       " WHERE superctg_name = '%s'"
       " AND type = '%s'"
       " AND chr.chromosome_id = a.chromosome_id"
       " GROUP by superctg_name",
       supercontigName, assemblyType);
 
  sth = sa->prepare((BaseAdaptor *)sa,qStr,strlen(qStr));
  sth->execute(sth);
  
  if ((row = sth->fetchRow(sth))) {
    slice = Slice_new(row->getStringAt(row, 0),
                      row->getIntAt(row,2),
                      row->getIntAt(row,3),
                      row->getIntAt(row,1),
                      assemblyType, sa, 0, FALSE);
  }

  sth->finish(sth);

  if (!slice) {
    fprintf(stderr,"Error: Supercontig %s is not on the golden path. Cannot build Slice\n",
            supercontigName);
    exit(1);
  }
  
  return slice;
}

Vector *SliceAdaptor_listOverlappingSupercontigs(SliceAdaptor *sa, Slice *slice) {
  StatementHandle *sth;
  ResultRow *row;
  char qStr[1024];
  Vector *overlapping = Vector_new();

  sprintf(qStr,
      "SELECT DISTINCT superctg_name"
      "  FROM assembly a, chromosome c"
      " WHERE c.chromosome_id = a.chromosome_id"
      "   AND c.name = '%s'"
      "   AND a.type = '%s'"
      "   AND a.chr_end >= %d"
      "   AND a.chr_start <= %d",
      Slice_getChrName(slice),
      Slice_getAssemblyType(slice),
      Slice_getChrStart(slice),
      Slice_getChrEnd(slice));
     
  sth = sa->prepare((BaseAdaptor *)sa,qStr,strlen(qStr));
  sth->execute(sth);

  while ((row = sth->fetchRow(sth))) {
    Vector_addElement(overlapping,row->getStringAt(row,0));
  }
  sth->finish(sth);

  return overlapping;
}

Slice *SliceAdaptor_fetchByCloneAccession(SliceAdaptor *sa, char *clone, int *sizeP) {
  Slice *slice;
  GenomicRange *gr = NULL;
  char *type;
  StatementHandle *sth;
  ResultRow *row;
  char qStr[1024];
  
  if (!clone) {
    fprintf(stderr,"Error: Must have clone to fetch Slice of clone\n");
    exit(1);
  }  

  type = DBAdaptor_getAssemblyType(sa->dba);

  if (!type) {
    fprintf(stderr,"Error: No assembly type defined\n");
    exit(1);
  }

  sprintf(qStr,"SELECT  c.name,"
               "        a.chr_start,"
               "        a.chr_end,"
               "        chr.name "
               "    FROM    assembly a, "
               "        contig c, "
               "        clone  cl,"
               "        chromosome chr"
               "    WHERE c.clone_id = cl.clone_id"
               "    AND cl.name = '%s'  "
               "    AND c.contig_id = a.contig_id "
               "    AND a.type = '%s' "
               "    AND chr.chromosome_id = a.chromosome_id"
               "    ORDER BY a.chr_start", clone, type);

  sth = sa->prepare((BaseAdaptor *)sa,qStr,strlen(qStr));
  sth->execute(sth);

  while ((row = sth->fetchRow(sth))) {
    if (!gr) {
      gr = GenomicRange_new();
      GenomicRange_setChrStart(gr, row->getIntAt(row,1));
      GenomicRange_setChrName(gr, row->getStringAt(row,3));
    }
    GenomicRange_setChrEnd(gr, row->getIntAt(row,2));
  }
  sth->finish(sth);

  if (!gr) {
    fprintf(stderr,"Error: Clone %s is not on the golden path. Cannot build Slice\n",clone);
    exit(1);
  }
     
  if (sizeP) GenomicRange_expand(gr, *sizeP);

  slice = SliceAdaptor_fetchByChrStartEnd(sa,
                                          GenomicRange_getChrName(gr),
                                          GenomicRange_getChrStart(gr),
                                          GenomicRange_getChrEnd(gr));
  GenomicRange_free(gr);

  return slice;
}

/* NOTE NOTE NOTE Different implementation to perl */
Slice *SliceAdaptor_fetchByTranscriptStableId(SliceAdaptor *sa, char *transcriptStableId, int *sizeP) {
  char whereClause[256];
  GenomicRange *gr = NULL;
  Slice *slice;

  if (!transcriptStableId) {
    fprintf(stderr, "Error: Must have gene id to fetch Slice of gene\n");
    exit(1);
  }

  sprintf(whereClause,"tr.transcript_id = tsi.transcript_id"
                      " AND tsi.stable_id = '%s'", transcriptStableId);

  gr = SliceAdaptor_getChrStartEndOfTranscript(sa, "transcript_stable_id tsi,", whereClause);

  if (!gr) {
    char *type = DBAdaptor_getAssemblyType(sa->dba);

    if (!type) {
      fprintf(stderr,"Error: No assembly type defined\n");
      exit(1);
    }

    fprintf(stderr,"Error: Transcript [%s] is not on the golden path '%s'. " 
		   "Cannot build Slice.",transcriptStableId, type);
  }
     
  if (sizeP) GenomicRange_expand(gr, *sizeP);

  slice = SliceAdaptor_fetchByChrStartEnd(sa,
                                          GenomicRange_getChrName(gr),
                                          GenomicRange_getChrStart(gr),
                                          GenomicRange_getChrEnd(gr));
  GenomicRange_free(gr);

  return slice;
}

/* NOTE NOTE NOTE Different implementation to perl */
Slice *SliceAdaptor_fetchByTranscriptId(SliceAdaptor *sa, IDType transcriptId, int *sizeP) {
  char whereClause[256];
  GenomicRange *gr = NULL;
  Slice *slice;

  if (!transcriptId) {
    fprintf(stderr, "Error: Must have transcriptId id to fetch Slice of transcript\n");
    exit(1);
  }

  sprintf(whereClause,"tr.transcript_id = " IDFMTSTR, transcriptId);

  gr = SliceAdaptor_getChrStartEndOfTranscript(sa, "", whereClause);

  if (!gr) {
    char *type = DBAdaptor_getAssemblyType(sa->dba);

    if (!type) {
      fprintf(stderr,"Error: No assembly type defined\n");
      exit(1);
    }

    fprintf(stderr,"Error: Transcript [dbID = " IDFMTSTR "] is not on the golden path '%s'. " 
		   "Cannot build Slice.",transcriptId, type);
  }
     
  if (sizeP) GenomicRange_expand(gr, *sizeP);

  slice = SliceAdaptor_fetchByChrStartEnd(sa,
                                          GenomicRange_getChrName(gr),
                                          GenomicRange_getChrStart(gr),
                                          GenomicRange_getChrEnd(gr));
  GenomicRange_free(gr);

  return slice;
}

Slice *SliceAdaptor_fetchByGeneStableId(SliceAdaptor *sa, char *geneStableId, int *sizeP) {
  GenomicRange *gr = NULL;
  Slice *slice;

  if (!geneStableId) {
    fprintf(stderr, "Error: Must have gene id to fetch Slice of gene\n");
  }

  gr = SliceAdaptor_getChrStartEndOfGene(sa, geneStableId);

  if (!gr) {
    char *type = DBAdaptor_getAssemblyType(sa->dba);

    if (!type) {
      fprintf(stderr,"Error: No assembly type defined\n");
      exit(1);
    }

    fprintf(stderr,"Error: Gene [%s] is not on the golden path '%s'. " 
		   "Cannot build Slice.",geneStableId, type);
  }
     
  if (sizeP) GenomicRange_expand(gr, *sizeP);

  slice = SliceAdaptor_fetchByChrStartEnd(sa,
                                          GenomicRange_getChrName(gr),
                                          GenomicRange_getChrStart(gr),
                                          GenomicRange_getChrEnd(gr));
  GenomicRange_free(gr);

  return slice;
}

GenomicRange *SliceAdaptor_getChrStartEndOfContig(SliceAdaptor *sa, char *contigName) {
  char qStr[512];
  char *type;
  ResultRow *row;
  StatementHandle *sth;
  GenomicRange *gr = NULL;

  if (!contigName) {
    fprintf(stderr,"Error: Must have contigName to get contig extents\n");
    exit(1);
  }

  type = DBAdaptor_getAssemblyType(sa->dba);

  if (!type) {
    fprintf(stderr,"Error: No assembly type defined\n");
    exit(1);
  }

  sprintf(qStr,"SELECT  c.name,"
                  "     a.chr_start,"
                  "     a.chr_end,"
                  "     chr.name"
                  " FROM assembly a, contig c, chromosome chr"
                  " WHERE c.name = '%s'"
                  " AND c.contig_id = a.contig_id"
                  " AND a.type = '%s'"
                  " AND chr.chromosome_id = a.chromosome_id",
                   contigName, type);
  sth = sa->prepare((BaseAdaptor *)sa,qStr,strlen(qStr));
  sth->execute(sth);

  if ((row = sth->fetchRow(sth))) {
    gr = GenomicRange_new();
    GenomicRange_setChrName(gr, row->getStringAt(row,2));
    GenomicRange_setChrStart(gr, row->getIntAt(row,0));
    GenomicRange_setChrStart(gr, row->getIntAt(row,1));
  }
  
  sth->finish(sth);

  if (!gr) {
    fprintf(stderr, "Error: Contig %s is not on the golden path of type %s", contigName, type);
    exit(1);
  }

  return gr;
}

GenomicRange *SliceAdaptor_getChrStartEndOfGene(SliceAdaptor *sa, char *geneStableId) {
  StatementHandle *sth;
  ResultRow *row;
  char *type;
  GenomicRange *gr = NULL;
  char qStr[512];
  
  type = DBAdaptor_getAssemblyType(sa->dba);

  if (!type) {
    fprintf(stderr,"Error: No assembly type defined\n");
    exit(1);
  }

  sprintf(qStr,"SELECT "  
   "if(a.contig_ori=1,(e.contig_start-a.contig_start+a.chr_start),"
                     " (a.chr_start+a.contig_end-e.contig_end)),"
   " if(a.contig_ori=1,(e.contig_end-a.contig_start+a.chr_start),"
                      " (a.chr_start+a.contig_end-e.contig_start)),"
     " chr.name"
                   " FROM    exon e,"
                   "     transcript tr,"
                   "     exon_transcript et,"
                   "     assembly a,"
                   "     gene_stable_id gsi,"
                   "     chromosome chr"
                   " WHERE e.exon_id=et.exon_id "
                   " AND et.transcript_id =tr.transcript_id "
                   " AND a.contig_id=e.contig_id "
                   " AND a.type = '%s'"
                   " AND tr.gene_id = gsi.gene_id"
                   " AND gsi.stable_id = '%s'"
                   " AND a.chromosome_id = chr.chromosome_id",
                  type, geneStableId); 
  
  sth = sa->prepare((BaseAdaptor *)sa,qStr,strlen(qStr));
  sth->execute(sth);

  while ((row = sth->fetchRow(sth))) {
    int start = row->getIntAt(row,0);
    int end   = row->getIntAt(row,1);

    if (start > end) { int tmp = start; start = end; end = tmp; } // swap start and end - paranoia

    if (!gr) {
      gr = GenomicRange_new();
      GenomicRange_setChrName(gr, row->getStringAt(row,2));
      GenomicRange_setChrStart(gr, start);
      GenomicRange_setChrStart(gr, end);
    }

    if (start < GenomicRange_getChrStart(gr)) GenomicRange_setChrStart(gr, start); 
    if (end   > GenomicRange_getChrEnd(gr)) GenomicRange_setChrEnd(gr, end); 
  }   

  sth->finish(sth);

  return gr;
}

GenomicRange *SliceAdaptor_getChrStartEndOfTranscript(SliceAdaptor *sa, char *fromClause, char *whereClause) {
  StatementHandle *sth;
  ResultRow *row;
  char *type;
  GenomicRange *gr = NULL;
  char qStr[1024];
  
  type = DBAdaptor_getAssemblyType(sa->dba);

  if (!type) {
    fprintf(stderr,"Error: No assembly type defined\n");
    exit(1);
  }

  sprintf(qStr,"SELECT "  
   "if(a.contig_ori=1,(e.contig_start-a.contig_start+a.chr_start),"
                     " (a.chr_start+a.contig_end-e.contig_end)),"
   " if(a.contig_ori=1,(e.contig_end-a.contig_start+a.chr_start),"
                      " (a.chr_start+a.contig_end-e.contig_start)),"
     " chr.name"
                   " FROM    exon e,"
                   "     transcript tr,"
                   "     exon_transcript et,"
                   "     assembly a,"
                   "     %s"
                   "     chromosome chr"
                   " WHERE e.exon_id=et.exon_id "
                   " AND et.transcript_id =tr.transcript_id "
                   " AND a.contig_id=e.contig_id "
                   " AND %s"
                   " AND a.type = '%s'"
                   " AND a.chromosome_id = chr.chromosome_id",
                  fromClause, whereClause, type); 
  
  sth = sa->prepare((BaseAdaptor *)sa,qStr,strlen(qStr));
  sth->execute(sth);

  while ((row = sth->fetchRow(sth))) {
    int start = row->getIntAt(row,0);
    int end   = row->getIntAt(row,1);

    if (start > end) { int tmp = start; start = end; end = tmp; } // swap start and end - paranoia

    if (!gr) {
      gr = GenomicRange_new();
      GenomicRange_setChrName(gr, row->getStringAt(row,2));
      GenomicRange_setChrStart(gr, start);
      GenomicRange_setChrStart(gr, end);
    }

    if (start < GenomicRange_getChrStart(gr)) GenomicRange_setChrStart(gr, start); 
    if (end   > GenomicRange_getChrEnd(gr)) GenomicRange_setChrEnd(gr, end); 
  }   

  sth->finish(sth);

  return gr;
}
