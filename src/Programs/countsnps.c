#include "DBAdaptor.h"
#include "Vector.h"
#include "Slice.h"
#include "StringHash.h"
#include "Chromosome.h"
#include "DNAAlignFeature.h"
#include "ChromosomeAdaptor.h"
#include "SliceAdaptor.h"

#include <stdlib.h>

typedef enum CodingTypesEnum {
  UNDEFINED,
  INTRONIC,
  CODING,
  NONCODING,
  INTERGENIC,
  NCODINGTYPE
} CodingTypes;

char *codingTypeStrings[] = { "Undefined", "Intronic", "Coding", "NonCoding", "Intergenic" };
  


int main(int argc, char **argv) {
  char *host = "ecs2b";
  char *user = "ensro";
  char *dbname = "steve_chr6_new_typed_snp";
  char *path = "SANGER_13";
  int port = 3306;
  char *geneTypes[] = {"Known","Novel_CDS",NULL};
  char *chromosomes[] = {"6",NULL};
  char **chrName;
  char *snpTypes[] = {NULL};

  DBAdaptor *dba;
  SliceAdaptor *sa;
  GeneAdaptor *ga;
  DNAAlignFeatureAdaptor *dafa;
  ChromosomeAdaptor *ca;


  initEnsC();

  dba = DBAdaptor_new(host,user,NULL,dbname,port,NULL);

  printf("Default assembly type %s\n",DBAdaptor_getAssemblyType(dba));

  DBAdaptor_setAssemblyType(dba, path);

  sa   = DBAdaptor_getSliceAdaptor(dba);
  ga   = DBAdaptor_getGeneAdaptor(dba);
  dafa = DBAdaptor_getDNAAlignFeatureAdaptor(dba);
  ca   = DBAdaptor_getChromosomeAdaptor(dba);

  chrName = chromosomes;
  while (*chrName) {
    char **geneType;
    Vector *genes = Vector_new();
    Vector *snps = Vector_new();
    Chromosome *chrom = ChromosomeAdaptor_fetchByChrName(ca, *chrName);
    int chrStart = 1;
    int chrEnd = Chromosome_getLength(chrom);
    Slice *slice = SliceAdaptor_fetchByChrStartEnd(sa, *chrName,chrStart,chrEnd);
    int i;
    StringHash *snpCodingType = StringHash_new(STRINGHASH_SMALL);
    char **snpType;
    
    fprintf(stderr, "Chr %s from %d to %d\n", *chrName, chrStart, chrEnd);

    printf("Fetching genes\n");
    
    geneType = geneTypes;

    while (*geneType) {
      Vector *genesOfType = Slice_getAllGenesByType(slice, *geneType);
      printf("Got %d %s genes\n",Vector_getNumElement(genesOfType),*geneType);

      Vector_append(genes,genesOfType);

      Vector_free(genesOfType,NULL);
      geneType++;
    }

    printf("Done fetching genes (fetched %d)\n",Vector_getNumElement(genes));

    Vector_sort(genes, SeqFeature_startCompFunc);

    printf("Fetching SNPs\n");

    snpType = snpTypes;
    if (!(*snpType)) {
      snps = Slice_getAllDNAAlignFeatures(slice,"",NULL);
    } else {
      while (*snpType) {
        Vector *snpsOfType = Slice_getAllDNAAlignFeatures(slice,*snpType,NULL);
        printf("Got %d %s SNPs\n",Vector_getNumElement(snpsOfType),*snpType);

        Vector_append(snps,snpsOfType);

        Vector_free(snpsOfType,NULL);
        snpType++;
      }
    }

    printf("Done fetching SNPs (fetched %d)\n",Vector_getNumElement(snps));

    printf("Starting sorting SNPs\n");
    Vector_sort(snps, SeqFeature_startCompFunc);
    printf("Done sorting SNPs\n");

    printf("Starting transcript sorting\n");
    for (i=0;i<Vector_getNumElement(genes); i++) {
      Gene *gene = Vector_getElementAt(genes,i); 
      int j;
      for (j=0; j<Gene_getTranscriptCount(gene); j++) {
        Transcript *trans = Gene_getTranscriptAt(gene,j);
        Transcript_sort(trans);
      }
    }
    printf("Done transcript sorting\n");

    printf("Starting coding type calcs\n");

    for (i=0;i<Vector_getNumElement(snps); i++) {
      DNAAlignFeature *snp = Vector_getElementAt(snps,i);
      int j;
      char *logicName;
      CodingTypes codingType = UNDEFINED;
      int done=0;
      int *counts;
  
      if (!(i%1000)) {
        printf(".");
        fflush(stdout);
      }
  
      for (j=0; j < Vector_getNumElement(genes) && !done; j++) {
        Gene *gene = Vector_getElementAt(genes,j); 
        if (!gene) {
          continue;
        }
        if (DNAAlignFeature_getEnd(snp) >= Gene_getStart(gene) && 
            DNAAlignFeature_getStart(snp) <= Gene_getEnd(gene)) {
          int k;
  
          if (codingType == UNDEFINED) codingType = INTRONIC;
  
          for (k=0; k<Gene_getTranscriptCount(gene) && !done; k++) {
            Transcript *trans = Gene_getTranscriptAt(gene,k);
            /* 
            printf("Translation %p Coding region start = %d coding region end = %d\n",
                   Transcript_getTranslation(trans),
                   Transcript_getCodingRegionStart(trans),
                   Transcript_getCodingRegionEnd(trans)); 
            */
  
            if (DNAAlignFeature_getEnd(snp) >= Transcript_getStart(trans) && 
                DNAAlignFeature_getStart(snp) <= Transcript_getEnd(trans)) {
              int m;
       
              for (m=0; m<Transcript_getExonCount(trans) && !done; m++) {
                Exon *exon = Transcript_getExonAt(trans,m);
  
                if (DNAAlignFeature_getEnd(snp) >= Exon_getStart(exon) && 
                    DNAAlignFeature_getStart(snp) <= Exon_getEnd(exon)) {
                  if (Transcript_getTranslation(trans) && 
                      DNAAlignFeature_getEnd(snp) >= Transcript_getCodingRegionStart(trans) && 
                      DNAAlignFeature_getStart(snp) <= Transcript_getCodingRegionEnd(trans)) {
                    codingType = CODING;
                    done = 1;
                  } else {
                    codingType = NONCODING;
                  }
                }
              }
            }
          }
        } else if (Gene_getStart(gene) > DNAAlignFeature_getEnd(snp)) {
          done = 1;
        } else if (Gene_getEnd(gene) < DNAAlignFeature_getStart(snp)) {
          /* 
          printf("Removing gene %s with extent %d to %d\n", 
                 Gene_getStableId(gene), 
                 Gene_getStart(gene),
                 Gene_getEnd(gene));
          */
          Vector_setElementAt(genes,j,NULL);
        }
      }

      logicName = Analysis_getLogicName(DNAAlignFeature_getAnalysis(snp));
      if (codingType == UNDEFINED) {
        codingType = INTERGENIC;
      }
      if (!StringHash_contains(snpCodingType, logicName)) {
        if ((counts = (int *)calloc(NCODINGTYPE,sizeof(int))) == NULL) {
          fprintf(stderr,"Error: Unable to allocate space for counts\n");
          exit(1);
        }
        StringHash_add(snpCodingType, logicName, counts);
      }
      counts = StringHash_getValue(snpCodingType,logicName);
      counts[codingType]++;
    }

    { // This brace is just me being lazy 
      char **keys = StringHash_getKeys(snpCodingType);
      int j;

      printf("\n\nCODING TYPE COUNTS\n\n");

      for (j=0; j<StringHash_getNumValues(snpCodingType); j++) {
        int *counts = StringHash_getValue(snpCodingType, keys[j]);
        int k;

        printf(" SNP Type %s\n", keys[j]);

        for (k=0; k<NCODINGTYPE; k++) {
          printf("  Type %s count = %d\n",codingTypeStrings[k], counts[k]);
        }
        printf("\n");
      }
      free(keys);
    }
    printf("Done coding type calcs\n");
    Vector_free(snps,NULL);
    Vector_free(genes,NULL);
    StringHash_free(snpCodingType,free);
    chrName++;
  }
  printf("Done\n");

  return 0;
}
