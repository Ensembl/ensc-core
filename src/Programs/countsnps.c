#include "DBAdaptor.h"

typedef enum CodingTypesEnum {
  UNDEFINED,
  INTRONIC,
  CODING,
  NONCODING,
  INTERGENIC,
  NCODINGTYPE
} CodingTypes;

char **codingTypeStrings = { "Undefined", "Intronic", "Coding", "NonCoding", "Intergenic" };
  


int main(int argc, char *argv) {
  char *host = "ecs2b";
  char *user = "ensro";
  char *dbname = "steve_chr6_new_typed_snp";
  char *path = "SANGER_13";
  char **geneTypes = {"Known","Novel_CDS",NULL};
  char **chromosomes = {"6",NULL}
  char *chrName;
  char **snpTypes = NULL;

  DBAdaptor *dba;
  SliceAdaptor *sa;
  GeneAdaptor *ga;
  DNAAlignFeatureAdaptor *dafa;


  initEnsC();

  dba = DBAdaptor_new(host,user,NULL,dbname,3306,NULL);

  printf("Default assembly type %s\n",DBAdaptor_getAssemblyType(dba));



if (scalar(@snptypes)) {
  @snptypes = split(/,/,join(',',@snptypes));
}

  DBAdaptor_setAssemblyType(path);

  sa   = DBAdaptor_getSliceAdaptor(dba);
  ga   = DBAdaptor_getGeneAdaptor(dba);
  dafa = DBAdaptor_getDNAAlignFeatureAdaptor(dba);
  ca   = DBAdaptor_getChromosomeAdaptor(dba);

  chrName = chromosomes;
  while (chrName) {
    char *geneType;
    Vector *genes = Vector_new();
    Vector *snps = Vector_new();
    Chromosome *chrom = ChromosomeAdaptor_fetchByChrName(ca, chrName);
    int chrStart = 1;
    int chrEnd = Chromosome_getLength(chrom);
    Slice *slice = SliceAdaptor_fetchByChrStartEnd(sa, chrName,chrStart,chrEnd);
    int count=0;
    StringHash *snpCodingType = StringHash_new(STRINGHASH_SMALL);
    
    fprintf(stderr, "Chr %s from %d to %d\n", chr, chrStart, chrEnd);

    printf("Fetching genes\n");
    
    geneType = geneTypes;

    while (geneType) {
      Vector *genesOfType = Slice_getAllGenes(slice, geneType);
      printf("Got %d %s genes\n",Vector_getNumElement(genesOfType),geneType);

      Vector_append(genes,genesOfType);

      Vector_free(genesOfType,NULL);
      geneType++;
    }

    printf("Done fetching genes (fetched %d)\n",Vector_getNumElement(genes));

    Vector_sort(genes, SeqFeature_startCompFunc);

    printf("Fetching SNPs\n");

    snpType = snpTypes;
    if (!snpType) {
      snps = Slice_getAllDNAAlignFeatures(slice,"",NULL);
    } else {
      while (snpType) {
        Vector *snpsOfType = Slice_getAllDNAAlignFeatures(slice,snpType,NULL);
        printf("Got %d %s SNPs\n",Vector_getNumElement(snpsOfType),snpType);

        Vector_append(snps,snpsOfType);

        Vector_free(snpsOfType,NULL);
        snpType++;
      }
    }

    printf("Done fetching SNPs (fetched %d)\n",Vector_getNumElement(snps));

    Vector_sort(genes, SeqFeature_startCompFunc);


    printf("Starting coding type calcs\n");

    for (i=0;i<Vector_getNumElement(snps); i++) {
      DNADNAAlignFeature *snp = Vector_getElementAt(snps,i);
      int j;
      char *logicName;
      CodingTypes codingType = UNDEFINED;
      int done=0;
      int *counts;
  
      if (!(i%1000)) {
        print ".";
      }
  
      for (j=0; j < Vector_getNumElement(genes) && !done; j++) {
        Gene *gene = genes[j]; 
        if (!gene) {
          continue;
        }
        if (DNADNAAlignFeature_getEnd(snp) >= Gene_getStart(gene) && 
            DNADNAAlignFeature_getStart(snp) <= Gene_getEnd(gene)) {
          int k;
  
          if (codingType == UNDEFINED) codingType = INTRONIC;
  
          for (k=0; k<Gene_getTranscriptCount(gene) && !done; k++) {
            Transcript *trans = Gene_getTranscriptAt(gene,k);
  
            if (DNADNAAlignFeature_getEnd(snp) >= Transcript_getStart(trans) && 
                DNADNAAlignFeature_getStart(snp) <= Transcript_getEnd(trans)) {
              int m;
       
              for (m=0; m<Transcript_getExonCount(trans); m++) {
                Exon *exon = Transcript_getExonAt(trans,m);
  
                if (DNADNAAlignFeature_getEnd(snp) >= Exon_getStart(exon) && 
                    DNADNAAlignFeature_getStart(snp) <= Exon_getEnd(exon)) {
                  if (Transcript_getTranslation(trans) && 
                      DNADNAAlignFeature_getEnd(snp) >= getCodingRegionStart(trans) && 
                      DNADNAAlignFeature_getStart(snp) <= getCodingRegionEnd(trans))
                    codingType = CODING;
                    done = 1;
                  } else {
                    codingType = NONCODING;
                  }
                }
              }
            }
          }
        } elsif (Gene_getStart(gene) > DNADNAAlignFeature_getEnd(snp)) {
          done = 1;
        } elsif (Gene_getEnd(gene) < DNADNAAlignFeature_getStart(snp)) {
          printf("Removing gene %s\n", Gene_getStableId(gene));
          Vector_setElementAt(genes,j,NULL);
        }
      }

      logicName = Analysis_getLogicName(DNADNAAlignFeature_getAnalysis(snp));
      if (codingType == UNDEFINED)) {
        codingType = INTERGENIC;
      }
      if (!StringHash_contains(logicName)) {
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
      char *key;
      char **keys = StringHash_getKeys(snpCodingType);
      int j;

      printf("\n\nCODING TYPE COUNTS\n\n");

      for (j=0; j<StringHash_getNumValue(snpCodingType); j++) {
        int *counts = StringHash_getValue(snpCodingType, keys[j]);
        int k;

        printf(" SNP Type %s\n", keys[j])

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
  }
  printf("Done\n");
}

int getCodingRegionStart(Transcript *trans) {
  return 500000000;
}

int getCodingRegionEnd(Transcript *trans) {
  return 500000000;
}
