# EnsC makefile

#CC = /software/gcc-4.8.1/bin/gcc
CC = gcc

# use libs_r for thread safe version
MYSQLCFLAGS=${shell mysql_config --cflags}
MYSQLLIBS=${shell mysql_config --libs_r}

# Removed LIBCONFIG stuff because /software/ensembl/central does not exist
#LIBCONFIGCFLAGS=${shell ( export PKG_CONFIG_PATH=/software/ensembl/central/lib/pkgconfig; pkg-config --cflags libconfig ) }
#LIBCONFIGLIBS=${shell ( export PKG_CONFIG_PATH=/software/ensembl/central/lib/pkgconfig; pkg-config --libs libconfig ) }

# Define HAVE_SAMTOOLS to compile programs that use samtools_parallel
ifdef $(HAVE_SAMTOOLS)
  #SAMDIR=/ensemblweb/searle/ensdistprep_opt/longdir/longlonglonglonglonglong/long/bin/samtools-0.1.18
  #SAMDIR=/Users/searle/progs/ensc-dependencies/bin/samtools_parallel
  #SAMDIR=/nfs/users/nfs_s/searle/src/samtools_parallel
  #SAMDIR=/ensemblweb/searle/testcvs_pre70/ensc-core/src/samtools-0.1.18
  #SAMDIR=/ensemblweb/searle/testcvs_pre70/ensc-core/src/samtools_parallel
  #ZLIBDIR=/ensemblweb/searle/ensdistprep_opt/longdir/longlonglonglonglonglong/long/include
  #SAMDIR=/software/ensembl/central/bin/samtools-0.1.18
  SAMDIR=/nfs/ensembl/searle/progs/production_code/ensembl-trunk_0713/samtools_parallel
endif

# Removed ZLIBDIR because /software/ensembl/central does not exist
#ZLIBDIR=/software/ensembl/central/include

# *edgrif* -Ofast is not recognised...
#
CFLAGS = -Wall -Wno-format-security -Wno-unused-but-set-variable -IDataModel -IDataAdaptors -IUtil -IDataModel/Mapper -IDataModel/Basic -ICompara/DataAdaptors -ICompara/DataModel -IOBDA $(MYSQLCFLAGS) $(LIBCONFIGCFLAGS) -g -O0
#CFLAGS =  -Ofast -fprofile-arcs -ftest-coverage -Wno-format-security -IDataModel -IDataAdaptors -IUtil -IDataModel/Mapper -IDataModel/Basic -ICompara/DataAdaptors -ICompara/DataModel -IOBDA $(MYSQLCFLAGS)
#CFLAGS =  -Ofast -fbranch-probabilities -Wno-format-security -IDataModel -IDataAdaptors -IUtil -IDataModel/Mapper -IDataModel/Basic -ICompara/DataAdaptors -ICompara/DataModel -IOBDA $(MYSQLCFLAGS)

# *edgrif* -lunwind is not recognised....
#
#LDFLAGS = $(MYSQLLIBS) -L/usr/lib/system -lSystem -L/Users/searle/progs/ensc-dependencies/lib -ltcmalloc
#LDFLAGS = $(MYSQLLIBS) $(LIBCONFIGLIBS) -lunwind
LDFLAGS = $(MYSQLLIBS) $(LIBCONFIGLIBS) 



#CPPFLAGS = -I/Users/searle/progs/ensc-dependencies/include

BINDIR=../bin

BUILT_PROGRAMS = \
          directories \
          libensc.a \
          Tests/AssemblyMapperTest \
          Tests/CacheTest \
          Tests/ChainedAssemblyMapperTest \
          Tests/CigarStrUtilTest \
          Tests/ClassTest \
          Tests/ComparaDBAdaptorTest \
          Tests/ComparaDNAAlignFeatureAdaptorTest \
          Tests/CoordSystemTest \
          Tests/DBAdaptorTest \
          Tests/DNAPepAlignFeatureTest \
          Tests/DNAPepAlignFeatureWriteTest \
          Tests/EcoStringTest \
          Tests/HomologyTest \
          Tests/MapperTest \
          Tests/PredictionTranscriptTest \
          Tests/RepeatFeatureTest \
          Tests/RepeatFeatureWriteTest \
          Tests/SequenceAdaptorTest \
          Tests/SeqUtilTest \
          Tests/SimpleFeatureTest \
          Tests/SliceAdaptorTest \
          Tests/StreamTest \
          Tests/StrUtilTest \
          Tests/SyntenyTest \
          Tests/TopLevelAssemblyMapperTest \
          Tests/TranslateTest \
          Tests/VectorTest

# These Programs require samtools_parallel
ifdef $(HAVE_SAMTOOLS)
  BUILT_PROGRAMS += \
          Programs/testdbc \
          Programs/bamcount \
          Programs/bamcount_exon \
          Programs/bamcov \
          Programs/bammap \
          Programs/countsnps \
          Programs/indicate \
          Programs/fastasplit \
          Programs/moveperl
endif

# These require libtcmalloc and libconfig
ifdef $(HAVE_LIBTCMALLOC)
  ifdef $(HAVE_LIBCONFIG)
    BUILT_PROGRAMS += \
          Programs/RefineSolexaGenes \
          Tests/DNAAlignFeatureTest \
          Tests/DNAAlignFeatureWriteTest \
          Tests/FreeTest \
          Tests/GeneTest \
          Tests/GeneWriteTest \
          Tests/SimpleFeatureWriteTest
  endif
endif

all	: $(BUILT_PROGRAMS)

install	: $(BUILT_PROGRAMS)

directories	:
	mkdir -p objectfiles
	mkdir -p objectfiles/DataAdaptors
	mkdir -p objectfiles/DataModel
	mkdir -p objectfiles/OBDA
	mkdir -p objectfiles/Util
	mkdir -p objectfiles/Programs
	mkdir -p objectfiles/Tests
	mkdir -p objectfiles/DataModel/Basic
	mkdir -p objectfiles/DataModel/Mapper
	mkdir -p objectfiles/Compara/DataModel
	mkdir -p objectfiles/Compara/DataAdaptors

objectfiles/%.o: %.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS) -o $@ $<

libensc.a	: \
	objectfiles/Compara/DataAdaptors/BaseComparaAdaptor.o \
	objectfiles/Compara/DataAdaptors/ComparaDBAdaptor.o \
	objectfiles/Compara/DataAdaptors/ComparaDNAAlignFeatureAdaptor.o \
	objectfiles/Compara/DataAdaptors/DNAFragAdaptor.o \
	objectfiles/Compara/DataAdaptors/GenomeDBAdaptor.o \
	objectfiles/Compara/DataAdaptors/GenomicAlignAdaptor.o \
	objectfiles/Compara/DataAdaptors/HomologyAdaptor.o \
	objectfiles/Compara/DataAdaptors/SyntenyAdaptor.o \
	objectfiles/Compara/DataModel/DNAFrag.o \
	objectfiles/Compara/DataModel/GenomeDB.o \
	objectfiles/Compara/DataModel/GenomicAlign.o \
	objectfiles/Compara/DataModel/Homology.o \
	objectfiles/Compara/DataModel/SyntenyRegion.o \
	objectfiles/DataAdaptors/AnalysisAdaptor.o \
	objectfiles/DataAdaptors/AssemblyMapperAdaptor.o \
	objectfiles/DataAdaptors/AttributeAdaptor.o \
	objectfiles/DataAdaptors/BaseAdaptor.o \
	objectfiles/DataAdaptors/BaseFeatureAdaptor.o \
	objectfiles/DataAdaptors/CachingSequenceAdaptor.o \
	objectfiles/DataAdaptors/ChromosomeAdaptor.o \
	objectfiles/DataAdaptors/CoordSystemAdaptor.o \
	objectfiles/DataAdaptors/CloneAdaptor.o \
	objectfiles/DataAdaptors/DBAdaptor.o \
	objectfiles/DataAdaptors/DBConnection.o \
	objectfiles/DataAdaptors/DBEntryAdaptor.o \
	objectfiles/DataAdaptors/DNAAlignFeatureAdaptor.o \
	objectfiles/DataAdaptors/ExonAdaptor.o \
	objectfiles/DataAdaptors/GeneAdaptor.o \
	objectfiles/DataAdaptors/IntronSupportingEvidenceAdaptor.o \
	objectfiles/DataAdaptors/MetaContainer.o \
	objectfiles/DataAdaptors/MetaCoordContainer.o \
	objectfiles/DataAdaptors/MysqlResultRow.o \
	objectfiles/DataAdaptors/MysqlStatementHandle.o \
	objectfiles/DataAdaptors/PredictionTranscriptAdaptor.o \
	objectfiles/DataAdaptors/PredictionExonAdaptor.o \
	objectfiles/DataAdaptors/ProteinAlignFeatureAdaptor.o \
	objectfiles/DataAdaptors/RawContigAdaptor.o \
	objectfiles/DataAdaptors/RepeatConsensusAdaptor.o \
	objectfiles/DataAdaptors/RepeatFeatureAdaptor.o \
	objectfiles/DataAdaptors/SequenceAdaptor.o \
	objectfiles/DataAdaptors/SeqRegionCacheEntry.o \
	objectfiles/DataAdaptors/SimpleFeatureAdaptor.o \
	objectfiles/DataAdaptors/SliceAdaptor.o \
	objectfiles/DataAdaptors/SupportingFeatureAdaptor.o \
	objectfiles/DataAdaptors/TranscriptAdaptor.o \
	objectfiles/DataAdaptors/TranscriptSupportingFeatureAdaptor.o \
	objectfiles/DataAdaptors/TranslationAdaptor.o \
	objectfiles/DataModel/Analysis.o \
	objectfiles/DataModel/AnnotatedSeqFeature.o \
	objectfiles/DataModel/AssemblyMapper.o \
	objectfiles/DataModel/Attribute.o \
	objectfiles/DataModel/BaseAlignFeature.o \
	objectfiles/DataModel/BaseAssemblyMapper.o \
	objectfiles/DataModel/BaseContig.o \
	objectfiles/DataModel/Basic/Vector.o \
	objectfiles/DataModel/ChainedAssemblyMapper.o \
	objectfiles/DataModel/Chromosome.o \
	objectfiles/DataModel/CoordSystem.o \
	objectfiles/DataModel/CigarStrUtil.o \
	objectfiles/DataModel/Class.o \
	objectfiles/DataModel/Clone.o \
	objectfiles/DataModel/DBEntry.o \
	objectfiles/DataModel/DNAAlignFeature.o \
	objectfiles/DataModel/DNAPepAlignFeature.o \
	objectfiles/DataModel/EnsC.o \
	objectfiles/DataModel/Exon.o \
	objectfiles/DataModel/FeaturePair.o \
	objectfiles/DataModel/FeatureSet.o \
	objectfiles/DataModel/Gene.o \
	objectfiles/DataModel/GenomicRange.o \
	objectfiles/DataModel/IdentityXref.o \
	objectfiles/DataModel/Intron.o \
	objectfiles/DataModel/IntronSupportingEvidence.o \
	objectfiles/DataModel/Mapper.o \
	objectfiles/DataModel/Mapper/CoordPair.o \
	objectfiles/DataModel/Mapper/IndelCoordinate.o \
	objectfiles/DataModel/Mapper/MapperCoordinate.o \
	objectfiles/DataModel/Mapper/MapperGap.o \
	objectfiles/DataModel/Mapper/MapperPair.o \
	objectfiles/DataModel/Mapper/MapperPairSet.o \
	objectfiles/DataModel/Mapper/MapperRangeSet.o \
	objectfiles/DataModel/Mapper/MapperUnit.o \
	objectfiles/DataModel/Mapper/RangeRegistry.o \
	objectfiles/DataModel/Object.o \
	objectfiles/DataModel/PredictionTranscript.o \
	objectfiles/DataModel/PredictionExon.o \
	objectfiles/DataModel/ProjectionSegment.o \
	objectfiles/DataModel/RawContig.o \
	objectfiles/DataModel/RepeatConsensus.o \
	objectfiles/DataModel/RepeatFeature.o \
	objectfiles/DataModel/SeqFeature.o \
	objectfiles/DataModel/SeqFeatureFactory.o \
	objectfiles/DataModel/SeqRegionRange.o \
	objectfiles/DataModel/Sequence.o \
	objectfiles/DataModel/SeqEdit.o \
	objectfiles/DataModel/SimpleFeature.o \
	objectfiles/DataModel/Slice.o \
	objectfiles/DataModel/Species.o \
	objectfiles/DataModel/StableIdInfo.o \
	objectfiles/DataModel/Storable.o \
	objectfiles/DataModel/TopLevelAssemblyMapper.o \
	objectfiles/DataModel/Transcript.o \
	objectfiles/DataModel/Translation.o \
	objectfiles/OBDA/BioIndex.o \
	objectfiles/OBDA/parsers.o \
	objectfiles/Util/CHash.o \
	objectfiles/Util/Cache.o \
	objectfiles/Util/EcoString.o \
	objectfiles/Util/Error.o \
	objectfiles/Util/FileUtil.o \
	objectfiles/Util/IDHash.o \
	objectfiles/Util/LRUCache.o \
	objectfiles/Util/MysqlUtil.o \
	objectfiles/Util/ProcUtil.o \
	objectfiles/Util/SeqUtil.o \
	objectfiles/Util/StrUtil.o \
	objectfiles/Util/Stream.o \
	objectfiles/Util/StringHash.o \
	objectfiles/Util/translate.o
	ar rcv libensc.a objectfiles/DataModel/*.o objectfiles/DataAdaptors/*.o objectfiles/Util/*.o \
               objectfiles/DataModel/Basic/*.o objectfiles/DataModel/Mapper/*.o \
               objectfiles/Compara/DataAdaptors/*.o objectfiles/Compara/DataModel/*.o \
               objectfiles/OBDA/*.o
	ranlib libensc.a

Programs/fastasplit: \
	libensc.a \
	Programs/fastasplit.c
	$(CC) $(CFLAGS) -o Programs/fastasplit Programs/fastasplit.c -L. -lensc $(LDFLAGS)

Programs/indicate: \
	libensc.a \
	Programs/indicate.c
	$(CC) $(CFLAGS) -o Programs/indicate Programs/indicate.c -L. -lensc $(LDFLAGS)

Programs/testdbc: \
	libensc.a \
	Programs/testdbc.c
	$(CC) $(CFLAGS) -o Programs/testdbc Programs/testdbc.c -L. -lensc $(LDFLAGS)

Programs/bamcount: \
	libensc.a \
	Programs/bamcount.c
	$(CC) $(CFLAGS) -I$(SAMDIR) -I$(ZLIBDIR) -D_PBGZF_USE -o Programs/bamcount Programs/bamcount.c -L. -lensc  -L$(SAMDIR) -lbam $(LDFLAGS)

Programs/bamcount_exon: \
	libensc.a \
	Programs/bamcount_exon.c
	$(CC) $(CFLAGS) -I$(SAMDIR) -I$(ZLIBDIR) -D_PBGZF_USE -o Programs/bamcount_exon Programs/bamcount_exon.c -L. -lensc  -L$(SAMDIR) -lbam $(LDFLAGS)

Programs/bamcov: \
	libensc.a \
	Programs/bamcov.c
	$(CC) $(CFLAGS) -I$(SAMDIR) -I$(ZLIBDIR) -D_PBGZF_USE -o Programs/bamcov Programs/bamcov.c -L. -lensc  -L$(SAMDIR) -lbam $(LDFLAGS)

Programs/bammap: \
	libensc.a \
	Programs/bammap.c
	$(CC) $(CFLAGS) -I$(SAMDIR) -I$(ZLIBDIR) -D_PBGZF_USE -o Programs/bammap Programs/bammap.c -L. -lensc  -L$(SAMDIR) -lbam $(LDFLAGS)
#$(CC) $(CFLAGS) -I$(SAMDIR) -I$(ZLIBDIR) -o Programs/bammap Programs/bammap.c -L. -lensc  -L$(SAMDIR) -lbam $(LDFLAGS) -lrt

Programs/countsnps: \
	libensc.a \
	Programs/countsnps.c
	$(CC) $(CFLAGS) -o Programs/countsnps Programs/countsnps.c -L. -lensc $(LDFLAGS)

Programs/moveperl: \
	Programs/moveperl.c
	$(CC) $(CFLAGS) -o Programs/moveperl Programs/moveperl.c $(LDFLAGS)

Programs/RefineSolexaGenes: \
	libensc.a \
	Programs/RefineSolexaGenes.c
	$(CC) $(CFLAGS) -I$(SAMDIR) -I$(ZLIBDIR) -D_PBGZF_USE -o Programs/RefineSolexaGenes Programs/RefineSolexaGenes.c Util/translate.c -L. -lensc  -L$(SAMDIR) -lbam $(LDFLAGS)

Tests/AssemblyMapperTest: \
	libensc.a \
	Tests/AssemblyMapperTest.c
	$(CC) $(CFLAGS) -o Tests/AssemblyMapperTest Tests/AssemblyMapperTest.c -L. -lensc $(LDFLAGS)

Tests/CacheTest: \
	libensc.a \
	Tests/CacheTest.c
	$(CC) $(CFLAGS) -o Tests/CacheTest Tests/CacheTest.c -L. -lensc $(LDFLAGS)

Tests/ChainedAssemblyMapperTest: \
	libensc.a \
	Tests/ChainedAssemblyMapperTest.c
	$(CC) $(CFLAGS) -o Tests/ChainedAssemblyMapperTest Tests/ChainedAssemblyMapperTest.c -L. -lensc $(LDFLAGS)

Tests/ClassTest: \
	libensc.a \
	Tests/ClassTest.c
	$(CC) $(CFLAGS) -o Tests/ClassTest Tests/ClassTest.c -L. -lensc $(LDFLAGS)

Tests/CigarStrUtilTest: \
	libensc.a \
	Tests/CigarStrUtilTest.c
	$(CC) $(CFLAGS) -o Tests/CigarStrUtilTest Tests/CigarStrUtilTest.c -L. -lensc $(LDFLAGS)

Tests/ComparaDBAdaptorTest: \
	libensc.a \
	Tests/ComparaDBAdaptorTest.c
	$(CC) $(CFLAGS) -o Tests/ComparaDBAdaptorTest Tests/ComparaDBAdaptorTest.c -L. -lensc $(LDFLAGS)

Tests/ComparaDNAAlignFeatureAdaptorTest: \
	libensc.a \
	Tests/ComparaDNAAlignFeatureAdaptorTest.c
	$(CC) $(CFLAGS) -o Tests/ComparaDNAAlignFeatureAdaptorTest Tests/ComparaDNAAlignFeatureAdaptorTest.c -L. -lensc $(LDFLAGS)

Tests/CoordSystemTest: \
	libensc.a \
	Tests/CoordSystemTest.c
	$(CC) $(CFLAGS) -o Tests/CoordSystemTest Tests/CoordSystemTest.c -L. -lensc $(LDFLAGS)

Tests/DBAdaptorTest: \
	libensc.a \
	Tests/DBAdaptorTest.c
	$(CC) $(CFLAGS) -o Tests/DBAdaptorTest Tests/DBAdaptorTest.c -L. -lensc $(LDFLAGS)

Tests/DNAAlignFeatureTest: \
	libensc.a \
	Tests/BaseRODBTest.h \
	Tests/DNAAlignFeatureTest.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -o Tests/DNAAlignFeatureTest Tests/DNAAlignFeatureTest.c -L. -lensc $(LDFLAGS)

Tests/DNAAlignFeatureWriteTest: \
	libensc.a \
	Tests/BaseRODBTest.h \
	Tests/BaseRWDBTest.h \
	Tests/DNAAlignFeatureWriteTest.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -o Tests/DNAAlignFeatureWriteTest Tests/DNAAlignFeatureWriteTest.c -L. -lensc $(LDFLAGS)

Tests/DNAPepAlignFeatureTest: \
	libensc.a \
	Tests/BaseRODBTest.h \
	Tests/DNAPepAlignFeatureTest.c
	$(CC) $(CFLAGS) -o Tests/DNAPepAlignFeatureTest Tests/DNAPepAlignFeatureTest.c -L. -lensc $(LDFLAGS)

Tests/DNAPepAlignFeatureWriteTest: \
	libensc.a \
	Tests/BaseRODBTest.h \
	Tests/BaseRWDBTest.h \
	Tests/DNAPepAlignFeatureWriteTest.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -o Tests/DNAPepAlignFeatureWriteTest Tests/DNAPepAlignFeatureWriteTest.c -L. -lensc $(LDFLAGS)

Tests/EcoStringTest: \
	libensc.a \
	Tests/EcoStringTest.c
	$(CC) $(CFLAGS) -o Tests/EcoStringTest Tests/EcoStringTest.c -L. -lensc $(LDFLAGS)

Tests/FreeTest: \
	libensc.a \
	Tests/FreeTest.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -o Tests/FreeTest Tests/FreeTest.c -L. -lensc $(LDFLAGS)

Tests/GeneTest: \
	libensc.a \
	Tests/BaseRODBTest.h \
	Tests/GeneTest.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -o Tests/GeneTest Tests/GeneTest.c -L. -lensc $(LDFLAGS)

Tests/GeneWriteTest: \
	libensc.a \
	Tests/BaseRODBTest.h \
	Tests/BaseRWDBTest.h \
	Tests/GeneWriteTest.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -o Tests/GeneWriteTest Tests/GeneWriteTest.c -L. -lensc $(LDFLAGS)

Tests/HomologyTest: \
	libensc.a \
	Tests/HomologyTest.c
	$(CC) $(CFLAGS) -o Tests/HomologyTest Tests/HomologyTest.c -L. -lensc $(LDFLAGS)

Tests/LRUCacheTest: \
	libensc.a \
	Tests/LRUCacheTest.c
	$(CC) $(CFLAGS) -o Tests/LRUCacheTest Tests/LRUCacheTest.c -L. -lensc $(LDFLAGS)

Tests/MapperTest: \
	libensc.a \
	Tests/MapperTest.c
	$(CC) $(CFLAGS) -o Tests/MapperTest Tests/MapperTest.c -L. -lensc $(LDFLAGS)

Tests/PredictionTranscriptTest: \
	libensc.a \
	Tests/BaseRODBTest.h \
	Tests/PredictionTranscriptTest.c
	$(CC) $(CFLAGS) -o Tests/PredictionTranscriptTest Tests/PredictionTranscriptTest.c -L. -lensc $(LDFLAGS)

Tests/RepeatFeatureTest: \
	libensc.a \
	Tests/BaseRODBTest.h \
	Tests/RepeatFeatureTest.c
	$(CC) $(CFLAGS) -o Tests/RepeatFeatureTest Tests/RepeatFeatureTest.c -L. -lensc $(LDFLAGS)

Tests/RepeatFeatureWriteTest: \
	libensc.a \
	Tests/BaseRODBTest.h \
	Tests/BaseRWDBTest.h \
	Tests/RepeatFeatureWriteTest.c
	$(CC) $(CFLAGS) -o Tests/RepeatFeatureWriteTest Tests/RepeatFeatureWriteTest.c -L. -lensc $(LDFLAGS)

Tests/SequenceAdaptorTest: \
	libensc.a \
	Tests/BaseRODBTest.h \
	Tests/SequenceAdaptorTest.c
	$(CC) $(CFLAGS) -o Tests/SequenceAdaptorTest Tests/SequenceAdaptorTest.c -L. -lensc $(LDFLAGS)

Tests/SeqUtilTest: \
	libensc.a \
	Tests/SeqUtilTest.c
	$(CC) $(CFLAGS) -o Tests/SeqUtilTest Tests/SeqUtilTest.c -L. -lensc $(LDFLAGS)

Tests/SimpleFeatureTest: \
	libensc.a \
	Tests/BaseRODBTest.h \
	Tests/SimpleFeatureTest.c
	$(CC) $(CFLAGS) -o Tests/SimpleFeatureTest Tests/SimpleFeatureTest.c -L. -lensc $(LDFLAGS)

Tests/SimpleFeatureWriteTest: \
	libensc.a \
	Tests/BaseRODBTest.h \
	Tests/BaseRWDBTest.h \
	Tests/SimpleFeatureWriteTest.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o Tests/SimpleFeatureWriteTest Tests/SimpleFeatureWriteTest.c -L. -lensc $(LDFLAGS)

Tests/SliceAdaptorTest: \
	libensc.a \
	Tests/BaseRODBTest.h \
	Tests/SliceAdaptorTest.c
	$(CC) $(CFLAGS) -o Tests/SliceAdaptorTest Tests/SliceAdaptorTest.c -L. -lensc $(LDFLAGS)

Tests/StreamTest: \
	libensc.a \
	Tests/StreamTest.c
	$(CC) $(CFLAGS) -o Tests/StreamTest Tests/StreamTest.c -L. -lensc $(LDFLAGS)

Tests/StrUtilTest: \
	libensc.a \
	Tests/StrUtilTest.c
	$(CC) $(CFLAGS) -o Tests/StrUtilTest Tests/StrUtilTest.c -L. -lensc $(LDFLAGS)

Tests/SyntenyTest: \
	libensc.a \
	Tests/SyntenyTest.c
	$(CC) $(CFLAGS) -o Tests/SyntenyTest Tests/SyntenyTest.c -L. -lensc $(LDFLAGS)

Tests/TopLevelAssemblyMapperTest: \
	libensc.a \
	Tests/TopLevelAssemblyMapperTest.c
	$(CC) $(CFLAGS) -o Tests/TopLevelAssemblyMapperTest Tests/TopLevelAssemblyMapperTest.c -L. -lensc $(LDFLAGS)

Tests/TranslateTest: \
	libensc.a \
	Tests/TranslateTest.c
	$(CC) $(CFLAGS) -o Tests/TranslateTest Tests/TranslateTest.c -L. -lensc $(LDFLAGS)

Tests/VectorTest: \
	libensc.a \
	Tests/VectorTest.c
	$(CC) $(CFLAGS) -o Tests/VectorTest Tests/VectorTest.c -L. -lensc $(LDFLAGS)

# deletes compiled class files
clean	:	
	rm -f libensc.a
	rm -f Programs/testdbc
	rm -f Programs/bammap
	rm -f Programs/bamcount
	rm -f Programs/countsnps
	rm -f Programs/indicate
	find ./ -name "*.o" -exec rm {} \;
	find ./Tests ! -name "*.c" ! -name "*.h" -type f -perm -0700 -exec rm {} \;

distclean: clean
	-rm -f TAGS

# I haven't worked out what needs doing here, this is just a placeholder for now
# so our recursive make via the zmap src tree will find this target.
maintainer-clean:

#DO NOT REMOVE THIS LINE OR ADD ANYTHING AFTER IT - ADDMAKE DEPENDS ON IT FOR MAKING DEPENDENCIES

