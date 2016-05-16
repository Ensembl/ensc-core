
include $(top_srcdir)/Makefile.am.common

SUBDIRS = .

noinst_bindir = .

noinst_bin_PROGRAMS = \
AssemblyMapperTest \
CacheTest \
ChainedAssemblyMapperTest \
CigarStrUtilTest \
ClassTest \
ComparaDBAdaptorTest \
ComparaDNAAlignFeatureAdaptorTest \
CoordSystemTest \
DBAdaptorTest \
DNAPepAlignFeatureTest \
DNAPepAlignFeatureWriteTest \
EcoStringTest \
HomologyTest \
MapperTest \
PredictionTranscriptTest \
RepeatFeatureTest \
RepeatFeatureWriteTest \
SeqUtilTest \
SequenceAdaptorTest \
SimpleFeatureTest \
SliceAdaptorTest \
StrUtilTest \
StreamTest \
SyntenyTest \
TopLevelAssemblyMapperTest \
TranslateTest \
VectorTest

TEST_LIBS = \
$(top_builddir)/Compara/DataAdaptors/.libs/libComparaDataAdaptors.a \
$(top_builddir)/Compara/DataModel/.libs/libComparaDataModel.a \
$(top_builddir)/DataAdaptors/.libs/libDataAdaptors.a \
$(top_builddir)/DataModel/.libs/libDataModel.a \
$(top_builddir)/DataModel/Basic/.libs/libDataModelBasic.a \
$(top_builddir)/DataModel/Mapper/.libs/libDataModelMapper.a \
$(top_builddir)/OBDA/.libs/libOBDA.a \
$(top_builddir)/Util/.libs/libUtil.a \
$(NULL)

AssemblyMapperTest_SOURCES = AssemblyMapperTest.c BaseRODBTest.h BaseTest.h
CacheTest_SOURCES = CacheTest.c BaseTest.h
ChainedAssemblyMapperTest_SOURCES = ChainedAssemblyMapperTest.c BaseRODBTest.h BaseTest.h
CigarStrUtilTest_SOURCES = CigarStrUtilTest.c BaseTest.h
ClassTest_SOURCES = ClassTest.c BaseTest.h
ComparaDBAdaptorTest_SOURCES = ComparaDBAdaptorTest.c BaseTest.h
ComparaDNAAlignFeatureAdaptorTest_SOURCES = ComparaDNAAlignFeatureAdaptorTest.c
CoordSystemTest_SOURCES = CoordSystemTest.c BaseRODBTest.h BaseTest.h
DBAdaptorTest_SOURCES = DBAdaptorTest.c BaseTest.h
DNAPepAlignFeatureTest_SOURCES = DNAPepAlignFeatureTest.c BaseRODBTest.h BaseTest.h
DNAPepAlignFeatureWriteTest_SOURCES = DNAPepAlignFeatureWriteTest.c BaseRODBTest.h BaseRWDBTest.h BaseTest.h
EcoStringTest_SOURCES = EcoStringTest.c BaseTest.h
HomologyTest_SOURCES = HomologyTest.c BaseComparaDBTest.h BaseTest.h
MapperTest_SOURCES = MapperTest.c BaseRODBTest.h BaseTest.h
PredictionTranscriptTest_SOURCES = PredictionTranscriptTest.c BaseRODBTest.h BaseTest.h
RepeatFeatureTest_SOURCES = RepeatFeatureTest.c BaseRODBTest.h BaseTest.h
RepeatFeatureWriteTest_SOURCES = RepeatFeatureWriteTest.c BaseRODBTest.h BaseRWDBTest.h BaseTest.h
SeqUtilTest_SOURCES = SeqUtilTest.c BaseTest.h
SequenceAdaptorTest_SOURCES = SequenceAdaptorTest.c BaseTest.h
SimpleFeatureTest_SOURCES = SimpleFeatureTest.c BaseRODBTest.h BaseTest.h
SliceAdaptorTest_SOURCES = SliceAdaptorTest.c BaseTest.h
StrUtilTest_SOURCES = StrUtilTest.c BaseTest.h
StreamTest_SOURCES = StreamTest.c BaseTest.h
SyntenyTest_SOURCES = SyntenyTest.c BaseComparaDBTest.h BaseTest.h
TopLevelAssemblyMapperTest_SOURCES = TopLevelAssemblyMapperTest.c BaseRODBTest.h BaseTest.h
TranslateTest_SOURCES = TranslateTest.c BaseTest.h
VectorTest_SOURCES = VectorTest.c BaseTest.h                                                        

if HAVE_LIBCONFIG
  if HAVE_LIBTCMALLOC
    noinst_bin_PROGRAMS += \
          DNAAlignFeatureTest \
          DNAAlignFeatureWriteTest \
          FreeTest \
          GeneTest \
          GeneWriteTest \
          SimpleFeatureWriteTest

    TEST_LIBS += $(DEPS_LIBCONFIG) $(DEPS_LIBTCMALLOC)

    DNAAlignFeatureTest_SOURCES = DNAAlignFeatureTest.c BaseRODBTest.h BaseTest.h
    DNAAlignFeatureWriteTest_SOURCES = DNAAlignFeatureWriteTest.c BaseRODBTest.h BaseRWDBTest.h BaseTest.h
    FreeTest_SOURCES = FreeTest.c BaseRODBTest.h BaseTest.h
    GeneTest_SOURCES = GeneTest.c BaseRODBTest.h BaseTest.h
    GeneWriteTest_SOURCES = GeneWriteTest.c BaseRODBTest.h BaseRWDBTest.h BaseTest.h
    SimpleFeatureWriteTest_SOURCES = SimpleFeatureWriteTest.c BaseRODBTest.h BaseRWDBTest.h BaseTest.h
  endif
endif 


AssemblyMapperTest_LDADD = $(TEST_LIBS)
CacheTest_LDADD = $(TEST_LIBS)
ChainedAssemblyMapperTest_LDADD = $(TEST_LIBS)
CigarStrUtilTest_LDADD = $(TEST_LIBS)
ClassTest_LDADD = $(TEST_LIBS)
ComparaDBAdaptorTest_LDADD = $(TEST_LIBS)
ComparaDNAAlignFeatureAdaptorTest_LDADD = $(TEST_LIBS)
CoordSystemTest_LDADD = $(TEST_LIBS)
DBAdaptorTest_LDADD = $(TEST_LIBS)
DNAPepAlignFeatureTest_LDADD = $(TEST_LIBS)
DNAPepAlignFeatureWriteTest_LDADD = $(TEST_LIBS)
EcoStringTest_LDADD = $(TEST_LIBS)
HomologyTest_LDADD = $(TEST_LIBS)
MapperTest_LDADD = $(TEST_LIBS)
PredictionTranscriptTest_LDADD = $(TEST_LIBS)
RepeatFeatureTest_LDADD = $(TEST_LIBS)
RepeatFeatureWriteTest_LDADD = $(TEST_LIBS)
SeqUtilTest_LDADD = $(TEST_LIBS)
SequenceAdaptorTest_LDADD = $(TEST_LIBS)
SimpleFeatureTest_LDADD = $(TEST_LIBS)
SliceAdaptorTest_LDADD = $(TEST_LIBS)
StrUtilTest_LDADD = $(TEST_LIBS)
StreamTest_LDADD = $(TEST_LIBS)
SyntenyTest_LDADD = $(TEST_LIBS)
TopLevelAssemblyMapperTest_LDADD = $(TEST_LIBS)
TranslateTest_LDADD = $(TEST_LIBS)
VectorTest_LDADD = $(TEST_LIBS)

TEST_LIBS += $(MYSQL_LDFLAGS)


# Extra files to remove for the maintainer-clean target.
#
MAINTAINERCLEANFILES = $(top_srcdir)/Tests/Makefile.in

