
include $(top_srcdir)/Makefile.am.common

SUBDIRS = .

bin_PROGRAMS = countsnps fastasplit indicate moveperl testdbc

PROG_LIBS = \
$(top_builddir)/Util/libUtil.la \
$(top_builddir)/DataModel/Basic/libDataModelBasic.la \
$(top_builddir)/DataModel/Mapper/libDataModelMapper.la \
$(top_builddir)/DataModel/libDataModel.la \
$(top_builddir)/DataAdaptors/libDataAdaptors.la \
$(top_builddir)/Compara/DataAdaptors/libComparaDataAdaptors.la \
$(top_builddir)/Compara/DataModel/libComparaDataModel.la \
$(top_builddir)/OBDA/libOBDA.la \
$(NULL)

countsnps_SOURCES = countsnps.c
fastasplit_SOURCES = fastasplit.c
indicate_SOURCES = indicate.c indicate.h
moveperl_SOURCES = moveperl.c
testdbc_SOURCES = testdbc.c

countsnps_LDADD = $(PROG_LIBS)
fastasplit_LDADD = $(PROG_LIBS)
indicate_LDADD = $(PROG_LIBS)
moveperl_LDADD = $(PROG_LIBS)
testdbc_LDADD = $(PROG_LIBS)


if HAVE_SAMTOOLS
  bin_PROGRAMS += bamcount bamcount_exon bamcov bammap

  PROG_LIBS += -lhts -lz

  bamcount_SOURCES = bamcount.c bamhelper.h
  bamcount_exon_SOURCES = bamcount_exon.c bamhelper.h
  bamcov_SOURCES = bamcov.c
  bammap_SOURCES = bammap.c bamhelper.h

  bamcount_LDADD = $(PROG_LIBS)
  bamcount_exon_LDADD = $(PROG_LIBS)
  bamcov_LDADD = $(PROG_LIBS)
  bammap_LDADD = $(PROG_LIBS)

if HAVE_LIBTCMALLOC
  bin_PROGRAMS += RefineSolexaGenes
  RefineSolexaGenes_SOURCES = RefineSolexaGenes.c RefineSolexaGenes.h 

  PROG_LIBS += -lconfig -ltcmalloc
  RefineSolexaGenes_LDADD = $(PROG_LIBS)
endif
endif



# Extra files to remove for the maintainer-clean target.
#
MAINTAINERCLEANFILES = $(top_srcdir)/Programs/Makefile.in

