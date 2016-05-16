
include $(top_srcdir)/Makefile.am.common

SUBDIRS = .

noinst_LTLIBRARIES = libComparaDataAdaptors.la

include_HEADERS = \
BaseComparaAdaptor.h \
ComparaAdaptorTypes.h \
ComparaDBAdaptor.h \
ComparaDNAAlignFeatureAdaptor.h \
DNAFragAdaptor.h \
GenomeDBAdaptor.h \
GenomicAlignAdaptor.h \
HomologyAdaptor.h \
SyntenyAdaptor.h \
$(NULL)

libComparaDataAdaptors_la_SOURCES = \
BaseComparaAdaptor.c \
ComparaDBAdaptor.c \
ComparaDNAAlignFeatureAdaptor.c \
DNAFragAdaptor.c \
GenomeDBAdaptor.c \
GenomicAlignAdaptor.c \
HomologyAdaptor.c \
SyntenyAdaptor.c \
$(NULL)

# Extra files to remove for the maintainer-clean target.
#
MAINTAINERCLEANFILES = $(top_srcdir)/Compara/DataAdaptors/Makefile.in

