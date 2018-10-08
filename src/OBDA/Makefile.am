
include $(top_srcdir)/Makefile.am.common

SUBDIRS = .

noinst_LTLIBRARIES = libOBDA.la

pkginclude_HEADERS = \
parsers.h \
BioIndex.h \
$(NULL)

libOBDA_la_SOURCES = \
parsers.c \
BioIndex.c \
$(NULL)


libOBDA_la_LIBADD = \
$(top_builddir)/Util/libUtil.la \
$(top_builddir)/DataModel/Basic/libDataModelBasic.la \
$(top_builddir)/DataModel/libDataModel.la \
$(NULL)


# Extra files to remove for the maintainer-clean target.
#
MAINTAINERCLEANFILES = $(top_srcdir)/OBDA/Makefile.in

