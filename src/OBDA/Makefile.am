
include $(top_srcdir)/Makefile.am.common

SUBDIRS = .

noinst_LTLIBRARIES = libOBDA.la

include_HEADERS = \
parsers.h \
BioIndex.h \
$(NULL)

libOBDA_la_SOURCES = \
parsers.c \
BioIndex.c \
$(NULL)


# Extra files to remove for the maintainer-clean target.
#
MAINTAINERCLEANFILES = $(top_srcdir)/OBDA/Makefile.in

