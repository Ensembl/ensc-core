
include $(top_srcdir)/Makefile.am.common

SUBDIRS = .

noinst_LTLIBRARIES = libUtil.la

pkginclude_HEADERS = \
CHash.h \
Cache.h \
EcoString.h \
EnsC.h \
Error.h \
FileUtil.h \
LRUCache.h \
Message.h \
MysqlUtil.h \
ProcUtil.h \
SeqUtil.h \
StrUtil.h \
Stream.h \
StringHash.h \
tplib.h \
translate.h \
$(NULL)

libUtil_la_SOURCES = \
CHash.c \
Cache.c \
EcoString.c \
EnsC.c \
Error.c \
FileUtil.c \
LRUCache.c \
MysqlUtil.c \
ProcUtil.c \
SeqUtil.c \
StrUtil.c \
Stream.c \
StringHash.c \
translate.c \
$(NULL)


# Extra files to remove for the maintainer-clean target.
#
MAINTAINERCLEANFILES = $(top_srcdir)/Util/Makefile.in

