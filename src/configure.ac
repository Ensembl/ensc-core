#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

m4_include([./ensc_version.m4])
m4_include([./ax_lib_mysql.m4])

AC_PREREQ([2.61])
AC_INIT([EnsC], [VERSION_NUMBER], [annosoft@sanger.ac.uk])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_SRCDIR([Compara/DataAdaptors/BaseComparaAdaptor.c])
AM_INIT_AUTOMAKE([1.9 foreign])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([limits.h malloc.h stdlib.h string.h strings.h sys/param.h unistd.h])

AC_CONFIG_FILES([Makefile
Compara/Makefile
Compara/DataAdaptors/Makefile
Compara/DataModel/Makefile
DataAdaptors/Makefile
DataModel/Makefile
DataModel/Basic/Makefile
DataModel/Mapper/Makefile
OBDA/Makefile
Util/Makefile
Programs/Makefile
Tests/Makefile
])

# Check if this is a Git repository and, if so, include info about the current
# status in the version info. Otherwise use PACKAGE_VERSION (which will get
# set from VERSION_NUMBER in version.m4).
if test -d "$srcdir/../.git" 
then AC_DEFINE_UNQUOTED([ENSC_VERSION], ["`git describe --always --abbrev=1`"], ["EnsC package version"])
else AC_DEFINE([ENSC_VERSION], [PACKAGE_VERSION], ["EnsC package version"])
fi

AC_CONFIG_MACRO_DIR([m4])

LT_INIT([disable-shared])

LT_LANG([C++])

# keep libtool up to date automatically.
AC_SUBST([LIBTOOL_DEPS])


# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_INT32_T
AC_TYPE_INT8_T
AC_TYPE_OFF_T
AC_TYPE_SIZE_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([bzero getcwd memmove memset mkdir realpath stpcpy strcasecmp strchr strpbrk strrchr strstr strtol])


# Check for mysql
AX_LIB_MYSQL()           # Standard macro to test for mysql

AM_CONDITIONAL([HAVE_MYSQL_LIB], [ test -n "$MYSQL_VERSION" ])
AM_COND_IF([HAVE_MYSQL_LIB], [AC_SUBST([HAVE_MYSQL_LIB], ['mysql lib is available'])])
AM_COND_IF([HAVE_MYSQL_LIB], [AC_DEFINE([HAVE_MYSQL_LIB], [1], [define to 1 if mysql lib is available])])
AM_COND_IF([HAVE_MYSQL_LIB], [[ CFLAGS="$CFLAGS $MYSQL_CFLAGS" ]])

# Also need to check for headers - need mysql include directory in cflags or AC_CHECK_HEADERS
# fails
AC_CHECK_HEADERS([mysql.h], [HAVE_MYSQL_HEADER=1])
AM_CONDITIONAL([HAVE_MYSQL], [ test -n "$HAVE_MYSQL_HEADER" && test -n "$HAVE_MYSQL_LIB" ])
AM_COND_IF([HAVE_MYSQL], [ echo 'mysql is available' ])
AM_COND_IF([HAVE_MYSQL], [AC_SUBST(HAVE_MYSQL, 'mysql is available')])
AM_COND_IF([HAVE_MYSQL], [AC_DEFINE(HAVE_MYSQL, [1], [define to 1 if mysql is available])])

AC_CHECK_HEADERS([htslib/sam.h htslib/hts.h], [HAVE_SAMTOOLS_HEADER=1])
AM_CONDITIONAL([HAVE_SAMTOOLS], [ test -n "$HAVE_SAMTOOLS_HEADER" ])
AM_COND_IF([HAVE_SAMTOOLS], [ echo 'htslib is available' ])
AM_COND_IF([HAVE_SAMTOOLS], [AC_SUBST([HAVE_SAMTOOLS], ['htslib is available'])])
AM_COND_IF([HAVE_SAMTOOLS], [AC_DEFINE([HAVE_SAMTOOLS], [1], [define to 1 if htslib is available])])

AC_CHECK_HEADERS([libconfig.h], [HAVE_LIBCONFIG_HEADER=1])
AM_CONDITIONAL([HAVE_LIBCONFIG], [ test -n "$HAVE_LIBCONFIG_HEADER" ])
AM_COND_IF([HAVE_LIBCONFIG], [ echo 'libconfig is available' ])
AM_COND_IF([HAVE_LIBCONFIG], [AC_SUBST([HAVE_LIBCONFIG], ['libconfig is available'])])
AM_COND_IF([HAVE_LIBCONFIG], [AC_DEFINE([HAVE_LIBCONFIG], [1], [define to 1 if libconfig is available])])

AC_CHECK_HEADERS([gperftools/tcmalloc.h], [HAVE_LIBTCMALLOC_HEADER=1])
AM_CONDITIONAL([HAVE_LIBTCMALLOC], [ test -n "$HAVE_LIBTCMALLOC_HEADER" ])
AM_COND_IF([HAVE_LIBTCMALLOC], [ echo 'libtcmalloc is available' ])
AM_COND_IF([HAVE_LIBTCMALLOC], [AC_SUBST([HAVE_LIBTCMALLOC], ['libtcmalloc is available'])])
AM_COND_IF([HAVE_LIBTCMALLOC], [AC_DEFINE([HAVE_LIBTCMALLOC], [1], [define to 1 if libtcmalloc is available])])



AC_SUBST(PKG_CONFIG)
AC_SUBST(CFLAGS)
AC_SUBST(CPPFLAGS)
AC_SUBST(LDFLAGS)


AC_OUTPUT


echo \
"-------------------------------------------------

 ${PACKAGE_NAME} Version ${PACKAGE_VERSION}

 Prefix: '${prefix}'
 C++ Compiler: '${CXX} ${CXXFLAGS} ${CPPFLAGS}'
   C Compiler: '${CC} ${CFLAGS} ${CPPFLAGS}'

 Now type 'make @<:@<target>@:>@'
   where the optional @<:@<target>@:>@ is:
     all       - build all libraries (default)
     install   - build and install everything
"

if test -z "$MYSQL_VERSION"
then
  echo " Warning: MySql is not installed so ensc-core access is disabled."
  echo
fi

if test -z "$HAVE_SAMTOOLS"
then
  echo " samtools not available: bam* programs and RefineSolexaGenes will not be compiled"; echo
fi

if test -z "$HAVE_LIBCONFIG"
then
  echo " libconfig not available: RefineSolexaGenes will not be compiled"; echo
fi

if test -z "$HAVE_LIBTCMALLOC"
then
  echo " libtcmalloc not available: RefineSolexaGenes will be compiled using malloc"; echo
fi


echo "---------------------------------------------------------------------
"

