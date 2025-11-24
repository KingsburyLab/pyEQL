#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.61])
AC_INIT([IPhreeqc], [3.8.6-17100], [charlton@usgs.gov])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])
AC_CONFIG_SRCDIR([src/IPhreeqc.cpp])

# enable silent rules when available (automake 1.11 or later).
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

AC_DEFINE(SWIG_SHARED_OBJ)
AC_DEFINE(USE_PHRQ_ALLOC)
AC_DEFINE(NDEBUG)

# Checks for programs.
AC_PROG_CXX
AC_PROG_AWK
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
LT_INIT

# Check if the iphreeqc fortran module is enabled
AC_MSG_CHECKING([whether to enable the IPhreeqc Fortran module])
AC_ARG_ENABLE([fortran-module],
     [AS_HELP_STRING([--enable-fortran-module],[enable the IPhreeqc Fortran module @<:@default=yes@:>@])],
     [IPQ_FORTRAN_MODULE=$enableval],
     [IPQ_FORTRAN_MODULE=yes])

if test "X$IPQ_FORTRAN_MODULE" = "Xyes"; then
  AC_MSG_RESULT(yes)
else
  AC_MSG_RESULT(no)
  AC_DEFINE(IPHREEQC_NO_FORTRAN_MODULE)  
fi
AM_CONDITIONAL([FORTRAN_MODULE], [test "X$IPQ_FORTRAN_MODULE" = "Xyes"])

# Check if the fortran test should be included
AC_MSG_CHECKING([whether to enable the IPhreeqc Fortran test])
AC_ARG_ENABLE([fortran-test],
     [AS_HELP_STRING([--enable-fortran-test],[enable the IPhreeqc Fortran test @<:@default=no@:>@])],
     [IPQ_FORTRAN=$enableval],
     [IPQ_FORTRAN=no])

if test "X$IPQ_FORTRAN" = "Xyes"; then
  AC_MSG_RESULT(yes)
else
  AC_MSG_RESULT(no)
fi
AM_CONDITIONAL([BUILD_FORTRAN], [test "X$IPQ_FORTRAN" = "Xyes"])

if test "X$IPQ_FORTRAN_MODULE" = "Xno" || test "X$IPQ_FORTRAN" = "Xyes"; then
    AC_PROG_FC
    AC_FC_LIBRARY_LDFLAGS
    AC_FC_WRAPPERS 
    AC_FC_MAIN
    m4_ifdef([AC_FC_MODULE_EXTENSION], [AC_FC_MODULE_EXTENSION], [
      AX_F90_MODULE_EXTENSION
      FC_MODEXT=mod
      if test x$ax_cv_f90_modext != xunknown ; then
        FC_MODEXT=$ax_cv_f90_modext
      fi
      AC_SUBST(FC_MODEXT)
    ])
    AC_SUBST(F77, $FC)
    AC_SUBST(FFLAGS, $FCFLAGS)
fi

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([float.h limits.h memory.h stddef.h stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STRTOD

# libc functions
AC_CHECK_FUNCS([memmove], , AC_MSG_ERROR(memmove not found in libc))
AC_CHECK_FUNCS([memset],  , AC_MSG_ERROR(memset not found in libc))
AC_CHECK_FUNCS([strchr],  , AC_MSG_ERROR(strchr not found in libc))
AC_CHECK_FUNCS([strcspn], , AC_MSG_ERROR(strcspn not found in libc))
AC_CHECK_FUNCS([strtol],  , AC_MSG_ERROR(strtol not found in libc))

# libm functions
AC_CHECK_FUNCS([floor],   , AC_CHECK_LIB(m, floor, , AC_MSG_ERROR(cannot find floor)) )
AC_CHECK_FUNCS([pow],     , AC_CHECK_LIB(m, pow,   , AC_MSG_ERROR(cannot find pow)) )
AC_CHECK_FUNCS([sqrt],    , AC_CHECK_LIB(m, sqrt,  , AC_MSG_ERROR(cannot find sqrt)) )

# isfinite
AC_LANG_PUSH(C++)
AC_CACHE_CHECK([for isfinite], ac_cv_isfinite,
 [AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <math.h>]], [[double x; int y; y = isfinite(x);]])],[ac_cv_isfinite=yes],[ac_cv_isfinite=no
])])
if test x"$ac_cv_isfinite" = x"yes"; then
  AC_DEFINE(HAVE_ISFINITE, [], [Has isfinite])
fi
AC_LANG_POP(C++)

# finite
AC_LANG_PUSH(C++)
AC_CACHE_CHECK([for finite], ac_cv_finite,
 [AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <math.h>]], [[double x; int y; y = finite(x);]])],[ac_cv_finite=yes],[ac_cv_finite=no
])])
if test x"$ac_cv_finite" = x"yes"; then
  AC_DEFINE(HAVE_FINITE, [], [Has finite])
fi
AC_LANG_POP(C++)

# isnan
AC_CHECK_FUNCS([isnan],    , AC_CHECK_LIB(m, isnan))


AC_CONFIG_FILES([
  Makefile
  doc/Makefile
  examples/Makefile
  gtest/Makefile
  src/Makefile
  tests/Makefile
  tests/phreeqc.dat
  tests/ex2
  database/Makefile
])
AC_OUTPUT
