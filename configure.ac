# switch off caching, since that leads to insurmountable trouble if
# you call ./configure subsequently on different systems without
# clearing the cache
define([AC_CACHE_LOAD], )
define([AC_CACHE_SAVE], )
AC_REVISION($Revision: 1.1 $)


AC_INIT
AC_CONFIG_FILES(Make.common)

# Compute the canonical target system type variable "target", and its
# three individual parts: "target_cpu", "target_vendor", and "target_os"
AC_CANONICAL_TARGET

# The AC_SUBST(foo) command creates an output variable from a shell variable.
# This command also causes AC_OUTPUT to replace all instances of @VAR@ with
# the value of the shell variable VAR in the output files (e.g. one or more Makefiles).
AC_SUBST(host)

# This works because we always run configure from the directory where
# our script installs GMP and MPFR, but it is very surprising to me
# that things like $abs_srcdir and $ac_abs_srcdir don't work at all
# for this purpose...
MY_TOP_SRCDIR=`pwd`

# Location of the GMP installation.  Defaults to `pwd`/gmp, which is
# where the build_gmp_mpfr.sh script installs it.
AC_ARG_WITH(gmp-include,
            AC_HELP_STRING([--with-gmp-include=PATH],[Specify the path for GMP header files]),
            withgmpinc=$withval,
            withgmpinc=$MY_TOP_SRCDIR/gmp/include)
AC_SUBST(withgmpinc)

AC_ARG_WITH(gmp-lib,
            AC_HELP_STRING([--with-gmp-lib=PATH],[Specify the path for GMP libraries]),
            withgmplib=$withval,
            withgmplib=$MY_TOP_SRCDIR/gmp/lib)
AC_SUBST(withgmplib)


# Location of the MPFR installation.  Defaults to `pwd`/mpfr, which is
# where the build_gmp_mpfr.sh script installs it.
AC_ARG_WITH(mpfr-include,
            AC_HELP_STRING([--with-mpfr-include=PATH],[Specify the path for MPFR header files]),
            withmpfrinc=$withval,
            withmpfrinc=$MY_TOP_SRCDIR/mpfr/include)
AC_SUBST(withmpfrinc)

AC_ARG_WITH(mpfr-lib,
            AC_HELP_STRING([--with-mpfr-lib=PATH],[Specify the path for MPFR libraries]),
            withmpfrlib=$withval,
            withmpfrlib=$MY_TOP_SRCDIR/mpfr/lib)
AC_SUBST(withmpfrlib)

# Initialize libtool
LT_INIT

# Set the shell variable LIBTOOL_DEPS, which can be used to automatically update the libtool script.
# See the special rule for this purpose in Make.common.in
AC_SUBST([LIBTOOL_DEPS])

# Makes it so that libtool supports the --tag=CXX syntax
LT_LANG([C++])

#
# Testing GMP
#

# FIXME: what is the "right" way to mess with the CFLAGS/LDFLAGS/LIBS
# variables in configure scripts?  I've seen some people that save the
# original values and restore them after testing... maybe I should do
# that?
CFLAGS="$CFLAGS -I$withgmpinc"

# Test compiling a non-trivial program that includes gmp.h.  Extra
# paths are passed to the compiler via the CFLAGS variable.  If GMP is
# not found, we error out during configure by calling AC_MSG_ERROR.
AC_TRY_COMPILE(
[
#include <gmp.h>
],
[
],
[AC_MSG_RESULT([gmp.h found in $withgmpinc])],
[AC_MSG_ERROR([GNU MP not found, see http://gmplib.org/])]
)

# Also check for linking with GMP libraries...  Need to add paths to
# LDFLAGS for this.
LDFLAGS="$LDFLAGS -L$withgmplib"

AC_CHECK_LIB(
gmp,
__gmpz_init,
[AC_MSG_RESULT([__gmpz_init found within libgmp])],
[AC_MSG_ERROR([Could not find __gmpz_init within libgmp, see http://gmplib.org/])]
)

# Try actually linking a non-trivial program which requires -lgmp.
# To do this, we need to set up the right LIBS.
LIBS="$LIBS -lgmpxx -lgmp"

AC_TRY_LINK(
[
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
],
[
mpz_t x;
mpz_t y;
mpz_t result;

mpz_init(x);
mpz_init(y);
mpz_init(result);

mpz_set_str(x, "7612058254738945", 10);
mpz_set_str(y, "9263591128439081", 10);

mpz_mul(result, x, y);
gmp_printf("\n    %Zd\n*\n    %Zd\n--------------------\n%Zd\n\n", x, y, result);

/* free used memory */
mpz_clear(x);
mpz_clear(y);
mpz_clear(result);
return EXIT_SUCCESS;
],
[AC_MSG_RESULT([Non-trivial GMP program successfully linked!])],
[AC_MSG_ERROR([Trouble linking non-trivial GMP program!])]
)



#
# Testing MPFR
#

# More configure testing, this time for MPFR.  Note that these test
# programs will be compiled with GMP include paths as well...
CFLAGS="$CFLAGS -I$withmpfrinc"

# Test compiling a non-trivial program that includes mpfr.h
AC_TRY_COMPILE(
[
#include <mpfr.h>
],
[
],
[AC_MSG_RESULT([mpfr.h found in $withmpfrinc])],
[AC_MSG_ERROR([GNU MPFR not found, see http://www.mpfr.org/])]
)


# Also check for linking with MPFR libraries...  Need to add paths to
# LDFLAGS for this.
LDFLAGS="$LDFLAGS -L$withmpfrlib"

AC_CHECK_LIB(
mpfr,
mpfr_init,
[AC_MSG_RESULT([mpfr_init found within libmpfr])],
[AC_MSG_ERROR([Could not find mpfr_init within libmpfr, see http://www.mpfr.org/])]
)


# Try actually linking a non-trivial program which requires -lmpfr.
# To do this, we need to set up the right LIBS.
LIBS="$LIBS -lmpfr"

AC_TRY_LINK(
[
#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>
],
[
unsigned int i;
mpfr_t s, t, u;

mpfr_init2 (t, 200);
mpfr_set_d (t, 1.0, MPFR_RNDD);
mpfr_init2 (s, 200);
mpfr_set_d (s, 1.0, MPFR_RNDD);
mpfr_init2 (u, 200);
for (i = 1; i <= 100; i++)
  {
    mpfr_mul_ui (t, t, i, MPFR_RNDU);
    mpfr_set_d (u, 1.0, MPFR_RNDD);
    mpfr_div (u, u, t, MPFR_RNDD);
    mpfr_add (s, s, u, MPFR_RNDD);
  }
printf ("Sum is ");
mpfr_out_str (stdout, 10, 0, s, MPFR_RNDD);
putchar ('\n');
mpfr_clear (s);
mpfr_clear (t);
mpfr_clear (u);
return 0;
],
[AC_MSG_RESULT([Non-trivial MPFR program successfully linked!])],
[AC_MSG_ERROR([Trouble linking non-trivial MPFR program!])]
)


AC_OUTPUT
rm -f conftest* config.cache

# Local Variables:
# mode: m4
# End: