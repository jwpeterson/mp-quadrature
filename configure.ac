# If you change any of the autotools files, run "autoreconf -if"
# to regenerate everyhing.

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

# Specify location of nlopt headers (optional)
AC_ARG_WITH(nlopt-include,
            AS_HELP_STRING([--with-nlopt-include=PATH],[Specify the path for nlopt header files]),
            withnloptinc=$withval,
            withnloptinc=)
AC_SUBST(withnloptinc)

# Specify location of nlopt libs (optional)
AC_ARG_WITH(nlopt-lib,
            AS_HELP_STRING([--with-nlopt-lib=PATH],[Specify the path for nlopt libraries]),
            withnloptlib=$withval,
            withnloptlib=)
AC_SUBST(withnloptlib)

# Location of the GMP installation.  Defaults to `pwd`/gmp, which is
# where the build_gmp_mpfr.sh script installs it.
AC_ARG_WITH(gmp-include,
            AS_HELP_STRING([--with-gmp-include=PATH],[Specify the path for GMP header files]),
            withgmpinc=$withval,
            withgmpinc=$MY_TOP_SRCDIR/gmp/include)
AC_SUBST(withgmpinc)

AC_ARG_WITH(gmp-lib,
            AS_HELP_STRING([--with-gmp-lib=PATH],[Specify the path for GMP libraries]),
            withgmplib=$withval,
            withgmplib=$MY_TOP_SRCDIR/gmp/lib)
AC_SUBST(withgmplib)


# Location of the MPFR installation.  Defaults to `pwd`/mpfr, which is
# where the build_gmp_mpfr.sh script installs it.
AC_ARG_WITH(mpfr-include,
            AS_HELP_STRING([--with-mpfr-include=PATH],[Specify the path for MPFR header files]),
            withmpfrinc=$withval,
            withmpfrinc=$MY_TOP_SRCDIR/mpfr/include)
AC_SUBST(withmpfrinc)

AC_ARG_WITH(mpfr-lib,
            AS_HELP_STRING([--with-mpfr-lib=PATH],[Specify the path for MPFR libraries]),
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
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
@%:@include <gmp.h>
]], [[
]])],[AC_MSG_RESULT([gmp.h found in $withgmpinc])],[AC_MSG_ERROR([GNU MP not found, see http://gmplib.org/])
])

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

AC_LINK_IFELSE([AC_LANG_PROGRAM([[
@%:@include <stdio.h>
@%:@include <stdlib.h>
@%:@include <gmp.h>
]], [[
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
]])],[AC_MSG_RESULT([Non-trivial GMP program successfully linked!])],[AC_MSG_ERROR([Trouble linking non-trivial GMP program!])
])



#
# Testing MPFR
#

# More configure testing, this time for MPFR.  Note that these test
# programs will be compiled with GMP include paths as well...
CFLAGS="$CFLAGS -I$withmpfrinc"

# Test compiling a non-trivial program that includes mpfr.h
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
@%:@include <mpfr.h>
]], [[
]])],[AC_MSG_RESULT([mpfr.h found in $withmpfrinc])],[AC_MSG_ERROR([GNU MPFR not found, see http://www.mpfr.org/])
])


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

AC_LINK_IFELSE([AC_LANG_PROGRAM([[
@%:@include <stdio.h>
@%:@include <gmp.h>
@%:@include <mpfr.h>
]], [[
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
]])],[AC_MSG_RESULT([Non-trivial MPFR program successfully linked!])],[AC_MSG_ERROR([Trouble linking non-trivial MPFR program!])
])


# Test for working nlopt installation
have_nlopt=no
NLOPT_INCLUDE=
NLOPT_LIB=

dnl Don't enable nlopt unless the user specified both include and lib locations.
AS_IF([test "x$withnloptinc" != "x" && test "x$withnloptlib" != "x"],
      [
      dnl Set potential NLOPT_INCLUDE and NLOPT_LIB variables.
      NLOPT_INCLUDE="-I$withnloptinc"
      NLOPT_LIB="-L$withnloptlib -lnlopt -lm"

      dnl Save original flags/libs so we can restore them afterward.
      CFLAGS_SAVED="$CFLAGS"
      LIBS_SAVED="$LIBS"
      CFLAGS="$CFLAGS $NLOPT_INCLUDE"
      LIBS="$LIBS $NLOPT_LIB"

      dnl try linking simple test program using C compiler
      AC_LINK_IFELSE(
      [
        AC_LANG_PROGRAM(
        [[
          @%:@include <math.h>
          @%:@include <nlopt.h>
          @%:@include <stdio.h>

          double myfunc(unsigned n, const double *x, double *grad, void *my_func_data)
          {
          if (grad) {
            grad[0] = 0.0;
            grad[1] = 0.5 / sqrt(x[1]);
           }
           return sqrt(x[1]);
          }

          typedef struct
          {
            double a, b;
          } my_constraint_data;

          double myconstraint(unsigned n, const double *x, double *grad, void *data)
          {
            my_constraint_data *d = (my_constraint_data *) data;
            double a = d->a, b = d->b;
            if (grad) {
              grad[0] = 3 * a * (a*x[0] + b) * (a*x[0] + b);
              grad[1] = -1.0;
            }
            return ((a*x[0] + b) * (a*x[0] + b) * (a*x[0] + b) - x[1]);
          }
        ]],
        [[
          double lb[2] = { -HUGE_VAL, 0 }; /* lower bounds */
          nlopt_opt opt;

          opt = nlopt_create(NLOPT_LD_MMA, 2); /* algorithm and dimensionality */
          nlopt_set_lower_bounds(opt, lb);
          nlopt_set_min_objective(opt, myfunc, NULL);

          my_constraint_data data[2] = { {2,0}, {-1,1} };
          nlopt_add_inequality_constraint(opt, myconstraint, &data[0], 1e-8);
          nlopt_add_inequality_constraint(opt, myconstraint, &data[1], 1e-8);
          nlopt_set_xtol_rel(opt, 1e-4);
          double x[2] = { 1.234, 5.678 };  /* some initial guess */
          double minf; /* the minimum objective value, upon return */

          if (nlopt_optimize(opt, x, &minf) < 0) {
            printf("nlopt failed!\n");
          }
          else {
            printf("found minimum at f(%g,%g) = %0.10g\n", x[0], x[1], minf);
          }
          nlopt_destroy(opt);
          return 0;
        ]])
      ],
      [
        AC_MSG_RESULT([nlopt test program successfully linked!])
        have_nlopt=yes
      ],
      [
        AC_MSG_RESULT([Trouble compiling/linking nlopt program, nlopt disabled])
        have_nlopt=no
        NLOPT_INCLUDE=
        NLOPT_LIB=
      ])

      dnl Restore original flags
      CFLAGS="$CFLAGS_SAVED"
      LIBS="$LIBS_SAVED"
      ])

AC_SUBST(have_nlopt)
AC_SUBST(NLOPT_INCLUDE)
AC_SUBST(NLOPT_LIB)



AC_OUTPUT
rm -f conftest* config.cache

# Local Variables:
# mode: autoconf
# End:
