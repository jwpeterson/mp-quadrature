# 1. INTRODUCTION

   The accurate calculation and tabulation of high-order quadrature
   rules (eg Gauss-Legendre, Gauss-Jacobi, Gauss-Lobatto, etc.) is
   essential in many areas of numerical analysis. Standard
   double-precision arithmetic is typically only sufficient to obtain
   14 (or fewer) digits of accuracy in the points and weights, and
   therefore multiple-precision algebra libraries are required to
   improve this scenario. Furthermore, while the standard techniques
   for computing quadrature rules have been known for some time, some
   methods are better than others for computing arbitrary precision
   rules. Here we have collected a (hopefully growing) number of
   algorithms based on the freely-available GMP, MPFR, and GMPFRXX
   libraries for generating quadrature rules. This code was used to
   tabulate some of the 1D quadrature rules in the [LibMesh](https://github.com/libMesh/libmesh) finite element library.

# 2. INSTALLATION

   To build the library, type

   ```
   ./configure
   make
   ```

   You must have both the GMP and MPFR libraries installed in order to
   build the mp-quadrature library.  There are at least two options:

   1. Run the included `build_gmp_mpfr.sh` script.  This will download, build,
      and install GMP and MPFR from source into the `./gmp` and `./mpfr` directories.
      The configure script will then automatically find those.

   1. Specify the locations of your system's GMP and MPFR
      installations using the following options to configure:
      `--with-gmp-include=/path/to/gmp/include`, `--with-gmp-lib=/path/to/gmp/lib`, `--with-mpfr-include=/path/to/mpfr/include`, `--with-mpfr-lib=/path/to/mpfr/lib`.

      Configure will print an error message if some simple test codes
      involving these libraries fail to compile.


# 3. DRIVER PROGRAMS

   There are several driver programs (drivers/*.C) which make use of
   the library which gets built in lib/. run_tests.sh is a script
   which runs several of the driver programs to verify that the
   installation is working.

   1. `drivers/print_gauss.C`
      Computes Legendre polynomial roots and weights for 1D Gaussian
      quadrature rules.  To compute values for the rule with 10 points, run: `./drivers/print_gauss 10`.
   1. `drivers/jacobi_rule.C`
      Computes and prints Jacobi quadrature rules for alpha=1, beta=0,
      rescaled to the interval [0,1], for rules with 2 through 22 points
      (i.e. through order 43).
   1. `drivers/conical_product_2D.C` and `drivers/conical_product_3D.C`
      Computes and prints conical product rule points and weights for
      rules having n^2 (2D) or n^3 (3D) points.  Run e.g.

      `./drivers/conical_product_2D 3`

      or

      `./drivers/conical_product_3D 3`

      to get points and weights for a rule having 3\*3=9 (2D) or 3\*3\*3=27
      (3D) points.  The sum of weights is also printed for verification,
      it should be 0.5, the area of the reference triangle (2D) or
      0.1666..., the volume of the reference tetrahedron (3D).


# 4. NOTES ON OBTAINING SUPPORTING SOFTWARE:

The simplest approach (if it works!) is to run the included
`build_gmp_mpfr.sh` script.  This will download and install the GMP
and MPFR libraries from source.  If this does not work for some reason,
follow the directions below to build from source...

1. [GNU GMP](https://gmplib.org/)

   To build GMP from source:
   ```
   cd /where/you/want/to/build
   curl -O https://ftp.gnu.org/gnu/gmp/gmp-5.1.3.tar.bz2
   tar jxvf gmp-5.1.3.tar.bz2
   ./configure --prefix=/location/to/install/gmp --enable-cxx
   make -j4
   make -j4 check  # this worked just fine for me on Snow Leopard,
                   # and much later using the clang compiler on Mavericks
   sudo make install
   ```

1. [GNU MPFR](http://www.mpfr.org/)

   To build MPFR from source:
   ```
   cd /where/you/want/to/build
   curl -O http://www.mpfr.org/mpfr-current/mpfr-3.1.2.tar.bz2
   tar jxvf mpfr-3.1.2.tar.bz2
   cd mpfr-3.1.2
   ./configure --with-gmp-include=/location/to/install/gmp/include \
               --with-gmp-lib=/location/to/install/gmp/lib \
               --prefix=/location/to/install/mpfr
   make -j4
   make -j4 check
   sudo make install
   ```

# 5. MISCELLANEOUS

   There are also some more-or-less OK Matlab/Octave implementations in
   the matlab/ directory, though these are strictly double-precision
   implementations.  The Gauss Matlab implementation in particular is
   very simplistic and should not be relied on for accurate results at
   high orders!
