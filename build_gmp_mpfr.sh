#!/bin/bash

# This script downloads, builds, and installs the GMP and MPFR
# libraries from source into the current directory.  You can then
# configure the invariant quadrature software to link against them by
# simply running:
#
# ./configure
#
# which will automatically search for GMP and MPFR in the ./gmp and
# ./mpfr directories, respectively.

# Top-level directory where we will install gmp and mpfr
install_dir=`pwd`

# 1.) Build GMP

# Start from a clean slate... don't try to deal with failed builds, etc.
rm -rf gmp
mkdir gmp
pushd gmp

# Download the tarball
curl -O https://ftp.gnu.org/gnu/gmp/gmp-5.1.3.tar.bz2

# *or* Debugging: copy over a previously-downloaded tarball
# cp /Users/petejw/software/gmp/5.1.3/gmp-5.1.3.tar.bz2 .

# Unpack the tarball
tar jxvf gmp-5.1.3.tar.bz2
pushd gmp-5.1.3

# Configure GMP
./configure --prefix=${install_dir}/gmp --enable-cxx

# Build and install GMP
make -j4
make install

# Go back
popd
popd



# 2.) Build MPFR

# Start from a clean slate... don't try to deal with failed builds, etc.
rm -rf mpfr
mkdir mpfr
pushd mpfr

# Download the tarball
curl -O http://www.mpfr.org/mpfr-current/mpfr-3.1.2.tar.bz2

# Unpack the tarball
tar jxvf mpfr-3.1.2.tar.bz2
pushd mpfr-3.1.2

# Configure using the GMP we just built
./configure --prefix=${install_dir}/mpfr --with-gmp-include=${install_dir}/gmp/include --with-gmp-lib=${install_dir}/gmp/lib

# Build and install MPFR
make -j4
make install

# Go back
popd
popd
