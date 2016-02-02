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

if [ ! -d gmp ]; then
  mkdir gmp

  # Change directories
  pushd gmp

  # Download the tarball if necessary
  if [ ! -f gmp-5.1.3.tar.bz2 ]; then
    curl -O https://ftp.gnu.org/gnu/gmp/gmp-5.1.3.tar.bz2
  fi

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

# If the gmp/ directory already exists, notify the user we're skipping building it.
else
    printf '\n'
    printf 'Warning! gmp directory already exists, skipping gmp build.\n'
    printf '(If you really want to rebuild gmp, please mv or rm the gmp directory.)\n'
    printf '\n'
fi



# 2.) Build MPFR

# If the mpfr/ directory doesn't exist, create it.
if [ ! -d mpfr ]; then
  mkdir mpfr

  # Change directories
  pushd mpfr

  # Download the tarball
  if [ ! -f mpfr-3.1.2.tar.bz2 ]; then
    curl -O http://www.mpfr.org/mpfr-3.1.2/mpfr-3.1.2.tar.bz2
  fi

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

# If the mpfr/ directory already exists, notify the user we're skipping building it.
else
    printf '\n'
    printf 'Warning! mpfr directory already exists, skipping mpfr build.\n'
    printf '(If you really want to rebuild mpfr, please mv or rm the mpfr directory.)\n'
    printf '\n'
fi
