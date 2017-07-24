#!/bin/bash

# This script downloads, builds, and installs the GMP and MPFR
# libraries from source into the current directory.  You can then
# configure the mp-quadrature software to link against them by
# simply running:
#
# ./configure
#
# which will automatically search for GMP and MPFR in the ./gmp and
# ./mpfr directories, respectively.

# Top-level directory where we will install gmp and mpfr
install_dir=`pwd`

# We will build/install everything unless there is already something there.
go=1

# If the libraries are already there, warn before deleting them.
if ([ -f gmp/lib/libgmp.a ] || [ -f mpfr/lib/libmpfr.a ]); then
    echo "Warning: this will delete existing local GMP or MPFR libraries."
    read -p "Do you want to continue? " -r
    if [[ ! $REPLY =~ ^[Yy]$ ]]
    then
        go=0
    fi
fi

# If the user didn't say yes, quit.
if [ $go -eq 0 ] ; then
    echo "Skipping build."
    exit 0
fi

# 1.) Build GMP

# If the gmp/ directory doesn't exist, create it.
if [ ! -d gmp ]; then
  mkdir gmp
fi

# Change directories
pushd gmp

# Download the tarball if necessary
if [ ! -f gmp-5.1.3.tar.bz2 ]; then
    curl -O https://ftp.gnu.org/gnu/gmp/gmp-5.1.3.tar.bz2
fi

# Clean up anything left over from a previous installation.
rm -rf gmp-5.1.3 include lib share

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

# If the mpfr/ directory doesn't exist, create it.
if [ ! -d mpfr ]; then
  mkdir mpfr
fi

# Change directories
pushd mpfr

# Download the tarball
if [ ! -f mpfr-3.1.2.tar.bz2 ]; then
    curl -O http://www.mpfr.org/mpfr-3.1.2/mpfr-3.1.2.tar.bz2
fi

# Clean up anything left over from a previous installation.
rm -rf mpfr-3.1.2 include lib share

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
