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

# Versions of the libraries we will use. These will be updated periodically.
gmp_version=6.1.2
mpfr_version=4.0.1

# Check whether we already have the required tarball versions.
have_gmp_tarball=no
have_mpfr_tarball=no
if test -f gmp/gmp-${gmp_version}.tar.bz2; then
    have_gmp_tarball=yes
fi
if test -f mpfr/mpfr-${mpfr_version}.tar.bz2; then
    have_mpfr_tarball=yes
fi

# If the required tarballs are not already available, we need to
# choose a command line downloader. By default we use curl, falling
# back on wget if it's not available, and throwing an error if neither
# is available.
if test $have_gmp_tarball = no || test $have_mpfr_tarball = no; then
    if type -P "curl" &>/dev/null; then
        echo "Using curl to download tarballs"
        cld=curl
        cld_options=-O
    elif type -P "wget" &>/dev/null; then
        echo "Using wget to download tarballs"
        cld=wget
        cld_options=
    else
        echo "Either curl or wget is required."
        exit 0
    fi
fi

# 1.) Build GMP

# If the gmp/ directory doesn't exist, create it.
if [ ! -d gmp ]; then
  mkdir gmp
fi

# Change directories
pushd gmp

# Download the tarball if necessary
if test $have_gmp_tarball = no; then
    $cld $cld_options https://ftp.gnu.org/gnu/gmp/gmp-${gmp_version}.tar.bz2
fi

# Clean up anything left over from a previous installation.
rm -rf gmp-${gmp_version} include lib share

# Unpack the tarball
tar jxvf gmp-${gmp_version}.tar.bz2
pushd gmp-${gmp_version}

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
if test $have_mpfr_tarball = no; then
    $cld $cld_options https://www.mpfr.org/mpfr-${mpfr_version}/mpfr-${mpfr_version}.tar.bz2
fi

# Clean up anything left over from a previous installation.
rm -rf mpfr-${mpfr_version} include lib share

# Unpack the tarball
tar jxvf mpfr-${mpfr_version}.tar.bz2
pushd mpfr-${mpfr_version}

# Configure using the GMP we just built
./configure --prefix=${install_dir}/mpfr --with-gmp-include=${install_dir}/gmp/include --with-gmp-lib=${install_dir}/gmp/lib

# Build and install MPFR
make -j4
make install

# Go back
popd
popd
