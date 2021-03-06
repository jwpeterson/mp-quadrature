#!/bin/bash

# Commands I ran to install all support software and build mp-quadrature.
# This assumes you have already cloned the mp-quadrature package into
# ~/projects/mp-quadrature and are sitting in that directory.
sudo apt-get update
sudo apt-get -y install libmpfr-dev
sudo apt-get -y install g++
sudo apt-get -y install libnlopt-dev
sudo apt-get -y install libnlopt0
sudo apt-get -y install make
sudo apt-get -y install emacs

# mkdir projects
# cd projects
# git clone https://github.com/jwpeterson/mp-quadrature.git
# cd mp-quadrature
./configure --with-nlopt-include=/usr/include --with-nlopt-lib=/usr/lib/x86_64-linux-gnu
make
