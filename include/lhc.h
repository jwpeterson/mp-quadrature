#ifndef LHC_H
#define LHC_H

// C++ includes
#include <vector>

// Generates an n_dim * n_tests array of random numbers in [0,1]
// using the Latin Hypercube (LHC) sampling method. Note that double
// precision numbers are used since it is assumed that arbitrary
// precision is not important when setting initial conditions.
std::vector<std::vector<double>>
lhc(unsigned int n_dim, unsigned int n_tests);

#endif
