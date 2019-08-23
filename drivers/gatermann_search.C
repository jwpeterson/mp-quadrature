// Functions for computing exact integrals of monomials
// on the reference triangle.
#include "exact.h"
#include "matrix.h"
#include "vect.h"
#include "solver_data.h"
#include "test_ro3.h"

// C++ includes
#include <vector>
#include <utility>
#include <iostream>
#include <algorithm> // std::sort
#include <numeric> // std::accumulate
#include <stdlib.h> // random
#include <time.h> // time()
#include <assert.h>

// Our goal here is to determine whether there might be more than one
// solution to the SEVENTH-order Ro3-invariant quadrature rule
// originally reported in Gatermann, 1988.  We are going to try
// several different sets of starting data and see if any of them
// converge to a different solution than the one reported in the
// paper.
//
// For more information, see gatermann.C
int main()
{
  std::cout.precision(32);
  std::cout.setf(std::ios_base::scientific);

  // # of binary digits
  // 53 binary digits is about what you normally get with a double.
  mpfr_set_default_prec(256);

  // Varies the sequence of pseudorandom numbers returned by random().
  // srandom(1566449865);

  // Use the current time since epoch as a seed.
  time_t seed = time(nullptr);
  std::cout << "seed=" << seed << std::endl;
  srandom(seed);

  // The max degree of polynomials the rule is supposed to integrate
  // exactly (user input).
  // Invalid degrees are: 2, 4, 11, 13, ...
  unsigned int d=15;

  // Number of tests to run per test set (user input).
  unsigned int n_tests = 1;

  // Set parameters to be used by solvers
  SolverData solver_data;
  solver_data.verbose = false;
  solver_data.maxits = 15;
  solver_data.residual_and_jacobian = ResidualAndJacobian(d);

  unsigned int testset_counter = 0;
  while (true)
    {
      ++testset_counter;
      std::cout << "Running test set " << testset_counter << std::endl;
      test_ro3(d, n_tests, solver_data);
      // Comment out this break to run forever.
      break;
    }

  return 0;
}
