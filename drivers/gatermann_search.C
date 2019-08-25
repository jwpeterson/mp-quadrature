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
  srandom(1566877491);

  // Use the current time since epoch as a seed.
  // time_t seed = time(nullptr);
  // std::cout << "seed=" << seed << std::endl;
  // srandom(seed);

  // The max degree of polynomials the rule is supposed to integrate
  // exactly (user input).
  // Invalid degrees are: 2, 4, 11, 13, ...
  unsigned int d=7;

  // Testing
  // std::cout << mpfr_class(exact_tri(0,0)) << std::endl; // =0.5
  // std::cout << mpfr_class(exact_tri(2,0)) << std::endl; // =0.0833333333333
  // std::cout << mpfr_class(exact_tri(3,0)) << std::endl; // =0.05
  // std::cout << mpfr_class(exact_tri(2,1)) << std::endl; // =0.0166666666667
  // std::cout << mpfr_class(exact_tri(4,0)) << std::endl; // =0.0333333333333
  // std::cout << mpfr_class(exact_tri(5,0)) << std::endl; // =0.0238095238095
  // std::cout << mpfr_class(exact_tri(4,1)) << std::endl; // =0.0047619047619
  // std::cout << mpfr_class(exact_tri(6,0)) << std::endl; // =0.0178571428571
  // std::cout << mpfr_class(exact_tri(5,1)) << std::endl; // =0.00297619047619
  // std::cout << mpfr_class(exact_tri(4,2)) << std::endl; // =0.00119047619048
  // std::cout << mpfr_class(exact_tri(7,0)) << std::endl; // =0.0138888888889
  // std::cout << mpfr_class(exact_tri(6,1)) << std::endl; // =0.00198412698413

  // Number of tests to run per test set (user input).
  unsigned int n_tests = 1;

  // Set parameters to be used by solvers
  SolverData solver_data;
  solver_data.verbose = true;
  solver_data.maxits = 150;
  solver_data.residual_and_jacobian = ResidualAndJacobian(d);
  solver_data.check_feasibility = CheckFeasibility(d);
  solver_data.do_backtracking = true;

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
