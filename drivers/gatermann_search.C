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

// Function used to (possibly) compute the residual and Jacobian at
// the input u.  If either "r" or "jac" is nullptr, their computation
// is skipped.
void residual_and_jacobian(std::vector<mpfr_class> * r,
                           Matrix<mpfr_class> * jac,
                           const std::vector<mpfr_class> & u);

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

  // Set parameters to be used by solvers
  SolverData solver_data;
  solver_data.verbose = false;
  solver_data.residual_and_jacobian = residual_and_jacobian;

  // The max degree of polynomials the rule is supposed to integrate
  // exactly (input).
  unsigned int d=7;

  unsigned int testset_counter = 0;
  while (true)
    {
      ++testset_counter;
      std::cout << "Running test set " << testset_counter << std::endl;
      test_ro3(d, solver_data);
      break;
    }

  return 0;
}



void
residual_and_jacobian(std::vector<mpfr_class> * r,
                      Matrix<mpfr_class> * jac,
                      const std::vector<mpfr_class> & u)
{
  // If there's nothing to do, then there's nothing to do.
  if (!r && !jac)
    return;

  // The list of monomial exponents that we are going to test.
  std::vector<std::pair<unsigned int, unsigned int>> polys =
    {
      {0,0}, {2,0}, {3,0}, {2,1}, {4,0}, {5,0},
      {4,1}, {6,0}, {5,1}, {7,0}, {6,1}, {4,2}
    };

  // The number of equations (and the number of unknowns).
  unsigned int n = polys.size();

  // Zero any previous values and allocate space for residual and Jacobian, if required.
  if (r)
    {
      r->clear();
      r->resize(n);
    }
  if (jac)
    {
      jac->clear();
      jac->resize(n, n);
    }

  // For each basis function i, we loop over all (w,x,y) triples in u and
  // compute the residual and Jacobian contributions, as requested.
  for (unsigned int i=0; i<n; i++)
    {
      // The powers of x and y in the monomial we are currently integrating.
      unsigned int xpower = polys[i].first;
      unsigned int ypower = polys[i].second;

      // Accumlate quadrature rule contribution for r(i).
      for (unsigned int q=0; q<u.size(); q+=3)
        {
          // The unknowns are ordered in terms of (w,x,y) triples.
          mpfr_class w = u[q];
          mpfr_class x = u[q+1];
          mpfr_class y = u[q+2];

          // The implied third barycentric coordinate
          mpfr_class z = mpfr_class(1) - x - y;

          // The "spatial" part is needed by both the residual and Jacobian,
          // so we can always compute it.
          mpfr_class spatial =
            pow(x, xpower) * pow(y, ypower) +
            pow(z, xpower) * pow(x, ypower) +
            pow(y, xpower) * pow(z, ypower);

          // Compute residual contribution, if required.
          if (r)
            (*r)[i] += w * spatial;

          // Compute Jacobian contribution, if required.
          if (jac)
            {
              // Derivative wrt w
              (*jac)(i, q) += spatial;

              // We are differentiating polynomials, so if xpower or
              // ypower is 0, the derivative of the corresponding term
              // will be zero. In that case we don't want to subtract
              // 1 from the unsigned variable which represents the
              // exponent, since that would cause it to wrap, possibly
              // leading to other issues. Therefore we will just
              // define the power to still be zero in that case...
              // Note that the -1 terms come from differentiating "z"
              // wrt either x or y.
              unsigned int xpm1 = xpower > 0 ? xpower-1 : 0;
              unsigned int ypm1 = ypower > 0 ? ypower-1 : 0;

              // Derivative wrt x
              (*jac)(i,q+1) += w * (
                xpower * pow(x, xpm1) * pow(y, ypower) +
                (mpfr_class(-1) * xpower * pow(z, xpm1) * pow(x, ypower) + pow(z, xpower) * ypower * pow(x, ypm1)) +
                pow(y, xpower) * mpfr_class(-1) * ypower * pow(z, ypm1));

              // Derivative wrt y
              (*jac)(i,q+2) += w * (
                pow(x, xpower) * ypower * pow(y, ypm1) +
                mpfr_class(-1) * xpower * pow(z, xpm1) * pow(x, ypower) +
                (xpower * pow(y, xpm1) * pow(z, ypower) + pow(y, xpower) * mpfr_class(-1) * ypower * pow(z, ypm1)));
            }
        }

      // Subtract off the true integral value, I(p_i)
      if (r)
        (*r)[i] -= exact_tri(xpower, ypower);
    }
}
