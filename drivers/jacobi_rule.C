#include <iostream>
#include <limits>
#include <cmath>
#include <cstdio>
#include <vector>
#include <iomanip>
#include <sstream>

// The GNU multi-precision library
#include "gmp.h"

// The MPFR library defines special functions
// like sin, cos, exp, etc.
#include "mpfr.h"

// The GMPFRXX library defines a C++ interface for
// the mpfr library.  You cannot include the gmpxx.h
// header *and* the gmpfrxx.h headers at the same time
// though!
#include "gmpfrxx.h"

// Header files for this project
#include "common_definitions.h"
#include "jacobi.h"

int main(int argc, char** argv)
{
  std::cout.precision(32);
  std::cout.setf(std::ios_base::scientific);

  // # of binary digits
  // 53 binary digits is about what you normally get with a double.
  mpfr_set_default_prec(256);

  // Case 1: weights sum to 1/2
  const Real alpha=1.0, beta=0.0;

  // Case 2: weights sum to 1/3
  // const Real alpha=2.0, beta=0.0;

  Jacobi jacobi_rule(alpha, beta);

  // Read number of points in rule from command line
  unsigned int n=6;
  if (argc > 1)
    n = atoi(argv[1]);

  // Make sure n is valid
  if (n==0)
    {
      std::cout << "Warning, could not determine valid rule order from command line." << std::endl;
      std::cout << "Running with default 6-point rule." << std::endl;
      n = 6;
    }

      std::cout << "Jacobi rule with alpha=" << alpha << ", beta=" << beta << ", "
                << n << " points, order=" << 2*n-1 << std::endl;

      // Valid for polynomials (not including the weighting function) of order = 2*n-1
      jacobi_rule.rule(n);

      // Scale Jacobi weights so they sum to 1/3 (alpha==2) or 1/2 (alpha==1)
      if (alpha == 2.0)
        jacobi_rule.scale_weights(mpfr_class(1.0)/mpfr_class(3.0));

      else if (alpha == 1.0)
        jacobi_rule.scale_weights(0.5);

      else
        std::cout << "Warning: weights unscaled!" << std::endl;

      // Scale Jacobi points so they lie on [0, 1]
      jacobi_rule.scale_points(mpfr_class(0.0), mpfr_class(1.0));

      // Print the result

      // The points and weights arrays are 1-based
      const std::vector<mpfr_class>& x = jacobi_rule.get_points();
      const std::vector<mpfr_class>& w = jacobi_rule.get_weights();

      for (unsigned i=1; i<x.size(); ++i)
        std::cout << "_points[" << std::setw(2) << i-1 << "](0) = " << fix_string(x[i]) << ";\n";

      // blank line
      std::cout << std::endl;

      for (unsigned i=1; i<w.size(); ++i)
        std::cout << "_weights[" << std::setw(2) << i-1 << "]   = " << fix_string(w[i]) << ";\n";

      std::cout << "\nVerifying rule:" << std::endl;

      mpfr_class sumweights = 0.;
      for (unsigned i=1; i<w.size(); ++i)
        sumweights += w[i];
      std::cout << "Sum of weights is: " << sumweights << std::endl;

      unsigned max_order = 2*n - 1;
      for (unsigned order=0; order <= max_order; ++order)
        {
          mpfr_class sum = 0.;
          for (unsigned n_qp=1; n_qp<x.size(); ++n_qp)
            sum += w[n_qp] * pow(x[n_qp], order);

          // Exact solutions for alpha=1:
          //
          // We use the scaled interval [0,1] here, so the exact
          // solution for alpha=1 to int((1-x)*x^p, x=0..1) is given
          // by:
          // 1 / (p^2 + 3p + 2)
          //
          // If one instead uses the interval [-1,1]:
          // int((1-x)*x^p, x=-1..1), the exact solution is given by:
          //        p           p
          //  2 (-1)  p + 3 (-1)  + 1
          //  -----------------------
          //      (p + 2) (p + 1)
          //
          // = { ( 2*p + 4) / (p + 2) / (p + 1), p even
          //   { (-2*p - 2) / (p + 2) / (p + 1), p odd
          mpfr_class exact = mpfr_class(1.0)/mpfr_class(order*order + 3*order + 2);

          // std::cout << "quadrature = " << sum << std::endl;
          // std::cout << "exact      = " << exact << std::endl;

          // Compute the absolute error:
          mpfr_class abs_err = abs(sum - exact);

          // Print message
          std::cout << "Computing int((1-x)^" << static_cast<unsigned>(alpha) << " * x^" << order << ", x=0..1)"
                    << ", abs_err = " << abs_err << std::endl;

          // Abort if error is too large.  Most of the results are
          // accurate to a tighter tolerance than 1.e-30, but 1.e-30 at
          // least guarantees that all the tests pass.
          if (abs_err > mpfr_class(1.e-30))
            {
              std::cerr << "Quadrature error too large, possible problem with points and weights!" << std::endl;
              std::abort();
            }
        }

      std::cout << "\n";

  return 0;
}
