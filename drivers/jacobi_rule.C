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

int main()
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

  for (unsigned int j=2; j<23; ++j)
    {
      std::cout << "================================================================================" << std::endl;
      std::cout << "Jacobi rule with alpha=" << alpha << ", beta=" << beta << ", "
                << j << " points, order=" << 2*j-1 << std::endl;
      jacobi_rule.rule(j); // order = 2*j-1

      // Scale Jacobi weights so they sum to 1/3 (alpha==2) or 1/2 (alpha==1)
      if (alpha==2.0)
        jacobi_rule.scale_weights(mpfr_class(1.0)/mpfr_class(3.0));

      else if (alpha==1.0)
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

      std::cout << "\n";
    }

  // jacobi_rule.rule(12); // order 23
  // jacobi_rule.rule(21); // Happened to need a slightly tighter tolerance...
  // jacobi_rule.rule(22); // order 43


  return 0;
}
