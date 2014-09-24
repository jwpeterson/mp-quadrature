#include <iostream>
#include <limits>
#include <cmath>
#include <cstdio>
#include <vector>
#include <iomanip>
#include <sstream>

// The GNU multi-precision library
#include "gmp.h"
// #include "gmpxx.h"

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
// #include "jacobi_constants.h"
#include "jacobi.h"

// Returns a string which contains the 32 decimal digits of x in the
// form: -9.4489927222288222340758013830322e-01L.  This should really
// be moved to its own header file, since both this file and the
// Gauss-Lobatto codes use it.
std::string fix_string(mpfr_class x);

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
        {
          mpfr_class one_third(1.0);
          one_third /= 3.0;
          jacobi_rule.scale_weights(one_third);
        }
      else if (alpha==1.0)
        {
          jacobi_rule.scale_weights(0.5);
        }
      else
        {
          std::cout << "Warning: weights unscaled!" << std::endl;
        }

      // Scale Jacobi points so they lie on [0, 1]
      mpfr_class zero(0.0), one(1.0);
      jacobi_rule.scale_points(zero, one);

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

      std::cout << "\n";
    }

  // jacobi_rule.rule(12); // order 23
  // jacobi_rule.rule(21); // Happened to need a slightly tighter tolerance...
  // jacobi_rule.rule(22); // order 43


  return 0;
}



std::string fix_string(mpfr_class x)
{
  // Go to some extra trouble to print e-01L instead of e-1
  std::ostringstream number_stream;
  number_stream.precision(32);
  number_stream.setf(std::ios_base::scientific);
  number_stream << x;
  std::string number = number_stream.str();

  // We almost always want to append an L (long double literal) but not always
  bool append_L = true;

  // Find the exponent 'e'.
  size_t e_pos = number.find('e');

  // If 'e' is not found...
  if (e_pos == std::string::npos)
    {
      if (number == "0")
        {
          // if number == "0", append a decimal point
          number.append(".");
          append_L = false;
        }
      else
        {
          // otherwise, append "e+00"
          number.append("e+00");
        }
    }
  else
    {
      // If 'e' is found, insert a leading zero to be consistent
      // with the Gauss rules in libmesh.  The number of digits
      // used in the exponent is implementation-dependent.  To be
      // totally general, we should check to see how many zeros
      // there are, and only insert exactly as many as we need...
      number.insert(e_pos+2, "0");
    }

  // The L suffix in C/C++ makes the compiler treat a number literal as 'long double'
  if (append_L)
    number.append("L");

  return number;
}
