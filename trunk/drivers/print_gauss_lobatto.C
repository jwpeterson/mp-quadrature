#include <stdlib.h> // atoi
#include <iostream>
#include <limits>
#include <cmath>
#include <cstdio>
#include <vector>
#include <iomanip>
#include <sstream>

// The GNU multi-precision library
#include "gmp.h"
#include "mpfr.h"
#include "gmpfrxx.h"

// Header files for this project
#include "common_definitions.h"
#include "gauss.h"

// Returns a string which contains the 32 decimal digits of x in the form:
// -9.4489927222288222340758013830322e-01L
std::string fix_string(mpfr_class x);

// This program prints the points and weights for 1D Gauss-Lobatto rules.
int main(int argc, char** argv)
{
  std::cout.precision(32);
  std::cout.setf(std::ios_base::scientific);

  // # of binary digits
  // 53 binary digits is about what you normally get with a double.
  mpfr_set_default_prec(256);

  unsigned int n=6;
  // Read number of points in rule from command line
  if (argc > 1)
    n = atoi(argv[1]);

  // Don't allow n==0, 1
  if (n<2)
    {
      std::cout << "n>=2 is required, running with default order 6 rule." << std::endl;
      n = 6;
    }

  std::cout << "\nComputing Gauss-Lobatto rule for n=" << n << std::endl;
  std::vector<mpfr_class> x, w;
  gauss_lobatto_rule(n, x, w);

  // For certain rules we have analytical expressions for the points
  // and weights.
  bool analytical_check = ((n>=2) && (n<=5));

  if (analytical_check)
    {
      std::vector<mpfr_class> x_analytical, w_analytical;

      switch (n)
        {
        case 2:
          {
            // x= +/-1
            x_analytical.resize(2);
            x_analytical[0] = -1.0;
            x_analytical[1] =  1.0;

            // w=1 for this case
            w_analytical.resize(2);
            w_analytical[0] = 1.0;
            w_analytical[1] = 1.0;
            break;
          }

        case 3:
          {
            // x = 0, +/-1
            x_analytical.resize(3);
            x_analytical[0] = -1.0;
            x_analytical[1] =  0.0;
            x_analytical[2] =  1.0;

            // w = 1/3, 4/3, 1/3
            w_analytical.resize(3);
            w_analytical[0] = mpfr_class(1.0)/mpfr_class(3.0);
            w_analytical[1] = mpfr_class(4.0)/mpfr_class(3.0);
            w_analytical[2] = w_analytical[0];

            break;
          }

        case 4:
          {
            // x = +/-1, +/-1/sqrt(5)
            x_analytical.resize(4);
            x_analytical[0] = -1.0;
            x_analytical[1] = -1.0 / sqrt(mpfr_class(5.0));
            x_analytical[2] = -1.0 * x_analytical[1];
            x_analytical[3] = -1.0 * x_analytical[0];

            // w = 1/6, 5/6
            w_analytical.resize(4);
            w_analytical[0] = mpfr_class(1.0)/mpfr_class(6.0);
            w_analytical[1] = mpfr_class(5.0)/mpfr_class(6.0);
            w_analytical[2] = w_analytical[1];
            w_analytical[3] = w_analytical[0];

            break;
          }

        case 5:
          {
            // x = +/-1, +/-sqrt(3/7), 0
            x_analytical.resize(5);
            x_analytical[0] = -1.0;
            x_analytical[1] = -sqrt(mpfr_class(3.0)/mpfr_class(7.0));
            x_analytical[2] = 0.0;
            x_analytical[3] = -1.0*x_analytical[1];
            x_analytical[4] = -1.0*x_analytical[0];

            // w = 1/10, 49/90, 32/45
            w_analytical.resize(5);
            w_analytical[0] = mpfr_class(1.0)/mpfr_class(10.0);
            w_analytical[1] = mpfr_class(49.0)/mpfr_class(90.0);
            w_analytical[2] = mpfr_class(32.0)/mpfr_class(45.0);
            w_analytical[3] = w_analytical[1];
            w_analytical[4] = w_analytical[0];

            break;
          }

        default:
          {
            // No analytical form known
          }
        }

      // Compute error in x, w.  Note, the values returned in x, w are 1-based.
      if ((x.size()-1 == x_analytical.size()) && (w.size()-1 == w_analytical.size()))
        {
          std::cout << "\nError in points and weights: " << std::endl;

          for (unsigned int i=0; i<x_analytical.size(); ++i)
            {
              mpfr_class delta_x = x[i+1] - x_analytical[i];
              mpfr_class delta_w = w[i+1] - w_analytical[i];

              std::cout << "delta_x[" << i << "]=" << delta_x << ", "
                        << "delta_w[" << i << "]=" << delta_w << std::endl;
            }
          std::cout << std::endl;
        }
      else
        {
          std::cout << "Different number of points computed for numerical and analytical rules!" << std::endl;
          std::cout << "x.size()=" << x.size() << ", x_analytical.size()=" << x_analytical.size() << std::endl;
        }

      // Print analytical points/weights values
      // for (unsigned int i=0; i<x_analytical.size(); ++i)
      //   std::cout << "x_analytical["<<i<<"]=" << x_analytical[i] << std::endl;
      //
      // std::cout << std::endl;
      //
      // for (unsigned int i=0; i<w_analytical.size(); ++i)
      //   std::cout << "w_analytical["<<i<<"]=" << w_analytical[i] << std::endl;
    } // end if analytical_check

  // Print points and weights using standard C-style indexing and
  // by setting the second half of the points and weights equal
  // to the mirror of the first half.
  const unsigned int m=(n+1)/2;

  // Points first, then weights.  This is how libmesh displays the rules.
  std::cout << std::endl;

  for (unsigned int i=1; i<=m; ++i)
    std::cout << "_points[" << std::setw(2) << i-1 << "](0) = " << fix_string(x[i]) << ";\n";

  for (unsigned int i=m+1,d=(n%2)?m-1:m; i<=n; ++i,--d)
    std::cout << "_points[" << std::setw(2) << i-1 << "]    = -_points[" << d-1 << "];\n";

  std::cout << std::endl;

  for (unsigned int i=1; i<=m; ++i)
    std::cout << "_weights[" << std::setw(2) << i-1 << "]   = " << fix_string(w[i]) << ";\n";

  for (unsigned int i=m+1,d=(n%2)?m-1:m; i<=n; ++i,--d)
    std::cout << "_weights[" << std::setw(2) << i-1 << "]   = _weights[" << d-1 << "];" << std::endl;

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
