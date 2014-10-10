#include <cstdlib> // atoi, std::abort
#include <iostream>
#include <limits>
#include <cmath>
#include <cstdio>
#include <vector>
#include <iomanip>

// The GNU multi-precision library
#include "gmp.h"
#include "mpfr.h"
#include "gmpfrxx.h"

// Header files for this project
#include "common_definitions.h"
#include "gauss_lobatto.h"

// This program prints the points and weights for 1D Gauss-Lobatto rules.
int main(int argc, char** argv)
{
  std::cout.precision(32);
  std::cout.setf(std::ios_base::scientific);

  // # of binary digits
  // 53 binary digits is about what you normally get with a double.
  mpfr_set_default_prec(256);

  // Read number of points in rule from command line
  unsigned int n=6;
  if (argc > 1)
    n = atoi(argv[1]);

  // Don't allow n==0, 1
  if (n<2)
    {
      std::cout << "n>=2 is required, running with default 6 point rule." << std::endl;
      n = 6;
    }

  std::cout << "\nComputing Gauss-Lobatto rule with " << n << " points" <<  std::endl;
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

  // Verify that polynomials up to the supported order are integrated
  // exactly (to within the arbitrary precision we are using).  Be
  // careful with the 1-based indexing we are using...
  std::cout << "\nVerifying rule:" << std::endl;
  unsigned max_order = 2*n - 3;
  for (unsigned order=0; order <= max_order; ++order)
    {
      mpfr_class sum = 0.;
      for (unsigned n_qp=1; n_qp<x.size(); ++n_qp)
        {
          sum += w[n_qp] * pow(x[n_qp], order);
        }

      // The exact solution to int(x^p, x=-1, 1) is:
      // x^{p+1} / (p+1) |_{-1}^{+1}
      // = (1 - (-1)^{p+1}) / (p+1)
      // = { 2/(p+1), p even
      //   { 0,       p odd
      mpfr_class exact = (order%2==0) ? mpfr_class(2.)/mpfr_class(order+1) : mpfr_class(0.);

      // std::cout << "quadrature = " << sum << std::endl;
      // std::cout << "exact      = " << exact << std::endl;

      // Compute the absolute error:
      mpfr_class abs_err = abs(sum-exact);

      // Print message
      std::cout << "Computing int(x^" << order << ", x=-1..1)"
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

  return 0;
}
