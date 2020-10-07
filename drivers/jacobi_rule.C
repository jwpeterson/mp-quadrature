// mp-quadrature includes
#include "common_definitions.h"
#include "jacobi.h"

// The GMPFRXX library defines a C++ interface for
// the mpfr library.  You cannot include the gmpxx.h
// header *and* the gmpfrxx.h headers at the same time
// though!
#include "gmpfrxx.h"

// The MPFR library defines special functions
// like sin, cos, exp, etc.
#include "mpfr.h"

// The GNU multi-precision library
#include "gmp.h"

// C++ includes
#include <cstdlib> // std::abort
#include <iostream>
#include <limits>
#include <cmath>
#include <cstdio>
#include <vector>
#include <iomanip>
#include <sstream>
#include <getopt.h> // getopt_long()

// Print a simple usage message when unrecognized command line arguments are encountered
void usage();

int main(int argc, char** argv)
{
  std::cout.precision(32);
  std::cout.setf(std::ios_base::scientific);

  // # of binary digits
  // 53 binary digits is about what you normally get with a double.
  mpfr_set_default_prec(256);

  // Case 1: alpha=1, beta=0, weights sum to 1/2
  // Case 2: alpha=2, beta=0, weights sum to 1/3
  unsigned
    alpha = 1,
    beta = 0;

  // Number of points in rule, to be read from command line
  unsigned int n = 6;

  // options descriptor - this associates the following "long" and "short" command line options
  // --alpha, -a
  // --beta, -b
  // --npoints, -n
  static struct option longopts[] =
    {
      {"alpha",   required_argument, NULL, 'a'},
      {"beta",    required_argument, NULL, 'b'},
      {"npoints", required_argument, NULL, 'n'},
      {"help",    no_argument,       NULL, 'h'},
      { NULL,     0,                 NULL,  0 }
    };

  // Parse command line options using getopt_long()
  int ch = -1;
  while ((ch = getopt_long(argc, argv, "ha:b:n:", longopts, NULL)) != -1)
    {
      switch (ch)
        {
        case 'a':
          {
            alpha = atoi(optarg);
            if (!(alpha == 1 || alpha == 2))
              {
                std::cerr << "Error, only alpha=1, 2 are supported, but alpha = " << alpha << " was provided!" << std::endl;
                std::abort();
              }

            break;
          }

        case 'b':
          {
            beta = atoi(optarg);

            if (beta != 0)
              {
                std::cerr << "Error, only beta=0 is supported, but beta = " << beta << " was provided!" << std::endl;
                std::abort();
              }
            break;
          }

        case 'n':
          {
            n = atoi(optarg);

            // Make sure n is valid
            if (n==0)
              {
                std::cout << "Warning, could not determine valid rule order from command line." << std::endl;
                std::cout << "Running with default 6-point rule." << std::endl;
                n = 6;
              }
            break;
          }

        case 'h':
          {
            usage();
            return 0;
          }

        default:
          // We could error here, print a usage command, or just ignore unrecognized arguments...
          usage();
        }
    } // end while

  Jacobi jacobi_rule(alpha, beta);

  std::cout << "Jacobi rule with alpha=" << alpha << ", beta=" << beta << ", "
            << n << " points, order=" << 2*n-1 << std::endl;

  // Valid for polynomials (not including the weighting function) of order = 2*n-1
  jacobi_rule.rule(n);

  // Get a reference to the 1-based points and weights arrays
  const std::vector<mpfr_class>& x = jacobi_rule.get_points();
  const std::vector<mpfr_class>& w = jacobi_rule.get_weights();

  // Scale Jacobi weights so they sum to 1/3 (alpha==2) or 1/2 (alpha==1)
  if (alpha == 2)
    jacobi_rule.scale_weights(mpfr_class(1.0)/mpfr_class(3.0));

  else if (alpha == 1)
    jacobi_rule.scale_weights(0.5);

  else
    std::cout << "Warning: weights unscaled!" << std::endl;

  // Scale Jacobi points so they lie on [0, 1]
  jacobi_rule.scale_points(mpfr_class(0.0), mpfr_class(1.0));

  // Print the result
  for (unsigned i=1; i<x.size(); ++i)
    std::cout << "_points[" << std::setw(2) << i-1 << "](0) = " << fix_string(x[i]) << ";\n";

  // blank line
  std::cout << std::endl;

  for (unsigned i=1; i<w.size(); ++i)
    std::cout << "_weights[" << std::setw(2) << i-1 << "]   = " << fix_string(w[i]) << ";\n";


  // Do numerical verification
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
      // int((1-x)*x^p, x=0..1) = 1 / (p^2 + 3p + 2)

      // Exact solution for alpha=2:
      //
      // int((1-x)^2*x^p, x=0..1) = 2 / (p^3 + 6*p^2 + 11*p + 6)
      mpfr_class exact = (alpha==1) ?
        mpfr_class(1.0)/mpfr_class(order*order + 3*order + 2) :
        mpfr_class(2.0)/mpfr_class(order*order*order + 6*order*order + 11*order + 6) ;

      // std::cout << "quadrature = " << sum << std::endl;
      // std::cout << "exact      = " << exact << std::endl;

      // Compute the absolute error:
      mpfr_class abs_err = my_abs(mpfr_class(sum - exact));

      // Print message
      std::cout << "Computing int((1-x)^" << alpha << " * x^" << order << ", x=0..1)"
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

// We used the scaled interval [0,1] here, if one instead uses
// the interval [-1,1]: int((1-x)*x^p, x=-1..1), the exact
// solution is given by:
//        p           p
//  2 (-1)  p + 3 (-1)  + 1
//  -----------------------
//      (p + 2) (p + 1)
//
// = { ( 2*p + 4) / (p + 2) / (p + 1), p even
//   { (-2*p - 2) / (p + 2) / (p + 1), p odd

void usage()
{
  std::cout << "\n";
  std::cout << "This program generates the points and weights for a Jacobi quadrature rule.\n";
  std::cout << "\n";
  std::cout << "Valid command line options are:\n";
  std::cout << "--npoints, -n # = Build Jacobi rule with # points.\n";
  std::cout << "--alpha, -a #   = Build Jacobi rule with alpha = #.  Only alpha=1, 2 are supported.\n";
  std::cout << "--beta, -b #    = Build Jacobi rule with beta = #.  Only beta=0 is supported.\n";
  std::cout << "--help, -h      = Print this message.\n";
  std::cout << "\n";
}
