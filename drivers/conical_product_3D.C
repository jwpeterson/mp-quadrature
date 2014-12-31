#include <cstdlib> // std::abort, atoi
#include <iomanip>
#include <algorithm> // std::sort
#include "exact.h"
#include "conical.h"

// In this program we compute multi-precision
// points and weights to be used in 3D Conical Product
// formulae for tetrahedra.
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
    n=atoi(argv[1]);

  const unsigned int n_total = n*n*n;
  std::cout << "\nComputing 3D conical product rule with n^3=" << n_total << " points." << std::endl;

  // Construct a conical product rule of order 2*n-1, where n is the requested number of points
  Conical conical;
  conical.rule3D(2*n-1);

  // Get a reference to the conical product rule points and weights
  const std::vector<Point<mpfr_class> > & rule_points = conical.get_points();
  const std::vector<mpfr_class> & rule_weights = conical.get_weights();

  // Print out the conical product points and weights in a form we can use
  // (These are now in 0-based storage)
  for (unsigned int i=0; i<n_total; ++i)
    {
      std::cout << "_points[" << std::setw(2) << i << "](0)=" << fix_string(rule_points[i](0)) << ";\n";
      std::cout << "_points[" << std::setw(2) << i << "](1)=" << fix_string(rule_points[i](1)) << ";\n";
      std::cout << "_points[" << std::setw(2) << i << "](2)=" << fix_string(rule_points[i](2)) << ";\n";
    }

  std::cout << "\n";
  for (unsigned int i=0; i<n_total; ++i)
    std::cout << "_weights[" << std::setw(2) << i << "]=" << fix_string(rule_weights[i]) << ";\n";

  // As a check of the method, compute the sum of the weights
  mpfr_class sumweights=0.0;
  for (unsigned int i=0; i<n_total; ++i)
    sumweights += rule_weights[i];

  std::cout << "Sum of weights = " << fix_string(sumweights) << std::endl;

  // The rule with n^3 points should be able to integrate a polynomial of _total_ degree 2*n-1
  std::cout << "\nVerifying rule..." << std::endl;
  unsigned max_order = 2*n-1;
  for (unsigned x_power=0; x_power<max_order; ++x_power)
    for (unsigned y_power=0; y_power<max_order; ++y_power)
      for (unsigned z_power=0; z_power<max_order; ++z_power)
        {
          // Only try to integrate polynomials we can integrate exactly
          if (x_power + y_power + z_power > max_order)
            continue;

          mpfr_class sum = 0.;
          for (unsigned n_qp=0; n_qp<rule_points.size(); ++n_qp)
            sum += rule_weights[n_qp]
              * pow(rule_points[n_qp](0), x_power)
              * pow(rule_points[n_qp](1), y_power)
              * pow(rule_points[n_qp](2), z_power);

          // std::cout << "quadrature = " << sum << std::endl;

          // Compute the analytical integral value.
          mpfr_class analytical = exact_tet(x_power, y_power, z_power);

          // std::cout << "analytical = " << analytical << std::endl;

          // Compute the absolute error:
          mpfr_class abs_err = abs(sum-analytical);

          // Print message.  In 3D, this is just way too much output.
          //std::cout << "Computing integral of: x^" << x_power << " y^" << y_power << " z^" << z_power
          //          << ", abs_err = " << abs_err << std::endl;

          // Abort if error is too large.
          if (abs_err > mpfr_class(1.e-30))
            {
              std::cerr << "Quadrature error too large, possible problem with points and weights!" << std::endl;
              std::abort();
            }
        }

  std::cout << "... Verification complete!" << std::endl;

  return 0;
}
