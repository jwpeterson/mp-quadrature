#include <cstdlib> // std::abort, atoi
#include <iomanip>
#include "gauss.h"
#include "jacobi.h"
#include "exact.h"

// In this program we compute multi-precision
// points and weights to be used in Conical Product
// formulae.
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
  std::cout << "\nComputing conical product rule with n^2=" << n*n << " points." << std::endl;

  // Compute the Gauss rule points/weights
  std::vector<mpfr_class> gauss_x, gauss_w;
  gauss_rule(n, gauss_x, gauss_w);

  // Scale the Gauss points so they lie in [0,1].  We should probably make
  // a Gauss class and have this be a member. See eg: Jacobi
  mpfr_class zero(0.0), one(1.0);
  {
    mpfr_class a ( 0.5*(one-zero) );
    mpfr_class b ( 0.5*(zero+one) );
    for (unsigned int j=1; j<gauss_x.size(); ++j)
      {
        gauss_x[j] = a*gauss_x[j] + b;
        gauss_w[j] *= 0.5;
      }
  }

  // Compute the Jacobi rule points/weights
  const Real alpha=1.0, beta=0.0;
  Jacobi p(alpha,beta);
  p.rule(n);

  // Scale Jacobi weights so they sum to 1/3 (alpha==2) or 1/2 (alpha==1)
  if (alpha==2.0)
    {
      mpfr_class one_third(1.0);
      one_third /= 3.0;
      p.scale_weights(one_third);
    }
  else if (alpha==1.0)
    {
      p.scale_weights(0.5);
    }
  else
    {
      std::cout << "Warning: weights unscaled!" << std::endl;
    }

  // Scale Jacobi points so they lie on [0, 1]
  p.scale_points(zero, one);

  // Print the result
  // p.printxw();

  // Get const references to the jacobi points and weights vectors.
  const std::vector<mpfr_class>& jacobi_x = p.get_points();
  const std::vector<mpfr_class>& jacobi_w = p.get_weights();

  // Compute the conical product rule, with space for n^2 entries.
  // See also the code in LibMesh's src/quadrature/quadrature.C
  std::vector<mpfr_class> conical_x(n*n), conical_y(n*n), conical_w(n*n);
  unsigned int gp = 0;
  for (unsigned int i=0; i<n; i++)
    for (unsigned int j=0; j<n; j++)
      {
        // Note: Access the 1D arrays from [1] ... [n]
        conical_x[gp] = jacobi_x[j+1];                     //s[j];
        conical_y[gp] = gauss_x[i+1] * (1.-jacobi_x[j+1]); //r[i]*(1.-s[j]);
        conical_w[gp] = gauss_w[i+1] * jacobi_w[j+1];      //A[i]*B[j];
        gp++;
      }

  // Print out the conical product points and weights in a form we can use
  for (unsigned int i=0; i<n*n; ++i)
    {
      std::cout << "_points[" << std::setw(2) << i << "](0)=" << fix_string(conical_x[i]) << ";\n";
      std::cout << "_points[" << std::setw(2) << i << "](1)=" << fix_string(conical_y[i]) << ";\n";
    }

  std::cout << "\n";
  for (unsigned int i=0; i<n*n; ++i)
    std::cout << "_weights[" << std::setw(2) << i << "]=" << fix_string(conical_w[i]) << ";\n";

  // As a check of the method, compute the sum of the weights
  mpfr_class sumweights=0.0;
  for (unsigned int i=0; i<n*n; ++i)
    sumweights += conical_w[i];

  std::cout << "Sum of weights = " << fix_string(sumweights) << std::endl;

  // Maple commands:
  // assume(a, integer, a, nonnegative);
  // assume(b, integer, b, nonnegative);
  // simplify(int(int(x^a * y^b, y=0..1-x), x=0..1));
  // Result is:
  // 1/(b+1)*Beta(a+1,b+2)
  // where Beta(p,q) := (p-1)! (q-1)! / (p+q-1)! is the beta function, as defined for positive integers.
  // After simplification, the result is:
  // a! b! / (a + b + 2)!

  // The rule with n^2 points should be able to integrate a polynomial of _total_ degree 2*n-1
  std::cout << "\nVerifying rule:" << std::endl;
  unsigned max_order = 2*n-1;
  for (unsigned x_power=0; x_power<max_order; ++x_power)
    for (unsigned y_power=0; y_power<max_order; ++y_power)
      {
        // Only try to integrate polynomials we can integrate exactly
        if (x_power + y_power > max_order)
          continue;

        // conical_w, conical_x, and conical_y use 0-based array access!
        mpfr_class sum = 0.;
        for (unsigned n_qp=0; n_qp<conical_x.size(); ++n_qp)
          sum += conical_w[n_qp] * pow(conical_x[n_qp], x_power) * pow(conical_y[n_qp], y_power);

        // std::cout << "quadrature = " << sum << std::endl;

        // Compute the analytical integral value.
        mpfr_class analytical = exact_tri(x_power, y_power);

        // std::cout << "analytical = " << analytical << std::endl;

        // Compute the absolute error:
        mpfr_class abs_err = abs(sum-analytical);

        // Print message
        std::cout << "Computing integral of: x^" << x_power << " y^" << y_power
                  << ", abs_err = " << abs_err << std::endl;

        // Abort if error is too large.
        if (abs_err > mpfr_class(1.e-30))
          {
            std::cerr << "Quadrature error too large, possible problem with points and weights!" << std::endl;
            std::abort();
          }
      }

  return 0;
}
