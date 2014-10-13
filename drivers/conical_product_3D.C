#include <cstdlib> // std::abort, atoi
#include <iomanip>
#include <algorithm> // std::sort
#include "gauss.h"
#include "jacobi.h"

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

  // Compute the Jacobi rules' points/weights:
  //
  // p1 := alpha=1, beta=0
  const Real alpha1=1.0, beta1=0.0;
  Jacobi p1(alpha1,beta1);
  p1.rule(n);
  p1.scale_weights(0.5);
  p1.scale_points(zero, one);
  // p1.printxw();  // Print the result


  // p2 := alpha=2, beta=0
  const Real alpha2=2.0, beta2=0.0;
  Jacobi p2(alpha2,beta2);
  p2.rule(n);
  {
    mpfr_class one_third(1.0);
    one_third /= 3.0;
    p2.scale_weights(one_third);
  }
  p2.scale_points(zero, one);
  // p2.printxw(); // Print the result

  // Get const references to the jacobi points and weights vectors.
  const std::vector<mpfr_class>& jacobi1_x = p1.get_points();
  const std::vector<mpfr_class>& jacobi1_w = p1.get_weights();
  const std::vector<mpfr_class>& jacobi2_x = p2.get_points();
  const std::vector<mpfr_class>& jacobi2_w = p2.get_weights();

  // Compute the conical product rule, with space for n^3 entries.
  // See also the code in LibMesh's src/quadrature/quadrature.C
  std::vector<mpfr_class> conical_x(n_total), conical_y(n_total), conical_z(n_total), conical_w(n_total);
  unsigned int gp = 0;
  for (unsigned int i=0; i<n; i++)
    for (unsigned int j=0; j<n; j++)
      for (unsigned int k=0; k<n; k++)
        {
          // Note: Access the 1D arrays from [1] ... [n]
          conical_x[gp] = jacobi2_x[k+1]; // jacB1D.qp(k)(0);
          conical_y[gp] = jacobi1_x[j+1] * (1.-jacobi2_x[k+1]); // jacA1D.qp(j)(0) * (1.-jacB1D.qp(k)(0));
          conical_z[gp] = gauss_x[i+1] * (1.-jacobi1_x[j+1]) * (1.-jacobi2_x[k+1]); // gauss1D.qp(i)(0) * (1.-jacA1D.qp(j)(0)) * (1.-jacB1D.qp(k)(0));
          conical_w[gp] = gauss_w[i+1] * jacobi1_w[j+1] * jacobi2_w[k+1]; //gauss1D.w(i) * jacA1D.w(j) * jacB1D.w(k);
          gp++;
        }

  // Print out the conical product points and weights in a form we can use
  // (These are now in 0-based storage)
  for (unsigned int i=0; i<n_total; ++i)
    {
      std::cout << "_points[" << std::setw(2) << i << "](0)=" << fix_string(conical_x[i]) << ";\n";
      std::cout << "_points[" << std::setw(2) << i << "](1)=" << fix_string(conical_y[i]) << ";\n";
      std::cout << "_points[" << std::setw(2) << i << "](2)=" << fix_string(conical_z[i]) << ";\n";
    }

  std::cout << "\n";
  for (unsigned int i=0; i<n_total; ++i)
    std::cout << "_weights[" << std::setw(2) << i << "]=" << fix_string(conical_w[i]) << ";\n";

  // As a check of the method, compute the sum of the weights
  mpfr_class sumweights=0.0;
  for (unsigned int i=0; i<n_total; ++i)
    sumweights += conical_w[i];

  std::cout << "Sum of weights = " << fix_string(sumweights) << std::endl;

  // Maple has a lot of trouble computing this integral exactly. Try for yourself:
  // assume(a, integer, a, nonnegative);
  // assume(b, integer, b, nonnegative);
  // assume(c, integer, c, nonnegative);
  // int(int(int(x^a * y^b * z^c, z=0..1-x-y), y=0..1-x), x=0..1);
  //
  // The main trick is to perform a change of variables to get the
  // second integral into the form of the beta function integral...
  //
  // 1.) Integrate in z: result is I = int(int((x^a * y^b * (1-x-y)^(c+1))/(c+1), y=0..1-x), x=0..1)
  // 2.) Introduce change of variables eta = y/(1-x)
  //     Then I = Beta(b+1, c+2) / (c+1) * int(x^a * (1-x)^(b+c+2), x=0..1)
  // 3.) I can now be computed using the Beta function integral again to get:
  //     I = Beta(b+1, c+2) / (c+1) * Beta(a+1, b+c+3)
  // 4.) Substituting in the Beta function definition,
  //     Beta(p,q) := (p-1)! (q-1)! / (p+q-1)!,
  //     and canceling terms, we finally obtain:
  //     I = a! b! c! / (a + b + c + 3)!

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

          // conical_w, conical_x, and conical_y use 0-based array access!
          mpfr_class sum = 0.;
          for (unsigned n_qp=0; n_qp<conical_x.size(); ++n_qp)
            sum += conical_w[n_qp] * pow(conical_x[n_qp], x_power) * pow(conical_y[n_qp], y_power) * pow(conical_z[n_qp], z_power);

          // std::cout << "quadrature = " << sum << std::endl;

          // Compute the exact solution as described above, divide
          // term-by-term to avoid huge numbers/overflows (even though
          // this is GMP).
          mpfr_class analytical = 1.0;
          {
            // Sort the a, b, c values
            unsigned sorted_powers[3] = {x_power, y_power, z_power};
            std::sort(sorted_powers, sorted_powers+3);

            // Cancel the largest power with the denominator, fill in the
            // entries for the remaining numerator terms and the denominator.
            std::vector<unsigned>
              numerator_1(sorted_powers[0] > 1 ? sorted_powers[0]-1 : 0),
              numerator_2(sorted_powers[1] > 1 ? sorted_powers[1]-1 : 0),
              denominator(3 + sorted_powers[0] + sorted_powers[1]);

            // Fill up the vectors with sequences starting at the right values.
            iota(numerator_1.begin(), numerator_1.end(), 2);
            iota(numerator_2.begin(), numerator_2.end(), 2);
            iota(denominator.begin(), denominator.end(), sorted_powers[2]+1);

            // The denominator is guaranteed to have the most terms...
            for (unsigned i=0; i<denominator.size(); ++i)
              {
                if (i < numerator_1.size())
                  analytical *= numerator_1[i];

                if (i < numerator_2.size())
                  analytical *= numerator_2[i];

                analytical /= denominator[i];
              }
          }

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
