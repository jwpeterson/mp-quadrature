#include "gauss.h"
#include "jacobi.h"

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
  std::cout << "\nComputing conical product rule with n=" << n << " points." << std::endl;

  // Compute the Gauss rule points/weights
  std::vector<mpfr_class> gauss_x, gauss_w;
  gauss_rule(n, gauss_x, gauss_w);

  // Compute the Jacobi rule points/weights
  Jacobi p(/*alpha=*/2.0, /*beta=*/0.0);
  p.rule(n);
  
  // Scale Jacobi weights so they sum to 1/3.  This is the one for
  // triangles.
  mpfr_class one_third(1.0);
  one_third /= 3.0;
  p.scale_weights(one_third);

  // Scale Jacobi points so they lie on [0, 1]
  mpfr_class zero(0.0), one(1.0);
  p.scale_points(zero, one);

  // Print the result
  p.printxw();
  
  return 0;
}
