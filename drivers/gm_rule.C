#include <stdlib.h> // atoi
#include "grundmann_moller.h"

int main(int argc, char** argv)
{
  std::cout.precision(32);
  std::cout.setf(std::ios_base::scientific);

  // # of binary digits
  // 53 binary digits is about what you normally get with a double.
  mpfr_set_default_prec(256);

  // Read rule index, s, from command line.  Any integer s>=0 is valid.
  unsigned s = 0;
  if (argc > 1)
    s = atoi(argv[1]);

  GrundmannMoller gm;
  gm.rule(s);

  // Get a reference to the points and weights arrays
  const std::vector<Point>& x = gm.get_points();
  const std::vector<mpfr_class>& w = gm.get_weights();

  std::cout << "Index " << s
            << " rule has degree " << 2*s+1
            << ", and " << x.size() << " point(s)." << std::endl;

  mpfr_class sum_weights = 0.;

  // Print the points and weights
  for (unsigned i=0; i<x.size(); ++i)
    {
      std::cout << fix_string(x[i](0)) << ", "
                << fix_string(x[i](1)) << ", "
                << fix_string(x[i](2)) << ", "
                << fix_string(w[i]) << std::endl;
      sum_weights += w[i];
    }

  // std::cout << "Sum of weights = " << sum_weights << std::endl;
  mpfr_class error_weights = abs(sum_weights - mpfr_class(1.)/mpfr_class(6.));
  std::cout << "Error in sum of weights = " << error_weights << std::endl;

  return 0;
}
