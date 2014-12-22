#include "common_definitions.h"
#include "dubiner.h"

// This test driver verifies the orthogonality of the Dubiner
// polynomials, both symbolic and numeric versions.
int main(int argc, char** argv)
{
  std::cout.precision(32);
  std::cout.setf(std::ios_base::scientific);
  // 53 binary digits is about what you normally get with a double.
  mpfr_set_default_prec(256);

  Dubiner dubiner;

  const unsigned max_dubiner_degree = 10;

  for (unsigned dubiner_degree=0; dubiner_degree <= max_dubiner_degree; ++dubiner_degree)
    {
      // The points and weights from a conical product quadrature rule
      // of high enough order to integrate a Dubiner mass matrix
      // exactly.
      std::vector<Point<mpfr_class> > rule_points;
      std::vector<mpfr_class> rule_weights;

      std::vector<mpfr_class> current_vals;

      for (unsigned i=0; i<rule_points.size(); ++i)
        {
          // Evaluate all the Dubiner polynomials at the current qp
          dubiner.p_numeric(dubiner_degree,
                            /*xi=*/ rule_points[i](0),
                            /*eta=*/ rule_points[i](1),
                            current_vals);
        }
    }

  return 0;
}
