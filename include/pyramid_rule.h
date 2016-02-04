#ifndef PYRAMID_RULE_H
#define PYRAMID_RULE_H

#include "jacobi.h"
#include "gauss.h"

/**
 * This class combines the Gauss and Jacobi 3D rules for Pyramid regions.
 */
class PyramidRule
{
public:
  // Generate the pyramid quadrature rule that can integrate polynomials of
  // total order 'order' exactly.
  void generate(unsigned int order);

  // Constant access to the points and weights vectors
  const std::vector<Point<mpfr_class> >& get_points()  { return x; }
  const std::vector<mpfr_class>& get_weights() { return w; }

private:
  // The points of the quadrature rule.
  std::vector<Point<mpfr_class> > x;

  // Weights for a Jacobi quadrature rule of degree n.
  std::vector<mpfr_class> w;
};

#endif

