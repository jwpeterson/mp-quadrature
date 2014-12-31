#ifndef CONICAL_H
#define CONICAL_H

#include "jacobi.h"
#include "gauss.h"

/**
 * This class combines the Gauss and Jacobi rules to make 2D and 3D
 * conical product rules.
 */
class Conical
{
public:
  // Compute the 2D Conical quadrature rule which can integrate polynomials of
  // total order 'order' exactly.
  void rule2D(unsigned int order);

  // Compute the 3D Conical quadrature rule which can integrate polynomials of
  // total order 'order' exactly.
  void rule3D(unsigned int order);

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

