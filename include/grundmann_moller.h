#ifndef __grundmann_moller_h__
#define __grundmann_moller_h__

#include <vector>

#include "gmpfrxx.h"
#include "common_definitions.h"

class GrundmannMoller
{
public:
  // Constructor.
  GrundmannMoller();

  // Destructor.
  ~GrundmannMoller() {}

  // Compute and print the dim=3 Grundmann-Moller quadrature rule capable of
  // exactly integrating polynomials of degree <= d.  Such a rule will have
  // (dim+1+s)! / (dim+1)! / s!
  // points, where s := (d-1)/2.
  void rule(unsigned int d);

  // Constant access to the points and weights vectors
  const std::vector<Point>& get_points()  { return x; }
  const std::vector<mpfr_class>& get_weights() { return w; }

private:
  // Roots of the Jacobi polynomial of degree n.
  std::vector<Point> x;

  // Weights for a Jacobi quadrature rule of degree n.
  std::vector<mpfr_class> w;

  // Copy constructor. Unimplemented.
  GrundmannMoller(const GrundmannMoller&);

  // Assignment operator.  Unimplemented
  GrundmannMoller& operator=(const GrundmannMoller&);
};

#endif
