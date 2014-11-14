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

  // Compute and print the dim=3 Grundmann-Moller quadrature rule with
  // index s.  Such a rule will be capable of exactly integrating
  // polynomials of degree <= 2*s+1, and will have
  // (dim+1+s)! / (dim+1)! / s!
  // points.  Note that it is not valid to ask for an even degree rule,
  // only odd-order GM rules exist.
  void rule(unsigned s);

  // The points and weights of the GM rules are rational numbers, so
  // it is possible to use the mpq_class type here.
  const std::vector<Point<mpq_class> >& get_points()  { return x; }
  const std::vector<mpq_class>& get_weights() { return w; }

private:
  // Roots of the Jacobi polynomial of degree n.
  std::vector<Point<mpq_class> > x;

  // Weights for a Jacobi quadrature rule of degree n.
  std::vector<mpq_class> w;

  // Copy constructor. Unimplemented.
  GrundmannMoller(const GrundmannMoller&);

  // Assignment operator.  Unimplemented
  GrundmannMoller& operator=(const GrundmannMoller&);
};

#endif
