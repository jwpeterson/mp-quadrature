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

  // Constant access to the points and weights vectors
  const std::vector<Point>& get_points()  { return x; }
  const std::vector<mpfr_class>& get_weights() { return w; }

private:
  // Roots of the Jacobi polynomial of degree n.
  std::vector<Point> x;

  // Weights for a Jacobi quadrature rule of degree n.
  std::vector<mpfr_class> w;

  /**
   * Helper function for computing compositions
   * This routine for computing compositions and their permutations is
   * originally due to:
   *
   * Albert Nijenhuis, Herbert Wilf,
   * Combinatorial Algorithms for Computers and Calculators,
   * Second Edition,
   * Academic Press, 1978,
   * ISBN: 0-12-519260-6,
   * LC: QA164.N54.
   *
   * I still don't understand exactly why it works, but it does.
   * s = number to be compositioned
   * p = number of partitions
   */
  void compose_all(unsigned s,
                   unsigned p,
                   std::vector<std::vector<unsigned> >& result);

  // Copy constructor. Unimplemented.
  GrundmannMoller(const GrundmannMoller&);

  // Assignment operator.  Unimplemented
  GrundmannMoller& operator=(const GrundmannMoller&);
};

#endif
