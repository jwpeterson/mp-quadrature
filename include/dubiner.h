#ifndef DUBINER_H
#define DUBINER_H

#include <vector>

#include "gmpfrxx.h"
#include "matrix.h"
#include "common_definitions.h"

class Dubiner
{
public:
  Dubiner() {}

  // Evaluate the Dubiner polynomials up to and including degree d at
  // the point (xi, eta).  Fill up the vector vals with the polynomial
  // values in the canonical ordering.  This routine evaluates the
  // Dubiner polynomials "numerically" using the method in Burgers'
  // paper, instead of using generated code.  Once this is working, it
  // should be a lot less code than the generated code route...
  void p_numeric(unsigned d,
                 const mpfr_class & xi,
                 const mpfr_class & eta,
                 std::vector<mpfr_class> & vals);

  // Same as above, but returns grad(p)
  void dp(unsigned d,
          const mpfr_class & xi,
          const mpfr_class & eta,
          std::vector<Point<mpfr_class> > & vals);

  // Builds the coefficient matrix K_{ij} = int(phi_i*phi_j +
  // grad(phi_i)*grad(phi_j)) for the Dubiner polynomials of degree d.
  // Currently only works for d=10...
  void build_H1_projection_matrix(unsigned d,
                                  Matrix<mpfr_class> & matrix);

private:
  // Compute the value of the nth Jacobi polynomial, P_n^{alpha,beta}
  // at x.  Note that this implementation has a slightly different
  // normalization than the one in the Jacobi class, which is designed
  // specifically to create conical product quadrture rules...  for
  // P_1, the Jacobi class produces:
  // x + (alpha-beta)/(alpha+beta+2)
  // while this class produces:
  // (1/2)*((alpha+beta+2)*x + (alpha-beta))
  mpfr_class jacobi(unsigned n, unsigned alpha, unsigned beta, mpfr_class x);

  // Same as above, but the first derivative
  mpfr_class djacobi(unsigned n, unsigned alpha, unsigned beta, mpfr_class x);
};

#endif
