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
  // the point (xi, eta).  Fill up the vector 'vals' with the
  // polynomial values, and the vector 'gradients' with the gradients.
  // This routine evaluates the Dubiner polynomials "numerically"
  // using the method described in Burgers' paper.
  // R. B\"{u}rger, M. Sepulveda, and T. Voitovich, "On the
  // Proriol-Koornwinder-Dubiner hierarchical orthogonal polynomial
  // basis for the DG-FEM".
  // http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.582.3100&rep=rep1&type=pdf
  void p(unsigned d,
         const mpfr_class & xi,
         const mpfr_class & eta,
         std::vector<mpfr_class> & vals,
         std::vector<Point<mpfr_class> > & gradients);

  // Builds the coefficient matrix K_{ij} = int(phi_i*phi_j +
  // grad(phi_i)*grad(phi_j)) for the Dubiner polynomials of degree d.
  void build_H1_projection_matrix(unsigned d,
                                  Matrix<mpfr_class> & matrix);

  // Compare multiple precision and double precision Jacobi polynomial
  // values.  This was used in developing/verifying a double-precision
  // version of the Jacobi polynomial recurrence relation.
  void compare_jacobi();

private:
  // Compute the value and first derivative of the nth Jacobi
  // polynomial, P_n^{alpha,beta}, at x.  Note that this
  // implementation has a slightly different normalization than the
  // one in the Jacobi class, which is designed specifically to create
  // conical product quadrture rules...  for P_1, the Jacobi class
  // produces:
  // x + (alpha-beta)/(alpha+beta+2)
  // while this class produces:
  // (1/2)*((alpha+beta+2)*x + (alpha-beta))
  std::pair<mpfr_class, mpfr_class> jacobi(unsigned n, unsigned alpha, unsigned beta, mpfr_class x);

  // Compute Jacobi polynomial value and first derivative:
  // P^{(\alpha,\beta)}_n(x)
  // d/dx P^{(\alpha,\beta)}_n(x)
  // using double precision values. This is here so we can compare it
  // to the multiprecision implementation above. It therefore uses the
  // same normalization as the implementation above. The recurrence
  // relation we are using can be found on Wikipedia, and is similar
  // in structure to the recurrence relation for Legendre polynomials,
  // so we use a similar approach to compute it. The page also describes
  // recurrence relations for computing the first derivative.
  std::pair<double, double>
  jacobi(unsigned n, unsigned alpha, unsigned beta, double x);
};

#endif
