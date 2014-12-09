#ifndef LEGENDRE_H
#define LEGENDRE_H

#include <vector>

#include "gmpfrxx.h"
#include "matrix.h"

class Legendre
{
public:
  Legendre() {}

  // Evaluate the Legendre polynomials up to and including degree d at
  // the point x.  Fill up the vector vals with the polynomial values
  // in the canonical ordering.  These values were computed and
  // verified by the legendre.py script.
  void p(unsigned d,
         const mpfr_class & x,
         std::vector<mpfr_class> & vals);

  // Return a vector of the orthogonality coefficients, these are the
  // integrals of phi_i^2 for each phi_i
  void orthogonality_coeffs(unsigned d,
                            std::vector<mpq_class> & coeffs);

  // Builds the coefficient matrix K_{ij} = int(phi_i*phi_j +
  // grad(phi_i)*grad(phi_j)) for Legendre polynomials of degree d.
  void build_H1_projection_matrix(unsigned d,
                                  Matrix<mpfr_class> & matrix);

private:
  // Used to fill in the H1 projection matrix
  static const unsigned laplace_matrix[31][31];
};

#endif
