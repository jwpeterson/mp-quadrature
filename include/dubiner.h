#ifndef DUBINER_H
#define DUBINER_H

#include <vector>

#include "gmpfrxx.h"
#include "matrix.h"

class Dubiner
{
public:
  Dubiner() {}

  // Evaluate the Dubiner polynomials up to and including degree d at
  // the point (xi, eta).  Fill up the vector vals with the polynomial
  // values in the canonical ordering.  These values were computed and
  // verified by the dubiner.py script.
  void p(unsigned d,
         const mpfr_class & xi,
         const mpfr_class & eta,
         std::vector<mpfr_class> & vals);

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

  // Return a vector of the orthogonality coefficients, these are the
  // integrals of phi_i^2 for each phi_i
  void orthogonality_coeffs(unsigned d,
                            std::vector<mpq_class> & coeffs);

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

  // I split the guts of the "p" function across several files so we
  // could compile it in parallel more easily.
  void dubiner_5th(const mpfr_class & zeta0,
                   const mpfr_class & zeta1,
                   const mpfr_class & zeta2,
                   std::vector<mpfr_class> & vals);

  void dubiner_6th(const mpfr_class & zeta0,
                   const mpfr_class & zeta1,
                   const mpfr_class & zeta2,
                   std::vector<mpfr_class> & vals);

  void dubiner_7th(const mpfr_class & zeta0,
                   const mpfr_class & zeta1,
                   const mpfr_class & zeta2,
                   std::vector<mpfr_class> & vals);

  void dubiner_8th(const mpfr_class & zeta0,
                   const mpfr_class & zeta1,
                   const mpfr_class & zeta2,
                   std::vector<mpfr_class> & vals);

  void dubiner_9th(const mpfr_class & zeta0,
                   const mpfr_class & zeta1,
                   const mpfr_class & zeta2,
                   std::vector<mpfr_class> & vals);

  void dubiner_10th(const mpfr_class & zeta0,
                    const mpfr_class & zeta1,
                    const mpfr_class & zeta2,
                    std::vector<mpfr_class> & vals);

  void dubiner_11th(const mpfr_class & zeta0,
                    const mpfr_class & zeta1,
                    const mpfr_class & zeta2,
                    std::vector<mpfr_class> & vals);

  void dubiner_12th(const mpfr_class & zeta0,
                    const mpfr_class & zeta1,
                    const mpfr_class & zeta2,
                    std::vector<mpfr_class> & vals);

  void dubiner_13th(const mpfr_class & zeta0,
                    const mpfr_class & zeta1,
                    const mpfr_class & zeta2,
                    std::vector<mpfr_class> & vals);

  void dubiner_14th(const mpfr_class & zeta0,
                    const mpfr_class & zeta1,
                    const mpfr_class & zeta2,
                    std::vector<mpfr_class> & vals);

  void dubiner_15th(const mpfr_class & zeta0,
                    const mpfr_class & zeta1,
                    const mpfr_class & zeta2,
                    std::vector<mpfr_class> & vals);

  void dubiner_16th(const mpfr_class & zeta0,
                    const mpfr_class & zeta1,
                    const mpfr_class & zeta2,
                    std::vector<mpfr_class> & vals);

  void dubiner_17th(const mpfr_class & zeta0,
                    const mpfr_class & zeta1,
                    const mpfr_class & zeta2,
                    std::vector<mpfr_class> & vals);

  void dubiner_18th(const mpfr_class & zeta0,
                    const mpfr_class & zeta1,
                    const mpfr_class & zeta2,
                    std::vector<mpfr_class> & vals);

  void dubiner_19th(const mpfr_class & zeta0,
                    const mpfr_class & zeta1,
                    const mpfr_class & zeta2,
                    std::vector<mpfr_class> & vals);

  void dubiner_20th(const mpfr_class & zeta0,
                    const mpfr_class & zeta1,
                    const mpfr_class & zeta2,
                    std::vector<mpfr_class> & vals);

  void dubiner_21st(const mpfr_class & zeta0,
                    const mpfr_class & zeta1,
                    const mpfr_class & zeta2,
                    std::vector<mpfr_class> & vals);

  void dubiner_22nd(const mpfr_class & zeta0,
                    const mpfr_class & zeta1,
                    const mpfr_class & zeta2,
                    std::vector<mpfr_class> & vals);

  void dubiner_23rd(const mpfr_class & zeta0,
                    const mpfr_class & zeta1,
                    const mpfr_class & zeta2,
                    std::vector<mpfr_class> & vals);

  void dubiner_24th(const mpfr_class & zeta0,
                    const mpfr_class & zeta1,
                    const mpfr_class & zeta2,
                    std::vector<mpfr_class> & vals);

  void dubiner_25th(const mpfr_class & zeta0,
                    const mpfr_class & zeta1,
                    const mpfr_class & zeta2,
                    std::vector<mpfr_class> & vals);

  void dubiner_26th(const mpfr_class & zeta0,
                    const mpfr_class & zeta1,
                    const mpfr_class & zeta2,
                    std::vector<mpfr_class> & vals);

  void dubiner_27th(const mpfr_class & zeta0,
                    const mpfr_class & zeta1,
                    const mpfr_class & zeta2,
                    std::vector<mpfr_class> & vals);

  void dubiner_28th(const mpfr_class & zeta0,
                    const mpfr_class & zeta1,
                    const mpfr_class & zeta2,
                    std::vector<mpfr_class> & vals);

  void dubiner_29th(const mpfr_class & zeta0,
                    const mpfr_class & zeta1,
                    const mpfr_class & zeta2,
                    std::vector<mpfr_class> & vals);

  void dubiner_30th(const mpfr_class & zeta0,
                    const mpfr_class & zeta1,
                    const mpfr_class & zeta2,
                    std::vector<mpfr_class> & vals);

  // Used to fill in the H1 projection matrix - this is only valid for d=10!
  static const mpfr_class laplace_matrix[66][66];
};

#endif
