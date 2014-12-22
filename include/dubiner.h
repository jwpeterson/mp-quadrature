#ifndef DUBINER_H
#define DUBINER_H

#include <vector>

#include "gmp.h"
#include "mpfr.h"
#include "gmpfrxx.h"

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

  // Return a vector of the orthogonality coefficients, these are the
  // integrals of phi_i^2 for each phi_i
  void orthogonality_coeffs(unsigned d,
                            std::vector<mpq_class> & coeffs);

private:
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
};

#endif