#ifndef __jacobi_h__
#define __jacobi_h__

#include <vector>

#include "gmp.h"
#include "mpfr.h"
#include "gmpfrxx.h"
#include "common_definitions.h"


/**
 * This class encapsulates the different components
 * required for computing multiprecision Jacobi polynomials,
 * finding their roots via Newton's Method, and computing
 * weights for quadrature.
 */
class Jacobi
{
public:
  // Constructor.  Initializes the alpha, beta
  // parameters.
  Jacobi(Real alpha, Real beta);

  // Destructor.
  ~Jacobi() {}

  // Set or re-set the alpha parameter.
  void set_alpha(Real alpha);

  // Set or re-set the beta parameter.
  void set_beta(Real beta);
  
  // Compute and print the Jacobi quadrature rule of
  // degree n.  If previous rules have been computed,
  // this routine may try to re-use some data.
  void rule(unsigned int n);

  // Scale the weights by multiplying them all by 'scale_factor'.
  // (By default, the weights sum to 1.0)
  void scale_weights(const mpfr_class& scale_factor);

  // Scale the points to lie in the inverval [x1, x2]
  // (By default, the points lie in [-1 1])
  void scale_points(const mpfr_class& x1, const mpfr_class& x2);

  // Constant access to the points and weights vectors
  const std::vector<mpfr_class>& get_points()  { return x; }
  const std::vector<mpfr_class>& get_weights() { return w; }
  
private:
  // Multi-precision versions of the parameters
  mpfr_class mp_alpha, mp_beta;

  // Order of Jacobi polynomial we're computing.
  //unsigned int n;

  // Multi-precision object to hold P_n(x)
  mpfr_class p;

  // Multi-precision object to hold P'_n(x)
  mpfr_class dp;

  // Multi-precision object to hold P_{n-1}(x)
  mpfr_class pnm1;
  
  // Jacobi recurrence relation coefficients.
  std::vector<mpfr_class> b, c;

  // Roots of the Jacobi polynomial of degree n.
  std::vector<mpfr_class> x;

  // Weights for a Jacobi quadrature rule of degree n.
  std::vector<mpfr_class> w;
  
  // Function which computes Jacobi polynomial recurrence relation
  // constants.  Only computes those which have not already been
  // computed.
  void constants(unsigned int n);

  // Resets any previously computed values.  To be used when
  // alpha or beta changes.
  void reset();

  // Computes p, dp, and pnm1 for the Jacobi polynomial of
  // degree n at the point x.  x must lie in the
  // domain [-1 1].  A version of this could also
  // be available to be called publicly... 
  // The input value x is a multi-precision object.  It would
  // be nice to have a version of this function where x is a
  // vector, but not essential.
  void value(const mpfr_class& xval, unsigned int n);

  // Computes the roots of the Jacobi polynomial of degree n,
  // using the initial guesses published by Stroud, and Newton's
  // method to converge the root to the required tolerance.  Fills
  // up the vector x located in this class.
  void points(unsigned int n);

  // Computes the root of P_n(x) closest to the guess value x.  The
  // root is returned in the same multi-precision input object as the
  // guess.
  void newton(unsigned int n, mpfr_class& xroot);



  // Sum up the entries in the w vector.
  void sumweights();

  // Print out the recurrence relation constants (the b_i and c_i)
  void print_constants();
  
  // Copy constructor. Unimplemented.
  Jacobi(const Jacobi&);

  // Assignment operator.  Unimplemented
  Jacobi& operator=(const Jacobi&);
};

#endif
