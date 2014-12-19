#ifndef INTEGRAND_H
#define INTEGRAND_H

#include "gmpfrxx.h"

// These are some 1D integrands that I used to investigate the a
// priori error estimates under varying conditions.

// Class representing a function, its exact integral, L2, and H1 norms.
// We assume that the integrand has at least one parameter named x0.
struct Integrand
{
  // Constructor
  Integrand(const mpfr_class & x0_in);

  // Destructor
  virtual ~Integrand();

  // The parameter
  mpfr_class x0;

  // Returns the value of f at the point x
  virtual mpfr_class f(const mpfr_class & x) const = 0;

  // Returns the exact integral f
  virtual mpfr_class exact() const = 0;

  // Returns the L2-norm, squared
  virtual mpfr_class L2_norm_squared() const = 0;

  // Returns the H1-semi norm, squared
  virtual mpfr_class H1_semi_norm_squared() const = 0;

  // Returns the H1-norm.
  mpfr_class H1_norm() const;
};



// Class representing the function |x-x_0|^{alpha}. The function is in
// H1 if alpha>1/2.
struct AlphaSingularity : Integrand
{
  AlphaSingularity(const mpfr_class & x0_in,
                   const mpfr_class & alpha_in);

  mpfr_class alpha;

  // f = abs(x-x_0)^{\alpha}, \alpha > \frac{1}{2}
  virtual mpfr_class f(const mpfr_class & x) const;

  // I(f) = 1/(\alpha+1) [ (1+x_0)^{\alpha+1} + (1-x_0)^{\alpha+1} ]
  virtual mpfr_class exact() const;

  // Compute the L2-norm, squared
  virtual mpfr_class L2_norm_squared() const;

  // Compute the H1-semi norm, squared
  virtual mpfr_class H1_semi_norm_squared() const;
};



// Class representing the piecewise integrand
// { 1 + (x-x0)/(1+x0)
// { 1 - (x-x0)/(1-x0)
struct Piecewise : Integrand
{
  Piecewise(const mpfr_class & x0_in);

  // Return the function value
  virtual mpfr_class f(const mpfr_class & x) const;

  // I(f) = 1
  virtual mpfr_class exact() const;

  // Compute the L2-norm, squared
  virtual mpfr_class L2_norm_squared() const;

  // Compute the H1-semi norm, squared
  virtual mpfr_class H1_semi_norm_squared() const;
};




// Class representing the exponential integrand
// exp(-alpha*|x-x0|)
struct Exponential : Integrand
{
  Exponential(const mpfr_class & x0_in,
              const mpfr_class & alpha_in);

  mpfr_class alpha;

  // Return the function value
  virtual mpfr_class f(const mpfr_class & x) const;

  // True integral
  virtual mpfr_class exact() const;

  // Compute the L2-norm, squared
  virtual mpfr_class L2_norm_squared() const;

  // Compute the H1-semi norm, squared
  virtual mpfr_class H1_semi_norm_squared() const;
};



// Class representing the sinusoidal integrand
// sin(alpha*(x-x0))
struct Sinusoidal : Integrand
{
  Sinusoidal(const mpfr_class & x0_in,
             const mpfr_class & alpha_in);

  mpfr_class alpha;

  // Return the function value
  virtual mpfr_class f(const mpfr_class & x) const;

  // True integral
  virtual mpfr_class exact() const;

  // Compute the L2-norm, squared
  virtual mpfr_class L2_norm_squared() const;

  // Compute the H1-semi norm, squared
  virtual mpfr_class H1_semi_norm_squared() const;
};

#endif
