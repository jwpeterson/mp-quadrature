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
  Integrand(const mpfr_class & x0_in) :
    x0(x0_in)
  {}

  // Destructor
  virtual ~Integrand() {}

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
  mpfr_class H1_norm() const
  {
    return sqrt(this->L2_norm_squared() + this->H1_semi_norm_squared());
  }
};



// Class representing the function |x-x_0|^{alpha}. The function is in
// H1 if alpha>1/2.
struct AlphaSingularity : Integrand
{
  AlphaSingularity(const mpfr_class & x0_in,
                   const mpfr_class & alpha_in) :
    Integrand(x0_in),
    alpha(alpha_in)
  {
    if (alpha <= mpfr_class(1.)/2.)
      {
        std::cerr << "alpha must be greater than 1/2." << std::endl;
        std::abort();
      }
  }

  mpfr_class alpha;

  // f = abs(x-x_0)^{\alpha}, \alpha > \frac{1}{2}
  virtual mpfr_class f(const mpfr_class & x) const
  {
    return pow(abs(x - x0), alpha);
  }

  // I(f) = 1/(\alpha+1) [ (1+x_0)^{\alpha+1} + (1-x_0)^{\alpha+1} ]
  virtual mpfr_class exact() const
  {
    return mpfr_class(1.)/(alpha+1.) * (pow(1+x0, alpha+1) + pow(1-x0, alpha+1));
  }

  // Compute the L2-norm, squared
  virtual mpfr_class L2_norm_squared() const
  {
    return mpfr_class(1.)/(2.*alpha+1.) * (pow(1+x0, 2.*alpha+1) + pow(1-x0, 2.*alpha+1));
  }

  // Compute the H1-semi norm, squared
  virtual mpfr_class H1_semi_norm_squared() const
  {
    return alpha*alpha/(2.*alpha - 1.) * (pow(1+x0, 2.*alpha-1) + pow(1-x0, 2.*alpha-1));
  }
};



// Class representing the piecewise integrand
// { 1 + (x-x0)/(1+x0)
// { 1 - (x-x0)/(1-x0)
struct Piecewise : Integrand
{
  Piecewise(const mpfr_class & x0_in) :
    Integrand(x0_in)
  {}

  // Return the function value
  virtual mpfr_class f(const mpfr_class & x) const
  {
    // Handle edge cases without dividing by zero.
    if (x0 == -1.)
      return 0.5*(1.-x);
    if (x0 == 1.)
      return 0.5*(1.+x);

    // Handle in-between cases
    if (x <= x0)
      return mpfr_class(1.) + (x-x0)/(1.+x0);
    else
      return mpfr_class(1.) - (x-x0)/(1.-x0);
  }

  // I(f) = 1
  virtual mpfr_class exact() const
  {
    return 1.;
  }

  // Compute the L2-norm, squared
  virtual mpfr_class L2_norm_squared() const
  {
    return mpfr_class(2.)/3.;
  }

  // Compute the H1-semi norm, squared
  virtual mpfr_class H1_semi_norm_squared() const
  {
    // Handle edge cases without dividing by zero.
    if (x0 == -1. || x0 == 1.)
      return 0.5;

    return mpfr_class(1.)/(1. + x0) + mpfr_class(1.)/(1. - x0);
  }
};




// Class representing the exponential integrand
// exp(-alpha*|x-x0|)
struct Exponential : Integrand
{
  Exponential(const mpfr_class & x0_in,
              const mpfr_class & alpha_in) :
    Integrand(x0_in),
    alpha(alpha_in)
  {}

  mpfr_class alpha;

  // Return the function value
  virtual mpfr_class f(const mpfr_class & x) const
  {
    return exp(-alpha*abs(x-x0));
  }

  // True integral
  virtual mpfr_class exact() const
  {
    return 2./alpha * (1. - exp(-alpha)*cosh(alpha*x0));
  }

  // Compute the L2-norm, squared
  virtual mpfr_class L2_norm_squared() const
  {
    return 1./alpha * (1. - exp(-2.*alpha)*cosh(2.*alpha*x0));
  }

  // Compute the H1-semi norm, squared
  virtual mpfr_class H1_semi_norm_squared() const
  {
    return alpha * (1. - exp(-2.*alpha)*cosh(2.*alpha*x0));
  }
};



// Class representing the sinusoidal integrand
// sin(alpha*(x-x0))
struct Sinusoidal : Integrand
{
  Sinusoidal(const mpfr_class & x0_in,
             const mpfr_class & alpha_in) :
    Integrand(x0_in),
    alpha(alpha_in)
  {}

  mpfr_class alpha;

  // Return the function value
  virtual mpfr_class f(const mpfr_class & x) const
  {
    return sin(alpha*(x-x0));
  }

  // True integral
  virtual mpfr_class exact() const
  {
    return -2./alpha * sin(alpha) * sin(alpha * x0);
  }

  // Compute the L2-norm, squared
  virtual mpfr_class L2_norm_squared() const
  {
    return 1. - sin(2*alpha) * cos(2*alpha*x0) / (2.*alpha);
  }

  // Compute the H1-semi norm, squared
  virtual mpfr_class H1_semi_norm_squared() const
  {
    return alpha*alpha * (1. + sin(2*alpha) * cos(2*alpha*x0) / (2.*alpha));
  }
};

#endif
