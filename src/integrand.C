#include <cstdlib> // std::abort()
#include "integrand.h"

Integrand::Integrand(const mpfr_class & x0_in) :
  x0(x0_in)
{}


Integrand::~Integrand()
{}


mpfr_class Integrand::H1_norm() const
{
  return sqrt(this->L2_norm_squared() + this->H1_semi_norm_squared());
}



AlphaSingularity::AlphaSingularity(const mpfr_class & x0_in,
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



mpfr_class AlphaSingularity::f(const mpfr_class & x) const
{
  return pow(abs(x - x0), alpha);
}



mpfr_class AlphaSingularity::exact() const
{
  return mpfr_class(1.)/(alpha+1.) * (pow(1+x0, alpha+1) + pow(1-x0, alpha+1));
}



mpfr_class AlphaSingularity::L2_norm_squared() const
{
  return mpfr_class(1.)/(2.*alpha+1.) * (pow(1+x0, 2.*alpha+1) + pow(1-x0, 2.*alpha+1));
}



mpfr_class AlphaSingularity::H1_semi_norm_squared() const
{
  return alpha*alpha/(2.*alpha - 1.) * (pow(1+x0, 2.*alpha-1) + pow(1-x0, 2.*alpha-1));
}



Piecewise::Piecewise(const mpfr_class & x0_in) :
  Integrand(x0_in)
{}



mpfr_class Piecewise::f(const mpfr_class & x) const
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



mpfr_class Piecewise::exact() const
{
  return 1.;
}



mpfr_class Piecewise::L2_norm_squared() const
{
  return mpfr_class(2.)/3.;
}



mpfr_class Piecewise::H1_semi_norm_squared() const
{
  // Handle edge cases without dividing by zero.
  if (x0 == -1. || x0 == 1.)
    return 0.5;

  return mpfr_class(1.)/(1. + x0) + mpfr_class(1.)/(1. - x0);
}



Exponential::Exponential(const mpfr_class & x0_in,
                         const mpfr_class & alpha_in) :
  Integrand(x0_in),
  alpha(alpha_in)
{}



mpfr_class Exponential::f(const mpfr_class & x) const
{
  return exp(-alpha*abs(x-x0));
}



mpfr_class Exponential::exact() const
{
  return 2./alpha * (1. - exp(-alpha)*cosh(alpha*x0));
}



mpfr_class Exponential::L2_norm_squared() const
{
  return 1./alpha * (1. - exp(-2.*alpha)*cosh(2.*alpha*x0));
}



mpfr_class Exponential::H1_semi_norm_squared() const
{
  return alpha * (1. - exp(-2.*alpha)*cosh(2.*alpha*x0));
}



Sinusoidal::Sinusoidal(const mpfr_class & x0_in,
                       const mpfr_class & alpha_in) :
  Integrand(x0_in),
  alpha(alpha_in)
{}



mpfr_class Sinusoidal::f(const mpfr_class & x) const
{
  return sin(alpha*(x-x0));
}



mpfr_class Sinusoidal::exact() const
{
  return -2./alpha * sin(alpha) * sin(alpha * x0);
}



mpfr_class Sinusoidal::L2_norm_squared() const
{
  return 1. - sin(2*alpha) * cos(2*alpha*x0) / (2.*alpha);
}



mpfr_class Sinusoidal::H1_semi_norm_squared() const
{
  return alpha*alpha * (1. + sin(2*alpha) * cos(2*alpha*x0) / (2.*alpha));
}
