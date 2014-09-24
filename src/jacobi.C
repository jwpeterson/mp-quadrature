#include <cassert>
#include <iomanip> // std::setw
#include <cmath> // fabs
#include <algorithm> // std::sort

#include "jacobi.h"


Jacobi::Jacobi(Real alpha, Real beta)
  : mp_alpha(alpha),
    mp_beta(beta),
    //n(0),
    p(0.0),
    dp(0.0),
    pnm1(0.0)
{
}




void Jacobi::rule(unsigned int n)
{
  // Allocate space for the weights.  As with everything
  // else, allocate one extra space and skip the [0] index.
  w.resize(n+1);
  
  // Compute the quadrature points for a degree n rule.
  this->points(n);

  // Sort the points in increasing order.  This is just
  // to make the printed out data match what other programs
  // will expect.  Be careful when sorting not to sort the
  // [0] entry!
  std::sort(++(x.begin()), x.end());
  
  // Compute the product of the c_j, from j=2 to n.  The
  // constants have already been computed by the call to
  // points, so this should be a no-op...
  this->constants(n);
  //this->print_constants();
  
  mpfr_class cprod(1.0);
  for (unsigned int j=2; j<=n; ++j)
    cprod *= c[j];
  
  // For each point, compute the weight
  for (unsigned int j=1; j<=n; ++j)
    {
      this->value(x[j], n);
      
      // Need to check dp, pnm1 against NaN?
      // std::cout << "dp=" << dp << std::endl;
      // std::cout << "pnm1=" << pnm1 << std::endl;
      
      w[j] = cprod / dp / pnm1;
    }

  // Rescale the rule to lie on [0,1] and have weights
  // which sum to 1.  This is to be handled by the user
  // through the public this->scale_weights() and
  // this->scale_points() interfaces.
  // this->scale();

  // Compute the sum of the weights as a check on the algorithm.
  this->sumweights();
}




void Jacobi::sumweights()
{
  mpfr_class sumweights(0.0);
  for (unsigned int j=1; j<w.size(); ++j)
    sumweights += w[j];
  
  //std::cout << "Sum of weights=" << sumweights << std::endl;
}



void Jacobi::scale_weights(const mpfr_class& scale_factor)
{
  for (unsigned int j=1; j<w.size(); ++j)
    w[j] *= scale_factor;
}


void Jacobi::scale_points(const mpfr_class& x1, const mpfr_class& x2)
{
  //{
  // x[j] = half*(x[j] + one);
  // w[j] /= three;
  // x[j] = 0.5*(x[j] + 1.0);
  // w[j] *= 0.5; // alpha=1 case
  //}

  mpfr_class a ( 0.5*(x2-x1) );
  mpfr_class b ( 0.5*(x1+x2) );
  for (unsigned int j=1; j<x.size(); ++j)
    x[j] = a*x[j] + b;
}



// Reset any computed data which depends on alpha and beta.
void Jacobi::set_alpha(Real alpha)
{
  mp_alpha = alpha;
  this->reset();
}



// Reset any computed data which depends on alpha and beta.
void Jacobi::set_beta(Real beta)
{
  mp_beta = beta;
  this->reset();
}



void Jacobi::constants(unsigned int n)
{
  assert (n>=1);
  
  // Check for easy return: we may have already computed and stored
  // enough b and c coefficients.
  if ((b.size() >= n+1) && (c.size() >= n+1))
    return;

  // Reserve space in b,c.  Note: request one
  // additional storage space because we skip the [0] entry and
  // the user can still refer to b[i] without worrying about indexing.
  b.reserve(n+1);
  c.reserve(n+1);

  //std::cout << "b.size()=" << b.size() << std::endl;
    
  // Otherwise, compute entries through n
  for (unsigned int i=b.size(); i<n+1; ++i)
    {
      // debugging
      //std::cout << "Computing b[i], c[i]  for i=" << i << std::endl;
      
      // No entry for i==0: b[i], c[i] defined for i=1,2,3,...
      // This is easier than trying to rescale all the formulas.
      if (i==0)
	{
	  b.push_back( 0.0 );
	  c.push_back( 0.0 );
	  continue;
	}
      
      const Real ireal = static_cast<Real>(i);

      b.push_back( (mp_alpha+mp_beta) / (mp_alpha+mp_beta+2.*ireal) *
		   (mp_beta-mp_alpha) / (mp_alpha+mp_beta+2.*ireal-2.) );
      
      c.push_back( 4.*(ireal-1.) / (mp_alpha+mp_beta+2.*ireal-1.) * (mp_alpha+ireal-1.) /
		   (mp_alpha+mp_beta+2.*ireal-2.)/(mp_alpha+mp_beta+2.*ireal-2.) *
		   (mp_beta+ireal-1.) / (mp_alpha+mp_beta+2.*ireal-3.) * (mp_alpha+mp_beta+ireal-1.) );
    }

  // Debugging
  //this->print_constants();
}



// This function computes the values of the Jacobi
// polynomial P(alpha,beta,n) and its derivative at x
// using the recursion given by Stroud:
//
// P_n = (x-b_n)*P_{n-1} - c_n*P_{n-2} , n >= 2
//
// The polynomials are defined on [-1 1].  The routine is initialized by
// P_0 = 1
// P_1 = x + (a-b)/(a+b+2)
//
// And the coefficients b_n and c_n take care of the rest of the scaling.
// * Uses for this function include: Newton's Method for root finding.
// * It may be more accurate than pre-computing all n+1 polynomial coeffs
//   and evaluating the polynomial that way.
// * This works with x scalar or an array.
// See: Stroud and Secrest, "Gaussian Quadrature Formulas," 1966, Prentice Hall
void Jacobi::value(const mpfr_class& xval, unsigned int n)
{
  // Special return for n=0 and 1
  if (n==0)
    {
      p=1.;
      dp=0.;
      pnm1=0.; // Not well-defined for n==0
      return;
    }


  if (n==1)
    {
      p=0.5*(mp_alpha+mp_beta+2.)*xval + 0.5*(mp_alpha-mp_beta);
      dp=0.5*(mp_alpha+mp_beta+2.);
      pnm1=1.;
      return; 
    }

#ifdef DEBUG
  std::cout << "================================================================================" << std::endl;
  std::cout << "Jacobi::value(n=" << n << ") - Pre-initialized Values" << std::endl;
  
  // Recurrence relation polynomial values.
  std::cout << "p=" << p << std::endl;
  std::cout << "pnm1=" << pnm1 << std::endl;
  
  // Recurrence relation derivative values
  std::cout << "dp=" << dp << std::endl;
  std::cout << "x=" << xval << std::endl;
#endif

      
  // Initialize for the recurrence relation:
  // P_{n-2} and P'_{n-2}
  mpfr_class pnm2 = 1.0, dpnm2 = 0.0;

  // P_{n-1} and P'_{n-1}
  pnm1 = xval + (mp_alpha-mp_beta)/(mp_alpha+mp_beta+2.); // Underflows ?!
  mpfr_class dpnm1 = 1.0;

#ifdef DEBUG
      std::cout << "Jacobi::value(n=" << n << ") - Initial Values" << std::endl;

      // Recurrence relation polynomial values.
      std::cout << "p=" << p << std::endl;
      std::cout << "pnm1=" << pnm1 << std::endl;
      std::cout << "pnm2=" << pnm2 << std::endl;

      // Recurrence relation derivative values
      std::cout << "dp=" << dp << std::endl;
      std::cout << "dpnm1=" << dpnm1 << std::endl;
      std::cout << "dpnm2=" << dpnm2 << std::endl;
#endif
  
  // Compute all recurrence relation constants for this n
  this->constants(n);

  // Use the recurrence formula to construct the desired
  // polynomial value.
  for (unsigned int j=2; j<=n; ++j)
    {
      p  = (xval-b[j])*pnm1 - c[j]*pnm2;
      dp = (xval-b[j])*dpnm1 + pnm1 - c[j]*dpnm2;

#ifdef DEBUG
      std::cout << "Jacobi::value(n=" << n << ", j=" << j << ")" << std::endl;

      // Recurrence relation polynomial values.
      std::cout << "p=" << p << std::endl;
      std::cout << "pnm1=" << pnm1 << std::endl;
      std::cout << "pnm2=" << pnm2 << std::endl;

      // Recurrence relation derivative values
      std::cout << "dp=" << dp << std::endl;
      std::cout << "dpnm1=" << dpnm1 << std::endl;
      std::cout << "dpnm2=" << dpnm2 << std::endl;
#endif
      
      // If we're at the last step, return before
      // updating values.
      if (j==n)
	return; // or break; to print return values

      // Update polynomial values
      pnm2 = pnm1;
      pnm1 = p;
  
      // Update polynomial derivs
      dpnm2 = dpnm1;
      dpnm1 = dp;
    }

  
  // Debugging
  // std::cout << "p(x)=" << p << std::endl;
  // std::cout << "dp(x)=" << dp << std::endl;
  // std::cout << "pnm1(x)=" << pnm1 << std::endl;
}



void Jacobi::reset()
{
  b.clear();
  c.clear();
}


void Jacobi::points(unsigned int n)
{
  // Make enough room to store all the roots.  (We skip the [0] index)
  // The default initialization for these values is @NaN@.
  x.resize(n+1);

  // Real versions of a, b, n
  const Real nreal = static_cast<Real>(n);
  const Real a = mp_alpha.get_d();
  const Real b = mp_beta.get_d();
  const Real an = a / nreal;
  const Real bn = b / nreal;
  
  // We probably don't need to compute the guess in multi-precision, but
  // the returned root will be an MP object.
  //mpfr_class guess (0.0);
  
  for (unsigned int j=1; j<=n; ++j)
    {
      // Stroud's "largest zero," the root closest to x=1
      if (j==1)
	{
	  const Real R1 = (1.+a)*(2.78/(4.+nreal*nreal) + .768*an/nreal);
	  const Real R2 = 1. + 1.48*an + .96*bn + .452*an*an + .83*an*bn;
	  x[j] = 1.- R1/R2;
	}

      // Stroud's "second zero"
      if (j==2)
	{
	  const Real R1 = (4.1+a)/((1.+a)*(1.+.156*a));
	  const Real R2 = 1. + .06*(nreal-8.)*(1.+.12*a)/nreal;
	  const Real R3 = 1. + .012*b*(1. + .25*fabs(a))/nreal;
	  x[j] = x[1] - R1*R2*R3*(1.-x[1]);
	}

      // Stroud's "third zero"
      if (j==3)
	{
	  const Real R1 = (1.67 + .28*a)/(1. + .37*a);
	  const Real R2 = 1. + .22*(nreal-8.)/nreal;
	  const Real R3 = 1. + 8.*b/((6.28 + b)*nreal*nreal);
	  x[j] = x[2] - R1*R2*R3*(x[1] - x[2]);
	}

      // Stroud's "middle zeros"
      if ((j > 3) && (j<n-1))
	x[j] = 3.*x[j-1] - 3.*x[j-2] + x[j-3];

      // Stroud's "second last zero"
      if ((j==n-1) && (n>3))
	{
	  const Real R1 = (1. + .235*b)/(.766 + .119*b);
	  const Real R2 = 1. / (1. + .639*(nreal-4.)/(1. + .71*(nreal-4.)));  
	  const Real R3 = 1. / (1. + 20.*a/((7.5+a)*nreal*nreal) );
	  x[j] = x[j-1] + R1*R2*R3*(x[j-1] - x[j-2]);
	}

      // Stroud's "last zero"
      if ((j==n) && (n>3))
	{
	  const Real R1 = (1. + .37*b)/(1.67 + .28*b);
	  const Real R2 = 1. / (1. + .22*(nreal-8.)/nreal);
	  const Real R3 = 1. / (1. + 8.*a/((6.28+a)*nreal*nreal));
	  x[j] = x[j-1] + R1*R2*R3*(x[j-1] - x[j-2]);
	}

      
      // Compute root using x[j] as initial guess.
      this->newton(n, x[j]);

      // Debugging
      // std::cout << "x[" << j << "]=" << x[j] << std::endl;

      // What is the best way to check for NaN?
//       if ( isnan(x[j]) )
// 	{
// 	  std::cerr << "Polynomial root is NaN!" << std::endl;
// 	  abort();
// 	}
    }
}




void Jacobi::newton(unsigned int n, mpfr_class& xroot)
{
  // Not sure if the residual and tolerance need to be multi-precision.
  mpfr_class tol = 1.e-25;
  mpfr_class residual = 1.;

  // Iteration counters
  unsigned int n_its=0;
  const unsigned int max_its=30;

  while ((residual > tol) && (n_its<max_its))
    {
      // Compute Jacobi values at the current x-value
      this->value(xroot, n);

      // Compute the Newton update
      xroot = xroot - (p/dp);

      // Compute the residual (the distance of p(x) from zero) 
      residual = abs(p);

      // Debugging: check quadratic Newton convergence
      // std::cout << "residual[" << n_its << "]=" << residual << std::endl;
      
      // Update the iteration count
      n_its++;
    }

  if ((n_its==max_its) && (residual > tol))
    {
      std::cerr << "Iterations did not converge!" << std::endl;
      abort();
    }

  // Debugging
  // std::cout << "Newton converged in: " << n_its << " iterations." << std::endl;
  // std::cout << "With residual: " << residual << std::endl;
}



void Jacobi::print_constants()
{
  for (unsigned int j=1; j<b.size(); ++j)
    {
      std::cout << "b["<<j<<"]=" << b[j] << "\t";
      std::cout << "c["<<j<<"]=" << c[j] << std::endl;
    }
}
