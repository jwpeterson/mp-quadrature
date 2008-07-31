#include "gauss.h"

void gauss_rule(unsigned int n, // Number of points (not order!)
		std::vector<mpfr_class>& x, // Quadrature points
		std::vector<mpfr_class>& w) // Quadrature weights
{
  // Allocate space for points and weights.  As is the
  // case with a lot of these numerical codes, we skip the
  // [0] entry of the vectors...
  x.resize(n+1);
  w.resize(n+1);
  
  // Find only half the roots because of symmetry
  const unsigned int m=(n+1)/2; 

  // Maximum number of Newton iterations allowed
  const unsigned int max_its = 30;

  // Tolerance to be use in Newton's method
  const mpfr_class tol = 1.e-30;
  
  // Three recurrence relation values and one derivative value.
  mpfr_class
    pn(0.0),  // P_{n}
    pnm1(0.0),// P_{n-1}
    pnm2(0.0),// P_{n-2}
    dpn(0.0); // P'_{n}

  for (unsigned int i=1; i<=m; ++i)
    {
      // If n is odd, the rule always contains the point x==0.
      // In this case we don't need Newton iterations.
      const bool skip_newton = (n%2) && (i==m);
      if ( skip_newton )
	{
	  x[i] = 0.0;
	}
	  
      else
	{
	  // Remarkably, this simple relation provides a very
	  // good initial guess for x_i.  See, for example,
	  // F. G. Lether and P. R. Wenston
	  // Journal of Computational and Applied Mathematics  
	  // Minimax approximations to the zeros of Pn(x) and
	  // Gauss-Legendre quadrature
	  // Volume 59 ,  Issue 2  (May 1995) table of contents
	  // Pages: 245 - 252  
	  // 1995 
	  x[i] = cos(3.141592654*(i-0.25)/(n+0.5));
	}
	  
      // Newton loop iteration counter
      unsigned int n_its = 0;

      // do loop stopping boolean
      bool continue_while=true;
      
      // Begin Newton iterations
      do
	{
	  // Initialize recurrence relation
	  pn=1.0; pnm1=0.0;

	  // Use recurrence relation to compute p(x[i])
	  for (unsigned int j=1; j<=n; ++j)
	    {
	      pnm2=pnm1; pnm1=pn;
	      pn=((2.0*j-1.0)*x[i]*pnm1 - (j-1.0)*pnm2)/static_cast<Real>(j); 
	    }
	  
	  // A recurrence relation also gives the derivative.
	  dpn=n*(x[i]*pn-pnm1)/(x[i]*x[i]-1.0);

	  // Compute Newton update
	  if (!skip_newton)
	    x[i] -= pn/dpn;

	  // Increment counter
	  n_its++;

	  // Determine while loop exit status
	  continue_while = (abs(pn) > tol) && (n_its < max_its);

	  // Don't do any more iterations if we skipped Newton!
	  if (skip_newton)
	    continue_while = false;
	  
	} while ( continue_while );

      // Test for convergence failure
      if (n_its>=max_its)
	{
	  std::cerr << "Error! Max iterations reached!" << std::endl;
	  abort();
	}

      // Debugging/status info
      //       if (!skip_newton)
      // 	{
      // 	  std::cout << "Newton converged in " << n_its << " iterations, ";
      // 	  std::cout << "with tolerance=" << abs(pn) << std::endl;
      // 	}

      // Set x[i] and its mirror image.  We set these in increasing order.
      x[i]     = -x[i];
      x[n+1-i] = -x[i];

      // Compute the weight w[i], its mirror is the same value
      w[i]     = 2.0 / ( (1.0-x[i]*x[i])*dpn*dpn );
      w[n+1-i] = w[i];
    } // end for

  // Debugging/verification info
  mpfr_class sumweights (0.0);
  for (unsigned int j=1; j<w.size(); ++j)
    sumweights += w[j];
  
  //std::cout << "Sum of weights=" << sumweights << std::endl;
  
}
