#include <stdlib.h> // abort()
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
      // {
      //   std::cout << "Newton converged in " << n_its << " iterations, ";
      //   std::cout << "with tolerance=" << abs(pn) << std::endl;
      // }

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



void gauss_lobatto_rule(unsigned int n,
                        std::vector<mpfr_class>& x,
                        std::vector<mpfr_class>& w)
{
  // TODO:
  //
  // From: http://keisan.casio.com/exec/system/1280801905
  //
  // Definitions:
  //
  // For an n-point rule:
  // nodes: P'_{n-1}(x) = 0
  //        -1, +1
  // weights: w_1 = w_n = 2/n/(n-1)
  //          w_i = 2 / n / (n-1) / P_{n-1}(x_i)^2, i=2..n-1
  //
  // Algorithm: for i=2..n-1:
  // 1.) Initial guess: x_i = (1 - 3*(n-2)/8/(n-1)^3)*cos( (4*i-3)/(4*(n-1)+1)*pi )
  // 2.) Solve P'_{n-1}(x_i) = 0 for x_i
  //     Halley's method: x <- x - 2*y*y'/(2*y'*y' - y*y'')
  //                      where: y    := P'_{n-1}(x_i)
  //                             y'   := P''_{n-1}(x_i)
  //                             y''  := P'''_{n-1}(x_i)
  //                      using the recurrence relation:
  //                      P'_n(x)   = n*(P_{n-1}(x) - x*P_n(x))/(1-x^2)
  //                      P''_n(x)  = (2*x*P'_n(x) - n*(n+1)*P_n(x))/(1-x^2)
  //                      P'''_n(x) = (2*x*P''_n(x) - (n*(n+1)-2)*P'_n(x))/(1-x^2)
  // 3.) Compute the weights:
  //     w_i = 2 / n / (n-1) / P_{n-1}(x_i)^2

  // Allocate space for points and weights.  As is the
  // case with a lot of these numerical codes, we skip the
  // [0] entry of the vectors...
  x.resize(n+1);
  w.resize(n+1);

  // Find only half the roots because of symmetry
  const unsigned int m=(n+1)/2;

  // Maximum number of Halley iterations allowed
  const unsigned int max_its = 30;

  // Tolerance to be use in Halley's method
  const mpfr_class tol = 1.e-30;

  // Three Legendre polynomial values
  mpfr_class
    pn(0.0),  // P_{n}
    pnm1(0.0),// P_{n-1}
    pnm2(0.0); // P_{n-2}

  // Three Legendre derivative values
  mpfr_class
    dpnm1(0.0), // P'_{n-1}
    d2pnm1(0.0), // P''_{n-1}
    d3pnm1(0.0); // P'''_{n-1}

  // The first/last points are always at x=-1/+1 and have weight=2/n/(n-1).
  // Note the explicit casts in the computation of w[1].  If we don't
  // do that, we don't retain the full precision of the mpfr_class
  // type, since the compiler first converts the 2.0 to a double, does
  // double-precision division, and finally converts the result to an
  // mpfr_class object.
  x[1] = -1.0;
  x[n] =  1.0;
  w[1] = w[n] = static_cast<mpfr_class>(2.0) / static_cast<mpfr_class>(n*(n-1));

  // Now compute the rest of the points/weights:
  for (unsigned int i=2; i<=m; ++i)
    {
      // Initial guess - only needs to be double precision.
      x[i] = (1. - 3.*(n-2)/8./(n-1)/(n-1)/(n-1))*cos(3.141592654*(4.*i - 3.)/(4.*(n-1.) + 1.));

      // Halley iteration counter
      unsigned int n_its = 0;

      // do loop stopping boolean
      bool keep_going=true;

      // Begin Halley iterations
      do
        {
          // Initialize Legendre polynomial recurrence relation
          pn = 1.0;
          pnm1 = 0.0;

          // Use recurrence relation to compute up to P_{n}.  Note
          // that we only need up to P_{n-1} for Gauss-Lobatto
          // quadrature, the final pn value will not be used.
          for (unsigned int j=1; j<=n; ++j)
            {
              pnm2 = pnm1;
              pnm1 = pn;
              pn = ((2.0*j-1.0)*x[i]*pnm1 - (j-1.0)*pnm2)/static_cast<mpfr_class>(j);
            }

          // Compute derivatives needed by Halley iteration.  FIXME:
          // This should be a function parameterized on 'n', rather
          // than hard-coded for n-1...

          // P'_{n-1}
          dpnm1 = (n-1)*(x[i]*pnm1 - pnm2)/(x[i]*x[i]-1.0);

          // P''_{n-1}
          d2pnm1 = ((n-1)*n*pnm1 - 2.*x[i]*dpnm1)/(x[i]*x[i]-1.0);

          // P'''_n
          d3pnm1 = (((n-1)*n - 2.)*dpnm1 - 2.*x[i]*d2pnm1)/(x[i]*x[i]-1.0);

          // Save old value of x, so we can compute convergence

          // Compute the Halley update
          x[i] -= 2.*dpnm1*d2pnm1 / (2.*d2pnm1*d2pnm1 - dpnm1*d3pnm1);

          // Increment counter
          n_its++;

          // Debugging: print out current point value
          // std::cout << "x_current = " << x[i] << std::endl;

          // x[i] should be a root of P'_{n-1}, so this determines the convergence check.
          keep_going = (abs(dpnm1) > tol) && (n_its < max_its);

        } while (keep_going);

      // Test for convergence failure
      if (n_its >= max_its)
        {
          std::cerr << "Error! Max iterations reached!" << std::endl;
          abort();
        }

      // Our initial guess makes us solve for roots > 0, but we want
      // the negative roots in the vector first, so swap it.
      x[i]     = -x[i];
      x[n+1-i] = -x[i];

      // Compute the weight w[i], its mirror is the same value
      w[i]     = static_cast<mpfr_class>(2.0) / static_cast<mpfr_class>(n*(n-1)) / pnm1 / pnm1;
      w[n+1-i] = w[i];

      // Debugging: print point and weight
      // std::cout << "x[" << i << "]=" << x[i] << std::endl;
      // std::cout << "w[" << i << "]=" << w[i] << std::endl;
    } // end for

  // Debugging/verification info
  mpfr_class sumweights (0.0);
  for (unsigned int j=1; j<w.size(); ++j)
    {
      // std::cout << "w[j]=" << w[j] << std::endl;
      sumweights += w[j];
    }

  std::cout << "Sum of weights=" << sumweights << std::endl;
}
