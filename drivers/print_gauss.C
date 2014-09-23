#include <iostream>
#include <limits>
#include <cmath>
#include <cstdio>
#include <vector>
#include <iomanip>

// The GNU multi-precision library
#include "gmp.h"
#include "mpfr.h"
#include "gmpfrxx.h"

// Header files for this project
#include "common_definitions.h"
#include "gauss.h"

// This program implements what is probably a Numerical
// recipes algorithm for computing Legendre polynomial
// zeros via Newton's method and computing weights for
// Gaussian quadrature formulae.  We're using a multi-precision
// library to achieve a larger number of significant
// digits.

int main(int argc, char** argv)
{
  std::cout.precision(32);
  std::cout.setf(std::ios_base::scientific);

  // # of binary digits
  // 53 binary digits is about what you normally get with a double.
  mpfr_set_default_prec(256); 

  unsigned int n=6;
  // Read number of points in rule from command line
  if (argc > 1)
    n=atoi(argv[1]);
  
  // Probably failed to read the number of points from the command line...
  if (n==0)
    {
      std::cout << "Warning, could not determine valid rule order from command line." << std::endl;
      std::cout << "Running with default order 6 rule." << std::endl;
      n=6;
    }

  std::cout << "\nComputing Gauss rule for n=" << n << std::endl;
  std::vector<mpfr_class> x, w;
  gauss_rule(n, x, w);

  // For certain rules we have analytical expressions for the points
  // and weights.
  bool analytical_check=false;
  
  if ((n>=2) && (n<=5))
    analytical_check = true;

  if (analytical_check)
    {
      std::vector<mpfr_class> x_analytical, w_analytical;

      switch (n)
	{
	case 2:
	  {
	    // x= +/-sqrt(3)/3
	    x_analytical.resize(2);
	    x_analytical[0] = 3.0;
	    x_analytical[0] = sqrt(x_analytical[0]);
	    x_analytical[0] /= 3.0;
	    x_analytical[0] *= -1.0;
	    
	    x_analytical[1] = -1.0*x_analytical[0];

	    // w=1 for this case
	    w_analytical.resize(2);
	    w_analytical[0] = 1.0;
	    w_analytical[1] = 1.0;
	    break;
	  }

	case 3:
	  {
	    // x= 0, +/-sqrt(15)/5
	    x_analytical.resize(3);
	    x_analytical[0] = 15.0;
	    x_analytical[0] = sqrt(x_analytical[0]);
	    x_analytical[0] /= 5.0;
	    x_analytical[0] *= -1.0;

	    x_analytical[1] = 0.0;
	    x_analytical[2] = -1.0 * x_analytical[0];

	    // w=8/9, 5/9
	    w_analytical.resize(3);
	    w_analytical[0] = 5.0;
	    w_analytical[0] /= 9.0;

	    w_analytical[1] = 8.0;
	    w_analytical[1] /= 9.0;
	    
	    w_analytical[2] = w_analytical[0];

	    break;
	  }

	case 4:
	  {
	    // x= +/- sqrt(525 - 70*sqrt(30))/35,
	    //    +/- sqrt(525 + 70*sqrt(30))/35,
	    x_analytical.resize(4);
	    mpfr_class seventy_root_thirty = 30.0;
	    seventy_root_thirty = sqrt(seventy_root_thirty);
	    seventy_root_thirty *= 70.0;
	  
	    x_analytical[0] = 525. + seventy_root_thirty;
	    x_analytical[0] = sqrt(x_analytical[0]);
	    x_analytical[0] /= 35.0;
	    x_analytical[0] *= -1.0;

	    x_analytical[1] = 525. - seventy_root_thirty;
	    x_analytical[1] = sqrt(x_analytical[1]);
	    x_analytical[1] /= 35.0;
	    x_analytical[1] *= -1.0;

	    x_analytical[2] = -1.0 * x_analytical[1];
	    x_analytical[3] = -1.0 * x_analytical[0];
	  
	    // w= 1/36 * (18 + sqrt(30)),
	    //    1/36 * (18 - sqrt(30))
	    w_analytical.resize(4);
	    mpfr_class sqrt_thirty = 30.0;
	    sqrt_thirty = sqrt(sqrt_thirty);
	  
	    w_analytical[0] = 18. - sqrt_thirty;
	    w_analytical[0] /= 36.0;

	    w_analytical[1] = 18. + sqrt_thirty;
	    w_analytical[1] /= 36.0;

	    w_analytical[2] = w_analytical[1];
	    w_analytical[3] = w_analytical[0];
	  
	    break;
	  }

	case 5:
	  {
	    // x = +/- sqrt(245 +/- 14*sqrt(70))/21
	    x_analytical.resize(5);
	    mpfr_class fourteen_root_seventy = 70.0;
	    fourteen_root_seventy = sqrt(fourteen_root_seventy);
	    fourteen_root_seventy *= 14.0;

	    // We'll print the negative roots, since that is what the
	    // algorithm computes
	    x_analytical[0] = 245.0 + fourteen_root_seventy;
	    x_analytical[0] = sqrt(x_analytical[0]);
	    x_analytical[0] /= 21.0;
	    x_analytical[0] *= -1.0;

	    x_analytical[1] = 245.0 - fourteen_root_seventy;
	    x_analytical[1] = sqrt(x_analytical[1]);
	    x_analytical[1] /= 21.0;
	    x_analytical[1] *= -1.0;

	    x_analytical[2] = 0.0;

	    x_analytical[3] = -1.0*x_analytical[1];
	    x_analytical[4] = -1.0*x_analytical[0];

	    // w = (322 +/- 13*sqrt(70))/900, 128/225
	    w_analytical.resize(5);
	    mpfr_class thirteen_root_seventy = 70.0;
	    thirteen_root_seventy = sqrt(thirteen_root_seventy);
	    thirteen_root_seventy *= 13.0;

	    w_analytical[0] = 322. - thirteen_root_seventy;
	    w_analytical[0] /= 900.;
	  
	    w_analytical[1] = 322. + thirteen_root_seventy;
	    w_analytical[1] /= 900.;

	    w_analytical[2] = 128.0;
	    w_analytical[2] /= 225.;

	    w_analytical[3] = w_analytical[1];
	    w_analytical[4] = w_analytical[0];

	    break;
	  }

	default:
	  {
	    // No analytical form known
	  }
	}

      // Compute error in x, w.  Note, the values returned in x, w are 1-based. 
      if ((x.size()-1 == x_analytical.size()) && (w.size()-1 == w_analytical.size()))
	{
	  std::cout << "\nError in points and weights: " << std::endl;
	  
	  for (unsigned int i=0; i<x_analytical.size(); ++i)
	    {
	      mpfr_class delta_x = x[i+1] - x_analytical[i];
	      mpfr_class delta_w = w[i+1] - w_analytical[i];

	      std::cout << "delta_x[" << i << "]=" << delta_x << ", "
			<< "delta_w[" << i << "]=" << delta_w << std::endl;
	    }
	  std::cout << std::endl;
	}
      else
	{
	  std::cout << "Different number of points computed for numerical and analytical rules!" << std::endl;
	  std::cout << "x.size()=" << x.size() << ", x_analytical.size()=" << x_analytical.size() << std::endl;
	}


      
      for (unsigned int i=0; i<x_analytical.size(); ++i)
	std::cout << "x_analytical["<<i<<"]=" << x_analytical[i] << std::endl;
      std::cout << std::endl;
      for (unsigned int i=0; i<w_analytical.size(); ++i)
	std::cout << "w_analytical["<<i<<"]=" << w_analytical[i] << std::endl;
    } // end if analytical_check
  
  // Print points and weights using standard C-style indexing and
  // by setting the second half of the points and weights equal
  // to the mirror of the first half.
  const unsigned int m=(n+1)/2;

  // Side-by-side output.  This is a little hard to get right...
//   for (unsigned int i=1; i<=m; ++i)
//     {
//       std::cout << "_points[" << i-1 << "](0)=" << x[i] << ";\t";
//       std::cout << "_weights[" << i-1 << "]=" << w[i] << std::endl;
//     }
//   for (unsigned int i=m+1,d=m-1; i<=n; ++i,--d)
//     {
//       std::cout << "_points[" << i-1 << "](0)=-_points[" << d-1 << "];\t";
//       std::cout << "_weights[" << i-1 << "]=_weights[" << d-1 << "];" << std::endl;
//     }

  
  // Points first, then weights.  This is how the library currently is...
  std::cout << std::endl;

  for (unsigned int i=1; i<=m; ++i)
    std::cout << "_points[" << std::setw(2) << i-1 << "](0) = " << x[i] << ";\n";
  for (unsigned int i=m+1,d=(n%2)?m-1:m; i<=n; ++i,--d)
    std::cout << "_points[" << std::setw(2) << i-1 << "]    = -_points[" << d-1 << "];\n";

  std::cout << std::endl;
  
  for (unsigned int i=1; i<=m; ++i)
    std::cout << "_weights[" << std::setw(2) << i-1 << "]   = " << w[i] << ";\n";
  for (unsigned int i=m+1,d=(n%2)?m-1:m; i<=n; ++i,--d)
    std::cout << "_weights[" << std::setw(2) << i-1 << "]   = _weights[" << d-1 << "];" << std::endl;
  
  return 0;
}
  
