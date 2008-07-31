#include <iostream>
#include <limits>
#include <cmath>
#include <cstdio>
#include <vector>

// The GNU multi-precision library
#include "gmp.h"
// #include "gmpxx.h"

// The MPFR library defines special functions
// like sin, cos, exp, etc.
#include "mpfr.h"

// The GMPFRXX library defines a C++ interface for
// the mpfr library.  You cannot include the gmpxx.h
// header *and* the gmpfrxx.h headers at the same time
// though!
#include "gmpfrxx.h"

// Header files for this project
#include "common_definitions.h"
// #include "jacobi_constants.h"
#include "jacobi.h"



int main()
{
  std::cout.precision(32);
  std::cout.setf(std::ios_base::scientific);

  // # of binary digits
  // 53 binary digits is about what you normally get with a double.
  mpfr_set_default_prec(256); 

  // Case 1: weights sum to 1/2
  // const Real alpha=1.0, beta=0.0;
  
  // Case 2: weights sum to 1/3
  const Real alpha=2.0, beta=0.0;

  Jacobi p(alpha, beta);

  for (unsigned int j=2; j<23; ++j)
    {
      std::cout << "================================================================================" << std::endl;
      std::cout << "Jacobi rule with alpha=" << alpha << ", beta=" << beta << ", "
		<< j << " points, order=" << 2*j-1 << std::endl;
      p.rule(j); // order = 2*j-1
      std::cout << "\n";
    }
  
  // p.rule(12); // order 23 
  // p.rule(21); // Happened to need a slightly tighter tolerance...
  // p.rule(22); // order 43


  return 0;
}
  
