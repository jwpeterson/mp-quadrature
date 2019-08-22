#ifndef NEWTON_H
#define NEWTON_H

// Local includes
#include "matrix.h"

// Support library headers
#include "mpfr.h"
#include "gmpfrxx.h"

// C++ includes
#include <vector>

// Create a typedef for working with function pointers
using ResidualAndJacobianFunctionPtr =
  void (*) (std::vector<mpfr_class> * r,
            Matrix<mpfr_class> * jac,
            const std::vector<mpfr_class> & u);

// Given an initial guess u, uses Newton's method to see if it can
// obtain a solution.  If the Newton iterations fail, this is reported
// back to the user.
bool newton(std::vector<mpfr_class> & u,
            ResidualAndJacobianFunctionPtr residual_and_jacobian);

// Use a quasi-Newton minimization algorithm to try and
// find a mininum of the scalar function f = (1/2) \vec{r} \cdot
// \vec{r}, where \vec{r} is the residual vector. The Hessian
// is approximated by finite differencing the Jacobian.
bool newton_min(std::vector<mpfr_class> & u,
                ResidualAndJacobianFunctionPtr residual_and_jacobian);

#endif

