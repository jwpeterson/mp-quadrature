#ifndef SOLVER_DATA_H
#define SOLVER_DATA_H

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

// Parameters that control the behavior of solvers.
struct SolverData
{
  SolverData() :
    tol(1.e-36),
    divtol(1.e16),
    maxits(20),
    alphamin(1.e-3),
    verbose(false),
    residual_and_jacobian(nullptr) {}

  // Initial guess/solution vector
  std::vector<mpfr_class> u;
  mpfr_class tol;
  mpfr_class divtol;
  unsigned int maxits;
  // smallest backtracking backtracking linesearch fraction
  mpfr_class alphamin;
  bool verbose;
  ResidualAndJacobianFunctionPtr residual_and_jacobian;
};

#endif

