#ifndef SOLVER_DATA_H
#define SOLVER_DATA_H

// Local includes
#include "matrix.h"
#include "exact.h"
#include "vect.h"
#include "ro3.h"

// Support library headers
#include "mpfr.h"
#include "gmpfrxx.h"

// C++ includes
#include <vector>

// Parameters that control the behavior of solvers.
struct SolverData
{
  SolverData(const Ro3 & ro3_in) :
    tol(1.e-36),
    divtol(1.e16),
    maxits(20),
    alphamin(1.e-3),
    do_backtracking(false),
    residual_reduction_required(false),
    verbose(false),
    ro3(ro3_in)
  {}

  // Initial guess/solution vector, residual, and Jacobian.
  std::vector<mpfr_class> u;
  std::vector<mpfr_class> r;
  Matrix<mpfr_class> jac;

  // Tolerances and other solver parameters.
  mpfr_class tol;
  mpfr_class divtol;
  unsigned int maxits;
  // smallest backtracking backtracking linesearch fraction
  mpfr_class alphamin;
  bool do_backtracking;
  bool residual_reduction_required;
  bool verbose;
  Ro3 ro3;
};

#endif
