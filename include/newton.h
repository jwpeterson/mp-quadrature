#ifndef NEWTON_H
#define NEWTON_H

// Forward declares
class SolverData;

// Given an initial guess u, uses Newton's method to see if it can
// obtain a solution.  If the Newton iterations fail, this is reported
// back to the user.
bool newton(SolverData & solver_data);

// Use a quasi-Newton minimization algorithm to try and
// find a mininum of the scalar function f = (1/2) \vec{r} \cdot
// \vec{r}, where \vec{r} is the residual vector. The Hessian
// is approximated by finite differencing the Jacobian.
bool newton_min(SolverData & solver_data);

#endif

