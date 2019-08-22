#ifndef NLCG_H
#define NLCG_H

// Forward declares
class SolverData;

// Use nonlinear conjugate gradient method to search for a minimum
// https://en.wikipedia.org/wiki/Nonlinear_conjugate_gradient_method
bool nlcg(SolverData & solver_data);

#endif
