#ifndef GRADIENT_DESCENT_H
#define GRADIENT_DESCENT_H

// Forward declares
class SolverData;

// Implements the gradient descent method for the scalar function
// f = (1/2) \vec{r} \cdot \vec{r},
// where \vec{r} is the residual vector.
bool gradient_descent(SolverData & solver_data);

#endif
