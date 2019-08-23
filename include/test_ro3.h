#ifndef TEST_RO3_H
#define TEST_RO3_H

// Forward declaration
class SolverData;

// Generates random initial guesses for an Ro3-invariant quadrature
// rule which would be exact for polynomials of degree d and tests to
// see whether any of them converge to a valid solution.
void test_ro3(unsigned int d,
              unsigned int n_tests,
              SolverData & solver_data);

#endif
