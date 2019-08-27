#ifndef SOLVER_DATA_H
#define SOLVER_DATA_H

// Local includes
#include "matrix.h"
#include "exact.h"
#include "vect.h"

// Support library headers
#include "mpfr.h"
#include "gmpfrxx.h"

// C++ includes
#include <vector>

// A functor class for computing the residual and Jacobian for
// Ro3-invariant rules which are exact for degree d polynomials.
struct ResidualAndJacobian
{
  ResidualAndJacobian(unsigned int d_in = 1) :
    d(d_in),
    N((d*d + 3*d + 6)/6),
    n_centroid(N % 3),
    n_ro3(N / 3),
    one_third(mpfr_class(1) / 3)
  {
    // Initialize the polynomial exponents array.
    polys =
      {
        {0,0},                                         // const
        {2,0},                                         // 2nd
        {3,0},  {2,1},                                 // 3rd
        {4,0},                                         // 4th
        {5,0},  {4,1},                                 // 5th
        {6,0},  {5,1},  {4,2},                         // 6th
        {7,0},  {6,1},                                 // 7th
        {8,0},  {7,1},  {6,2},                         // 8th
        {9,0},  {8,1},  {7,2},  {6,3},                 // 9th
        {10,0}, {9,1},  {8,2},                         // 10th
        {11,0}, {10,1}, {9,2},  {8,3},                 // 11th
        {12,0}, {11,1}, {10,2}, {9,3},  {8,4},         // 12th
        {13,0}, {12,1}, {11,2}, {10,3},                // 13th
        {14,0}, {13,1}, {12,2}, {11,3}, {10,4},        // 14th
        {15,0}, {14,1}, {13,2}, {12,3}, {11,4}, {10,5} // 15th
      };

    // Throw an error if we have not tabulated enough basis polynomials yet.
    if (polys.size() < N)
      {
        std::cout << "Not enough polynomials for d = " << d << std::endl;
        abort();
      }
  }

  // polynomial degree
  unsigned int d;

  // Number of unknowns in the residual vector
  unsigned int N;

  // The number of centroid points
  unsigned int n_centroid;

  // The number of Ro3-invariant points
  unsigned int n_ro3;

  // A handy multi-precision constant
  mpfr_class one_third;

  // The list of monomial exponents that we are going to test.
  // We must integrate the first N of these exactly
  std::vector<std::pair<unsigned int, unsigned int>> polys;

  // Compute the residual and Jacobian at u.
  void operator() (std::vector<mpfr_class> * r,
                   Matrix<mpfr_class> * jac,
                   const std::vector<mpfr_class> & u)
  {
    // If there's nothing to do, then there's nothing to do.
    if (!r && !jac)
      return;

    // Zero any previous values and allocate space for residual and Jacobian, if required.
    if (r)
      {
        r->clear();
        r->resize(N);
      }
    if (jac)
      {
        jac->clear();
        jac->resize(N, N);
      }

    // Compute residual and Jacobian contributions For each basis function i.
    for (unsigned int i=0; i<N; i++)
      {
        // The powers of x and y in the monomial we are currently integrating.
        unsigned int xpower = polys[i].first;
        unsigned int ypower = polys[i].second;

        // Residual contribution due to centroid point. The evaluation point is
        // fixed at (x,y) = 1/3.
        if (r && n_centroid)
          (*r)[i] += u[0] * pow(one_third, xpower) * pow(one_third, ypower);

        // Jacobian contribution due to centroid point. This is easy
        // because this term is linear in the weight parameter.
        if (jac && n_centroid)
          (*jac)(i,0) += pow(one_third, xpower) * pow(one_third, ypower);

        // Next loop Loop over all (w,x,y) triples in u and compute the residual
        // and Jacobian contributions.
        for (unsigned int q=n_centroid; q<u.size(); q+=3)
          {
            // The unknowns are ordered in terms of (w,x,y) triples.
            mpfr_class w = u[q];
            mpfr_class x = u[q+1];
            mpfr_class y = u[q+2];

            // The implied third barycentric coordinate
            mpfr_class z = mpfr_class(1) - x - y;

            // The "spatial" part is needed by both the residual and Jacobian,
            // so we can always compute it.
            mpfr_class spatial =
              pow(x, xpower) * pow(y, ypower) +
              pow(z, xpower) * pow(x, ypower) +
              pow(y, xpower) * pow(z, ypower);

            // Compute residual contribution, if required.
            if (r)
              (*r)[i] += w * spatial;

            // Compute Jacobian contribution, if required.
            if (jac)
              {
                // Derivative wrt w
                (*jac)(i, q) += spatial;

                // We are differentiating polynomials, so if xpower or
                // ypower is 0, the derivative of the corresponding term
                // will be zero. In that case we don't want to subtract
                // 1 from the unsigned variable which represents the
                // exponent, since that would cause it to wrap, possibly
                // leading to other issues. Therefore we will just
                // define the power to still be zero in that case...
                // Note that the -1 terms come from differentiating "z"
                // wrt either x or y.
                unsigned int xpm1 = xpower > 0 ? xpower-1 : 0;
                unsigned int ypm1 = ypower > 0 ? ypower-1 : 0;

                // Derivative wrt x
                (*jac)(i,q+1) += w * (
                                      xpower * pow(x, xpm1) * pow(y, ypower) +
                                      (mpfr_class(-1) * xpower * pow(z, xpm1) * pow(x, ypower) + pow(z, xpower) * ypower * pow(x, ypm1)) +
                                      pow(y, xpower) * mpfr_class(-1) * ypower * pow(z, ypm1));

                // Derivative wrt y
                (*jac)(i,q+2) += w * (
                                      pow(x, xpower) * ypower * pow(y, ypm1) +
                                      mpfr_class(-1) * xpower * pow(z, xpm1) * pow(x, ypower) +
                                      (xpower * pow(y, xpm1) * pow(z, ypower) + pow(y, xpower) * mpfr_class(-1) * ypower * pow(z, ypm1)));
              }
          }

        // Subtract off the true integral value, I(p_i)
        if (r)
          (*r)[i] -= exact_tri(xpower, ypower);
      }
  }
};



// Provides an operator() that returns true or false depending on
// whether the trial solution is "feasible", i.e. satisfies any
// required constraints/bounds which are separate from the nonlinear
// equations themselves.
struct CheckFeasibility
{
  CheckFeasibility(unsigned int d_in = 1) :
    d(d_in),
    N((d*d + 3*d + 6)/6),
    n_centroid(N % 3),
    n_ro3(N / 3)
  {
  }

  // polynomial degree
  unsigned int d;

  // Number of unknowns in the residual vector
  unsigned int N;

  // The number of centroid points
  unsigned int n_centroid;

  // The number of Ro3-invariant points
  unsigned int n_ro3;

  // Return true if the trial_u is feasible, false otherwise.
  bool operator() (const std::vector<mpfr_class> & trial_u)
  {
    if (N != trial_u.size())
      {
        std::cout << "Error, wrong size solution in CheckFeasibility." << std::endl;
        abort();
      }

    // 1.) No parameters can be negative (or exactly zero).
    // In fact, they probably shouldn't be too small, either!
    for (unsigned int q=0; q<trial_u.size(); ++q)
      if (trial_u[q] <= 1.e-4)
        {
          // Debugging:
          // std::cout << "trial_u is infeasible due to parameter <= 0!"
          //           << std::endl;
          // print(trial_u);

          return false;
        }

    // 2.) Check if 1-x-y <= 0
    for (unsigned int q=n_centroid; q<trial_u.size(); q+=3)
      {
        mpfr_class x = trial_u[q+1];
        mpfr_class y = trial_u[q+2];
        if (mpfr_class(1) - x - y <= 0)
          {
            // Debugging:
            // std::cout << "trial_u is infeasible due to 1-x-y <= 0!"
            //           << std::endl;
            // print(trial_u);

            return false;
          }
      }

    // If we made it here without returning, must be feasible!
    return true;
  }
};



// Parameters that control the behavior of solvers.
struct SolverData
{
  SolverData() :
    tol(1.e-36),
    divtol(1.e16),
    maxits(20),
    alphamin(1.e-3),
    do_backtracking(false),
    verbose(false) {}

  // Initial guess/solution vector
  std::vector<mpfr_class> u;
  mpfr_class tol;
  mpfr_class divtol;
  unsigned int maxits;
  // smallest backtracking backtracking linesearch fraction
  mpfr_class alphamin;
  bool do_backtracking;
  bool verbose;
  ResidualAndJacobian residual_and_jacobian;
  CheckFeasibility check_feasibility;
};

#endif
