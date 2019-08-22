// Functions for computing exact integrals of monomials
// on the reference triangle.
#include "exact.h"
#include "matrix.h"
#include "vect.h"
#include "newton.h"
#include "solver_data.h"

// C++ includes
#include <vector>
#include <utility>
#include <iostream>
#include <algorithm> // std::sort
#include <numeric> // std::accumulate
#include <stdlib.h> // random
#include <time.h> // time()

// Function used to (possibly) compute the residual and Jacobian at
// the input u.  If either "r" or "jac" is nullptr, their computation
// is skipped.
void residual_and_jacobian(std::vector<mpfr_class> * r,
                           Matrix<mpfr_class> * jac,
                           const std::vector<mpfr_class> & u);

// Generates a fixed-size, dense block of random samples using the
// Latin hypercube sampling method and checks whether any of them
// converge.
void generate_lhc_and_test();

// Our goal here is to determine whether there might be more than one
// solution to the SEVENTH-order Ro3-invariant quadrature rule
// originally reported in Gatermann, 1988.  We are going to try
// several different sets of starting data and see if any of them
// converge to a different solution than the one reported in the
// paper.
//
// For more information, see gatermann.C
int main()
{
  std::cout.precision(32);
  std::cout.setf(std::ios_base::scientific);

  // # of binary digits
  // 53 binary digits is about what you normally get with a double.
  mpfr_set_default_prec(256);

  // Varies the sequence of pseudorandom numbers returned by random().
  // srandom(1566142027);

  // Use the current time since epoch as a seed.
  time_t seed = time(nullptr);
  std::cout << "seed=" << seed << std::endl;
  srandom(seed);

  unsigned int testset_counter = 0;
  while (true)
    {
      ++testset_counter;
      std::cout << "Running test set " << testset_counter << std::endl;
      generate_lhc_and_test();
    }

  return 0;
}



void generate_lhc_and_test()
{
  // Set parameters to be used by solvers
  SolverData solver_data;
  solver_data.verbose = true;
  solver_data.residual_and_jacobian = residual_and_jacobian;

  // A "random" initial guess. Note that all parameters must be
  // positive, the sum of the weights cannot exceed the reference element
  // area. Also, each weight parameter corresponds to 3 quadrature
  // points, so actually it cannot exceed (1/3) * (1/2) = 1/6.
  // For the points, they must be chosen such that 1 - x - y > 0
  // and their relative order matters, i.e. reversing x and y will
  // give you a different solution.

  // The degree 6 Ro3-invariant rule has 10 unknowns:
  // 1 centroid weight
  // 3 coronas with 3 points each
  std::vector<mpfr_class> w_guess(4);
  std::vector<mpfr_class> x_guess(3);
  std::vector<mpfr_class> y_guess(3);

  // Latin hypercube sampling. n_dim is the number of random values that
  // are needed per initial guess.
  unsigned int n_dim = (w_guess.size() - 1) + x_guess.size() + y_guess.size();
  unsigned int n_tests = 100000;
  double bin_size = 1. / double(n_tests);

  // Status message
  std::cout << "Generating random initial guesses..." << std::endl;

  // random_numbers will be n_dim * n_tests in size. Each *column*
  // will therefore represent an n_dim-dimensional set of test
  // parameters. We will random_shuffle() the rows so that we get one
  // sample from each "bin" of a given random variable.
  std::vector<std::vector<double>> random_numbers(n_dim);
  for (unsigned int i=0; i<random_numbers.size(); ++i)
    {
      random_numbers[i].resize(n_tests);

      // Get random number in the current bin.
      for (unsigned int j=0; j<random_numbers[i].size(); ++j)
        {
          // Lower bound of current bin
          double lower = j * bin_size;

          random_numbers[i][j] = bin_size * double(random()) / RAND_MAX + lower;
        }
    }

  // random_shuffle all rows after the first one. This way we don't have
  // a sample from the *same* bin for each random variable.
  for (unsigned int i=0; i<random_numbers.size(); ++i)
    std::random_shuffle(random_numbers[i].begin(), random_numbers[i].end());

  // Print with fewer digits
  // auto flags = std::cout.flags();
  // std::cout.precision(8);
  // for (unsigned int i=0; i<random_numbers.size(); ++i)
  //   {
  //     for (unsigned int j=0; j<random_numbers[i].size(); ++j)
  //       std::cout << random_numbers[i][j] << " ";
  //     std::cout << std::endl;
  //   }
  // std::cout.flags(flags);

  // Status message
  std::cout << "Testing initial guesses..." << std::endl;

  // Counter to keep track of the number of initial guesses tested.
  unsigned int guess = 0;
  // unsigned int n_converged = 0;
  // unsigned int n_converged_last_message = 1;

  for (unsigned int t=0; t<n_tests; ++t)
    {
      ++guess;

      if (guess == n_tests/2)
        std::cout << "Half way done..." << std::endl;

      // if (n_converged > n_converged_last_message && (n_converged % 10 == 0))
      //   {
      //     std::cout << "Found " << n_converged << " duplicate converged solutions." << std::endl;
      //     n_converged_last_message = n_converged;
      //   }

      std::vector<double> current_random;
      current_random.reserve(n_dim);
      for (unsigned int j=0; j<n_dim; ++j)
        current_random.push_back(random_numbers[j][t]);

      // Every time we need a random number, post-increment "ctr".
      unsigned int ctr = 0;

      for (unsigned int i=0; i<w_guess.size(); ++i)
        {
          w_guess[i] = mpfr_class(1) / 6;

          for (unsigned int j=0; j<i; ++j)
            w_guess[i] -= w_guess[j];

          // The last weight doesn't get randomized, it's whatever is left
          // after choosing the first 3.
          if (i < w_guess.size() - 1)
            {
              // std::cout << "ctr=" << ctr << std::endl;
              w_guess[i] *= current_random[ctr++];
            }
        }


      // Set random initial guesses for x and y:
      for (auto & val : x_guess)
        {
          // std::cout << "ctr=" << ctr << std::endl;
          val = current_random[ctr++];
        }
      for (unsigned int i=0; i<y_guess.size(); ++i)
        {
          // std::cout << "ctr=" << ctr << std::endl;
          y_guess[i] = (1. - x_guess[i]) * current_random[ctr++];
        }

      // std::cout << "Done setting initial guesses." << std::endl;

      // We are done, so make sure we used all the random numbers.
      if (ctr != n_dim)
        {
          std::cout << "Incorrect number of random number used!" << std::endl;
          abort();
        }

      // Verify barycentric points are inside the element.
      for (unsigned int i=0; i<x_guess.size(); ++i)
        {
          // std::cout << "(x" << i << ", y" << i << ")="
          //           << "(" << x_guess[i] << ", " << y_guess[i] << ")" << std::endl;

          mpfr_class barycentric_point = mpfr_class(1) - x_guess[i] - y_guess[i];
          if (barycentric_point < 0)
            {
              std::cout << "Error, invalid initial guess!" << std::endl;
              std::abort();
            }
        }

      // Set initial guess vector
      std::vector<mpfr_class> & u = solver_data.u;
      u.clear();
      u.reserve(w_guess.size() + x_guess.size() + y_guess.size());
      u.push_back(w_guess[0]); // centroid weight is first dof
      for (unsigned int i=0; i<x_guess.size(); ++i)
        {
          u.push_back(w_guess[i+1]);
          u.push_back(x_guess[i]);
          u.push_back(y_guess[i]);
        }

      // Print initial guess
      // std::cout << "Initial guess=" << std::endl;
      // print(u);

      bool converged = false;

      // Attempt to improve initial guess with a Newton minimization iterations.
      converged = newton_min(solver_data);

      // We now use Newton iterations to obtain more digits of accuracy in the
      // points and weights.
      converged = newton(solver_data);

      if (converged)
        {
          std::cout << "Possible new solution u=" << std::endl;
          print(u);
        }
    } // end loop over n_tests
}



void
residual_and_jacobian(std::vector<mpfr_class> * r,
                      Matrix<mpfr_class> * jac,
                      const std::vector<mpfr_class> & u)
{
  // If there's nothing to do, then there's nothing to do.
  if (!r && !jac)
    return;

  // The list of the 10 monomials we must integrate exactly.
  std::vector<std::pair<unsigned int, unsigned int>> polys =
    {
      {0,0}, {2,0}, {3,0}, {2,1}, {4,0},
      {5,0}, {4,1}, {6,0}, {5,1}, {4,2}
    };

  // A handy constant
  const mpfr_class one_third = mpfr_class(1) / 3;

  // The number of equations (and the number of unknowns).
  unsigned int n = polys.size();

  // Zero any previous values and allocate space for residual and Jacobian, if required.
  if (r)
    {
      r->clear();
      r->resize(n);
    }
  if (jac)
    {
      jac->clear();
      jac->resize(n, n);
    }

  // Compute residual and Jacobian contributions For each basis function i.
  for (unsigned int i=0; i<n; i++)
    {
      // The powers of x and y in the monomial we are currently integrating.
      unsigned int xpower = polys[i].first;
      unsigned int ypower = polys[i].second;

      // Residual contribution due to centroid point. The evaluation point is
      // fixed at (x,y) = 1/3.
      if (r)
        (*r)[i] += u[0] * pow(one_third, xpower) * pow(one_third, ypower);

      // Jacobian contribution due to centroid point. This is easy
      // because this term is linear in the weight parameter.
      if (jac)
        (*jac)(i,0) += pow(one_third, xpower) * pow(one_third, ypower);

      // Next loop Loop over all (w,x,y) triples in u and compute the residual
      // and Jacobian contributions.
      for (unsigned int q=1; q<u.size(); q+=3)
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
