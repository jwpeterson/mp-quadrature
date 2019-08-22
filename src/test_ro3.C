#include "test_ro3.h"
#include "solver_data.h"
#include "lhc.h"
#include "gradient_descent.h"
#include "nlcg.h"
#include "newton.h"
#include "vect.h"

void test_ro3(unsigned int d,
              SolverData & solver_data)
{
  // The number of unknowns
  unsigned int N = (d*d + 3*d + 6)/6;
  // The number of centroid points (0 or 1)
  unsigned int n_centroid = N % 3;
  // The number of Ro3 points
  unsigned int n_ro3 = N / 3;

  std::cout << "N = " << N << std::endl;
  std::cout << "n_centroid = " << n_centroid << std::endl;
  std::cout << "n_ro3 = " << n_ro3 << std::endl;

  if (n_centroid > 1)
    {
      std::cout << "Error, can't have a rule with more than one centroid point." << std::endl;
      abort();
    }

  // A "random" initial guess. Note that all parameters must be
  // positive, the sum of the weights cannot exceed the reference element
  // area. Also, each weight parameter corresponds to 3 quadrature
  // points, so actually it cannot exceed (1/3) * (1/2) = 1/6.
  // For the points, they must be chosen such that 1 - x - y > 0
  // and their relative order matters, i.e. reversing x and y will
  // give you a different solution.
  std::vector<mpfr_class> w_guess(n_centroid + n_ro3);
  std::vector<mpfr_class> x_guess(n_ro3);
  std::vector<mpfr_class> y_guess(n_ro3);

  // The total number of random numbers we need is 1 less than the
  // number of parameters we are solving for, since we choose the
  // weights so that they sum up to the reference element area.
  unsigned int n_dim = n_centroid + (3 * n_ro3) - 1;
  unsigned int n_tests = 1;

  // Status message
  std::cout << "Generating random initial guesses..." << std::endl;

  // Latin hypercube sampling.
  // random_numbers will be n_dim * n_tests in size. Each *column*
  // will therefore represent an n_dim-dimensional set of test
  // parameters. We will random_shuffle() the rows so that we get one
  // sample from each "bin" of a given random variable.
  std::vector<std::vector<double>> random_numbers =
    lhc(n_dim, n_tests);

  // Status message
  std::cout << "Testing initial guesses..." << std::endl;

  // Counter to keep track of the number of initial guesses tested.
  unsigned int guess = 0;
  unsigned int n_converged = 0;
  unsigned int n_converged_last_message = 1;

  for (unsigned int t=0; t<n_tests; ++t)
    {
      ++guess;

      if (guess == n_tests/2)
        std::cout << "Half way done..." << std::endl;

      if (n_converged > n_converged_last_message && (n_converged % 10 == 0))
        {
          std::cout << "Found " << n_converged << " duplicate converged solutions." << std::endl;
          n_converged_last_message = n_converged;
        }

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
            w_guess[i] *= current_random[ctr++];
        }

      // Set random initial guesses for x and y:
      for (auto & val : x_guess)
        val = current_random[ctr++];
      for (unsigned int i=0; i<y_guess.size(); ++i)
        y_guess[i] = (1. - x_guess[i]) * current_random[ctr++];

      // We are done, so make sure we used all the random numbers.
      if (ctr != n_dim)
        {
          std::cout << "Incorrect number of random numbers used!" << std::endl;
          abort();
        }

      // // Verify initial guess is OK.
      // for (unsigned int i=0; i<w_guess.size(); ++i)
      //   {
      //     // std::cout << "(x" << i << ", y" << i << ")="
      //     //           << "(" << x_guess[i] << ", " << y_guess[i] << ")" << std::endl;
      //
      //     mpfr_class barycentric_point = mpfr_class(1) - x_guess[i] - y_guess[i];
      //     if (barycentric_point < 0)
      //       {
      //         std::cout << "Error, invalid initial guess!" << std::endl;
      //         std::abort();
      //       }
      //   }

      // Set initial guess vector
      std::vector<mpfr_class> & u = solver_data.u;
      u.clear();
      u.reserve(n_dim + 1);

      // If there's a centroid point, it's weight appears first in the
      // vector of unknowns.
      if (n_centroid)
        u.push_back(w_guess[0]);

      for (unsigned int i=0; i<x_guess.size(); ++i)
        {
          u.push_back(w_guess[i + n_centroid]);
          u.push_back(x_guess[i]);
          u.push_back(y_guess[i]);
        }

      // Debugging: set initial guess to known solution for 12 point, degree=7 case.
//      u =
//        {
//          2.65e-02, // w1
//          6.23e-02, // x1
//          6.75e-02, // y1
//
//          4.38e-02, // w2
//          5.52e-02, // x2
//          3.21e-01, // y2
//
//          2.87e-02, // w3
//          3.43e-02, // x3
//          6.60e-01, // y3
//
//          6.74e-02, // w4
//          5.15e-01, // x4
//          2.77e-01, // y4
//        };

      // Print initial guess
      // std::cout << "Initial guess=" << std::endl;
      // print(u);

      // Status flag returned by different methods.
      bool converged = false;

      // Try a few iterations of one of the minimization routines.
      // converged = newton_min(solver_data);
      converged = gradient_descent(solver_data);
      // converged = nlcg(solver_data);

      // Follow up initial guess minimization with Newton. For a
      // random initial guess, it might be worthwhile seeing if we can
      // improve it slightly with an optimization algorithm before
      // switching to Newton's method?
      converged = newton(solver_data);

      if (converged)
        {
          // Print converged solution
          std::cout << "Solution converged!" << std::endl;
          std::cout << "u=" << std::endl;
          print(u);

          // Keep track of the number of converged solutions.
          n_converged++;
        }
    } // end loop over n_tests
}
