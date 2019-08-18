// Functions for computing exact integrals of monomials
// on the reference triangle.
#include "exact.h"
#include "matrix.h"

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

// Given an initial guess u, uses Newton's method to see if it can obtain a solution.
// If the Newton iterations fail, this is reported back to the user.
bool newton(std::vector<mpfr_class> & u);

// Compute the l2-norm of vector r.
mpfr_class norm(const std::vector<mpfr_class> & r);

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
  // Seeds that produce initial guesses which fail to converge
  // (only the weights were chosen randomly).
  // 1566139255
  // 1566139317 // matrix is singular!
  // 1566139363
  // 1566139397
  // 1566139449
  // 1566139481
  // 1566139509
  // 1566139629
  // 1566139669
  // 1566139711
  // 1566139735
  // 1566139757
  // 1566139782
  // 1566139962
  // 1566142027 // backtracking seems like it *could* help this case, but doesn't.
  time_t seed = time(nullptr);
  std::cout << "seed=" << seed << std::endl;
  srandom(seed);

  // An initial guess which is known to converge to the solution in Gatermann's paper.
  // std::vector<mpfr_class> u =
  //   {
  //     2.65e-02, // w1
  //     6.23e-02, // x1
  //     6.75e-02, // y1
  //
  //     4.38e-02, // w2
  //     5.52e-02, // x2
  //     3.21e-01, // y2
  //
  //     2.87e-02, // w3
  //     3.43e-02, // x3
  //     6.60e-01, // y3
  //
  //     6.74e-02, // w4
  //     5.15e-01, // x4
  //     2.77e-01, // y4
  //   };

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
  // A "random" initial guess. Note that all parameters must be
  // positive, the sum of the weights cannot exceed the reference element
  // area. Also, each weight parameter corresponds to 3 quadrature
  // points, so actually it cannot exceed (1/3) * (1/2) = 1/6.
  // For the points, they must be chosen such that 1 - x - y > 0
  // and their relative order matters, i.e. reversing x and y will
  // give you a different solution.
  std::vector<mpfr_class> w_guess(4);
  std::vector<mpfr_class> x_guess(4);
  std::vector<mpfr_class> y_guess(4);

  // Latin hypercube sampling.
  // The total number of random numbers we need is n_dim = 11.  If the
  // total number of tests is n_tests, then we will need n_tests *
  // n_dim random numbers. We generate all the points at the same time
  // to ensure that they don't overlap.  1 / n_tests = "bin size", we
  // guarantee only one sample from each bin size over the entire set
  // of tests. Note also: n_tests = n_bins, so we only need one variable
  // for this. Note that we can't make the number of tests *too* big, since
  // we are building a dense matrix of random numbers here. So it may be
  // necessary to wrap this in an outer loop.
  unsigned int n_dim = 11;
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

      // Debugging: for a small n_tests, make sure that each random
      // number comes from a different bin.
      // std::cout << "Random numbers used for current test." << std::endl;
      // for (const auto & val : current_random)
      //   std::cout << val << " ";
      // std::cout << std::endl;

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
          std::cout << "Incorrect number of random number used!" << std::endl;
          abort();
        }

      // For the true solution, the sum of weight parameters is:
      // 2.65e-02 + 4.38e-02 + 2.87e-02 + 6.74e-02 = 0.1664,
      // so we should verify the same is true for our random initial
      // guess.
      // for (unsigned int i=0; i<w_guess.size(); ++i)
      //   std::cout << w_guess[i] << std::endl;
      // std::cout << "sum = "
      //           << std::accumulate(w_guess.begin(), w_guess.end(), mpfr_class(0))
      //           << std::endl;

      // Verify initial guess is OK.
      for (unsigned int i=0; i<w_guess.size(); ++i)
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
      std::vector<mpfr_class> u;
      u.reserve(12);
      for (unsigned int i=0; i<w_guess.size(); ++i)
        {
          u.push_back(w_guess[i]);
          u.push_back(x_guess[i]);
          u.push_back(y_guess[i]);
        }

      // Print initial guess
      // std::cout << "Initial guess=" << std::endl;
      // for (const auto & val : u)
      //   std::cout << val << std::endl;

      // We now use Newton iterations to obtain more digits of accuracy in the
      // points and weights.
      bool converged = newton(u);

      if (converged)
        {
          // We want to check whether this solution is a permutation
          // of the one that we already have. So we sort the weights
          // and compare them to the weights of the solution we
          // already have.
          std::vector<mpfr_class> weights = {u[0], u[3], u[6], u[9]};
          std::sort(weights.begin(), weights.end());

          const std::vector<mpfr_class> known_weights =
            {
              mpfr_class("2.6517028157436251428754180460739e-2"),
              mpfr_class("2.8775042784981585738445496900219e-2"),
              mpfr_class("4.3881408714446055036769903139288e-2"),
              mpfr_class("6.7493187009802774462697086166421e-2")
            };

          for (unsigned int i=0; i<weights.size(); ++i)
            weights[i] -= known_weights[i];

          mpfr_class norm_w = norm(weights);
          // std::cout << "norm_w=" << norm_w << std::endl;

          if (norm_w > mpfr_class(1.e-6))
            {
              std::cout << "Possible new solution u=" << std::endl;
              for (const auto & val : u)
                std::cout << val << std::endl;
            }
          else
            {
              // Keep track of the number of converged (but duplicate) solutions.
              n_converged++;
            }
        }
    } // end loop over n_tests
}



bool newton(std::vector<mpfr_class> & u)
{
  // Newton iteration parameters.
  const mpfr_class tol(1.e-36);
  const mpfr_class divtol(1.e16);
  const mpfr_class alphamin(1.e-3);
  const bool do_backtracking = false;
  unsigned iter = 0;
  const unsigned int maxits = 20;
  bool converged = false;

  // Storage for residual, Newton update, and Jacobian
  std::vector<mpfr_class> r, du;
  Matrix<mpfr_class> jac;

  // Store previous residual norm. Used for simplified line searching technique.
  mpfr_class old_residual_norm = 0.;

  while (true)
    {
      ++iter;
      if (iter > maxits)
        break;

      // Compute the residual (vector) at the current guess.
      residual_and_jacobian(&r, nullptr, u);

      // Check the norm of the residual vector to see if we are done.
      mpfr_class residual_norm = norm(r);

      // std::cout << "Iteration " << iter << ", residual_norm=" << residual_norm << std::endl;

      if (residual_norm < tol)
        {
          converged = true;
          break;
        }

      // If the residual is too large, just give up.
      if (residual_norm > divtol)
        break;

      // Compute Jacobian (only) at the current guess.
      residual_and_jacobian(nullptr, &jac, u);

      // Compute update: du = -jac^{-1} * r
      jac.lu_solve(du, r);

      // Compute next iterate, u -= alpha*du using simplified
      // backtracking, alpha_min < alpha <= 1.
      // This is only theoretically helpful at the moment, so far I have not encountered
      // an actual case where a failing case was able to converge by using backtracking.
      mpfr_class alpha(1);
      bool backtracking_converged = false;
      while (true)
        {
          if (alpha < alphamin)
            break;

          for (unsigned int i=0; i<u.size(); ++i)
            u[i] -= alpha * du[i];

          if (do_backtracking)
            {
              residual_and_jacobian(&r, nullptr, u);
              mpfr_class new_residual_norm = norm(r);

              std::cout << "Trying Newton step with alpha = " << alpha
                        << ", residual = " << new_residual_norm << std::endl;

              if (new_residual_norm < residual_norm)
                {
                  backtracking_converged = true;
                  break;
                }

              // Don't waste time backtracking if we already diverged
              if (new_residual_norm > divtol)
                break;

              alpha /= mpfr_class(2);
            }
          else
            {
              // If we aren't doing backtracking, then backtracking is "done".
              backtracking_converged = true;
              break;
            }
        }

      if (!backtracking_converged)
        break;
    } // end while

  // if (!converged)
  //   std::cout << "Newton iteration diverged, backtracking failed, or max iterations exceeded."
  //             << std::endl;

  return converged;
}


mpfr_class norm(const std::vector<mpfr_class> & r)
{
  mpfr_class residual_norm = 0.;
  for (unsigned int i=0; i<r.size(); ++i)
    residual_norm += r[i] * r[i];
  residual_norm = sqrt(residual_norm);
  return residual_norm;
}



void
residual_and_jacobian(std::vector<mpfr_class> * r,
                      Matrix<mpfr_class> * jac,
                      const std::vector<mpfr_class> & u)
{
  // If there's nothing to do, then there's nothing to do.
  if (!r && !jac)
    return;

  // The list of monomial exponents that we are going to test.
  std::vector<std::pair<unsigned int, unsigned int>> polys =
    {
      {0,0}, {2,0}, {3,0}, {2,1}, {4,0}, {5,0},
      {4,1}, {6,0}, {5,1}, {7,0}, {6,1}, {4,2}
    };

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

  // For each basis funciton i, we loop over all (w,x,y) triples in u and
  // compute the residual and Jacobian contributions, as requested.
  for (unsigned int i=0; i<n; i++)
    {
      // The powers of x and y in the monomial we are currently integrating.
      unsigned int xpower = polys[i].first;
      unsigned int ypower = polys[i].second;

      // Accumlate quadrature rule contribution for r(i).
      for (unsigned int q=0; q<u.size(); q+=3)
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
